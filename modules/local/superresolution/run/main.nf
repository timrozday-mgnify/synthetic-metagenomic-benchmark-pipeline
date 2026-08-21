// Run a superresolution pipeline (amplicon or shotgun) on the generated reads via
// a nested `nextflow run`, then normalise its composition CSV into the same
// truth-comparable profile sylph emits.
//
// Both pipelines share one interface — a YAML samplesheet of {id, reads,
// platform, references}, a single combined reference FASTA (no external DB), and
// `composition/<id>/<id>.inferred_composition.csv` as output — so `meta.profiler`
// only selects which repo to launch.
//
// Runs on the host (executor 'local', no container) so it reuses the host nextflow
// + container engine, exactly like RUN_AAP. These are lightweight nested-Nextflow
// launchers; the nested run manages its own images, executor, and heavy-job resources.
// The nested run does NOT inherit the outer -profile: set params.sr_profile (and
// params.sr_configs for extra -c files).
process BUILD_SUPERRESOLUTION_MISMAPPING {
    tag "${meta.reference_set} (${meta.profiler})"
    label 'process_single'
    executor 'local'

    input:
    tuple val(meta), val(read_paths), path(refs)

    output:
    tuple val(meta), path("${meta.id}.mismapping_matrix.csv"), emit: mismapping
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def repo = meta.profiler == 'sr_amplicon' ? params.sr_amplicon_repo : params.sr_shotgun_repo
    if (!repo) error "BUILD_SUPERRESOLUTION_MISMAPPING: params.${meta.profiler == 'sr_amplicon' ? 'sr_amplicon_repo' : 'sr_shotgun_repo'} is not set (reference set ${meta.reference_set})"
    def revArg = (repo.startsWith('/') || repo.startsWith('.')) ? '' : "-r ${params.sr_revision}"
    def profArg = meta.sr_profile ? "-profile ${meta.sr_profile}" : ''
    def extraCfg = (meta.sr_configs ?: []).collect { "-c ${file(it, checkIfExists: true)}" }.join(' ')
    def nestedArgs = [revArg, profArg, '--input sr_samplesheet.yml', '--outdir sr_out', extraCfg]
        .findAll { it }
        .join(' ')
    def platform = meta.platform ? "printf '  platform: %s\\n' '${meta.platform}' >> sr_samplesheet.yml" : 'true'
    def reads = (read_paths instanceof List ? read_paths : [read_paths]).collect { it.toString() }
    assert reads.every { it } : "BUILD_SUPERRESOLUTION_MISMAPPING: empty read path for ${meta.reference_set}"
    def readCmds = reads.collect { "printf '    - %s\\n' '${it}' >> sr_samplesheet.yml" }.join('\n    ')
    """
    export NXF_ASSETS="\$PWD/.nxf_assets"

    printf -- '- id: %s\\n' '${meta.id}' > sr_samplesheet.yml
    printf '  reads:\\n' >> sr_samplesheet.yml
    ${readCmds}
    ${platform}
    printf '  references: %s\\n' "\$(realpath ${refs})" >> sr_samplesheet.yml

    nextflow run ${repo} \\
        ${nestedArgs}

    # Current superresolution-amplicon and superresolution-shotgun revisions publish
    # generated matrices in an opaque-key bundle under mismapping/, rather than next
    # to a sample's composition. The representative nested run must yield one bundle.
    find sr_out/mismapping -type f -name mismapping_matrix.csv -print 2>/dev/null | sort > matrix_paths.txt
    matrix_count=\$(wc -l < matrix_paths.txt | tr -d ' ')
    if [ "\$matrix_count" -ne 1 ]; then
        echo "Expected exactly one nested superresolution mismapping matrix under sr_out/mismapping/, found \$matrix_count." >&2
        if [ "\$matrix_count" -gt 0 ]; then
            sed 's/^/Discovered matrix: /' matrix_paths.txt >&2
        else
            echo 'Files produced below sr_out/:' >&2
            find sr_out -type f -print 2>/dev/null | sort >&2 || true
        fi
        exit 1
    fi
    matrix_path=\$(sed -n '1p' matrix_paths.txt)
    cp "\$matrix_path" ${meta.id}.mismapping_matrix.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${repo}: ${params.sr_revision}
        nextflow: \$(nextflow -version 2>&1 | grep -oE 'version [0-9.]+' | sed 's/version //')
    END_VERSIONS
    """

    stub:
    def repo = meta.profiler == 'sr_amplicon' ? params.sr_amplicon_repo : params.sr_shotgun_repo
    """
    echo 'src,dst,prob' > ${meta.id}.mismapping_matrix.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${repo}: stub
        nextflow: stub
    END_VERSIONS
    """
}

process RUN_SUPERRESOLUTION {
    tag "${meta.id} (${meta.profiler})"
    label 'process_single'
    executor 'local'

    input:
    // Reads are NOT staged: they're passed as absolute host paths (val), matching
    // RUN_AAP — this process is executor 'local' so the nested run reads them
    // directly, and the superresolution main.nf resolves relative paths against its
    // OWN projectDir, which a staged basename would break.
    tuple val(meta), val(read_paths), path(refs), path(mismapping_matrix)

    output:
    tuple val(meta), path("${meta.id}.sr_profile.tsv"),          emit: profile
    tuple val(meta), path("sr_out/composition/${meta.id}/*"),    emit: composition
    path "versions.yml",                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def repo = meta.profiler == 'sr_amplicon' ? params.sr_amplicon_repo : params.sr_shotgun_repo
    if (!repo) error "RUN_SUPERRESOLUTION: params.${meta.profiler == 'sr_amplicon' ? 'sr_amplicon_repo' : 'sr_shotgun_repo'} is not set (sample ${meta.id})"
    // -r only applies to a Nextflow project name; a local checkout path takes none.
    def rev_arg   = (repo.startsWith('/') || repo.startsWith('.')) ? '' : "-r ${params.sr_revision}"
    def prof_arg  = meta.sr_profile ? "-profile ${meta.sr_profile}" : ''
    def extra_cfg = (meta.sr_configs ?: []).collect { "-c ${file(it, checkIfExists: true)}" }.join(' ')
    // Pass presence-gate settings explicitly rather than via sr_configs: the nested
    // amplicon and shotgun pipelines have different defaults and can be swept
    // independently in one benchmark invocation.
    def srKind = meta.profiler == 'sr_amplicon' ? 'amplicon' : 'shotgun'
    def presence = params["sr_${srKind}_infer_presence"]
    def presencePrior = params["sr_${srKind}_infer_presence_prior"]
    def presenceTemp = params["sr_${srKind}_infer_presence_temp"]
    def presenceArgs = []
    if (presence != null) presenceArgs << "--infer_presence ${presence}"
    if (presencePrior != null) presenceArgs << "--infer_presence_prior ${presencePrior}"
    if (presenceTemp != null) presenceArgs << "--infer_presence_temp ${presenceTemp}"
    def presenceArg = presenceArgs.join(' ')
    def nestedArgs = [rev_arg, prof_arg, '--input sr_samplesheet.yml', '--outdir sr_out', extra_cfg,
                      "--mismapping_matrix ${mismapping_matrix}", presenceArg]
        .findAll { it }
        .join(' ')
    def platform  = meta.platform ? "printf '  platform: %s\\n' '${meta.platform}' >> sr_samplesheet.yml" : 'true'
    def reads     = (read_paths instanceof List ? read_paths : [read_paths]).collect { it.toString() }
    assert reads.every { it } : "RUN_SUPERRESOLUTION: empty read path for ${meta.id}"
    def read_cmds = reads.collect { "printf '    - %s\\n' '${it}' >> sr_samplesheet.yml" }.join('\n    ')
    """
    # Nextflow pulls a project into a SHARED asset cache (\$NXF_HOME/assets). One task
    # per sample means several nested runs pull at once on a cold cache and one reads
    # another's half-written clone: "Repository may be corrupted" — which then sticks
    # until `nextflow drop`. A task-local asset dir removes the shared state, and as a
    # bonus makes each run fetch the requested revision rather than a stale cached one.
    # ponytail: costs one small clone per task; pre-warm a shared cache instead only if
    # that ever shows up in the runtime.
    export NXF_ASSETS="\$PWD/.nxf_assets"

    # One-sample YAML samplesheet; `references` must be absolute (the nested pipeline
    # resolves relative paths against its own projectDir).
    printf -- '- id: %s\\n' '${meta.id}' > sr_samplesheet.yml
    printf '  reads:\\n' >> sr_samplesheet.yml
    ${read_cmds}
    ${platform}
    printf '  references: %s\\n' "\$(realpath ${refs})" >> sr_samplesheet.yml

    nextflow run ${repo} \\
        ${nestedArgs}

    normalize_sr_profile.py \\
        --composition sr_out/composition/${meta.id}/${meta.id}.inferred_composition.csv \\
        --output ${meta.id}.sr_profile.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${repo}: ${params.sr_revision}
        nextflow: \$(nextflow -version 2>&1 | grep -oE 'version [0-9.]+' | sed 's/version //')
    END_VERSIONS
    """

    stub:
    def repo = meta.profiler == 'sr_amplicon' ? params.sr_amplicon_repo : params.sr_shotgun_repo
    """
    mkdir -p sr_out/composition/${meta.id}
    printf 'sample,genome_id,observed_rel_abundance,inferred_mean,inferred_lo,inferred_hi\\n' \\
        > sr_out/composition/${meta.id}/${meta.id}.inferred_composition.csv
    printf 'genome_id\\tpredicted_rel_abundance\\tpredicted_tax_rel_abundance\\n' > ${meta.id}.sr_profile.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${repo}: stub
        nextflow: stub
    END_VERSIONS
    """
}
