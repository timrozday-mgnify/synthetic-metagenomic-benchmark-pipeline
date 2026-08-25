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
// Each nested run keeps its work dir and cache OUTSIDE the task directory (under
// <workDir>/nested/sr/<key>) and is launched with -resume, so a retried or re-run task
// resumes the nested pipeline rather than repeating it. The key is the reference set
// (matrix build) or the reference set + a digest of the batch's sample ids (inference):
// stable across outer runs, and never shared by two tasks that could run at once.
// Wipe <workDir>/nested to start clean.
// Inference is BATCHED: every sample sharing a reference set goes into one nested run's
// multi-row samplesheet, exactly as RUN_AAP batches by DB config. The reference set is
// the largest safe batch — the matrix and primer pair are per-run CLI flags, not
// per-row samplesheet fields.
// Pull each nested pipeline ONCE per run and make its helpers noexec-safe there.
// Every SR task used to pull for itself (a task-local asset dir, to dodge the shared
// asset cache's concurrent-clone corruption); with any real fan-out that hammers the
// GitHub API into 504s. One producer task per repo keeps the isolation and drops the
// call count to one, and consumers stage the result rather than fetching it.
// The nested amplicon pipeline extracts its reference amplicons with its OWN in-silico
// PCR, defaulting to 515F/806R (V4). A sample amplified with any other pair must be told
// which region to cut, or its reads and the reference amplicons cover different parts of
// the 16S and NOT ONE READ hits a reference. meta.primer_sets/meta.primer carry the pair
// the reads were generated with (--step all); a profile-only row has neither and falls
// back to the nested default.
def srPrimerArgs(meta) {
    if (meta.profiler != 'sr_amplicon') return ''
    def pair = (meta.primer_sets ?: []).find { it[0] == meta.primer }
    pair ? "--fwd_primer ${pair[1]} --rev_primer ${pair[2]}" : ''
}

process SR_PULL_REPO {
    tag "${profiler}"
    label 'process_single'
    executor 'local'

    input:
    val profiler

    output:
    tuple val(profiler), path("sr_assets"), emit: assets
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def repo = profiler == 'sr_amplicon' ? params.sr_amplicon_repo : params.sr_shotgun_repo
    if (!repo) error "SR_PULL_REPO: params.${profiler == 'sr_amplicon' ? 'sr_amplicon_repo' : 'sr_shotgun_repo'} is not set"
    def remoteRepo = !(repo.startsWith('/') || repo.startsWith('.'))
    def revArg = remoteRepo ? "-r ${params.sr_revision}" : ''
    """
    # A supplied local checkout is launched in place and left untouched (it can carry
    # its own equivalent patch), so this stays an empty marker for the channel join.
    mkdir -p sr_assets
    export NXF_ASSETS="\$PWD/sr_assets"

    if [ '${remoteRepo}' = 'true' ]; then
        # GitHub's API intermittently 504s; retry rather than fail the whole run.
        for attempt in 1 2 3; do
            nextflow pull ${repo} ${revArg} && break
            [ "\$attempt" = 3 ] && exit 1
            sleep \$((attempt * 30))
        done

        python "\$(command -v patch_sr_helpers.py)" "\$NXF_ASSETS" --repo ${repo}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${repo}: ${remoteRepo ? params.sr_revision : 'local'}
        nextflow: \$(nextflow -version 2>&1 | grep -oE 'version [0-9.]+' | sed 's/version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p sr_assets

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repo: stub
    END_VERSIONS
    """
}

process BUILD_SUPERRESOLUTION_MISMAPPING {
    tag "${meta.reference_set} (${meta.profiler})"
    label 'process_single'
    executor 'local'

    input:
    tuple val(meta), val(read_paths), path(refs), path(sr_assets)

    output:
    tuple val(meta), path("${meta.id}.mismapping_matrix.csv"), emit: mismapping
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def repo = meta.profiler == 'sr_amplicon' ? params.sr_amplicon_repo : params.sr_shotgun_repo
    if (!repo) error "BUILD_SUPERRESOLUTION_MISMAPPING: params.${meta.profiler == 'sr_amplicon' ? 'sr_amplicon_repo' : 'sr_shotgun_repo'} is not set (reference set ${meta.reference_set})"
    def remoteRepo = !(repo.startsWith('/') || repo.startsWith('.'))
    def revArg = remoteRepo ? "-r ${params.sr_revision}" : ''
    // A named repository is pulled into this task's isolated asset directory before
    // launch, so the non-executable upstream helpers can be fixed there. A supplied
    // local checkout remains untouched and can carry its own equivalent patch.
    // The asset layout is Nextflow-version-dependent, so the patch helper resolves
    // and prints the checkout rather than this building the path.
    def launchRepo = remoteRepo ? '"\$launch_repo"' : repo
    def profArg = meta.sr_profile ? "-profile ${meta.sr_profile}" : ''
    def extraCfg = (meta.sr_configs ?: []).collect { "-c ${file(it, checkIfExists: true)}" }.join(' ')
    def nestedDir = "${workflow.workDir}/nested/sr/${meta.id.replaceAll(/[^A-Za-z0-9._-]+/, '_')}"
    def nestedArgs = [profArg, '--input sr_samplesheet.yml', '--outdir sr_out', extraCfg,
                      srPrimerArgs(meta), "-w '${nestedDir}/work'", '-resume']
        .findAll { it }
        .join(' ')
    def platform = meta.platform ? "printf '  platform: %s\\n' '${meta.platform}' >> sr_samplesheet.yml" : 'true'
    def reads = (read_paths instanceof List ? read_paths : [read_paths]).collect { it.toString() }
    assert reads.every { it } : "BUILD_SUPERRESOLUTION_MISMAPPING: empty read path for ${meta.reference_set}"
    def readCmds = reads.collect { "printf '    - %s\\n' '${it}' >> sr_samplesheet.yml" }.join('\n    ')
    """
    printf -- '- id: %s\\n' '${meta.id}' > sr_samplesheet.yml
    printf '  reads:\\n' >> sr_samplesheet.yml
    ${readCmds}
    ${platform}
    printf '  references: %s\\n' "\$(realpath ${refs})" >> sr_samplesheet.yml

    if [ '${remoteRepo}' = 'true' ]; then
        # SR_PULL_REPO already pulled and patched this checkout; the helper only
        # resolves its version-dependent location here (it rewrites nothing when the
        # patch is already applied, so concurrent consumers never write to the
        # shared clone).
        launch_repo=\$(python "\$(command -v patch_sr_helpers.py)" sr_assets --repo ${repo})
    fi

    # NXF_CACHE_DIR moves the nested run's .nextflow cache + history out of this task
    # directory, which -resume needs to find its previous session. The launch dir stays
    # the task dir, so the relative --input/--outdir paths above keep working.
    export NXF_CACHE_DIR='${nestedDir}/cache'
    mkdir -p "\$NXF_CACHE_DIR"

    # The nested run also infers composition for the representative sample, which can
    # fail for reasons that do not affect the matrix (e.g. a sample whose reads hit no
    # reference amplicon). The matrix is this task's only deliverable, so a nested
    # failure is fatal only if it left no matrix behind.
    set +e
    nextflow run ${launchRepo} \\
        ${nestedArgs}
    nested_status=\$?
    set -e

    # Current superresolution-amplicon and superresolution-shotgun revisions publish
    # generated matrices in an opaque-key bundle under mismapping/, rather than next
    # to a sample's composition. The representative nested run must yield one bundle.
    find sr_out/mismapping -type f -name mismapping_matrix.csv -print 2>/dev/null | sort > matrix_paths.txt
    matrix_count=\$(wc -l < matrix_paths.txt | tr -d ' ')

    # A nested failure downstream of the matrix can leave nothing published under
    # sr_out/. The matrix itself is still in the nested work dir — in several copies,
    # since Nextflow stages it into every consumer task — so fall back to those when
    # they are all byte-identical (a differing set means two producer tasks disagreed,
    # which is not something to guess at).
    if [ "\$matrix_count" -eq 0 ] && [ "\$nested_status" -ne 0 ]; then
        find '${nestedDir}/work' -type f -name mismapping_matrix.csv -print 2>/dev/null | sort > work_matrices.txt
        distinct=\$(while read -r m; do cksum "\$m"; done < work_matrices.txt | awk '{print \$1, \$2}' | sort -u | wc -l | tr -d ' ')
        if [ -s work_matrices.txt ] && [ "\$distinct" -eq 1 ]; then
            echo "WARN: nested run exited \$nested_status without publishing a matrix; recovering it from the nested work dir." >&2
            sed -n '1p' work_matrices.txt > matrix_paths.txt
            matrix_count=1
        fi
    fi

    if [ "\$matrix_count" -ne 1 ]; then
        echo "Expected exactly one nested superresolution mismapping matrix under sr_out/mismapping/, found \$matrix_count (nested exit status \$nested_status)." >&2
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
    if [ "\$nested_status" -ne 0 ]; then
        echo "WARN: nested superresolution run exited \$nested_status but produced the mismapping matrix; continuing." >&2
    fi

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
    tag "${metas[0].reference_set ?: metas[0].database} (${metas[0].profiler}, ${metas.size()})"
    label 'process_single'
    executor 'local'

    input:
    // One batch per reference set — the nested samplesheet is a multi-row YAML list, and
    // everything the command line fixes (matrix, primer pair, repo, presence params) is
    // constant across the set. layout is one [id, platform, [read paths]] per sample.
    // Reads are NOT staged: they're passed as absolute host paths (val), matching
    // RUN_AAP — this process is executor 'local' so the nested run reads them
    // directly, and the superresolution main.nf resolves relative paths against its
    // OWN projectDir, which a staged basename would break.
    tuple val(metas), val(layout), path(refs), path(mismapping_matrix), path(sr_assets)

    output:
    tuple val(metas), path("*.sr_profile.tsv"),      emit: profile
    tuple val(metas), path("sr_out/composition/**"), emit: composition
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Group-invariant: the reference set keys on profiler and primer, and sr_profile /
    // sr_configs / the presence params are run-global. So metas[0] speaks for the batch.
    def meta = metas[0]
    def repo = meta.profiler == 'sr_amplicon' ? params.sr_amplicon_repo : params.sr_shotgun_repo
    if (!repo) error "RUN_SUPERRESOLUTION: params.${meta.profiler == 'sr_amplicon' ? 'sr_amplicon_repo' : 'sr_shotgun_repo'} is not set (samples ${metas*.id})"
    // -r only applies to a Nextflow project name; a local checkout path takes none.
    def remoteRepo = !(repo.startsWith('/') || repo.startsWith('.'))
    def rev_arg   = remoteRepo ? "-r ${params.sr_revision}" : ''
    // As in BUILD_SUPERRESOLUTION_MISMAPPING: a named repo is pulled into this task's
    // own asset dir first so its helpers can be made noexec-safe there, and the run is
    // then launched from that resolved directory (which takes no -r).
    def launchRepo = remoteRepo ? '"\$launch_repo"' : repo
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
    // Persistent home for the nested run's cache + work dir, keyed by the batch identity
    // (reference set + the exact set of samples): stable across outer runs and distinct
    // between concurrently running batches — two runs must never share one nested cache.
    def batch_id   = java.security.MessageDigest.getInstance('MD5')
                         .digest(metas*.id.sort().join(',').bytes).encodeHex().toString()[0..7]
    def set_dir    = (meta.reference_set ?: meta.id).replaceAll(/[^A-Za-z0-9._-]+/, '_')
    def nestedDir  = "${workflow.workDir}/nested/sr/${set_dir}-${batch_id}"
    def nestedArgs = [prof_arg, '--input sr_samplesheet.yml', '--outdir sr_out', extra_cfg,
                      "--mismapping_matrix ${mismapping_matrix}", presenceArg, srPrimerArgs(meta),
                      "-w '${nestedDir}/work'", '-resume']
        .findAll { it }
        .join(' ')
    assert layout.every { it[2] } : "RUN_SUPERRESOLUTION: empty read paths in layout for ${metas*.id}"
    def sheet_cmds = layout.collect { id, platform, reads ->
        ([ "printf -- '- id: %s\\n' '${id}' >> sr_samplesheet.yml",
           "printf '  reads:\\n' >> sr_samplesheet.yml" ] +
         reads.collect { "printf '    - %s\\n' '${it}' >> sr_samplesheet.yml" } +
         (platform ? [ "printf '  platform: %s\\n' '${platform}' >> sr_samplesheet.yml" ] : []) +
         [ "printf '  references: %s\\n' \"\$refs_abs\" >> sr_samplesheet.yml" ]).join('\n    ')
    }.join('\n    ')
    def norm_cmds = metas*.id.collect { id ->
        """comp=sr_out/composition/${id}/${id}.inferred_composition.csv
    [ -f "\$comp" ] || { echo "RUN_SUPERRESOLUTION: nested run produced no composition for ${id}" >&2; exit 1; }
    python "\$(command -v normalize_sr_profile.py)" --composition "\$comp" --output ${id}.sr_profile.tsv"""
    }.join('\n    ')
    """
    # Multi-sample YAML samplesheet, one entry per batched sample; `references` must be
    # absolute (the nested pipeline resolves relative paths against its own projectDir)
    # and is shared by the whole reference set.
    refs_abs=\$(realpath ${refs})
    : > sr_samplesheet.yml
    ${sheet_cmds}
    # Fail loud if a row was dropped rather than silently profiling a subset.
    [ \$(grep -c '^- id:' sr_samplesheet.yml) -eq ${layout.size()} ] || { echo "RUN_SUPERRESOLUTION: samplesheet row count != ${layout.size()}" >&2; exit 1; }

    if [ '${remoteRepo}' = 'true' ]; then
        # SR_PULL_REPO already pulled and patched this checkout; the helper only
        # resolves its version-dependent location here (it rewrites nothing when the
        # patch is already applied, so concurrent consumers never write to the
        # shared clone).
        launch_repo=\$(python "\$(command -v patch_sr_helpers.py)" sr_assets --repo ${repo})
    fi

    # NXF_CACHE_DIR moves the nested run's .nextflow cache + history out of this task
    # directory, which -resume needs to find its previous session. The launch dir stays
    # the task dir, so the relative --input/--outdir paths above keep working.
    export NXF_CACHE_DIR='${nestedDir}/cache'
    mkdir -p "\$NXF_CACHE_DIR"

    # A sample whose reads hit no reference is a legitimate benchmark outcome (the
    # profiler found nothing), not a pipeline error — the nested pipelines report it as
    # an all-zero composition with status=no_reference_hits rather than aborting, so one
    # empty sample no longer takes down the rest of its batch. Any nested failure is real.
    nextflow run ${launchRepo} \\
        ${nestedArgs}

    ${norm_cmds}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${repo}: ${params.sr_revision}
        nextflow: \$(nextflow -version 2>&1 | grep -oE 'version [0-9.]+' | sed 's/version //')
    END_VERSIONS
    """

    stub:
    def repo = params["sr_${metas[0].profiler == 'sr_amplicon' ? 'amplicon' : 'shotgun'}_repo"]
    def stub_cmds = metas*.id.collect { id ->
        """mkdir -p sr_out/composition/${id}
    printf 'sample,genome_id,observed_rel_abundance,inferred_mean,inferred_lo,inferred_hi\\n' > sr_out/composition/${id}/${id}.inferred_composition.csv
    printf 'genome_id\\tpredicted_rel_abundance\\tpredicted_tax_rel_abundance\\n' > ${id}.sr_profile.tsv"""
    }.join('\n    ')
    """
    ${stub_cmds}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${repo}: stub
        nextflow: stub
    END_VERSIONS
    """
}
