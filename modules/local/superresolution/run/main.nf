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
// + container engine, exactly like RUN_AAP; the nested run manages its own images.
// The nested run does NOT inherit the outer -profile: set params.sr_profile (and
// params.sr_configs for extra -c files).
process RUN_SUPERRESOLUTION {
    tag "${meta.id} (${meta.profiler})"
    label 'process_single'
    executor 'local'

    input:
    // Reads are NOT staged: they're passed as absolute host paths (val), matching
    // RUN_AAP — this process is executor 'local' so the nested run reads them
    // directly, and the superresolution main.nf resolves relative paths against its
    // OWN projectDir, which a staged basename would break.
    tuple val(meta), val(read_paths), path(refs)

    output:
    tuple val(meta), path("${meta.id}.sr_profile.tsv"),          emit: profile
    tuple val(meta), path("sr_out/composition/${meta.id}/*"),    emit: composition
    path "versions.yml",                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def repo = meta.profiler == 'sr_amplicon' ? params.sr_amplicon_repo : params.sr_shotgun_repo
    if (!repo) error "RUN_SUPERRESOLUTION: params.${meta.profiler == 'sr_amplicon' ? 'sr_amplicon_repo' : 'sr_shotgun_repo'} is not set (sample ${meta.id})"
    // A local checkout path is run as-is; an `org/name` project is cloned here first
    // (see the clone command below for why we don't let Nextflow pull it).
    def is_path   = repo.startsWith('/') || repo.startsWith('.')
    def repo_url  = repo =~ /^[a-z]+:\/\// ? repo : "https://github.com/${repo}.git"
    def prof_arg  = meta.sr_profile ? "-profile ${meta.sr_profile}" : ''
    def extra_cfg = (meta.sr_configs ?: []).collect { "-c ${file(it, checkIfExists: true)}" }.join(' ')
    def platform  = meta.platform ? "printf '  platform: %s\\n' '${meta.platform}' >> sr_samplesheet.yml" : 'true'
    def reads     = (read_paths instanceof List ? read_paths : [read_paths]).collect { it.toString() }
    assert reads.every { it } : "RUN_SUPERRESOLUTION: empty read path for ${meta.id}"
    def read_cmds = reads.collect { "printf '    - %s\\n' '${it}' >> sr_samplesheet.yml" }.join('\n    ')
    """
    # Clone the pipeline ourselves rather than letting `nextflow run <org>/<name>` pull
    # it, for two reasons:
    #   1. Both superresolution repos carry a `vendor/skiver` submodule whose .gitmodules
    #      URL is `git@github.com:` (SSH). Nextflow's JGit initialises submodules on run
    #      and has no SSH session factory, so the pull dies with the misleading
    #      "Repository may be corrupted" — and leaves a broken clone behind. A plain
    #      clone skips submodules, which are only needed to BUILD the images anyway.
    #   2. Nextflow's asset cache is shared (\$NXF_HOME/assets); one task per sample
    #      means concurrent pulls of the same repo racing on it.
    # The chmod is the other thing Nextflow's own pull does for us: some of those
    # repos' bin/ scripts are committed 0644, and without it the nested run dies with
    # "Permission denied".
    # ponytail: costs one small clone per task. If the repos ever drop the submodule
    # (or point it at https), this collapses back to `nextflow run ${repo} -r <rev>`.
    ${is_path ? "SR_PIPELINE='${repo}'" : """git clone --quiet '${repo_url}' sr_pipeline
    git -C sr_pipeline checkout --quiet '${params.sr_revision}'
    chmod -R +x sr_pipeline/bin
    SR_PIPELINE=sr_pipeline"""}

    # One-sample YAML samplesheet; `references` must be absolute (the nested pipeline
    # resolves relative paths against its own projectDir).
    printf -- '- id: %s\\n' '${meta.id}' > sr_samplesheet.yml
    printf '  reads:\\n' >> sr_samplesheet.yml
    ${read_cmds}
    ${platform}
    printf '  references: %s\\n' "\$(realpath ${refs})" >> sr_samplesheet.yml

    nextflow run "\$SR_PIPELINE" \\
        ${prof_arg} \\
        --input sr_samplesheet.yml \\
        --outdir sr_out \\
        ${extra_cfg}

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
