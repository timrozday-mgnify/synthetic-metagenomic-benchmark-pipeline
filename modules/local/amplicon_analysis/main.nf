// Run the EBI-Metagenomics amplicon-analysis-pipeline (AAP) on the generated
// reads via a nested `nextflow run`. AAP is a full pipeline whose workflow is
// driven by params.input (no channel input), so we compose it as an ordinary
// process: write an AAP-format samplesheet from the staged reads and launch it.
//
// Runs on the host (executor 'local', no container) so it reuses the host
// nextflow + container engine — AAP manages its own containers/DBs internally.
// mapseq_databases and other AAP params come from the optional -c config.
process RUN_AAP {
    tag "$meta.id"
    label 'process_single'
    executor 'local'

    input:
    // Optional inputs are stageAs'd to distinct fixed names so the shared NO_FILE
    // placeholder doesn't collide across slots. `use_built` selects the DB source.
    // rfam_cm / rfam_claninfo are absolute host paths (val, not staged): the nested
    // AAP run reads them directly (this process is executor 'local', no container).
    tuple val(meta), path(reads), val(use_built), path(aap_config, stageAs: 'aap_config_in'), path(mapseq_fasta, stageAs: 'mapseq_db.fasta'), path(mapseq_tax, stageAs: 'mapseq_db.tax'), path(mapseq_otu, stageAs: 'mapseq_db.otu'), path(mapseq_mscluster, stageAs: 'mapseq_db.mscluster'), val(rfam_cm), val(rfam_claninfo)

    output:
    tuple val(meta), path("aap_out/**"), emit: results
    path "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def reads_list = reads instanceof List ? reads : [reads]
    def paired     = reads_list.size() > 1
    def fastq_1    = reads_list[0]
    def fastq_2    = paired ? reads_list[1] : ''
    def single_end = paired ? 'false' : 'true'
    // A pipeline-built mapseq DB wins over the pass-through params.aap_config.
    def dbname     = meta.database ?: 'community'
    def db_cfg     = use_built ? '-c aap.config' : (params.aap_config ? "-c ${aap_config}" : '')
    // Engine + extra -c files come from meta (samplesheet aap_configs/aap_profile or the
    // params fallback, resolved to absolute paths in main.nf). DB config first, so later
    // files override earlier and the engine config wins. -profile only when requested.
    def extra_cfg  = (meta.aap_configs ?: []).collect { "-c ${file(it, checkIfExists: true)}" }.join(' ')
    def prof_arg   = meta.aap_profile ? "-profile ${meta.aap_profile}" : ''
    """
    ${use_built ? "write_aap_config.py --name '${dbname}' --fasta ${mapseq_fasta} --tax ${mapseq_tax} --otu ${mapseq_otu} --mscluster ${mapseq_mscluster} --rfam-covariance-model '${rfam_cm}' --rfam-claninfo '${rfam_claninfo}' --output aap.config" : "true"}

    # AAP samplesheet: sample,fastq_1,fastq_2,single_end (absolute paths to staged reads).
    printf 'sample,fastq_1,fastq_2,single_end\\n' > aap_samplesheet.csv
    printf '%s,%s,%s,%s\\n' "${meta.id}" "\$(readlink -f ${fastq_1})" "${fastq_2 ? "\$(readlink -f ${fastq_2})" : ''}" "${single_end}" >> aap_samplesheet.csv

    nextflow run ebi-metagenomics/amplicon-analysis-pipeline \\
        -r ${params.aap_revision} \\
        ${prof_arg} \\
        --input aap_samplesheet.csv \\
        --outdir aap_out \\
        ${db_cfg} ${extra_cfg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amplicon-analysis-pipeline: ${params.aap_revision}
        nextflow: \$(nextflow -version 2>&1 | grep -oE 'version [0-9.]+' | sed 's/version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p aap_out/taxonomy-summary
    touch aap_out/taxonomy-summary/${meta.id}.krona.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amplicon-analysis-pipeline: stub
        nextflow: stub
    END_VERSIONS
    """
}
