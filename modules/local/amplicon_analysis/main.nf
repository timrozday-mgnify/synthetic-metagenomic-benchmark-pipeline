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
    tuple val(meta), path(reads), path(aap_config)

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
    def cfg_arg    = aap_config.name != 'NO_FILE' ? "-c ${aap_config}" : ''
    """
    # AAP samplesheet: sample,fastq_1,fastq_2,single_end (absolute paths to staged reads).
    printf 'sample,fastq_1,fastq_2,single_end\\n' > aap_samplesheet.csv
    printf '%s,%s,%s,%s\\n' "${meta.id}" "\$(readlink -f ${fastq_1})" "${fastq_2 ? "\$(readlink -f ${fastq_2})" : ''}" "${single_end}" >> aap_samplesheet.csv

    nextflow run ebi-metagenomics/amplicon-analysis-pipeline \\
        -r ${params.aap_revision} \\
        -profile ${params.aap_profile} \\
        --input aap_samplesheet.csv \\
        --outdir aap_out \\
        ${cfg_arg}

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
