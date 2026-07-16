// nf-core sylph/profile, adapted to this pipeline's versions.yml convention
// (upstream: https://github.com/nf-core/modules/tree/master/modules/nf-core/sylph/profile).
process SYLPH_PROFILE {
    tag "$meta.id"
    label 'process_high'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sylph:0.9.0--ha6fb395_0'
        : 'quay.io/biocontainers/sylph:0.9.0--ha6fb395_0'}"

    input:
    tuple val(meta), path(reads)
    path database

    output:
    tuple val(meta), path("${meta.id}.sylph.tsv"), emit: profile
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_reads = (reads instanceof List && reads.size() > 1)
        ? "-1 ${reads[0]} -2 ${reads[1]}"
        : "-r ${reads instanceof List ? reads[0] : reads}"
    """
    sylph profile \\
        ${database} \\
        ${input_reads} \\
        -t ${task.cpus} \\
        $args \\
        -o ${prefix}.sylph.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(sylph -V | sed 's/sylph //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'Sample_file\\tGenome_file\\tTaxonomic_abundance\\tSequence_abundance\\tContig_name\\n' > ${prefix}.sylph.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: stub
    END_VERSIONS
    """
}
