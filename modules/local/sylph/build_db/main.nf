// Build a sylph database (.syldb) from reference genome FASTAs.
// Used for the `database = self` case: the DB is the exact set of genomes the
// reads were generated from, so sylph's genome-level profile lines up with the
// genome-level ground truth (truth.tsv).
process SYLPH_BUILD_DB {
    tag "$meta.id"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sylph:0.9.0--ha6fb395_0'
        : 'quay.io/biocontainers/sylph:0.9.0--ha6fb395_0'}"

    input:
    tuple val(meta), path(fastas)

    output:
    tuple val(meta), path("${meta.id}.syldb"), emit: db
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sylph sketch \\
        -g ${fastas} \\
        -o ${prefix} \\
        -t ${task.cpus} \\
        $args
    # `sylph sketch -o PREFIX` writes PREFIX.syldb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(sylph -V | sed 's/sylph //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.syldb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: stub
    END_VERSIONS
    """
}
