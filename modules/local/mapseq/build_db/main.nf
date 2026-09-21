// Have mapseq build & cache a reference FASTA's clustering file (<fasta>.mscluster).
// We only keep that side-effect; the search results printed to stdout are discarded.
process MAPSEQ_CLUSTER {
    tag "$meta.id"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mapseq:2.1.1b--h3ab3c3b_0'
        : 'quay.io/biocontainers/mapseq:2.1.1b--h3ab3c3b_0'}"

    input:
    tuple val(meta), path(fasta), path(tax)

    output:
    tuple val(meta), path("${fasta}.mscluster"), emit: mscluster
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # mapseq clusters the whole reference set whatever the query, so a one-record query
    # builds the same .mscluster without the self-search (~15% of the run on SILVA NR99).
    awk '/^>/ { n++ } n <= 1' ${fasta} > probe.fasta
    mapseq probe.fasta ${fasta} ${tax} -nthreads ${task.cpus} $args > /dev/null

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mapseq: 2.1.1b
    END_VERSIONS
    """

    stub:
    """
    printf '0 0\\n' > ${fasta}.mscluster

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mapseq: stub
    END_VERSIONS
    """
}
