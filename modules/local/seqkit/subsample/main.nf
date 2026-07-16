// Subsample a (paired) read set to an exact read count for the generate stage.
// Emits the subsampled reads plus a names.txt of the kept read names so the
// ground-truth BAM can be filtered to the *same* reads (see GROUND_TRUTH).
//
// meta.subsample is the scalar depth for this run (null = passthrough: keep
// every read, still emit the full names list so downstream filtering is uniform).
process SEQKIT_SUBSAMPLE {
    tag "$meta.id"
    label 'process_low'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("sub_*"),                emit: reads
    tuple val(meta), path("${meta.id}.names.txt"), emit: names
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def reads_list = reads instanceof List ? reads : [reads]
    def r1 = reads_list[0]
    def r2 = reads_list.size() > 1 ? reads_list[1] : null
    def n  = meta.subsample
    // Names come from mate 1 with any /1,/2 suffix stripped so they match the
    // bare BAM query_name. ponytail: assumes both mates share record order + a
    // common base name (genome-blender convention) — that convention is the knob.
    def names_cmd = "seqkit seq -n -i sub_${r1} | sed 's,/[12]\$,,' > ${meta.id}.names.txt"
    if (n == null) {
        """
        ln -s ${r1} sub_${r1}
        ${r2 ? "ln -s ${r2} sub_${r2}" : ''}
        ${names_cmd}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(seqkit version | sed 's/seqkit v//')
        END_VERSIONS
        """
    } else {
        // Shuffle each mate with the same seed → identical permutation → head -n
        // keeps matching pairs, exact N, no cross-mate name lookup needed.
        // ponytail: seqkit shuffle is in-memory; swap for reservoir sampling if
        // draws ever get huge.
        """
        seqkit shuffle -s ${params.seed} ${r1} | seqkit head -n ${n} -o sub_${r1}
        ${r2 ? "seqkit shuffle -s ${params.seed} ${r2} | seqkit head -n ${n} -o sub_${r2}" : ''}
        ${names_cmd}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(seqkit version | sed 's/seqkit v//')
        END_VERSIONS
        """
    }

    stub:
    def reads_list = reads instanceof List ? reads : [reads]
    """
    ${reads_list.collect { "echo stub | gzip > sub_${it}" }.join('\n    ')}
    touch ${meta.id}.names.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: stub
    END_VERSIONS
    """
}
