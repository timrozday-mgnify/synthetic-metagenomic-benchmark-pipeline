process GROUND_TRUTH {
    tag "$meta.id"
    label 'process_single'

    // Reuse the genome-blender image: it already ships pysam + python, so no
    // separate samtools container is needed for sort/index/idxstats.
    container "ghcr.io/timrozday-mgnify/smb-genome-blender:${params.smb_genome_blender_tag}"

    input:
    tuple val(meta), path(bam), path(genomes_csv), path(names)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.sorted.bam.bai"), emit: bai
    tuple val(meta), path("${meta.id}.truth.tsv"),      emit: truth
    path "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    build_truth.py \\
        --bam ${bam} \\
        --genomes-csv ${genomes_csv} \\
        --keep-names ${names} \\
        --sorted-bam ${prefix}.sorted.bam \\
        --truth-tsv ${prefix}.truth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sorted.bam ${prefix}.sorted.bam.bai
    printf 'genome_id\\ttarget_rel_abundance\\trealized_n_reads\\trealized_rel_abundance\\n' > ${prefix}.truth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: stub
    END_VERSIONS
    """
}
