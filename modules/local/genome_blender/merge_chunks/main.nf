process MERGE_GENOME_BLENDER_CHUNKS {
    tag "$meta.id"
    label 'process_single'

    // Reuse the genome-blender image: it already ships pysam (wraps samtools
    // merge), so no separate samtools container is needed.
    container "ghcr.io/timrozday-mgnify/smb-genome-blender:${params.smb_genome_blender_tag}"

    input:
    tuple val(meta), path(reads), path(bams)

    output:
    tuple val(meta), path("${meta.id}*.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.bam"),       emit: bam
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = meta.id
    def r2      = meta.paired_end ? "cat *_R2.fastq.gz > ${prefix}_R2.fastq.gz" : ""
    def r1_glob = meta.paired_end ? "*_R1.fastq.gz" : "*.fastq.gz"
    def r1_out  = meta.paired_end ? "${prefix}_R1.fastq.gz" : "${prefix}.fastq.gz"
    """
    # Concatenated gzip members decompress fine with any gzip-aware reader.
    cat ${r1_glob} > ${r1_out}
    ${r2}
    python3 -c "import pysam,sys; pysam.merge('-f','-o','${prefix}.bam', *sys.argv[1:])" ${bams}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    def r2 = meta.paired_end ? "touch ${prefix}_R2.fastq.gz" : ""
    def r1 = meta.paired_end ? "${prefix}_R1.fastq.gz" : "${prefix}.fastq.gz"
    """
    touch ${r1}
    ${r2}
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: stub
    END_VERSIONS
    """
}
