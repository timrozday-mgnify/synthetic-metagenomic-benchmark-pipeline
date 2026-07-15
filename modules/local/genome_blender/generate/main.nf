process GENOME_BLENDER_GENERATE {
    tag "$meta.id"
    label 'process_medium'

    container "ghcr.io/timrozday-mgnify/smb-genome-blender:${params.smb_genome_blender_tag}"

    input:
    tuple val(meta), path(genomes_csv), path(fastas), path(model_pt), path(phred_cal)

    output:
    tuple val(meta), path("${meta.id}*.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.bam"),       emit: bam
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode_flag = [
        paired:   '--paired-end',
        single:   '--single-end',
        long:     '--long-read',
        amplicon: '--amplicon',
    ][meta.mode ?: 'single']
    if (mode_flag == null) {
        error "Unknown generation mode '${meta.mode}' for sample ${meta.id} (expected paired|single|long|amplicon)"
    }
    """
    # Point the genomes CSV at the locally staged FASTA basenames.
    rewrite_genomes_csv.py ${genomes_csv} genomes.local.csv

    generate-reads \\
        --input-csv genomes.local.csv \\
        --num-reads ${meta.num_reads} \\
        --output-prefix ${prefix} \\
        --skiver-model ${model_pt} \\
        --skiver-phred-calibration ${phred_cal} \\
        --seed ${params.seed} \\
        ${mode_flag} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genome-blender: \$(generate-reads --version 2>&1 | tail -n1 || echo "unknown")
        skiver-generate: \$(skiver-generate --version 2>&1 | tail -n1 || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def r2 = (meta.mode == 'paired') ? "echo stub | gzip > ${prefix}_R2.fastq.gz" : ""
    def r1 = (meta.mode == 'paired') ? "${prefix}_R1.fastq.gz" : "${prefix}.fastq.gz"
    """
    echo stub | gzip > ${r1}
    ${r2}
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genome-blender: stub
        skiver-generate: stub
    END_VERSIONS
    """
}
