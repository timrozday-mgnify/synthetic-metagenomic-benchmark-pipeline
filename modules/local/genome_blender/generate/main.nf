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
        shotgun:  '',
        amplicon: '--amplicon',
        long:     '--long-read',
    ][meta.mode ?: 'shotgun']
    if (mode_flag == null) {
        error "Unknown generation mode '${meta.mode}' for sample ${meta.id} (expected shotgun|amplicon|long)"
    }
    // ponytail: --long-read is single-end only (mutually exclusive with --paired-end in the CLI).
    def paired = (meta.mode == 'long') ? false : meta.paired_end
    def pair_flag = paired ? '--paired-end' : '--single-end'
    // Drop read SEQ/QUAL from the ground-truth BAM; truth tables don't read them.
    def slim_flag = params.slim_bam ? '--minimal-bam' : ''
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
        ${pair_flag} \\
        ${slim_flag} \\
        --read-length-mean ${meta.read_length_mean} \\
        --read-length-variance ${meta.read_length_variance} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genome-blender: \$(generate-reads --version 2>&1 | tail -n1 || echo "unknown")
        skiver-generate: \$(skiver-generate --version 2>&1 | tail -n1 || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def r2 = meta.paired_end ? "echo stub | gzip > ${prefix}_R2.fastq.gz" : ""
    def r1 = meta.paired_end ? "${prefix}_R1.fastq.gz" : "${prefix}.fastq.gz"
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
