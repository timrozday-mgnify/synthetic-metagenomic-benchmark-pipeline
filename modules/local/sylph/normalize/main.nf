// Map a `sylph profile` TSV onto genome_id and renormalise, so the predicted
// profile sits next to (and is directly comparable with) truth.tsv.
process NORMALIZE_SYLPH {
    tag "$meta.id"
    label 'process_single'

    // Pure-stdlib python; reuse the genome-blender image (ships python).
    container "ghcr.io/timrozday-mgnify/smb-genome-blender:${params.smb_genome_blender_tag}"

    input:
    tuple val(meta), path(sylph_tsv), path(genomes_csv)

    output:
    tuple val(meta), path("${meta.id}.sylph_profile.tsv"), emit: profile
    path "versions.yml",                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    // genomes_csv is optional: staged as NO_FILE for config DBs (no mapping).
    def map_arg  = genomes_csv.name != 'NO_FILE' ? "--genomes-csv ${genomes_csv}" : ''
    """
    normalize_sylph_profile.py \\
        --sylph-tsv ${sylph_tsv} \\
        ${map_arg} \\
        --output ${prefix}.sylph_profile.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'genome_id\\tpredicted_rel_abundance\\tpredicted_tax_rel_abundance\\n' > ${prefix}.sylph_profile.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
