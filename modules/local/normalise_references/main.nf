process NORMALISE_REFERENCES {
    tag "$meta.id"
    label 'process_single'

    // Reuse the genome-blender image: it ships the `normalise-fasta` CLI.
    container "ghcr.io/timrozday-mgnify/smb-genome-blender:${params.smb_genome_blender_tag}"

    input:
    tuple val(meta), path(genomes_csv), path(fastas)

    output:
    tuple val(meta), path(genomes_csv), path("normalised/*"), emit: references
    tuple val(meta), path("*.header_map.tsv"),                emit: mappings
    path "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Strip each reference header to its unique first token so genome-blender
    // emits one @SQ line per contig; duplicate/whitespace names otherwise
    // corrupt the merged BAM header. Basenames are preserved so the genomes
    // CSV (rewritten to local basenames downstream) still resolves.
    """
    mkdir normalised
    for fa in ${fastas}; do
        base=\$(basename "\$fa")
        normalise-fasta "\$fa" \\
            -o "normalised/\$base" \\
            -m "\$base.header_map.tsv"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genome-blender: \$(generate-reads --version 2>&1 | tail -n1 || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    mkdir normalised
    for fa in ${fastas}; do
        base=\$(basename "\$fa")
        cp "\$fa" "normalised/\$base"
        touch "\$base.header_map.tsv"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genome-blender: stub
    END_VERSIONS
    """
}
