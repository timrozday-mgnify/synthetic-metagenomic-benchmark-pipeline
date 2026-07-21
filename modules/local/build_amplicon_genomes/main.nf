// Split a batched AmpliconHunter amplicons.fa back into per-genome FASTAs and a
// genomes CSV, so the extracted amplicons feed the normal amplicon-generation path
// with each genome's original abundance + ground-truth identity preserved.
//
// Reuses the genome-blender image (python:3.11-slim base); the split is stdlib-only.

process BUILD_AMPLICON_GENOMES {
    tag "${meta.id}"
    label 'process_single'

    container "ghcr.io/timrozday-mgnify/smb-genome-blender:${params.smb_genome_blender_tag}"

    input:
    tuple val(meta), path(amplicons_fa), path(genomes_csv)

    output:
    tuple val(meta), path("amplicon_genomes.csv"), path("amplicons/*"), emit: references
    path "versions.yml",                                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    build_amplicon_genomes_csv.py \\
        ${amplicons_fa} \\
        ${genomes_csv} \\
        amplicon_genomes.csv \\
        amplicons/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    // Fabricate a valid amplicon-genomes set from the input genomes CSV (one
    // placeholder amplicon FASTA per genome), so downstream stub generation runs.
    """
    mkdir -p amplicons
    { head -n1 ${genomes_csv}; } > amplicon_genomes.csv
    tail -n +2 ${genomes_csv} | while IFS=, read -r gid fpath abund; do
        [ -z "\$gid" ] && continue
        printf '>%s_amplicon0\\nACGT\\n' "\$gid" > "amplicons/\$gid.fa"
        printf '%s,amplicons/%s.fa,%s\\n' "\$gid" "\$gid" "\$abund" >> amplicon_genomes.csv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
