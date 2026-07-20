// Stage a named sequence collection's 16S into a single mapseq reference FASTA +
// taxonomy file (and a headers.json sidecar for the OTU step). Pure-python
// (bin/build_mapseq_refs.py), so it runs in the smb-skiver image.
process MAPSEQ_PREP {
    tag "$meta.id"
    label 'process_single'
    container "ghcr.io/timrozday-mgnify/smb-skiver:${params.smb_skiver_tag}"

    input:
    tuple val(meta), path(ssu_fastas), val(manifest_text)

    output:
    tuple val(meta), path("${meta.id}.mapseq.fasta"), path("${meta.id}.mapseq.tax"), path("${meta.id}.headers.json"), emit: refs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // manifest_text is `id\tssu_basename\ttaxonomy` rows joined by literal \n (single
    // line, so Nextflow's script stripIndent stays well-defined).
    // ponytail: assumes distinct ssu basenames across the collection; stageAs if that breaks.
    """
    printf '${manifest_text}\\n' > manifest.tsv

    build_mapseq_refs.py \\
        --manifest manifest.tsv \\
        --out-fasta ${prefix}.mapseq.fasta \\
        --out-tax ${prefix}.mapseq.tax \\
        --out-headers ${prefix}.headers.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '>x|0|s\\nACGT\\n' > ${prefix}.mapseq.fasta
    printf '#cutoff: 0.00:0.08\\n#name: custom\\n#levels: Kingdom Genus Species\\nx|0|s\\tBacteria;Genus;species\\n' > ${prefix}.mapseq.tax
    printf '{"fasta_order":["x|0|s"],"header_to_id":{"x|0|s":"x"},"id_to_tax":{"x":"Bacteria;Genus;species"}}' > ${prefix}.headers.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
