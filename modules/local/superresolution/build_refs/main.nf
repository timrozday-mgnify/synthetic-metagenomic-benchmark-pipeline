// Build the combined reference FASTA the superresolution pipelines take
// (`--references`): one record per contig / 16S copy, headers `{genome_id}|{n}|{orig}`
// so both pipelines recover the genome id as the text before the first '|'.
//
// Driven by a genomes CSV (genome_id,fasta_path,...) whose FASTAs are staged
// alongside; used both for `database: self` (the sample's own genomes CSV) and for
// a named `databases:` collection (BUILD_DATABASES writes an equivalent CSV).
// A `taxonomy` column on every row also yields `<prefix>.sr_refs.tax`, the MAPseq
// taxonomy the nested run takes as --taxonomy. Rows with a `taxon` (a collection's taxon
// panel entries) go to `<prefix>.panel_taxa.tsv` instead, the nested --panel_taxa; a
// taxa-only collection writes no FASTA (meta.taxa_only / meta.has_taxa tell the stub).
process SR_BUILD_REFS {
    tag "$meta.id"
    label 'process_single'

    // Pure-stdlib python; reuse the genome-blender image (ships python).
    container "ghcr.io/timrozday-mgnify/smb-genome-blender:${params.smb_genome_blender_tag}"

    input:
    // Either a real genomes CSV (`database: self`) or, when it's the NO_FILE
    // placeholder, a manifest string written here, header line included (a named
    // `databases:` collection has no CSV of its own) — same trick as MAPSEQ_PREP.
    tuple val(meta), path(genomes_csv), path(fastas), val(manifest)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.sr_refs.fasta"), emit: refs, optional: true
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.panel_taxa.tsv"), emit: taxa, optional: true
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.sr_refs.tax"), emit: tax, optional: true
    path "versions.yml",                                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def csv_cmd  = genomes_csv.name != 'NO_FILE'
        ? "cp ${genomes_csv} sr_genomes.csv"
        : "printf '${manifest}\\n' > sr_genomes.csv"
    """
    ${csv_cmd}

    python "\$(command -v build_sr_refs.py)" \\
        --genomes-csv sr_genomes.csv \\
        --output ${prefix}.sr_refs.fasta \\
        --tax-output ${prefix}.sr_refs.tax \\
        --taxa-output ${prefix}.panel_taxa.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ${meta.taxa_only ? '' : "printf '>stub|0|stub\\nACGT\\n' > ${prefix}.sr_refs.fasta"}
    ${manifest.startsWith('genome_id,fasta_path,taxonomy') ? "printf 'stub|0|stub\\tstub\\n' > ${prefix}.sr_refs.tax" : ''}
    ${meta.has_taxa ? "printf 'id\\ttaxon\\nstub\\tstub\\n' > ${prefix}.panel_taxa.tsv" : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
