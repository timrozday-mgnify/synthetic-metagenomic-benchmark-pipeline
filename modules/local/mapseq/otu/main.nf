// Majority-vote each mapseq cluster to one taxonomy string -> the mapseq .otu file.
// Pure-python (bin/build_mapseq_otu.py), so it runs in the smb-skiver image.
process MAPSEQ_OTU {
    tag "$meta.id"
    label 'process_single'
    container "ghcr.io/timrozday-mgnify/smb-skiver:${params.smb_skiver_tag}"

    input:
    tuple val(meta), path(mscluster), path(headers_json)

    output:
    tuple val(meta), path("${meta.id}.mapseq.otu"), emit: otu
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    build_mapseq_otu.py \\
        --mscluster ${mscluster} \\
        --headers ${headers_json} \\
        --out-otu ${prefix}.mapseq.otu

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '0\\tBacteria;Genus;species\\t0\\n' > ${prefix}.mapseq.otu

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
