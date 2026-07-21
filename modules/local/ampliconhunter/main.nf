// AmpliconHunter — in-silico PCR across a collection of genome FASTAs
//
// AmpliconHunter screens ONE primer pair against many genomes per invocation.
// This module therefore receives all of a sample's genome FASTAs as a list
// (staged into genomes/) and one primer pair as value strings. The genome list
// file (absolute paths, as the tool requires) is generated inside the script.
//
// Reference: https://github.com/rhowardstone/AmpliconHunter
// Ported from mimicc-primer-investigations/amplicon-primer-screen.
//
// Primer file format (AmpliconHunter reads with .strip().split()):
//   <forward_seq>\t<reverse_seq>   — one pair per file; IUPAC degenerate bases OK
//
// Useful ext.args for the 'run' subcommand (set via params.ampliconhunter_args):
//   --mismatches <int>   mismatches allowed in either primer (default: 0)
//   --clamp <int>        3'-end bases that must match perfectly (default: 5)
//   --Lmin <int>         minimum amplicon length (default: 50)
//   --Lmax <int>         maximum amplicon length (default: 5000)
//   --Tm <float>         minimum melting temperature filter (°C)
//   --trim-primers       remove primer sequences from amplicon output
//   --include-offtarget  also report FF/RR off-target amplicons

process AMPLICONHUNTER {
    tag "${meta.id}"
    label 'process_medium'

    container "ghcr.io/timrozday-mgnify/smb-ampliconhunter:${params.smb_ampliconhunter_tag}"

    input:
    // meta.id = "<sample>.<pair>" (unique per sample × primer pair).
    // fastas: all of the sample's genome FASTAs, staged into genomes/.
    // forward / reverse: primer sequences (IUPAC degenerate bases allowed).
    tuple val(meta), path(fastas, stageAs: 'genomes/*'), val(forward), val(reverse)

    output:
    tuple val(meta), path("${prefix}/amplicons.fa"),        emit: amplicons, optional: true
    tuple val(meta), path("${prefix}/run_statistics.json"), emit: stats,     optional: true
    tuple val(meta), path("${prefix}/", type: 'dir'),       emit: results_dir
    path "versions.yml",                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    # Decompress any gzipped genomes — AmpliconHunter reads raw FASTA.
    mkdir -p decompressed
    for f in genomes/*; do
        base=\$(basename "\$f")
        if [[ "\$f" == *.gz ]]; then
            gunzip -c "\$f" > "decompressed/\${base%.gz}"
        else
            cp "\$f" "decompressed/\$base"
        fi
    done

    # AmpliconHunter requires absolute paths in the genome list.
    printf '%s\\n' "\$PWD/decompressed/"* > genome_list.txt

    # Primer file: two whitespace-separated sequences on one line.
    printf '%s\\t%s\\n' '${forward}' '${reverse}' > primers.txt

    ampliconhunter --base-dir .ah_cache run \\
        --threads ${task.cpus} \\
        $args \\
        genome_list.txt \\
        primers.txt \\
        ${prefix}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampliconhunter: \$(ampliconhunter --version 2>&1 | head -1 || echo "unknown")
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    printf '>ENA|STUB| source=genomes.coordinates=1-4.orientation=FR.Tm=60\\nACGT\\n' > ${prefix}/amplicons.fa
    echo '{}' > ${prefix}/run_statistics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampliconhunter: stub
        python: stub
    END_VERSIONS
    """
}
