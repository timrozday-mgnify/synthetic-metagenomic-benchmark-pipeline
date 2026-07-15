#!/usr/bin/env nextflow

include { SYNTHETIC_METAGENOMIC_BENCHMARK } from './workflows/synthetic_metagenomic_benchmark'

// Resolve a samplesheet path: absolute / URL passes through, relative resolves
// against the pipeline projectDir (matches the samplesheet convention).
def resolveFile(String p) {
    (p.startsWith('/') || p =~ /^[a-z]+:\/\//)
        ? file(p, checkIfExists: true)
        : file("${workflow.projectDir}/${p}", checkIfExists: true)
}

workflow {
    main:
    if (!params.input) {
        error "Provide a samplesheet with --input (columns: sample,train_id,train_fastq_1,train_fastq_2,platform,genomes_csv,num_reads,mode)"
    }

    ch_rows = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)

    // Synthetic-sample channel: [ meta, genomes_csv, [ fasta files ] ]
    ch_samples = ch_rows.map { row ->
        def meta = [
            id:        row.sample,
            train_id:  row.train_id,
            platform:  row.platform,
            mode:      (row.mode ?: 'paired'),
            num_reads: (row.num_reads as long),
        ]
        def genomesCsv = resolveFile(row.genomes_csv)
        // Resolve the FASTA files referenced by the genomes CSV so Nextflow stages them.
        def fastas = genomesCsv.readLines()
            .drop(1)
            .findAll { it?.trim() }
            .collect { line -> resolveFile(line.split(',')[1].trim()) }
        [ meta, genomesCsv, fastas ]
    }

    // Training channel, deduped per train_id: [ meta_train, [ reads ] ]
    ch_train = ch_rows
        .map { row ->
            def reads = [ resolveFile(row.train_fastq_1) ]
            if (row.train_fastq_2?.trim()) reads << resolveFile(row.train_fastq_2)
            def meta_train = [ id: row.train_id, platform: row.platform ]
            [ row.train_id, meta_train, reads ]
        }
        .unique { it[0] }
        .map { train_id, meta_train, reads -> [ meta_train, reads ] }

    SYNTHETIC_METAGENOMIC_BENCHMARK(ch_samples, ch_train)
}
