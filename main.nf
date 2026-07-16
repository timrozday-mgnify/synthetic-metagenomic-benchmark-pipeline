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
        error "Provide a samplesheet with --input"
    }
    if (!(params.step in ['all', 'generate', 'profile'])) {
        error "params.step must be one of: all | generate | profile (got '${params.step}')"
    }

    ch_rows = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)

    ch_samples    = Channel.empty()
    ch_train      = Channel.empty()
    ch_profile_in = Channel.empty()

    if (params.step in ['all', 'generate']) {
        // Generate samplesheet columns:
        //   sample,train_id,train_fastq_1,train_fastq_2,platform,genomes_csv,num_reads,mode,profiler,database
        ch_samples = ch_rows.map { row ->
            def meta = [
                id:        row.sample,
                train_id:  row.train_id,
                platform:  row.platform,
                mode:      (row.mode ?: 'shotgun').trim(),
                paired_end: (row.paired_end?.trim()) ? row.paired_end.toBoolean() : params.paired_end,
                read_length_mean:     (row.read_length_mean?.trim() ?: params.read_length_mean) as double,
                read_length_variance: (row.read_length_variance?.trim() ?: params.read_length_variance) as double,
                num_reads: (row.num_reads as long),
                profiler:  (row.profiler ?: '').trim(),
                database:  (row.database ?: '').trim(),
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
    }
    else {
        // Profile-only samplesheet columns: sample,profiler,benchmark_dir,database
        // Reads are discovered inside each benchmark_dir (the layout this pipeline
        // publishes to ${outdir}/${sample}). The predicted profile is published back
        // to ${outdir}/${sample}; point --outdir at the benchmark root to co-locate
        // it with the existing truth.tsv.
        ch_profile_in = ch_rows.map { row ->
            def dir = resolveFile(row.benchmark_dir)
            def reads = files("${dir}/*.fastq.gz").sort()
            if (!reads) error "No *.fastq.gz found in benchmark_dir '${dir}' for sample ${row.sample}"
            def meta = [
                id:       row.sample,
                mode:     (row.mode ?: 'paired'),
                profiler: (row.profiler ?: '').trim(),
                database: (row.database ?: '').trim(),
            ]
            [ meta, reads ]
        }
    }

    SYNTHETIC_METAGENOMIC_BENCHMARK(ch_samples, ch_train, ch_profile_in)
}
