#!/usr/bin/env nextflow

include { SYNTHETIC_METAGENOMIC_BENCHMARK } from './workflows/synthetic_metagenomic_benchmark'

// Resolve a samplesheet path: absolute / URL passes through, relative resolves
// against the pipeline projectDir (matches the samplesheet convention).
def resolveFile(String p) {
    (p.startsWith('/') || p =~ /^[a-z]+:\/\//)
        ? file(p, checkIfExists: true)
        : file("${workflow.projectDir}/${p}", checkIfExists: true)
}

// Requested subsample depths for a sample. Absolute read counts; none/null/
// empty/missing => a single passthrough run (keep all reads, scalar null).
def parseSubsamples(v) {
    def scalar = { x -> (x == null || x.toString().toLowerCase() in ['none', 'null', '']) ? null : (x as long) }
    def list = (v == null) ? [null] : ((v instanceof List) ? v : [v]).collect(scalar).unique()
    list ?: [null]
}

workflow {
    main:
    if (!params.input) {
        error "Provide a samplesheet with --input"
    }
    if (!(params.step in ['all', 'generate', 'profile'])) {
        error "params.step must be one of: all | generate | profile (got '${params.step}')"
    }

    // YAML samplesheet: a list of sample maps (allows nested fields like the
    // per-sample `subsample` list). See README for the schema.
    def rows = new org.yaml.snakeyaml.Yaml().load(file(params.input, checkIfExists: true).text)
    if (!(rows instanceof List)) {
        error "Samplesheet ${params.input} must be a YAML list of sample entries"
    }
    ch_rows = Channel.fromList(rows)

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
                mode:      (row.mode ?: 'shotgun'),
                paired_end: (row.paired_end != null) ? row.paired_end : params.paired_end,
                read_length_mean:     (row.read_length_mean     ?: params.read_length_mean) as double,
                read_length_variance: (row.read_length_variance ?: params.read_length_variance) as double,
                num_reads: (row.num_reads as long),
                profiler:  (row.profiler ?: ''),
                database:  (row.database ?: ''),
                subsamples: parseSubsamples(row.subsample),
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
                sample:   row.sample,   // profile-only: no subsampling, one run per sample
                publish_subdir: '',
                mode:     (row.mode ?: 'paired'),
                profiler: (row.profiler ?: ''),
                database: (row.database ?: ''),
            ]
            [ meta, reads ]
        }
    }

    SYNTHETIC_METAGENOMIC_BENCHMARK(ch_samples, ch_train, ch_profile_in)
}
