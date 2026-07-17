#!/usr/bin/env nextflow

include { SYNTHETIC_METAGENOMIC_BENCHMARK } from './workflows/synthetic_metagenomic_benchmark'

// Resolve a samplesheet path: absolute / URL passes through, relative resolves
// against the pipeline projectDir (matches the samplesheet convention).
def resolveFile(String p) {
    (p.startsWith('/') || p =~ /^[a-z]+:\/\//)
        ? file(p, checkIfExists: true)
        : file("${workflow.projectDir}/${p}", checkIfExists: true)
}

// Absolute read count, or null for none/null/empty/missing (no subsampling).
def parseSubsampleScalar(v) {
    (v == null || v.toString().toLowerCase() in ['none', 'null', '']) ? null : (v as long)
}

// Requested subsample depths for a sample. Absolute read counts; none/null/
// empty/missing => a single passthrough run (keep all reads, scalar null).
def parseSubsamples(v) {
    def list = (v == null) ? [null] : ((v instanceof List) ? v : [v]).collect { x -> parseSubsampleScalar(x) }.unique()
    list ?: [null]
}

workflow {
    main:
    if (!params.input) {
        error "Provide a samplesheet with --input"
    }
    if (!(params.step in ['all', 'generate', 'profile', 'train'])) {
        error "params.step must be one of: all | generate | profile | train (got '${params.step}')"
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
    ch_pretrained = Channel.empty()
    ch_profile_in = Channel.empty()
    def pretrainedIds = [] as Set

    if (params.step in ['all', 'generate']) {
        // Generate samplesheet columns:
        //   sample,train_id,train_fastq_1,train_fastq_2,train_subsample,platform,genomes_csv,num_reads,mode,profiler,database,chunks,error_model_dir
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
                chunks:    ((row.chunks ?: params.chunks) as int),
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

        // Rows pointing at an already-trained error-model dir (from a prior
        // `--step train` run): reach in for the model + calibration, keyed by
        // train_id, and skip training for those train_ids.
        def pretrainedRows = rows.findAll { it.error_model_dir }
        pretrainedIds = pretrainedRows.collect { it.train_id } as Set  // reassigns the outer local
        ch_pretrained = Channel
            .fromList(pretrainedRows.collect { row ->
                def dir = resolveFile(row.error_model_dir)
                def m = files("${dir}/*.model.pt")
                def c = files("${dir}/*.phred_calibration.json")
                if (m.size() != 1 || c.size() != 1) {
                    error "error_model_dir '${dir}' for train_id ${row.train_id} must contain exactly one *.model.pt and one *.phred_calibration.json"
                }
                [ row.train_id, m[0], c[0] ]
            })
            .unique { it[0] }
    }

    if (params.step in ['all', 'generate', 'train']) {
        // Training channel, deduped per train_id: [ meta_train, [ reads ] ].
        // Pretrained train_ids are excluded (their model comes from disk).
        ch_train = ch_rows
            .filter { row -> !(row.train_id in pretrainedIds) }
            .map { row ->
                def reads = [ resolveFile(row.train_fastq_1) ]
                if (row.train_fastq_2?.trim()) reads << resolveFile(row.train_fastq_2)
                def meta_train = [ id: row.train_id, platform: row.platform, subsample: parseSubsampleScalar(row.train_subsample) ]
                [ row.train_id, meta_train, reads ]
            }
            .unique { it[0] }
            .map { train_id, meta_train, reads -> [ meta_train, reads ] }
    }

    if (params.step == 'profile') {
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

    SYNTHETIC_METAGENOMIC_BENCHMARK(ch_samples, ch_train, ch_pretrained, ch_profile_in)
}
