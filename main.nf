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

// Parse a sample's `primers` into [ [pair_id, forward, reverse], ... ]. Accepts
// either an inline YAML list in the samplesheet (each entry a 3-item [id, fwd, rev]
// list or a map {pair_id/id, forward[_primer], reverse[_primer]}), or a path to a
// TSV with those columns (mirroring mimicc's data/primers/primer_pairs.tsv).
// null / empty => [] (the sample skips extraction: existing "amplicon references
// already provided" passthrough).
def parsePrimerPairs(v) {
    if (v == null) return []
    if (v instanceof List) {
        return v.findAll { it != null }.collect { e ->
            if (e instanceof Map) {
                def id = e.pair_id ?: e.id; def fwd = e.forward ?: e.forward_primer; def rev = e.reverse ?: e.reverse_primer
                if (!id || !fwd || !rev) error "primers entry ${e} needs pair_id, forward, reverse"
                [ id.toString().trim(), fwd.toString().trim(), rev.toString().trim() ]
            }
            else {
                def c = (e instanceof List) ? e : e.toString().split(/[\t,]/).toList()
                if (c.size() < 3) error "primers entry '${e}' needs 3 fields: pair_id, forward, reverse"
                [ c[0].toString().trim(), c[1].toString().trim(), c[2].toString().trim() ]
            }
        }
    }
    if (v.toString().trim() == '') return []
    def lines = resolveFile(v.toString().trim()).readLines().findAll { it?.trim() }
    if (!lines) return []
    def header = lines[0].split('\t').collect { it.trim() }
    def iId = header.findIndexOf { it in ['pair_id', 'id'] }
    def iF  = header.findIndexOf { it in ['forward', 'forward_primer'] }
    def iR  = header.findIndexOf { it in ['reverse', 'reverse_primer'] }
    if (iId < 0 || iF < 0 || iR < 0) {
        error "primers TSV '${v}' must have columns pair_id, forward[_primer], reverse[_primer]"
    }
    lines.drop(1).collect { line ->
        def c = line.split('\t')
        [ c[iId].trim(), c[iF].trim(), c[iR].trim() ]
    }
}

// Resolve nested-AAP `-c` config files to absolute path strings (a YAML list, or a
// comma-separated string from the CLI). Relative paths resolve against projectDir.
def parseAapConfigs(v) {
    def list = (v == null) ? [] : (v instanceof List ? v : v.toString().tokenize(','))
    list.findAll { it?.toString()?.trim() }.collect { p ->
        def s = p.toString().trim()
        (s.startsWith('/') || s =~ /^[a-z]+:\/\//) ? s : "${workflow.projectDir}/${s}".toString()
    }
}

workflow {
    main:
    if (!params.input) {
        error "Provide a samplesheet with --input"
    }
    if (!(params.step in ['all', 'generate', 'profile', 'train'])) {
        error "params.step must be one of: all | generate | profile | train (got '${params.step}')"
    }

    // YAML samplesheet: either a bare list of sample maps, or a map with
    // `samples:` (the list) and an optional `databases:` block of named sequence
    // collections used to build/select profiler DBs. See README for the schema.
    def loaded = new org.yaml.snakeyaml.Yaml().load(file(params.input, checkIfExists: true).text)
    def rows
    def dbDefs
    if (loaded instanceof List) {
        rows = loaded
        dbDefs = [:]
    }
    else if (loaded instanceof Map) {
        rows = loaded.samples
        dbDefs = (loaded.databases ?: [:])
        if (!(rows instanceof List)) {
            error "Samplesheet ${params.input}: 'samples:' must be a YAML list of sample entries"
        }
    }
    else {
        error "Samplesheet ${params.input} must be a YAML list, or a map with 'samples:' (and optional 'databases:')"
    }
    ch_rows = Channel.fromList(rows)

    // Optional samplesheet-level nested-AAP settings (map form); fall back to params.
    // Global for the run (the container engine is a site-wide setting, not per-sample).
    def effAapConfigs = parseAapConfigs((loaded instanceof Map && loaded.aap_configs != null) ? loaded.aap_configs : params.aap_configs)
    def effAapProfile = (loaded instanceof Map ? loaded.aap_profile : null) ?: params.aap_profile

    //
    // Named sequence collections -> profiler DBs. A collection is built (or its
    // pre-built dir consumed) only if some sample references it by `database` name
    // with a matching `profiler`. Names not defined under `databases:` fall back to
    // params.sylph_databases / params.aap_config (unchanged behaviour).
    //
    def dbProfilers = [:]
    rows.each { row ->
        def name = row.database
        def prof = row.profiler
        if (name && name != 'self' && prof in ['sylph', 'aap']) {
            dbProfilers.computeIfAbsent(name) { [] as Set } << prof
        }
    }
    def dbSpecs = []
    dbProfilers.each { name, profs ->
        def d = dbDefs[name]
        if (d == null) {
            return  // not a YAML-defined collection -> params fallback in PROFILE
        }
        if (d.path && d.sequences) {
            error "database '${name}': set either 'path' or 'sequences', not both"
        }
        // Rfam rRNA-detection DBs the nested AAP run needs (params.rrnas_rfam_*);
        // collection-level and required whenever the collection feeds 'aap'.
        def rfamCm   = d.rfam_covariance_model ? resolveFile(d.rfam_covariance_model) : null
        def rfamClan = d.rfam_claninfo         ? resolveFile(d.rfam_claninfo)         : null
        if ('aap' in profs && (!rfamCm || !rfamClan)) {
            error "database '${name}': profiler 'aap' requires 'rfam_covariance_model' and 'rfam_claninfo'"
        }
        if (d.path) {
            dbSpecs << [ name: name, profilers: profs, prebuilt_dir: resolveFile(d.path), sequences: null,
                         rfam_cm: rfamCm, rfam_claninfo: rfamClan ]
        }
        else if (d.sequences) {
            def seqs = d.sequences.collect { s ->
                [ id:       s.id,
                  genome:   s.genome ? resolveFile(s.genome) : null,
                  ssu:      s.ssu    ? resolveFile(s.ssu)    : null,
                  taxonomy: s.taxonomy ]
            }
            dbSpecs << [ name: name, profilers: profs, prebuilt_dir: null, sequences: seqs,
                         rfam_cm: rfamCm, rfam_claninfo: rfamClan ]
        }
        else {
            error "database '${name}': must define 'sequences:' or 'path:'"
        }
    }
    def builtSylphNames  = dbSpecs.findAll { 'sylph' in it.profilers }.collect { it.name } as Set
    def builtMapseqNames = dbSpecs.findAll { 'aap'   in it.profilers }.collect { it.name } as Set
    ch_db_specs = Channel.fromList(dbSpecs)

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
                aap_configs: effAapConfigs,
                aap_profile: effAapProfile,
                // Primer pairs for in-silico PCR (empty => no extraction, use the
                // genomes_csv FASTAs directly). Each pair is extracted + run separately.
                primer_sets: parsePrimerPairs(row.primers),
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
                aap_configs: effAapConfigs,
                aap_profile: effAapProfile,
            ]
            [ meta, reads ]
        }
    }

    SYNTHETIC_METAGENOMIC_BENCHMARK(ch_samples, ch_train, ch_pretrained, ch_profile_in,
        ch_db_specs, builtSylphNames, builtMapseqNames)
}
