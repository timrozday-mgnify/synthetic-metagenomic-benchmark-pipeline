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

// Resolve a nested run's `-c` config files (AAP, superresolution) to absolute path
// strings (a YAML list, or a comma-separated string from the CLI). Relative paths
// resolve against projectDir.
def parseNestedConfigs(v) {
    def list = (v == null) ? [] : (v instanceof List ? v : v.toString().tokenize(','))
    list.findAll { it?.toString()?.trim() }.collect { p ->
        def s = p.toString().trim()
        (s.startsWith('/') || s =~ /^[a-z]+:\/\//) ? s : "${workflow.projectDir}/${s}".toString()
    }
}

// Every profiler a sample's reads are benchmarked with: the row's `profilers`
// (a list, or comma-separated), else the samplesheet/params default. Deduped, order
// preserved. The reads are generated once and each profiler writes its own profile
// into the same benchmark dir.
def parseProfilers(row, defaultProfilers) {
    def v = (row.profilers != null) ? row.profilers : defaultProfilers
    def list = (v == null) ? [] : (v instanceof List ? v : v.toString().tokenize(','))
    list.collect { it?.toString()?.trim() }.findAll { it }.unique()
}

// Knobs a superresolution `sr_settings:` entry may set. The first four reach the
// matrix build, the rest the inference run; a key absent from an entry falls back to
// the matching sr_<amplicon|shotgun>_<key> param. Kept as a method, not a top-level
// `def`: a script-level variable is local to the run body and invisible in here.
def srSettingKeys() {
    ['mismapping_method', 'align_backend', 'align_tau', 'matrix_args',
     'infer_presence', 'infer_presence_prior', 'infer_presence_temp', 'inference_args']
}

// Superresolution parameter fan-out. A row's (or the samplesheet's) `sr_settings:` is
// a list of named knob maps; the sample's sr_* profilers run once per entry, each with
// its own mis-mapping matrix, its own nested run and its own `<id>.<name>.sr_profile.tsv`.
// Absent => one unnamed entry taking every knob from params — the pre-fan-out behaviour,
// which is why the fallback entry must have a null name (it must not rename anything).
def parseSrSettings(row, defaultSettings) {
    def v = (row.sr_settings != null) ? row.sr_settings : defaultSettings
    if (v == null) return [[ name: null, opts: [:] ]]
    if (!(v instanceof List)) {
        error "sr_settings must be a YAML list of named knob maps (got ${v.getClass().simpleName})"
    }
    def known = srSettingKeys()
    def out = v.findAll { it != null }.collect { e ->
        if (!(e instanceof Map)) error "sr_settings entry '${e}' must be a map with a 'name'"
        def name = (e.name ?: e.id)?.toString()?.trim()
        if (!name) error "sr_settings entry ${e} needs a 'name' (it names the output file)"
        if (!(name ==~ /[A-Za-z0-9._-]+/)) {
            error "sr_settings name '${name}' must be [A-Za-z0-9._-]+ (it becomes part of a filename)"
        }
        def unknown = e.keySet().findAll { !(it in known) && !(it in ['name', 'id']) }
        if (unknown) error "sr_settings '${name}': unknown key(s) ${unknown} (expected ${known})"
        [ name: name, opts: known.collectEntries { k -> [ (k): e[k] ] }.findAll { k, val -> val != null } ]
    }
    if (out*.name.unique().size() != out.size()) {
        error "sr_settings names must be unique (got ${out*.name})"
    }
    out ?: [[ name: null, opts: [:] ]]
}

workflow {
    main:
    if (!params.input) {
        error "Provide a samplesheet with --input"
    }
    if (!(params.step in ['all', 'generate', 'profile', 'train'])) {
        error "params.step must be one of: all | generate | profile | train (got '${params.step}')"
    }

    // superresolution mis-mapping mode. Checked here rather than only inside the nested
    // run, which is minutes of read generation away: a typo in a sweep should fail in
    // the first second, not after the reads are made.
    if (params.sr_amplicon_mismapping_method != null
        && !(params.sr_amplicon_mismapping_method in ['simulate', 'align'])) {
        error "sr_amplicon_mismapping_method must be 'simulate' or 'align' " +
              "(got '${params.sr_amplicon_mismapping_method}')"
    }
    if (params.sr_amplicon_align_backend != null) {
        if (!(params.sr_amplicon_align_backend in ['minimap2', 'exact-hash', 'kmer'])) {
            error "sr_amplicon_align_backend must be 'minimap2', 'exact-hash' or 'kmer' " +
                  "(got '${params.sr_amplicon_align_backend}')"
        }
        if (params.sr_amplicon_mismapping_method == 'simulate') {
            error "sr_amplicon_align_backend only applies to " +
                  "sr_amplicon_mismapping_method = 'align'"
        }
    }
    if (params.sr_amplicon_align_tau != null && (params.sr_amplicon_align_tau as int) < 0) {
        error "sr_amplicon_align_tau must be >= 0 (got '${params.sr_amplicon_align_tau}')"
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

    // Optional samplesheet-level settings for the nested runs (AAP, superresolution);
    // fall back to params. Global for the run (the container engine is a site-wide
    // setting, not per-sample).
    def effAapConfigs = parseNestedConfigs((loaded instanceof Map && loaded.aap_configs != null) ? loaded.aap_configs : params.aap_configs)
    def effAapProfile = (loaded instanceof Map ? loaded.aap_profile : null) ?: params.aap_profile
    def effSrConfigs  = parseNestedConfigs((loaded instanceof Map && loaded.sr_configs != null) ? loaded.sr_configs : params.sr_configs)
    def effSrProfile  = (loaded instanceof Map ? loaded.sr_profile : null) ?: params.sr_profile
    // Samplesheet-level default for the superresolution fan-out; a row's own
    // `sr_settings:` wins. No params equivalent — a sweep is a samplesheet, not a flag.
    def defaultSrSettings = (loaded instanceof Map) ? loaded.sr_settings : null

    //
    // Named sequence collections -> profiler DBs. A collection is built (or its
    // pre-built dir consumed) only if some sample references it by `database` name
    // with a matching profiler. Names not defined under
    // `databases:` fall back to params.sylph_databases / params.aap_config.
    //
    def knownProfilers = ['sylph', 'aap', 'sr_amplicon', 'sr_shotgun']
    def defaultProfilers = (loaded instanceof Map && loaded.profilers != null)
        ? loaded.profilers : params.profilers
    rows.each { row ->
        def unknown = parseProfilers(row, defaultProfilers).findAll { !(it in knownProfilers) }
        if (unknown) {
            error "Sample ${row.sample ?: row.id}: unknown profiler(s) ${unknown} (expected ${knownProfilers})"
        }
    }
    def dbProfilers = [:]
    rows.each { row ->
        def name = row.database
        if (name && name != 'self') {
            parseProfilers(row, defaultProfilers).each { prof ->
                dbProfilers.computeIfAbsent(name) { [] as Set } << prof
            }
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
        // The superresolution reference FASTA is built from whole genomes (shotgun)
        // or 16S (amplicon); fail here rather than deep in BUILD_DATABASES.
        if (d.sequences) {
            [ sr_shotgun: 'genome', sr_amplicon: 'ssu' ].each { prof, field ->
                if (prof in profs && d.sequences.any { !it[field] }) {
                    error "database '${name}': profiler '${prof}' requires '${field}' on every sequence"
                }
            }
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
    // Collection names resolved per profiler; PROFILE uses these to decide
    // built-vs-fallback. One map so adding a profiler doesn't add a workflow arg.
    def builtNames = knownProfilers.collectEntries { prof ->
        [ (prof): dbSpecs.findAll { prof in it.profilers }.collect { it.name } as Set ]
    }
    ch_db_specs = Channel.fromList(dbSpecs)

    ch_samples    = Channel.empty()
    ch_train      = Channel.empty()
    ch_pretrained = Channel.empty()
    ch_profile_in = Channel.empty()
    def pretrainedIds = [] as Set

    if (params.step in ['all', 'generate']) {
        // Generate samplesheet columns:
        //   sample,train_id,train_fastq_1,train_fastq_2,train_subsample,platform,genomes_csv,num_reads,mode,profilers,database,chunks,error_model_dir
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
                // Every profiler these reads are benchmarked with. PROFILE fans out
                // one run per entry, all sharing this benchmark dir.
                profilers: parseProfilers(row, defaultProfilers),
                database:  (row.database ?: ''),
                subsamples: parseSubsamples(row.subsample),
                aap_configs: effAapConfigs,
                aap_profile: effAapProfile,
                sr_configs:  effSrConfigs,
                sr_profile:  effSrProfile,
                // Primer pairs for in-silico PCR (empty => no extraction, use the
                // genomes_csv FASTAs directly). Each pair is extracted + run separately.
                primer_sets: parsePrimerPairs(row.primers),
                // superresolution knob sets this sample is benchmarked under (one
                // nested run + matrix each); [[name:null, opts:[:]]] = params only.
                sr_settings: parseSrSettings(row, defaultSrSettings),
                // Precomputed mapseq classification handed to sr_amplicon instead of
                // letting it map the reads again. An absolute path string, not a staged
                // file: RUN_SUPERRESOLUTION is executor 'local' and the nested run reads
                // it directly, exactly as it does the reads.
                mseq: (row.mseq ? resolveFile(row.mseq).toString() : null),
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
        // Profile-only samplesheet columns: sample,profilers,benchmark_dir,database
        // Reads are discovered inside each benchmark_dir (the layout this pipeline
        // publishes to ${outdir}/${sample}). The predicted profile is published back
        // to ${outdir}/${sample}; point --outdir at the benchmark root to co-locate
        // it with the existing truth.tsv.
        ch_profile_in = ch_rows.map { row ->
            def dir = resolveFile(row.benchmark_dir)
            def reads = files("${dir}/*.fastq.gz").sort()
            if (!reads) error "No *.fastq.gz found in benchmark_dir '${dir}' for sample ${row.sample}"
            // Amplicon reads: the pair they were amplified with, so sr_amplicon extracts
            // its reference amplicons from the same region (it defaults to V4 515F/806R,
            // which silently matches nothing against, say, V1-V3 reads). One row is one
            // set of already-generated reads, so exactly one pair.
            // One row is one already-generated benchmark dir. `subsample: N` says that dir
            // is the N-depth subsample of `sample`, and derives the same id/publish layout
            // the generate step used, so the profile lands next to that depth's truth.tsv.
            def sub = row.subsample
            if (sub instanceof List) {
                error "profile-only row ${row.sample}: `subsample` must be the single depth of this benchmark_dir, not a list — emit one row per depth"
            }
            def subN = (sub == null || sub.toString().trim() in ['', 'none']) ? null : sub
            def pairs = parsePrimerPairs(row.primers)
            if (pairs.size() > 1) {
                error "profile-only row ${row.sample}: `primers` must name the ONE pair these reads were amplified with, got ${pairs*.getAt(0)}"
            }
            def meta = [
                id:       subN != null ? "${row.sample}.sub${subN}" : row.sample,
                sample:   row.sample,
                publish_subdir: subN != null ? "subsample_${subN}" : '',
                mode:     (row.mode ?: 'paired'),
                platform: row.platform,
                profilers: parseProfilers(row, defaultProfilers),
                primer_sets: pairs,
                primer: pairs ? pairs[0][0] : null,
                database: (row.database ?: ''),
                aap_configs: effAapConfigs,
                aap_profile: effAapProfile,
                sr_configs:  effSrConfigs,
                sr_profile:  effSrProfile,
                sr_settings: parseSrSettings(row, defaultSrSettings),
                mseq: (row.mseq ? resolveFile(row.mseq).toString() : null),
            ]
            [ meta, reads ]
        }
    }

    SYNTHETIC_METAGENOMIC_BENCHMARK(ch_samples, ch_train, ch_pretrained, ch_profile_in,
        ch_db_specs, builtNames)
}
