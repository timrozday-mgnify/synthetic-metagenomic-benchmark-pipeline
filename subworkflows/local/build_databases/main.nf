//
// Build profiler databases from YAML-defined named sequence collections, or
// resolve a pre-built DB directory laid out like ${outdir}/databases/<name>/.
// Emits DBs keyed by collection name for PROFILE to select via a sample's
// `database` column.
//   - sylph  (WGS)      : SYLPH_BUILD_DB over the collection's genomes.
//   - mapseq (amplicon) : MAPSEQ_PREP -> MAPSEQ_CLUSTER -> MAPSEQ_OTU over 16S + taxonomy.
//   - superresolution   : one combined reference FASTA over the collection's genomes
//                         (sr_shotgun) or 16S (sr_amplicon); the nested pipeline builds
//                         its own index from it, so there's nothing else to build.
//
// Each spec (one per referenced collection) is a map:
//   [ name, profilers(Set of 'sylph'/'aap'/'sr_amplicon'/'sr_shotgun'), prebuilt_dir(file|null),
//     sequences([ {id, genome(file|null), ssu(file|null), taxonomy} ]|null) ]
//

include { SYLPH_BUILD_DB as SYLPH_BUILD_COLLECTION } from '../../../modules/local/sylph/build_db/main'
include { MAPSEQ_PREP    } from '../../../modules/local/mapseq/prep/main'
include { MAPSEQ_CLUSTER } from '../../../modules/local/mapseq/build_db/main'
include { MAPSEQ_OTU     } from '../../../modules/local/mapseq/otu/main'
include { SR_BUILD_REFS as SR_BUILD_COLLECTION_REFS } from '../../../modules/local/superresolution/build_refs/main'

// Panel FASTA field feeding each superresolution flavour's reference set.
def srSources() { [ sr_shotgun: 'genome', sr_amplicon: 'ssu' ] }

// Resolve exactly one file matching `pat` inside a pre-built DB directory.
def globOne(dir, pat, name) {
    def hits = files("${dir}/${pat}")
    if (hits.size() != 1) {
        error "prebuilt database '${name}': expected exactly one ${pat} in ${dir} (found ${hits.size()})"
    }
    hits[0]
}

workflow BUILD_DATABASES {
    take:
    ch_specs // channel of spec maps (see header)

    main:
    ch_versions = Channel.empty()
    // Placeholder for SR_BUILD_REFS' optional genomes-CSV slot (a collection has no
    // CSV of its own; the module writes one from the manifest instead).
    def no_file = file("${projectDir}/assets/NO_FILE")

    ch_specs
        .branch { spec ->
            prebuilt: spec.prebuilt_dir != null
            build:    true
        }
        .set { ch_b }

    //
    // Pre-built directories: resolve DB files by the published-name convention
    // (conf/modules.config publishes into ${outdir}/databases/<name>/).
    //
    ch_pre_sylph = ch_b.prebuilt
        .filter { 'sylph' in it.profilers }
        .map { spec -> [ spec.name, globOne(spec.prebuilt_dir, '*.syldb', spec.name) ] }

    ch_pre_mapseq = ch_b.prebuilt
        .filter { 'aap' in it.profilers }
        .map { spec ->
            [ spec.name,
              globOne(spec.prebuilt_dir, '*.mapseq.fasta', spec.name),
              globOne(spec.prebuilt_dir, '*.mapseq.tax',   spec.name),
              globOne(spec.prebuilt_dir, '*.mapseq.otu',   spec.name),
              globOne(spec.prebuilt_dir, '*.mscluster',    spec.name) ]
        }

    //
    // Build sylph DB from the collection's genomes.
    //
    ch_sylph_in = ch_b.build
        .filter { 'sylph' in it.profilers }
        .map { spec ->
            def genomes = spec.sequences.collect { it.genome }
            if (genomes.any { it == null }) {
                error "database '${spec.name}': every sequence needs a 'genome' for a sylph DB"
            }
            [ [ id: spec.name ], genomes ]
        }
    SYLPH_BUILD_COLLECTION(ch_sylph_in)
    ch_versions = ch_versions.mix(SYLPH_BUILD_COLLECTION.out.versions.first())
    ch_built_sylph = SYLPH_BUILD_COLLECTION.out.db.map { meta, db -> [ meta.id, db ] }

    //
    // Build mapseq DB from the collection's 16S + explicit taxonomy.
    //
    ch_prep_in = ch_b.build
        .filter { 'aap' in it.profilers }
        .map { spec ->
            if (spec.sequences.any { it.ssu == null || it.taxonomy == null }) {
                error "database '${spec.name}': every sequence needs 'ssu' and 'taxonomy' for a mapseq DB"
            }
            def ssu = spec.sequences.collect { it.ssu }
            // Single-line manifest (literal \t / \n) so the module's printf stays one line.
            def manifest = spec.sequences.collect { s -> "${s.id}\\t${s.ssu.name}\\t${s.taxonomy}" }.join('\\n')
            [ [ id: spec.name ], ssu, manifest ]
        }
    MAPSEQ_PREP(ch_prep_in)
    ch_versions = ch_versions.mix(MAPSEQ_PREP.out.versions.first())

    MAPSEQ_CLUSTER(MAPSEQ_PREP.out.refs.map { meta, fasta, tax, headers -> [ meta, fasta, tax ] })
    ch_versions = ch_versions.mix(MAPSEQ_CLUSTER.out.versions.first())

    ch_otu_in = MAPSEQ_CLUSTER.out.mscluster
        .map { meta, ms -> [ meta.id, meta, ms ] }
        .join(MAPSEQ_PREP.out.refs.map { meta, fasta, tax, headers -> [ meta.id, headers ] }, by: 0)
        .map { id, meta, ms, headers -> [ meta, ms, headers ] }
    MAPSEQ_OTU(ch_otu_in)
    ch_versions = ch_versions.mix(MAPSEQ_OTU.out.versions.first())

    // Assemble [ name, fasta, tax, otu, mscluster ].
    ch_built_mapseq = MAPSEQ_PREP.out.refs.map { meta, fasta, tax, headers -> [ meta.id, fasta, tax ] }
        .join(MAPSEQ_OTU.out.otu.map { meta, otu -> [ meta.id, otu ] }, by: 0)
        .join(MAPSEQ_CLUSTER.out.mscluster.map { meta, ms -> [ meta.id, ms ] }, by: 0)

    // Rfam DBs ride alongside as absolute-path strings (pass-through to the nested
    // AAP run; main.nf guarantees both are set for any 'aap' collection).
    ch_rfam = ch_b.build.filter { 'aap' in it.profilers }
        .mix(ch_b.prebuilt.filter { 'aap' in it.profilers })
        .map { spec -> [ spec.name, spec.rfam_cm.toString(), spec.rfam_claninfo.toString() ] }

    // [ name, fasta, tax, otu, mscluster, rfam_cm, rfam_claninfo ]
    ch_mapseq_dbs = ch_built_mapseq.mix(ch_pre_mapseq).join(ch_rfam, by: 0)

    //
    // Build the superresolution reference FASTA. A collection may feed both
    // flavours (different source FASTAs), so each is keyed "<name>:<source>".
    //
    ch_sr_in = ch_b.build
        .flatMap { spec ->
            srSources().findAll { prof, field -> prof in spec.profilers }.collect { prof, field ->
                def fastas = spec.sequences.collect { it[field] }
                // Single-line manifest (literal \n) so the module's printf stays one line.
                def manifest = spec.sequences.collect { s -> "${s.id},${s[field].name}" }.join('\\n')
                // meta.key is what PROFILE joins on; meta.id also names the published file.
                [ [ id: "${spec.name}_${field}", key: "${spec.name}:${field}" ], no_file, fastas, manifest ]
            }
        }
    SR_BUILD_COLLECTION_REFS(ch_sr_in)
    ch_versions = ch_versions.mix(SR_BUILD_COLLECTION_REFS.out.versions.first())
    ch_built_sr = SR_BUILD_COLLECTION_REFS.out.refs.map { meta, refs -> [ meta.key, refs ] }

    // Pre-built: the published layout is `<name>_<source>.sr_refs.fasta` (see
    // conf/modules.config), so each flavour resolves its own file.
    ch_pre_sr = ch_b.prebuilt
        .flatMap { spec ->
            srSources().findAll { prof, field -> prof in spec.profilers }.collect { prof, field ->
                [ "${spec.name}:${field}", globOne(spec.prebuilt_dir, "*_${field}.sr_refs.fasta", spec.name) ]
            }
        }

    emit:
    sylph_dbs  = ch_built_sylph.mix(ch_pre_sylph)
    mapseq_dbs = ch_mapseq_dbs
    sr_dbs     = ch_built_sr.mix(ch_pre_sr)   // [ "<name>:<genome|ssu>", refs_fasta ]
    versions   = ch_versions
}
