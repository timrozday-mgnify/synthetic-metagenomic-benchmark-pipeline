//
// Profile generated (or previously generated) reads and drop the predicted
// profile next to the ground truth. Per-sample profiler is chosen by meta.profiler:
//   - sylph : WGS. DB is either a config entry (params.sylph_databases[<database>])
//             or 'self' (built from the sample's own reference genomes).
//   - aap   : amplicon-analysis-pipeline, run via a nested `nextflow run`.
//

include { SYLPH_BUILD_DB  } from '../../../modules/local/sylph/build_db/main'
include { SYLPH_PROFILE   } from '../../../modules/nf-core/sylph/profile/main'
include { NORMALIZE_SYLPH } from '../../../modules/local/sylph/normalize/main'
include { RUN_AAP         } from '../../../modules/local/amplicon_analysis/main'

workflow PROFILE {
    take:
    ch_reads         // [ val(meta), reads ]                 meta: id, mode, profiler, database
    ch_aux           // [ id, genomes_csv, [ fasta ] ]       (empty in profile-only step)
    ch_sylph_dbs     // [ name, syldb ]                      built/prebuilt sylph DBs by name
    ch_mapseq_dbs    // [ name, fasta, tax, otu, mscluster ] built/prebuilt mapseq DBs by name
    builtSylphNames  // Set of collection names resolved for sylph
    builtMapseqNames // Set of collection names resolved for mapseq/aap

    main:
    ch_versions = Channel.empty()
    def no_file = file("${projectDir}/assets/NO_FILE")

    ch_reads
        .branch { meta, reads ->
            sylph: meta.profiler == 'sylph'
            aap:   meta.profiler == 'aap'
            other: true
        }
        .set { ch_by_prof }

    //
    // sylph
    //
    // self DB: build a .syldb from the sample's reference genomes. This needs the
    // genomes_csv + fastas from ch_aux (only present in the generate/all step), so
    // we inner-join by id — a 'self' row without reference genomes (profile-only)
    // simply has nothing to join and is dropped (see README).
    // ponytail: self DB is rebuilt per run (keyed by unique meta.id) even though
    // it depends only on the genomes; dedupe by meta.sample if it ever matters.
    ch_self = ch_by_prof.sylph
        .filter { it[0].database == 'self' }
        .map { meta, reads -> [ meta.sample ?: meta.id, meta, reads ] }
        .join(ch_aux, by: 0)
        .map { id, meta, reads, csv, fastas -> [ meta, reads, csv, fastas ] }
    SYLPH_BUILD_DB( ch_self.map { meta, reads, csv, fastas -> [ meta, fastas ] } )
    ch_versions = ch_versions.mix(SYLPH_BUILD_DB.out.versions.first())

    ch_self_in = ch_self
        .map { meta, reads, csv, fastas -> [ meta.id, meta, reads, csv ] }
        .join(SYLPH_BUILD_DB.out.db.map { meta, db -> [ meta.id, db ] }, by: 0)
        .map { id, meta, reads, csv, db -> [ meta, reads, db, csv ] }

    // named DB: a `database` defined under the samplesheet `databases:` block is
    // built (or its prebuilt dir resolved) by BUILD_DATABASES and joined here by
    // name; any other name falls back to params.sylph_databases.
    ch_by_prof.sylph
        .filter { it[0].database && it[0].database != 'self' }
        .branch { meta, reads ->
            built:  meta.database in builtSylphNames
            config: true
        }
        .set { ch_named }

    ch_built_in = ch_named.built
        .map { meta, reads -> [ meta.database, meta, reads ] }
        .combine(ch_sylph_dbs, by: 0)
        .map { name, meta, reads, db -> [ meta, reads, db, no_file ] }

    // config DB: resolve the .syldb path from params.sylph_databases (no ch_aux needed).
    ch_cfg_in = ch_named.config
        .map { meta, reads ->
            def entry = (params.sylph_databases ?: [:])[meta.database]
            if (!entry?.syldb) error "No sylph database '${meta.database}' in params.sylph_databases or samplesheet databases: block (sample ${meta.id})"
            [ meta, reads, file(entry.syldb, checkIfExists: true), no_file ]
        }

    ch_sylph_in = ch_self_in.mix(ch_built_in).mix(ch_cfg_in)

    // SYLPH_PROFILE takes reads + db as two inputs; multiMap keeps them aligned.
    ch_prof = ch_sylph_in.multiMap { meta, reads, db, csv ->
        reads: [ meta, reads ]
        db:    db
    }
    SYLPH_PROFILE(ch_prof.reads, ch_prof.db)
    ch_versions = ch_versions.mix(SYLPH_PROFILE.out.versions.first())

    // Normalise to a genome_id profile next to truth.tsv (attach genomes_csv / NO_FILE).
    ch_norm_in = SYLPH_PROFILE.out.profile
        .map { meta, tsv -> [ meta.id, meta, tsv ] }
        .join(ch_sylph_in.map { meta, reads, db, csv -> [ meta.id, csv ] }, by: 0)
        .map { id, meta, tsv, csv -> [ meta, tsv, csv ] }
    NORMALIZE_SYLPH(ch_norm_in)
    ch_versions = ch_versions.mix(NORMALIZE_SYLPH.out.versions.first())

    //
    // aap (nested nextflow run). A `database` naming a built/prebuilt mapseq
    // collection attaches its fasta/tax/otu/mscluster (RUN_AAP writes the aap.config);
    // any other name keeps the params.aap_config pass-through.
    //
    ch_by_prof.aap
        .branch { meta, reads ->
            built: meta.database in builtMapseqNames
            other: true
        }
        .set { ch_aap }

    ch_aap_built = ch_aap.built
        .map { meta, reads -> [ meta.database, meta, reads ] }
        .combine(ch_mapseq_dbs, by: 0)
        .map { name, meta, reads, fasta, tax, otu, mscluster ->
            [ meta, reads, true, no_file, fasta, tax, otu, mscluster ]
        }

    ch_aap_other = ch_aap.other.map { meta, reads ->
        [ meta, reads, false, (params.aap_config ? file(params.aap_config, checkIfExists: true) : no_file),
          no_file, no_file, no_file, no_file ]
    }

    RUN_AAP(ch_aap_built.mix(ch_aap_other))
    ch_versions = ch_versions.mix(RUN_AAP.out.versions.first())

    emit:
    sylph    = NORMALIZE_SYLPH.out.profile // [ meta, sylph_profile.tsv ]
    aap      = RUN_AAP.out.results         // [ meta, aap_out/** ]
    versions = ch_versions
}
