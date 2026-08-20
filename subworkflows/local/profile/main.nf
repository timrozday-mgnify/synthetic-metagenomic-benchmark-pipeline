//
// Profile generated (or previously generated) reads and drop the predicted
// profile next to the ground truth. Per-sample profiler is chosen by meta.profiler:
//   - sylph : WGS. DB is either a config entry (params.sylph_databases[<database>])
//             or 'self' (built from the sample's own reference genomes).
//   - aap   : amplicon-analysis-pipeline, run via a nested `nextflow run`.
//   - sr_amplicon / sr_shotgun : superresolution, also a nested `nextflow run`. Its
//             "DB" is one combined reference FASTA, built from the sample's own
//             genomes ('self') or from a named collection.
//

include { SYLPH_BUILD_DB      } from '../../../modules/local/sylph/build_db/main'
include { SYLPH_PROFILE       } from '../../../modules/nf-core/sylph/profile/main'
include { NORMALIZE_SYLPH     } from '../../../modules/local/sylph/normalize/main'
include { RUN_AAP             } from '../../../modules/local/amplicon_analysis/main'
include { SR_BUILD_REFS       } from '../../../modules/local/superresolution/build_refs/main'
include { BUILD_SUPERRESOLUTION_MISMAPPING } from '../../../modules/local/superresolution/run/main'
include { RUN_SUPERRESOLUTION               } from '../../../modules/local/superresolution/run/main'

// Panel FASTA field each superresolution flavour's references are built from
// (keep in sync with the same map in build_databases).
def srSources() { [ sr_shotgun: 'genome', sr_amplicon: 'ssu' ] }

workflow PROFILE {
    take:
    ch_reads         // [ val(meta), reads ]                 meta: id, mode, profiler, database
    ch_aux           // [ id, genomes_csv, [ fasta ] ]       (empty in profile-only step)
    ch_sylph_dbs     // [ name, syldb ]                      built/prebuilt sylph DBs by name
    ch_mapseq_dbs    // [ name, fasta, tax, otu, mscluster ] built/prebuilt mapseq DBs by name
    ch_sr_dbs        // [ "<name>:<genome|ssu>", refs_fasta ] built/prebuilt superresolution refs
    builtNames       // [ profiler: Set of collection names resolved for it ]

    main:
    ch_versions = Channel.empty()
    def no_file = file("${projectDir}/assets/NO_FILE")

    ch_reads
        .branch { meta, reads ->
            sylph: meta.profiler == 'sylph'
            aap:   meta.profiler == 'aap'
            sr:    meta.profiler in srSources().keySet()
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
    // combine, not join: join is 1:1 and consumes the key, so with extra_profilers
    // fanning a sample into several entries it would silently drop all but one.
    ch_self = ch_by_prof.sylph
        .filter { it[0].database == 'self' }
        .map { meta, reads -> [ meta.sample ?: meta.id, meta, reads ] }
        .combine(ch_aux, by: 0)
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
            built:  meta.database in builtNames.sylph
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
            built: meta.database in builtNames.aap
            other: true
        }
        .set { ch_aap }

    ch_aap_built = ch_aap.built
        .map { meta, reads -> [ meta.database, meta, reads ] }
        .combine(ch_mapseq_dbs, by: 0)
        .map { name, meta, reads, fasta, tax, otu, mscluster, rfam_cm, rfam_claninfo ->
            [ meta, reads, true, no_file, fasta, tax, otu, mscluster, rfam_cm, rfam_claninfo ]
        }

    // Passthrough branch supplies its own rfam DBs via params.aap_config, so the
    // rfam slots are empty (write_aap_config.py isn't run for this branch).
    ch_aap_other = ch_aap.other.map { meta, reads ->
        [ meta, reads, false, (params.aap_config ? file(params.aap_config, checkIfExists: true) : no_file),
          no_file, no_file, no_file, no_file, '', '' ]
    }

    // Batch samples that share a DB config into one nested AAP run (AAP's samplesheet
    // is multi-row). Group key = everything that fixes the nested command + its -c files:
    // the DB (use_built + database), the extra engine configs, and the -profile. The DB
    // slots are identical within a group, so we keep one copy ([0]); metas + layout carry
    // the per-sample rows. paired/single-end is NOT part of the key — it's per-row in the
    // AAP samplesheet.
    ch_aap_grouped = ch_aap_built.mix(ch_aap_other)
        .map { meta, reads, use_built, aap_config, fasta, tax, otu, mscluster, rfam_cm, rfam_claninfo ->
            def key = [ use_built, meta.database, meta.aap_configs, meta.aap_profile ]
            [ key, meta, reads, use_built, aap_config, fasta, tax, otu, mscluster, rfam_cm, rfam_claninfo ]
        }
        .groupTuple(by: 0)
        .map { key, metas, readsL, useBuilts, aapConfigs, fastas, taxs, otus, msclusters, rfamCms, rfamClaninfos ->
            // layout: one [id, single_end, fastq_1, fastq_2] per sample. Reads are passed as
            // absolute paths (val), not staged — see RUN_AAP (executor local, avoids sub_* collisions).
            def layout = [metas, readsL].transpose().collect { m, r ->
                def rl = r instanceof List ? r : [r]
                [ m.id, rl.size() > 1 ? 'false' : 'true', rl[0].toString(), rl.size() > 1 ? rl[1].toString() : '' ]
            }
            [ metas, layout, useBuilts[0], aapConfigs[0], fastas[0], taxs[0], otus[0], msclusters[0], rfamCms[0], rfamClaninfos[0] ]
        }

    RUN_AAP(ch_aap_grouped)
    ch_versions = ch_versions.mix(RUN_AAP.out.versions.first())

    // Demux the batched output back to per-sample [meta, files] for the emitted contract
    // (AAP namespaces files under aap_out/<id>/; match by that path segment). Publishing
    // itself is per-sample via the RUN_AAP publishDir in modules.config.
    ch_aap_out = RUN_AAP.out.results
        .flatMap { metas, files ->
            def fl = files instanceof List ? files : [files]
            metas.collect { m -> [ m, fl.findAll { it.toString().contains("/aap_out/${m.id}/") } ] }
        }

    //
    // superresolution (nested nextflow run). Its only "database" is a combined
    // reference FASTA: either built from the sample's own genomes ('self', needs
    // ch_aux) or from a named collection built by BUILD_DATABASES. There is no
    // params-configured fallback — an unknown name is an error.
    //
    ch_by_prof.sr
        .branch { meta, reads ->
            self:  meta.database == 'self' || !meta.database
            built: meta.database in builtNames[meta.profiler]
            other: true
        }
        .set { ch_sr }

    // Fail loudly rather than silently dropping the sample.
    ch_sr.other.map { meta, reads ->
        error "No superresolution reference collection '${meta.database}' in the samplesheet databases: block for profiler '${meta.profiler}' (sample ${meta.id}); use 'self' or define it"
    }

    // 'self': reuse the sample's genomes_csv + staged FASTAs (same inner-join as sylph;
    // a 'self' row without reference genomes — profile-only — has nothing to join).
    ch_sr_self = ch_sr.self
        .map { meta, reads -> [ meta.sample ?: meta.id, meta, reads ] }
        .combine(ch_aux, by: 0)
        .map { id, meta, reads, csv, fastas -> [ meta, reads, csv, fastas ] }
    SR_BUILD_REFS(ch_sr_self.map { meta, reads, csv, fastas -> [ meta, csv, fastas, '' ] })
    ch_versions = ch_versions.mix(SR_BUILD_REFS.out.versions.first())

    // Keyed by id+profiler: a sample may run both superresolution flavours, and each
    // gets its own reference set (genomes vs 16S).
    ch_sr_self_in = ch_sr_self
        .map { meta, reads, csv, fastas -> [ "${meta.id}:${meta.profiler}".toString(), meta, reads ] }
        .join(SR_BUILD_REFS.out.refs.map { meta, refs -> [ "${meta.id}:${meta.profiler}".toString(), refs ] }, by: 0)
        // A self reference set belongs to the source sample, so all of its
        // subsampling depths reuse one matrix while distinct samples remain isolated.
        .map { key, meta, reads, refs -> [ "self:${meta.sample ?: meta.id}:${meta.profiler}", meta, reads, refs ] }

    // Named collection: join by "<name>:<source>", the key BUILD_DATABASES emits.
    ch_sr_built_in = ch_sr.built
        .map { meta, reads -> [ "${meta.database}:${srSources()[meta.profiler]}".toString(), meta, reads ] }
        .combine(ch_sr_dbs, by: 0)
        // Named collections are reference sets shared by every matching sample.
        .map { key, meta, reads, refs -> [ "${meta.database}:${meta.profiler}", meta, reads, refs ] }

    ch_sr_runs = ch_sr_self_in.mix(ch_sr_built_in)

    // Materialise exactly one matrix per reference set. The nested pipelines need a
    // sample-shaped input to build the simulation matrix, so select one representative
    // run; all subsequent sample runs reuse its matrix through --mismapping_matrix.
    ch_sr_mismapping_in = ch_sr_runs
        .groupTuple(by: 0)
        .map { referenceSet, metas, readsList, refsList ->
            def representative = metas[0] + [
                id: "mismapping_${referenceSet.replaceAll(/[^A-Za-z0-9._-]+/, '_')}",
                reference_set: referenceSet,
                reference_set_dir: referenceSet.replaceAll(/[^A-Za-z0-9._-]+/, '_'),
            ]
            [ representative, (readsList[0] instanceof List ? readsList[0] : [readsList[0]])*.toString(), refsList[0] ]
        }
    BUILD_SUPERRESOLUTION_MISMAPPING(ch_sr_mismapping_in)
    ch_versions = ch_versions.mix(BUILD_SUPERRESOLUTION_MISMAPPING.out.versions.first())

    // Reads go through as absolute path strings (val) — see RUN_SUPERRESOLUTION.
    // Joining restores the one shared matrix to every run in its reference set.
    RUN_SUPERRESOLUTION(
        ch_sr_runs
            .map { referenceSet, meta, reads, refs -> [ referenceSet, meta, reads, refs ] }
            .join(BUILD_SUPERRESOLUTION_MISMAPPING.out.mismapping.map { meta, matrix -> [ meta.reference_set, matrix ] }, by: 0)
            .map { referenceSet, meta, reads, refs, matrix ->
                [ meta, (reads instanceof List ? reads : [reads])*.toString(), refs, matrix ]
            }
    )
    ch_versions = ch_versions.mix(RUN_SUPERRESOLUTION.out.versions.first())

    emit:
    sylph    = NORMALIZE_SYLPH.out.profile      // [ meta, sylph_profile.tsv ]
    aap      = ch_aap_out                       // [ meta, aap_out/<id> ]
    sr       = RUN_SUPERRESOLUTION.out.profile  // [ meta, sr_profile.tsv ]
    versions = ch_versions
}
