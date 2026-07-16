//
// Synthetic Metagenomic Benchmark: train an error profile from a natural
// metagenome, generate synthetic reads from reference genomes at chosen
// abundances/depth, publish reads + ground-truth BAM + ground-truth profile,
// and (optionally) taxonomically profile the reads next to that ground truth.
//
// params.step selects the phases:
//   generate : train + generate + ground truth (no profiling)
//   profile  : profile reads from existing benchmark dirs (no generation)
//   all      : generate then profile in one run (default)
//

include { TRAIN_ERROR_MODEL           } from '../subworkflows/local/train_error_model/main'
include { GENOME_BLENDER_GENERATE     } from '../modules/local/genome_blender/generate/main'
include { MERGE_GENOME_BLENDER_CHUNKS } from '../modules/local/genome_blender/merge_chunks/main'
include { SEQKIT_SUBSAMPLE            } from '../modules/local/seqkit/subsample/main'
include { GROUND_TRUTH                } from '../modules/local/ground_truth/main'
include { PROFILE                     } from '../subworkflows/local/profile/main'

// Derive a run-specific meta for one subsample depth (scalar null = full-depth
// passthrough). Keeps the original sample id in meta.sample for joins/publishing.
def runMeta(meta, n) {
    meta + [
        sample:         meta.id,
        subsample:      n,
        id:             n != null ? "${meta.id}.sub${n}" : meta.id,
        publish_subdir: n != null ? "subsample_${n}" : '',
    ]
}

// Round-robin remainder distribution so chunk read counts sum exactly to the total.
def chunkReadCounts(total, n) {
    def base = total.intdiv(n)
    def rem  = total % n
    (0..<n).collect { i -> base + (i < rem ? 1 : 0) }
}

// Per-chunk meta: unique id (for distinct output filenames), distinct seed
// (so chunks don't emit identical reads), orig_id to regroup on after generate.
def chunkMeta(meta, i, n_reads_i, n) {
    meta + [
        orig_id:   meta.id,
        id:        "${meta.id}.chunk${i}",
        num_reads: n_reads_i,
        seed:      params.seed + i,
    ]
}

workflow SYNTHETIC_METAGENOMIC_BENCHMARK {
    take:
    ch_samples    // [ val(meta), path(genomes_csv), [ path(fasta) ] ]  (generate/all)
    ch_train      // [ val(meta_train), [ reads ] ]                     (generate/all)
    ch_profile_in // [ val(meta), reads ]                               (profile)

    main:
    ch_versions = Channel.empty()

    ch_reads = Channel.empty()
    ch_aux   = Channel.empty()

    if (params.step != 'profile') {
        //
        // Train one error model + calibration per unique train_id.
        //
        TRAIN_ERROR_MODEL(ch_train)
        ch_versions = ch_versions.mix(TRAIN_ERROR_MODEL.out.versions)

        // Key trained artifacts by train_id for joining back to samples.
        ch_model = TRAIN_ERROR_MODEL.out.model.map       { meta, m -> [ meta.id, m ] }
        ch_cal   = TRAIN_ERROR_MODEL.out.calibration.map  { meta, c -> [ meta.id, c ] }

        //
        // Attach each sample's trained model + calibration by train_id, then generate.
        //
        ch_gen = ch_samples
            .map { meta, csv, fastas -> [ meta.train_id, meta, csv, fastas ] }
            .combine(ch_model, by: 0)
            .combine(ch_cal,   by: 0)
            .map { train_id, meta, csv, fastas, model, cal -> [ meta, csv, fastas, model, cal ] }
            .flatMap { meta, csv, fastas, model, cal ->
                def n = meta.chunks ?: 1
                if (n <= 1) return [[ meta, csv, fastas, model, cal ]]
                def counts = chunkReadCounts(meta.num_reads, n)
                (0..<n).collect { i -> [ chunkMeta(meta, i, counts[i], n), csv, fastas, model, cal ] }
            }

        GENOME_BLENDER_GENERATE(ch_gen)
        ch_versions = ch_versions.mix(GENOME_BLENDER_GENERATE.out.versions.first())

        //
        // Merge chunked generate calls back into one reads-set + BAM per sample
        // (skipped for samples with chunks<=1 — passed straight through).
        //
        ch_chunks  = GENOME_BLENDER_GENERATE.out.reads.join(GENOME_BLENDER_GENERATE.out.bam)
        ch_single  = ch_chunks.filter { meta, reads, bam -> (meta.chunks ?: 1) <= 1 }
        ch_grouped = ch_chunks
            .filter { meta, reads, bam -> (meta.chunks ?: 1) > 1 }
            .map { meta, reads, bam -> [ meta.orig_id, meta, reads, bam ] }
            .groupTuple(by: 0)
            .map { orig_id, metas, reads_lists, bams ->
                def base = metas[0].findAll { k, v -> !(k in ['orig_id', 'seed']) } +
                    [ id: orig_id, num_reads: metas.sum { it.num_reads } ]
                [ base, reads_lists.flatten(), bams ]
            }

        MERGE_GENOME_BLENDER_CHUNKS(ch_grouped)
        ch_versions = ch_versions.mix(MERGE_GENOME_BLENDER_CHUNKS.out.versions.first())

        ch_merged_pairs = MERGE_GENOME_BLENDER_CHUNKS.out.reads
            .join(MERGE_GENOME_BLENDER_CHUNKS.out.bam)
            .mix(ch_single)

        ch_gen_reads = ch_merged_pairs.map { meta, reads, bam -> [ meta, reads ] }
        ch_gen_bam   = ch_merged_pairs.map { meta, reads, bam -> [ meta, bam ] }

        //
        // Subsampling: fan each generated draw out over its requested depths.
        // Each run keeps the original sample id in meta.sample (for joins/publish)
        // and gets a unique meta.id. SEQKIT_SUBSAMPLE emits the (sub)reads plus a
        // names.txt so the ground-truth BAM is filtered to the same reads.
        //
        ch_sub_in = ch_gen_reads.flatMap { meta, reads ->
            meta.subsamples.collect { n -> [ runMeta(meta, n), reads ] }
        }
        SEQKIT_SUBSAMPLE(ch_sub_in)
        ch_versions = ch_versions.mix(SEQKIT_SUBSAMPLE.out.versions.first())

        //
        // Ground truth per run: filter the full BAM (keyed by sample) to this
        // run's reads (names.txt), then derive target + realized profiles.
        //
        ch_csv_by_id = ch_samples.map { meta, csv, fastas -> [ meta.id, csv ] }
        ch_full_bam  = ch_gen_bam.map { meta, bam -> [ meta.id, bam ] }
        ch_truth_in  = SEQKIT_SUBSAMPLE.out.names
            .map { meta, names -> [ meta.sample, meta, names ] }
            .combine(ch_full_bam,  by: 0)
            .combine(ch_csv_by_id, by: 0)
            .map { sample, meta, names, bam, csv -> [ meta, bam, csv, names ] }

        GROUND_TRUTH(ch_truth_in)
        ch_versions = ch_versions.mix(GROUND_TRUTH.out.versions.first())

        // Feed subsampled reads (+ reference genomes for database='self') to profiling.
        ch_reads = SEQKIT_SUBSAMPLE.out.reads
        ch_aux   = ch_samples.map { meta, csv, fastas -> [ meta.id, csv, fastas ] }
    }
    else {
        ch_reads = ch_profile_in
    }

    //
    // Profiling (all + profile steps).
    //
    if (params.step != 'generate') {
        PROFILE(ch_reads, ch_aux)
        ch_versions = ch_versions.mix(PROFILE.out.versions)
    }

    //
    // Collect software versions.
    //
    ch_versions
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info")

    emit:
    versions = ch_versions
}
