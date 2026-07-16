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

include { TRAIN_ERROR_MODEL        } from '../subworkflows/local/train_error_model/main'
include { GENOME_BLENDER_GENERATE  } from '../modules/local/genome_blender/generate/main'
include { GROUND_TRUTH             } from '../modules/local/ground_truth/main'
include { PROFILE                  } from '../subworkflows/local/profile/main'

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

        GENOME_BLENDER_GENERATE(ch_gen)
        ch_versions = ch_versions.mix(GENOME_BLENDER_GENERATE.out.versions.first())

        //
        // Ground truth: sort/index the BAM and derive target + realized profiles.
        //
        ch_csv_by_id = ch_samples.map { meta, csv, fastas -> [ meta.id, csv ] }
        ch_truth_in  = GENOME_BLENDER_GENERATE.out.bam
            .map { meta, bam -> [ meta.id, meta, bam ] }
            .combine(ch_csv_by_id, by: 0)
            .map { id, meta, bam, csv -> [ meta, bam, csv ] }

        GROUND_TRUTH(ch_truth_in)
        ch_versions = ch_versions.mix(GROUND_TRUTH.out.versions.first())

        // Feed generated reads (+ reference genomes for database='self') to profiling.
        ch_reads = GENOME_BLENDER_GENERATE.out.reads
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
