#!/usr/bin/env python3
"""Generate the sub-species V4 abundance-sweep benchmark input set.

Community: 20 background species (one genome each, fixed at equal abundance) plus
one ADDITIONAL genome of the same species as one background member. Across 20
samples only the split between that same-species pair changes (log-spaced); every
other species stays at equal abundance. Reads are V4 amplicon, error model trained
once from the real SC2200627 Illumina reads.

Fill in the PANEL paths (V4 amplicon FASTAs) and TRAIN_FASTQ_* below, then run:
    python generate_sweep.py
Writes samplesheet.yaml and genomes/sample_01..20.csv next to this script.
"""
import csv
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
# Must match run.sh's --outdir; error_model_dir points inside it (see generate_profile_samplesheet.py).
RESULTS_DIR = HERE.parent.parent / "results" / "subspecies_v4_sweep"

# --- Fill these in ---------------------------------------------------------
# Real reads the sequencing-error model is trained from (once, shared by all samples).
TRAIN_FASTQ_1 = "/Users/timrozday/Downloads/SC2200627-SC3/SC2200627-SC3_Batch_54_samples_26s000943-1-1_Finn_lane1C5_1_sequence.txt.gz"
TRAIN_FASTQ_2 = "/Users/timrozday/Downloads/SC2200627-SC3/SC2200627-SC3_Batch_54_samples_26s000943-1-1_Finn_lane1C5_2_sequence.txt.gz"

# Pre-trimmed 515-YF/806BR V4 amplicon FASTAs (mode=amplicon uses each record directly).
V4 = "/Users/timrozday/Documents/mimicc/mimicc-primer-investigations/amplicon-primer-screen_runs/20HM+/results/seqkit_amplicon/{s}/515-YF-806BR/{s}.amplicons.fasta"

# 21 genomes = 20 background species + 1 additional genome of an existing species.
# (genome_id, species, v4_fasta_path). Exactly one `species` appears twice: that pair
# is the sweep target. genome_id must be unique (it names the truth-table key / BAM prefix).
PANEL = [
    (s, s, V4.format(s=s)) for s in [
        "bacteroides_fragilis", "bacteroides_thetaiotaomicron", "bacteroides_uniformis",
        "bacteroides_vulgatus", "clostridium_bolteae", "clostridium_perfringens",
        "clostridium_ramosum", "clostridium_saccharolyticum", "collinsella_aerofaciens",
        "coprococcus_comes", "dorea_formicigenerans", "eggerthella_lenta",
        "eubacterium_rectale", "fusobacterium_nucleatum", "methanobrevibacter_smithii",
        "parabacteroides_merdae", "roseburia_intestinalis", "ruminococcus_gnavus",
        "streptococcus_parasanguinis", "streptococcus_salivarius",
    ]
] + [
    # The additional genome: same species as the first background member above.
    # Replace the path with the second strain's V4 amplicon FASTA.
    ("bacteroides_uniformis_strain2", "bacteroides_uniformis",
     "/PATH/TO/bacteroides_uniformis_strain2.v4.amplicons.fasta"),
]

PLATFORM = "hq-illumina"   # 2x300 MiSeq amplicon tails may fit 'lq-illumina' better
NUM_READS = 1_000_000      # 500k pairs x 2 (pipeline counts paired reads as pairs*2)
N_SAMPLES = 20
TRAIN_ID = "sc2200627"
SWEEP_STEEPNESS = 6.0   # higher -> more samples bunched near 0 and 1
# ---------------------------------------------------------------------------


def logistic_fracs(n, k):
    """n fractions in [0,1], symmetric, denser near the 0/1 extremes (logistic
    spacing). Endpoints are rescaled to land exactly on 0 and 1."""
    s = [1 / (1 + math.exp(-k * (2 * i / (n - 1) - 1))) for i in range(n)]
    lo, hi = s[0], s[-1]
    return [(v - lo) / (hi - lo) for v in s]


def main():
    # Identify the doubled species and its two genome_ids.
    species = [p[1] for p in PANEL]
    dup = next(sp for sp in species if species.count(sp) == 2)
    strain_ids = [p[0] for p in PANEL if p[1] == dup]
    singles = [p for p in PANEL if p[1] != dup]  # 19 background species
    assert len(singles) == 19, "expected exactly one duplicated species (20 species, 21 genomes)"

    strain_a, strain_b = strain_ids
    strain_fa = {gid: next(f for g, _sp, f in PANEL if g == gid) for gid in strain_ids}

    fracs = logistic_fracs(N_SAMPLES, SWEEP_STEEPNESS)

    (HERE / "genomes").mkdir(exist_ok=True)
    rows = []
    for i in range(1, N_SAMPLES + 1):
        # Sweep the intra-species split (logistic, denser at the extremes; endpoints
        # exactly 0/1): strain_a 0 -> 1, strain_b 1 -> 0, always summing to 1 (one
        # species' worth, like every other species).
        a = fracs[i - 1]
        b = 1 - a
        csv_path = HERE / "genomes" / f"sample_{i:02d}.csv"
        with open(csv_path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["genome_id", "fasta_path", "abundance"])
            for gid, _sp, fa in singles:
                w.writerow([gid, fa, 1])
            w.writerow([strain_a, strain_fa[strain_a], f"{a:.6g}"])
            w.writerow([strain_b, strain_fa[strain_b], f"{b:.6g}"])

        # ponytail: self-check — every species (incl. the doubled pair) has equal weight.
        weights = {}
        for gid, sp, _fa in singles:
            weights[sp] = weights.get(sp, 0) + 1
        weights[dup] = a + b
        total = sum(weights.values())
        per = [v / total for v in weights.values()]
        assert all(abs(x - per[0]) < 1e-9 for x in per), f"sample {i}: species not equal-weight"

        rows.append({
            "sample": f"S{i:02d}_a{a:.2f}", "train_id": TRAIN_ID,
            # run.sh trains once (--step train) into this dir; generate reuses it.
            "error_model_dir": str(RESULTS_DIR / "error_models" / TRAIN_ID),
            "platform": PLATFORM, "genomes_csv": str(csv_path),
            "num_reads": NUM_READS, "mode": "amplicon", "paired_end": "true",
            "read_length_mean": 300, "read_length_variance": 0,
        })

    # Train-only samplesheet (one row per train_id): consumed by `--step train`.
    with open(HERE / "train_samplesheet.yaml", "w") as fh:
        fh.write(f"- train_id: {TRAIN_ID}\n")
        fh.write(f"  train_fastq_1: {TRAIN_FASTQ_1}\n")
        fh.write(f"  train_fastq_2: {TRAIN_FASTQ_2}\n")
        fh.write(f"  platform: {PLATFORM}\n")

    # Hand-written YAML (no pyyaml dependency): a list of flat sample maps. Add a
    # `subsample: [none, N, ...]` line per sample to sweep read depths.
    fields = ["train_id", "error_model_dir", "platform",
              "genomes_csv", "num_reads", "mode", "paired_end",
              "read_length_mean", "read_length_variance"]
    with open(HERE / "samplesheet.yaml", "w") as fh:
        for r in rows:
            fh.write(f"- sample: {r['sample']}\n")
            for k in fields:
                fh.write(f"  {k}: {r[k]}\n")
            fh.write("\n")

    print(f"Wrote train_samplesheet.yaml, samplesheet.yaml ({len(rows)} samples) "
          f"and genomes/sample_01..{N_SAMPLES:02d}.csv")
    print(f"Swept intra-species split of {dup}: {strain_a} 0 -> 1, {strain_b} 1 -> 0 "
          f"(sum 1 in every sample); all other species fixed at 1.")


if __name__ == "__main__":
    main()
