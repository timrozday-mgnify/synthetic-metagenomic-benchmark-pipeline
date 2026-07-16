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
RATIO_LO, RATIO_HI = 1e-3, 1e3   # major:minor swept log-spaced 1:1000 -> 1000:1
# ---------------------------------------------------------------------------


def sweep_ratios(n):
    """n log-spaced major:minor ratios from RATIO_LO to RATIO_HI."""
    lo, hi = math.log10(RATIO_LO), math.log10(RATIO_HI)
    return [10 ** (lo + (hi - lo) * i / (n - 1)) for i in range(n)]


def main():
    # Identify the doubled species and its two genome_ids (major = first seen).
    species = [p[1] for p in PANEL]
    dup = next(sp for sp in species if species.count(sp) == 2)
    major_id, minor_id = [p[0] for p in PANEL if p[1] == dup]
    singles = [p for p in PANEL if p[1] != dup]  # 19 background species
    assert len(singles) == 19, "expected exactly one duplicated species (20 species, 21 genomes)"

    (HERE / "genomes").mkdir(exist_ok=True)
    rows = []
    for i, ratio in enumerate(sweep_ratios(N_SAMPLES), start=1):
        major, minor = ratio / (1 + ratio), 1 / (1 + ratio)  # pair sums to 1 = one species' worth
        csv_path = HERE / "genomes" / f"sample_{i:02d}.csv"
        with open(csv_path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["genome_id", "fasta_path", "abundance"])
            for gid, _sp, fa in singles:
                w.writerow([gid, fa, 1])
            major_fa = next(fa for gid, _sp, fa in PANEL if gid == major_id)
            minor_fa = next(fa for gid, _sp, fa in PANEL if gid == minor_id)
            w.writerow([major_id, major_fa, f"{major:.6g}"])
            w.writerow([minor_id, minor_fa, f"{minor:.6g}"])

        # ponytail: self-check — every species normalises to equal weight, pair sums to one.
        weights = {}
        for gid, sp, _fa in singles:
            weights[sp] = weights.get(sp, 0) + 1
        weights[dup] = major + minor
        total = sum(weights.values())
        per = [v / total for v in weights.values()]
        assert all(abs(x - per[0]) < 1e-9 for x in per), f"sample {i}: species not equal-weight"

        sample = f"S{i:02d}_minor{minor:.0e}".replace("-0", "-").replace("+0", "")
        rows.append({
            "sample": sample, "train_id": TRAIN_ID,
            "train_fastq_1": TRAIN_FASTQ_1, "train_fastq_2": TRAIN_FASTQ_2,
            "platform": PLATFORM, "genomes_csv": str(csv_path),
            "num_reads": NUM_READS, "mode": "amplicon", "paired_end": "true",
            "read_length_mean": 300, "read_length_variance": 0,
        })

    # Hand-written YAML (no pyyaml dependency): a list of flat sample maps. Add a
    # `subsample: [none, N, ...]` line per sample to sweep read depths.
    fields = ["train_id", "train_fastq_1", "train_fastq_2", "platform",
              "genomes_csv", "num_reads", "mode", "paired_end",
              "read_length_mean", "read_length_variance"]
    with open(HERE / "samplesheet.yaml", "w") as fh:
        for r in rows:
            fh.write(f"- sample: {r['sample']}\n")
            for k in fields:
                fh.write(f"  {k}: {r[k]}\n")
            fh.write("\n")

    print(f"Wrote samplesheet.yaml ({len(rows)} samples) and genomes/sample_01..{N_SAMPLES:02d}.csv")
    print(f"Swept species: {dup} ({major_id} : {minor_id}), minor frac "
          f"{1/(1+sweep_ratios(N_SAMPLES)[0]):.3g} -> {1/(1+sweep_ratios(N_SAMPLES)[-1]):.3g}")


if __name__ == "__main__":
    main()
