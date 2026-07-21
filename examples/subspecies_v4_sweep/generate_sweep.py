#!/usr/bin/env python3
"""Generate the sub-species V4 abundance-sweep benchmark input set from config.yaml.

Community: N-1 background species (one genome each, fixed at equal abundance) plus
one ADDITIONAL genome of the same species as one background member. Across the
samples only the split between that same-species pair changes (logistic spacing,
denser at the extremes); every other species stays at equal abundance. Reads are V4
amplicon; the error model is trained once from the real reads named in config.yaml.

Reads all paths and metadata from config.yaml (see sweep_config.py for the schema),
so nothing is hard-coded here. Run:
    python generate_sweep.py [config.yaml]
Writes samplesheet.yaml (a `databases:` block the pipeline builds the profiler DB
from, plus one `samples:` row per sweep sample) and genomes/sample_01..NN.csv.
"""
import csv
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))
import sweep_config as sc

HERE = Path(__file__).resolve().parent


def main():
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else HERE / "config.yaml"
    cfg = sc.load_config(cfg_path)

    panel = cfg["panel"]
    # With `primers:`, the genomes CSV points at each member's full `genome:` and the
    # pipeline extracts the amplicon in-silico (per primer pair); otherwise it uses the
    # pre-trimmed `amplicon:` FASTA directly.
    primers = cfg.get("primers")
    fa = {m["id"]: (m["genome"] if primers else m["amplicon"]) for m in panel}
    _dup, major_id, minor_id = sc.sweep_pair(cfg)
    singles = [m for m in panel if m["id"] not in (major_id, minor_id)]

    n = cfg["sweep"]["n_samples"]
    fracs = sc.logistic_fracs(n, cfg["sweep"]["steepness"])
    profiler = cfg["database"]["profilers"][0]  # primary profiler for the combined run
    db_name = cfg["database"]["name"]
    # Optional depth sweep: `reads.subsample` (scalar or list of absolute read counts,
    # or `none`) is written verbatim into every sample row; the pipeline runs one
    # (sub)sample per depth. Omitted => no field => full-depth passthrough.
    subsample = cfg["reads"].get("subsample")

    (HERE / "genomes").mkdir(exist_ok=True)
    rows = []
    for i in range(1, n + 1):
        # Sweep the intra-species split: major 0 -> 1, minor 1 -> 0, always summing
        # to 1 (one species' worth, like every other species).
        a = fracs[i - 1]
        b = 1 - a
        csv_path = HERE / "genomes" / f"sample_{i:02d}.csv"
        with open(csv_path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["genome_id", "fasta_path", "abundance"])
            for m in singles:
                w.writerow([m["id"], fa[m["id"]], 1])
            w.writerow([major_id, fa[major_id], f"{a:.6g}"])
            w.writerow([minor_id, fa[minor_id], f"{b:.6g}"])

        # ponytail: self-check - every species (incl. the doubled pair) has equal weight.
        weights = {m["species"]: 1 for m in singles}
        weights[_dup] = a + b
        per = [v / sum(weights.values()) for v in weights.values()]
        assert all(abs(x - per[0]) < 1e-9 for x in per), f"sample {i}: species not equal-weight"

        rows.append({
            "sample": f"S{i:02d}_a{a:.2f}",
            "train_id": cfg["train"]["id"],
            "train_fastq_1": cfg["train"]["fastq_1"],
            "train_fastq_2": cfg["train"]["fastq_2"],
            "platform": cfg["train"]["platform"],
            "genomes_csv": str(csv_path),
            "num_reads": cfg["reads"]["num_reads"],
            "mode": cfg["reads"]["mode"],
            "paired_end": cfg["reads"]["paired_end"],
            "read_length_mean": cfg["reads"]["read_length_mean"],
            "read_length_variance": cfg["reads"]["read_length_variance"],
            "profiler": profiler,
            "database": db_name,
            **({"subsample": subsample} if subsample is not None else {}),
            **({"primers": primers} if primers else {}),
        })

    # One combined samplesheet: the `databases:` block the pipeline builds the DB
    # from + the sweep samples (training is deduped by train_id, so `--step all`
    # trains once). Set `reads.subsample` in config.yaml to sweep depth per sample.
    doc = {"databases": sc.database_block(cfg), **sc.aap_settings(cfg), "samples": rows}
    with open(HERE / "samplesheet.yaml", "w") as fh:
        yaml.safe_dump(doc, fh, sort_keys=False, default_flow_style=False)

    print(f"Wrote samplesheet.yaml ({len(rows)} samples, database '{db_name}' with "
          f"{len(panel)} sequences, profiler '{profiler}') and "
          f"genomes/sample_01..{n:02d}.csv")
    print(f"Swept intra-species split of {_dup}: {major_id} 0 -> 1, {minor_id} 1 -> 0 "
          f"(sum 1 in every sample); all other species fixed at 1.")


if __name__ == "__main__":
    main()
