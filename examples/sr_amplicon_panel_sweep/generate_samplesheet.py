#!/usr/bin/env python3
"""Phase 1: config.yaml -> samplesheet.yaml + genomes/sample_NN.amplicon_16s.csv.

One `--step all` run off this samplesheet trains the error model, draws the
negative-binomial communities, generates their V4 amplicon reads at every subsample
depth, and profiles each depth ONCE against GTDB under `gtdb.map_setting`. It exists to
publish the two things phase 2 reuses:

  <benchmark_dir>/profiling/sr/<id>.<map>.obs.mseq.gz   the reads mapped against GTDB
  error_models/<train_id>/*.model.pt                    the model the reads were made with

The mapping setting is a row-level `sr_settings:` entry, not a benchmark.config pin, so
phase 2's panel entries inherit none of its align knobs.

    python generate_samplesheet.py [config.yaml]
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))
import panel_sweep as ps

HERE = Path(__file__).resolve().parent


def main():
    cfg = ps.load_config(sys.argv[1] if len(sys.argv) > 1 else HERE / "config.yaml")
    panel, train, gtdb = cfg["panel"], cfg["train"], cfg["gtdb"]
    abundances = ps.sample_abundances(cfg)          # (n_samples, n_genomes) integers

    (HERE / "genomes").mkdir(exist_ok=True)
    rows = []
    for gm in ps.generation_modes(cfg):
        reads = ps.mode_reads(cfg, gm)
        for i, abundance in enumerate(abundances, 1):
            csv_path = HERE / "genomes" / f"sample_{i:02d}.{gm['name']}.csv"
            with open(csv_path, "w", newline="") as fh:
                # lineterminator: csv.writer defaults to CRLF, which trips the
                # mixed-line-ending pre-commit hook on the committed CSVs.
                writer = csv.writer(fh, lineterminator="\n")
                writer.writerow(["genome_id", "fasta_path", "abundance"])
                for member, drawn in zip(panel, abundance):
                    if drawn > 0:                   # absent this sample
                        writer.writerow([member["id"], member[gm["source"]], int(drawn)])
            rows.append({
                "sample": f"S{i:02d}.{gm['name']}",
                "train_id": train["id"],
                "train_fastq_1": train["fastq_1"],
                "train_fastq_2": train["fastq_2"],
                "platform": train["platform"],
                "genomes_csv": str(csv_path),
                "num_reads": reads["num_reads"],
                "mode": gm["mode"],
                "paired_end": reads["paired_end"],
                "read_length_mean": reads["read_length_mean"],
                "read_length_variance": reads["read_length_variance"],
                "profilers": ["sr_amplicon"],
                "database": gtdb["name"],
                "sr_settings": [gtdb["map_setting"]],
                **({"subsample": reads["subsample"]} if reads.get("subsample") is not None else {}),
                **({"primers": gm["primers"]} if gm.get("primers") else {}),
            })

    doc = {"databases": ps.databases_block(cfg), "samples": rows}
    with open(HERE / "samplesheet.yaml", "w") as fh:
        ps.dump_yaml(doc, fh)

    print(f"Wrote samplesheet.yaml: {len(rows)} communit(ies) x {len(ps.depths(cfg))} "
          f"depth(s), mapped once against '{gtdb['name']}'; mean genomes present/sample: "
          f"{(abundances > 0).sum(axis=1).mean():.1f}/{len(panel)}")


if __name__ == "__main__":
    main()
