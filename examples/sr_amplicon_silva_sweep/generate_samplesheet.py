#!/usr/bin/env python3
"""Phase 1: config.yaml -> samplesheet.yaml + genomes/sample_NN.amplicon_16s.csv.

One `--step all` run off this samplesheet trains the error model, draws the
negative-binomial communities, generates their V4 amplicon reads at every subsample
depth, and profiles each depth with superresolution-amplicon ONCE against the pre-built
reference set. No `sr_settings:` here on purpose - phase 1 exists to produce the two
things phase 2 reuses: the reads (and their truth.tsv), and the mapseq classification of
those reads, published to `<benchmark_dir>/profiling/sr/<id>.obs.mseq.gz`.

Phase 1's own matrix knobs come from benchmark.config, which pins them to the grid's
cheap `exact` point - against a SILVA-sized reference set the nested pipeline's default
(`simulate`) is not an affordable baseline.

    python generate_samplesheet.py [config.yaml]
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))
import silva_sweep as gs

HERE = Path(__file__).resolve().parent


def main():
    cfg = gs.load_config(sys.argv[1] if len(sys.argv) > 1 else HERE / "config.yaml")
    panel = cfg["panel"]
    n = cfg["sampling"]["n_samples"]
    abundances = gs.sample_abundances(cfg)          # (n_samples, n_genomes) integers
    modes = gs.generation_modes(cfg)

    (HERE / "genomes").mkdir(exist_ok=True)
    rows = []
    for gm in modes:
        fasta = {m["id"]: m[gs.mode_source_field(gm)] for m in panel}
        reads = gs.mode_reads(cfg, gm)
        for i in range(1, n + 1):
            abundance = abundances[i - 1]
            assert int(abundance.sum()) > 0, f"sample {i}: empty community"
            csv_path = HERE / "genomes" / f"sample_{i:02d}.{gm['name']}.csv"
            with open(csv_path, "w", newline="") as fh:
                # lineterminator: csv.writer defaults to CRLF, which trips the
                # mixed-line-ending pre-commit hook on the committed CSVs.
                writer = csv.writer(fh, lineterminator="\n")
                writer.writerow(["genome_id", "fasta_path", "abundance"])
                for member, drawn in zip(panel, abundance):
                    if drawn > 0:                   # absent this sample
                        writer.writerow([member["id"], fasta[member["id"]], int(drawn)])
            rows.append(_row(cfg, gm, reads, f"S{i:02d}.{gm['name']}", csv_path))

    doc = {"databases": gs.database_block(cfg), "samples": rows}
    with open(HERE / "samplesheet.yaml", "w") as fh:
        gs.dump_yaml(doc, fh)

    print(f"Wrote samplesheet.yaml: {len(rows)} sample(s) x {len(gs.depths(cfg))} "
          f"depth(s) against '{cfg['database']['name']}'; mean genomes present/sample: "
          f"{(abundances > 0).sum(axis=1).mean():.1f}/{len(panel)}")


def _row(cfg, gm, reads, sample, csv_path):
    train = cfg["train"]
    subsample = reads.get("subsample")
    return {
        "sample": sample,
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
        "profilers": gs.mode_profilers(gm),
        "database": cfg["database"]["name"],
        **({"subsample": subsample} if subsample is not None else {}),
        **({"primers": gm["primers"]} if gm.get("primers") else {}),
    }


if __name__ == "__main__":
    main()
