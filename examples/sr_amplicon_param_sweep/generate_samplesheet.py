#!/usr/bin/env python3
"""Phase 1: config.yaml -> samplesheet.yaml + genomes/community.csv.

One `--step all` run off this samplesheet trains the error model, generates the
amplicon reads for every subsample depth, builds the `community_v4` reference set, and
profiles each depth with superresolution-amplicon ONCE, at the nested pipeline's own
default settings. That baseline run is what produces the two things phase 2 reuses:
the reads (and their truth.tsv), and the mapseq classification of those reads, published
to `<benchmark_dir>/profiling/sr/<id>.obs.mseq.gz`.

No `sr_settings:` here on purpose - phase 1 exists to map the reads, and the sweep is
what re-profiles them.

    python generate_samplesheet.py [config.yaml]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))
import sr_sweep as sw

HERE = Path(__file__).resolve().parent


def main():
    cfg = sw.load_config(sys.argv[1] if len(sys.argv) > 1 else HERE / "config.yaml")
    reads, train = cfg["reads"], cfg["train"]

    # One row per sweep sample: a flat equal-abundance `community` with no `sweep:`
    # block, otherwise `sweep.n_samples` samples varying the doubled species' split.
    (HERE / "genomes").mkdir(exist_ok=True)
    rows = []
    for sample, abundance in sw.samples(cfg):
        csv = HERE / "genomes" / f"{sample}.csv"
        with open(csv, "w") as fh:
            fh.write("genome_id,fasta_path,abundance\n")
            for m in cfg["panel"]:
                fh.write(f"{m['id']},{sw.resolve(cfg, m['ssu'])},{abundance[m['id']]:.6g}\n")
        rows.append(_row(cfg, sample, csv))

    doc = {"databases": sw.database_block(cfg), "samples": rows}
    with open(HERE / "samplesheet.yaml", "w") as fh:
        sw.dump_yaml(doc, fh)

    print(f"Wrote samplesheet.yaml: {len(rows)} sample(s) x {len(cfg['primers'])} "
          f"primer pair(s) x {len(sw.depths(cfg))} depth(s), {len(cfg['panel'])} "
          f"genomes in '{cfg['database']['name']}'")


def _row(cfg, sample, csv):
    reads, train = cfg["reads"], cfg["train"]
    return {
        "sample": sample,
        "train_id": train["id"],
        "train_fastq_1": sw.resolve(cfg, train["fastq_1"]),
        "train_fastq_2": sw.resolve(cfg, train["fastq_2"]),
        "platform": train["platform"],
        "genomes_csv": str(csv),
        "num_reads": reads["num_reads"],
        "mode": "amplicon",
        "paired_end": reads["paired_end"],
        "read_length_mean": reads["read_length_mean"],
        "read_length_variance": reads["read_length_variance"],
        "subsample": [d if d is not None else "none" for d in sw.depths(cfg)],
        "primers": cfg["primers"],
        "profilers": ["sr_amplicon"],
        "database": cfg["database"]["name"],
    }


if __name__ == "__main__":
    main()
