#!/usr/bin/env python3
"""Phase 1: config.yaml -> samplesheet.yaml + genomes/sample_NN.amplicon_16s.csv.

One `--step all` run off this samplesheet trains the error model(s), draws the
negative-binomial communities, generates their V4 amplicon reads at every subsample
depth, and profiles each depth with superresolution-amplicon ONCE against the pre-built
reference set. The communities are generated once per `error_models:` arm - the reads
themselves differ, so each arm is its own set of samples with its own train_id. Phase 1
exists to produce the two things phase 2 reuses: the reads (and their truth.tsv), and the
mapseq classification of those reads, published to
`<benchmark_dir>/profiling/sr/<id>.map.obs.mseq.gz`.

Its one `sr_settings:` entry, `map`, is `arm_point`'s knobs over the panel collection:
against SILVA a run without a panel would build its kernel over every SILVA V4 group.

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
    error_arms = gs.error_models(cfg)

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
            # Same community, one sample per error-model arm: the abundances are shared
            # so the arms differ only in how the reads were damaged.
            for em in error_arms:
                rows.append(_row(cfg, gm, reads, em,
                                 gs.sample_id(cfg, em, gm, i), csv_path))

    doc = {"databases": gs.database_block(cfg), "samples": rows}
    with open(HERE / "samplesheet.yaml", "w") as fh:
        gs.dump_yaml(doc, fh)

    arms = ", ".join(e["name"] or "default" for e in error_arms)
    print(f"Wrote samplesheet.yaml: {len(rows)} sample(s) x {len(gs.depths(cfg))} "
          f"depth(s) against '{cfg['database']['name']}'; error-model arms: {arms}; "
          f"mean genomes present/sample: "
          f"{(abundances > 0).sum(axis=1).mean():.1f}/{len(panel)}")


def _row(cfg, gm, reads, em, sample, csv_path):
    train = cfg["train"]
    subsample = reads.get("subsample")
    return {
        "sample": sample,
        # One train_id per error-model arm: the pipeline dedupes training by it, so two
        # arms sharing an id would silently share one model. A `preset:` arm trains
        # nothing at all and ignores the FASTQs below; they stay for the other arms.
        "train_id": em["train_id"],
        "train_fastq_1": train["fastq_1"],
        "train_fastq_2": train["fastq_2"],
        "platform": train["platform"],
        **({"error_model": em["preset"]} if em.get("preset") else {}),
        **({"error_model_components": em["components"]} if em.get("components") else {}),
        "genomes_csv": str(csv_path),
        "num_reads": reads["num_reads"],
        "mode": gm["mode"],
        "paired_end": reads["paired_end"],
        "read_length_mean": reads["read_length_mean"],
        "read_length_variance": reads["read_length_variance"],
        "profilers": gs.mode_profilers(gm),
        "database": cfg["database"]["name"],
        "sr_settings": [gs.map_setting(cfg)],
        **({"subsample": subsample} if subsample is not None else {}),
        **({"primers": gm["primers"]} if gm.get("primers") else {}),
    }


if __name__ == "__main__":
    main()
