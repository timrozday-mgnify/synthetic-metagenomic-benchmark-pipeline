#!/usr/bin/env python3
"""Generate the negative-binomial abundance-sample benchmark input set from config.yaml.

Community: every panel genome's abundance is drawn independently from a
negative-binomial distribution (mean + dispersion set in config.yaml, optionally
overridden per genome), so each of the `sampling.n_samples` samples is a different
random composition. Draws are seeded (`sampling.seed`) for reproducibility. Each
sample is emitted once per entry in config.yaml's `generation_modes:` (e.g. wgs +
amplicon), sharing the same drawn abundances; the error model is trained once from
the real reads named in config.yaml.

Reads all paths and metadata from config.yaml (see scripts/nb_config.py for the
schema), so nothing is hard-coded here. Run:
    python generate_samplesheet.py [config.yaml]
Writes samplesheet.yaml (a `databases:` block the pipeline builds the profiler DB
from, plus one `samples:` row per sample x generation-mode) and
genomes/sample_NN.<mode>.csv.
"""
import csv
import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))
import nb_config as nc

HERE = Path(__file__).resolve().parent


def main():
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else HERE / "config.yaml"
    cfg = nc.load_config(cfg_path)

    panel = cfg["panel"]
    n = cfg["sampling"]["n_samples"]
    abundances = nc.sample_abundances(cfg)  # (n_samples, n_genomes) integers
    db_name = cfg["database"]["name"]
    modes = nc.generation_modes(cfg)

    (HERE / "genomes").mkdir(exist_ok=True)
    rows = []
    # Every sample is emitted once per generation mode; a mode picks which panel
    # FASTA feeds read generation (genome/ssu/amplicon), the pipeline read mode, its
    # profiler, optional in-silico-PCR primers, and read-param overrides. The drawn
    # abundances are shared across modes (same biology, different read chemistry).
    for gm in modes:
        mname = gm["name"]
        fa = {m["id"]: m[nc.mode_source_field(gm)] for m in panel}
        reads = nc.mode_reads(cfg, gm)
        # Depth sweep: `reads.subsample` (scalar/list of absolute counts, or `none`)
        # written verbatim into each row; omitted => full-depth passthrough.
        subsample = reads.get("subsample")
        for i in range(1, n + 1):
            abund = abundances[i - 1]
            csv_path = HERE / "genomes" / f"sample_{i:02d}.{mname}.csv"
            with open(csv_path, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["genome_id", "fasta_path", "abundance"])
                for m, a in zip(panel, abund):
                    # Skip genomes the NB draw put at zero (absent this sample).
                    if a > 0:
                        w.writerow([m["id"], fa[m["id"]], int(a)])

            # ponytail: self-check - at least one genome present (guaranteed by the
            # sampler's redraw), CSV abundances match the drawn row.
            assert int(abund.sum()) > 0, f"sample {i}: empty community"

            rows.append({
                "sample": f"S{i:02d}.{mname}",
                "train_id": cfg["train"]["id"],
                "train_fastq_1": cfg["train"]["fastq_1"],
                "train_fastq_2": cfg["train"]["fastq_2"],
                "platform": cfg["train"]["platform"],
                "genomes_csv": str(csv_path),
                "num_reads": reads["num_reads"],
                "mode": gm["mode"],
                "paired_end": reads["paired_end"],
                "read_length_mean": reads["read_length_mean"],
                "read_length_variance": reads["read_length_variance"],
                "profiler": gm.get("profiler") or cfg["database"]["profilers"][0],
                "database": db_name,
                **({"subsample": subsample} if subsample is not None else {}),
                **({"primers": gm["primers"]} if gm.get("primers") else {}),
            })

    # One combined samplesheet: the `databases:` block the pipeline builds the DB
    # from + the samples (training is deduped by train_id, so `--step all` trains
    # once). Set `reads.subsample` in config.yaml to sweep depth per sample.
    doc = {"databases": nc.database_block(cfg), **nc.aap_settings(cfg), "samples": rows}
    with open(HERE / "samplesheet.yaml", "w") as fh:
        nc.dump_yaml(doc, fh)

    nb = cfg["sampling"]["negative_binomial"]
    print(f"Wrote samplesheet.yaml ({len(rows)} samples = {n} samples x "
          f"{len(modes)} mode(s) [{', '.join(gm['name'] for gm in modes)}], "
          f"database '{db_name}' with {len(panel)} sequences) and "
          f"genomes/sample_NN.<mode>.csv")
    print(f"Drew abundances ~ NegBinom(mean={nb['mean']}, dispersion={nb['dispersion']}) "
          f"per genome, seed {cfg['sampling'].get('seed', 0)}; "
          f"mean genomes present/sample: {(abundances > 0).sum(axis=1).mean():.1f}/{len(panel)}.")


if __name__ == "__main__":
    main()
