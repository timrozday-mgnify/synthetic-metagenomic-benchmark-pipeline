#!/usr/bin/env python3
"""Generate a profile-only samplesheet (`--step profile`) to re-profile reads
already written by an earlier run, without regenerating them.

Reuses the same config.yaml as generate_sweep.py: it emits the same `databases:`
block (so the pipeline builds the `community_v4` DB) and one `samples:` row per
already-generated sample, re-profiling each with the profiler its generation mode
used (read from samplesheet.yaml), so a wgs sample re-runs sylph and an amplicon
sample re-runs aap - without regenerating the reads.

    python generate_profile_samplesheet.py [results_dir] [config.yaml]
    nextflow run ../../main.nf -profile docker --step profile \\
        --input profile_samplesheet.yaml --outdir <results_dir>
"""
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))
import sweep_config as sc

HERE = Path(__file__).resolve().parent
DEFAULT_RESULTS_DIR = HERE.parent.parent / "results" / "subspecies_v4_sweep"


def generated_samples(samplesheet):
    """(sample_id, profiler) for each row of the generate samplesheet.yaml, so each
    sample re-profiles with the profiler its generation mode used."""
    if not samplesheet.exists():
        sys.exit(f"{samplesheet} not found - run generate_sweep.py first.")
    doc = yaml.safe_load(samplesheet.read_text())
    return [(s["sample"], s["profiler"]) for s in doc["samples"]]


def main():
    results_dir = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else DEFAULT_RESULTS_DIR
    cfg_path = sys.argv[2] if len(sys.argv) > 2 else HERE / "config.yaml"
    cfg = sc.load_config(cfg_path)

    samples = generated_samples(HERE / "samplesheet.yaml")
    db_name = cfg["database"]["name"]

    rows = [
        {"sample": sample, "profiler": profiler,
         "benchmark_dir": str(results_dir / sample), "database": db_name}
        for sample, profiler in samples
    ]
    doc = {"databases": sc.database_block(cfg), **sc.aap_settings(cfg), "samples": rows}
    with open(HERE / "profile_samplesheet.yaml", "w") as fh:
        yaml.safe_dump(doc, fh, sort_keys=False, default_flow_style=False)

    print(f"Wrote profile_samplesheet.yaml: {len(rows)} rows (one per generated "
          f"sample, each with its generation-mode profiler) against '{db_name}', "
          f"benchmark_dir root {results_dir}")


if __name__ == "__main__":
    main()
