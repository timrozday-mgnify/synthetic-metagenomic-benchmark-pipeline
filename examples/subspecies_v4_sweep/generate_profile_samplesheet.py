#!/usr/bin/env python3
"""Generate a profile-only samplesheet (`--step profile`) to re-profile reads
already written by an earlier run, without regenerating them.

Reuses the same config.yaml as generate_sweep.py: it emits the same `databases:`
block (so the pipeline builds the `community_v4` DB) and one `samples:` row per
(sample, profiler) - i.e. every profiler in `database.profilers` is run against the
already-generated reads, useful for comparing profilers without regeneration.

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


def sample_names(samplesheet):
    """Sample ids from the generate samplesheet.yaml (its `samples:` list)."""
    if not samplesheet.exists():
        sys.exit(f"{samplesheet} not found - run generate_sweep.py first.")
    doc = yaml.safe_load(samplesheet.read_text())
    return [s["sample"] for s in doc["samples"]]


def main():
    results_dir = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else DEFAULT_RESULTS_DIR
    cfg_path = sys.argv[2] if len(sys.argv) > 2 else HERE / "config.yaml"
    cfg = sc.load_config(cfg_path)

    samples = sample_names(HERE / "samplesheet.yaml")
    profilers = cfg["database"]["profilers"]
    db_name = cfg["database"]["name"]

    rows = [
        {"sample": sample, "profiler": profiler,
         "benchmark_dir": str(results_dir / sample), "database": db_name}
        for sample in samples
        for profiler in profilers
    ]
    doc = {"databases": sc.database_block(cfg), "samples": rows}
    with open(HERE / "profile_samplesheet.yaml", "w") as fh:
        yaml.safe_dump(doc, fh, sort_keys=False, default_flow_style=False)

    print(f"Wrote profile_samplesheet.yaml: {len(samples)} samples x {profilers} "
          f"against '{db_name}', benchmark_dir root {results_dir}")


if __name__ == "__main__":
    main()
