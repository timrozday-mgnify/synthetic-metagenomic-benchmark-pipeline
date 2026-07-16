#!/usr/bin/env python3
"""Generate a profile-only samplesheet from this example's samplesheet.yaml,
for re-profiling reads already written by an earlier `--step generate` run.

Fill in PROFILER/DATABASE below, then run:
    python generate_profile_samplesheet.py [results_dir]
Writes profile_samplesheet.yaml, then:
    nextflow run ../../main.nf -profile docker --step profile \\
        --input profile_samplesheet.yaml --outdir <results_dir>
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_RESULTS_DIR = HERE.parent.parent / "results" / "subspecies_v4_sweep"

# --- Fill these in ---------------------------------------------------------
PROFILER = "sylph"       # sylph | aap
DATABASE = "gtdb_r220"   # configured params.sylph_databases key ('self' not available here)
# ---------------------------------------------------------------------------


def main():
    results_dir = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else DEFAULT_RESULTS_DIR
    samples = [
        line.split(":", 1)[1].strip()
        for line in (HERE / "samplesheet.yaml").read_text().splitlines()
        if line.startswith("- sample:")
    ]

    with open(HERE / "profile_samplesheet.yaml", "w") as fh:
        for sample in samples:
            fh.write(f"- sample: {sample}\n")
            fh.write(f"  profiler: {PROFILER}\n")
            fh.write(f"  benchmark_dir: {results_dir / sample}\n")
            fh.write(f"  database: {DATABASE}\n\n")

    print(f"Wrote profile_samplesheet.yaml ({len(samples)} samples), benchmark_dir root {results_dir}")


if __name__ == "__main__":
    main()
