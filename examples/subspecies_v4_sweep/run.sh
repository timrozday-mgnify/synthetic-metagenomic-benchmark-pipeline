#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"

# Regenerate samplesheet.yaml + genomes/*.csv from the PANEL in generate_sweep.py.
python "$HERE/generate_sweep.py"

nextflow run "$REPO/main.nf" \
    -profile docker \
    -c "$HERE/benchmark.config" \
    --input "$HERE/samplesheet.yaml" \
    --outdir "$REPO/results/subspecies_v4_sweep" \
    --seed 42
