#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
OUTDIR="$REPO/results/subspecies_v4_sweep"

# Regenerate samplesheet.yaml + genomes/*.csv from config.yaml.
python "$HERE/generate_sweep.py"

# One combined run: --step all trains the error model once (deduped by train_id),
# generates reads for every sample, builds the community_v4 profiler DB from the
# samplesheet's `databases:` block, and profiles each sample against it.
nextflow run "$REPO/main.nf" \
    -profile docker \
    -c "$HERE/benchmark.config" \
    --step all \
    --input "$HERE/samplesheet.yaml" \
    --outdir "$OUTDIR" \
    --seed 42

# To re-profile already-generated reads without regenerating (e.g. to compare
# profilers), regenerate a profile-only samplesheet and run --step profile:
#   python "$HERE/generate_profile_samplesheet.py" "$OUTDIR"
#   nextflow run "$REPO/main.nf" -profile docker -c "$HERE/benchmark.config" \
#       --step profile --input "$HERE/profile_samplesheet.yaml" --outdir "$OUTDIR" --seed 42
