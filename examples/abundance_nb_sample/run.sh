#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
OUTDIR="$REPO/results/abundance_nb_sample"

# Regenerate samplesheet.yaml + genomes/*.csv from config.yaml (NB abundance draws).
python "$HERE/generate_samplesheet.py"

# One combined run: --step all trains the error model once (deduped by train_id),
# generates reads for every sample, builds the profiler DBs from the samplesheet's
# `databases:` block, and profiles each sample with EVERY profiler its row lists -
# the mode's `profiler:` plus its `extra_profilers:` (sylph + superresolution-shotgun
# for wgs, aap + superresolution-amplicon for amplicon_16s). Each method's profile
# lands in the same benchmark dir, so they compare directly against truth.tsv.
nextflow run "$REPO/main.nf" \
    -profile docker \
    -c "$HERE/benchmark.config" \
    --step all \
    --input "$HERE/samplesheet.yaml" \
    --outdir "$OUTDIR" \
    --seed 42

# To re-profile already-generated reads without regenerating them (e.g. after adding
# a profiler to config.yaml), build a profile-only samplesheet and run --step profile:
#   python "$HERE/generate_profile_samplesheet.py" "$OUTDIR"
#   nextflow run "$REPO/main.nf" -profile docker -resume -c "$HERE/benchmark.config" \
#       --step profile --input "$HERE/profile_samplesheet.yaml" --outdir "$OUTDIR" --seed 42
