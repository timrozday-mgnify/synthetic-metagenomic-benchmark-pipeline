#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
OUTDIR="$REPO/results/abundance_nb_sample"

# Regenerate samplesheet.yaml + genomes/*.csv from config.yaml (NB abundance draws).
python "$HERE/generate_samplesheet.py"

# One combined run: --step all trains the error model once (deduped by train_id),
# generates reads for every sample, builds the profiler DB from the samplesheet's
# `databases:` block, and profiles each sample against it.
nextflow run "$REPO/main.nf" \
    -profile docker \
    -c "$HERE/benchmark.config" \
    --step all \
    --input "$HERE/samplesheet.yaml" \
    --outdir "$OUTDIR" \
    --seed 42

# Second pass: benchmark the remaining profilers against the SAME reads. The
# profile samplesheet has one row per (sample, profiler) - each mode's primary
# `profiler:` plus its `extra_profilers:` - so the superresolution methods land
# their profiles next to the sylph/aap ones without regenerating anything.
# (The primary profiler re-runs here too; it's cached by -resume, and dropping it
# would mean maintaining two sample lists.)
python "$HERE/generate_profile_samplesheet.py" "$OUTDIR"

nextflow run "$REPO/main.nf" \
    -profile docker \
    -resume \
    -c "$HERE/benchmark.config" \
    --step profile \
    --input "$HERE/profile_samplesheet.yaml" \
    --outdir "$OUTDIR" \
    --seed 42
