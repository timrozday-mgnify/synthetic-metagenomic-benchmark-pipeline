#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
OUTDIR="$REPO/results/subspecies_v4_sweep"

run() {  # $1 = step, $2 = input samplesheet
    nextflow run "$REPO/main.nf" \
        -profile docker \
        -c "$HERE/benchmark.config" \
        --step "$1" \
        --input "$2" \
        --outdir "$OUTDIR" \
        --seed 42
}

# Regenerate train_samplesheet.yaml + samplesheet.yaml + genomes/*.csv from the PANEL.
python "$HERE/generate_sweep.py"

# 1. Train the error model once -> $OUTDIR/error_models/<train_id>/
run train "$HERE/train_samplesheet.yaml"

# 2. Generate reads for every sample, reusing that model (rows carry error_model_dir).
run generate "$HERE/samplesheet.yaml"

# 3. Profile the generated samples (fill in PROFILER/DATABASE in the script first).
python "$HERE/generate_profile_samplesheet.py" "$OUTDIR"
run profile "$HERE/profile_samplesheet.yaml"
