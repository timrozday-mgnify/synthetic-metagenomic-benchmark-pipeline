#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
OUTDIR="$REPO/results/sr_amplicon_aap"
cd "$REPO"   # the samplesheet's results/ paths are relative to the pipeline dir

run() {
    nextflow run "$REPO/main.nf" -profile docker -c "$HERE/benchmark.config" \
        --input "$HERE/samplesheet.yaml" --outdir "$OUTDIR" -resume "$@"
}

# 1. Train the error model once. The samplesheet points the next step at it twice: to
#    generate the reads, and as the pairs arm's per-mate model.
run --step train
# 2. Generate the sample, run AAP on it, then all four SR arms.
run --step all

python "$HERE/score.py" "$OUTDIR" > "$OUTDIR/scores.csv"
cat "$OUTDIR/scores.csv"
