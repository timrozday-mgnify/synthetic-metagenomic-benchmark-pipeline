#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
OUTDIR="$REPO/results/sr_amplicon_param_sweep"

# --- Phase 1: generate the reads and map them once -------------------------------
# `--step all` trains the error model, generates the amplicon reads at every subsample
# depth, builds the community_v4 reference set, and runs superresolution-amplicon once
# per depth at the nested pipeline's default settings. That run publishes each depth's
# mapseq classification to <benchmark_dir>/profiling/sr/<id>.obs.mseq.gz.
python "$HERE/generate_samplesheet.py"

nextflow run "$REPO/main.nf" \
    -profile docker \
    -c "$HERE/benchmark.config" \
    --step all \
    --input "$HERE/samplesheet.yaml" \
    --outdir "$OUTDIR" \
    --seed 42

# --- Phase 2: the parameter sweep -------------------------------------------------
# `--step profile` re-profiles those same reads once per grid point. Each row hands back
# the phase-1 mapseq output as `mseq:`, so no grid point re-maps anything; the run's only
# real cost is one mis-mapping matrix per distinct combination of the matrix knobs, plus
# one (cheap) inference per grid point.
python "$HERE/generate_sweep_samplesheet.py" "$OUTDIR"

nextflow run "$REPO/main.nf" \
    -profile docker \
    -c "$HERE/benchmark.config" \
    --step profile \
    --input "$HERE/sweep_samplesheet.yaml" \
    --outdir "$OUTDIR" \
    --seed 42

# Each grid point's profile lands next to the truth it is scored against:
#   $OUTDIR/community.515-YF-806BR/community.515-YF-806BR.<grid point>.sr_profile.tsv
# and each matrix under $OUTDIR/mismapping/<reference set>/.
