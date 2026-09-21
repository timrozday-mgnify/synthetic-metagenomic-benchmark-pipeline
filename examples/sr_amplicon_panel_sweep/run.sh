#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
OUTDIR="$REPO/results/sr_amplicon_panel_sweep"

# --- Phase 1: draw the communities, generate the reads, map them against SILVA once ----
# `--step all` trains the error model, draws the negative-binomial communities, generates
# their V4 amplicon reads at every subsample depth, and runs superresolution-amplicon once
# per depth against SILVA under config.yaml's `generic.map_setting` (align tau 0 over the panel: cheap
# at SILVA scale). That publishes each depth's SILVA classification to
# <benchmark_dir>/profiling/sr/<id>.map.obs.mseq.gz, and the trained model to
# $OUTDIR/error_models/<train_id>/.
python "$HERE/generate_samplesheet.py"

nextflow run "$REPO/main.nf" \
    -profile docker \
    -c "$HERE/benchmark.config" \
    --step all \
    --input "$HERE/samplesheet.yaml" \
    --outdir "$OUTDIR" \
    --seed 42

# --- Phase 2: all three arms, every grid point -----------------------------------------
# `--step profile` re-profiles those reads twice per benchmark dir: against the custom
# 20HM collection (its own grid), and against SILVA reusing phase 1's mapping, with every
# setting reinterpreting the labels through the 20HM panel or its species panel (the other
# two grids).
python "$HERE/generate_sweep_samplesheet.py" "$OUTDIR"

nextflow run "$REPO/main.nf" \
    -profile docker \
    -c "$HERE/benchmark.config" \
    --step profile \
    --input "$HERE/sweep_samplesheet.yaml" \
    --outdir "$OUTDIR" \
    --seed 42

# --- Score -----------------------------------------------------------------------------
# One row per benchmark dir x method: TV, strain split, false mass, background share.
python "$HERE/scripts/score_sweep.py" "$OUTDIR" > "$OUTDIR/panel_sweep_scores.csv"
echo "Scores: $OUTDIR/panel_sweep_scores.csv"
