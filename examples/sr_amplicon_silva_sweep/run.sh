#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
OUTDIR="$REPO/results/sr_amplicon_silva_sweep"

# --- Phase 1: draw the communities, generate the reads, map them once -------------
# `--step all` trains the error model(s), draws the negative-binomial communities,
# generates their V4 amplicon reads ONCE PER `error_models:` ARM at every subsample
# depth, and runs superresolution-amplicon once per depth against the pre-built reference
# set named in config.yaml's `database.path:`. The error-model arms are the one axis that
# cannot share work: different reads means a separate generation and a separate SILVA
# mapping pass each. Only `custom_database` is built; the pipeline resolves
# <path>/<name>_ssu.sr_refs.{fasta,tax} and hands them to the nested runs. That run publishes
# each depth's mapseq classification to <benchmark_dir>/profiling/sr/<id>.obs.mseq.gz -
# the expensive part, and the reason phase 2 is affordable.
python "$HERE/generate_samplesheet.py"

nextflow run "$REPO/main.nf" \
    -profile docker \
    -c "$HERE/benchmark.config" \
    --step all \
    --input "$HERE/samplesheet.yaml" \
    --outdir "$OUTDIR" \
    --seed 42

# --- Phase 2: the parameter sweep -------------------------------------------------
# `--step profile` re-profiles those same reads four ways: the `silva` grid sweep, the
# `custom` arm against the panel-built collection, the `silva_panel` reinterpretation, and
# `aap` (no superresolution at all). The `silva` and `silva_panel` rows hand back the
# phase-1 mapseq output as `mseq:`, so no grid point re-maps anything against SILVA; the
# run's real cost is one mis-mapping matrix per distinct combination of the matrix knobs,
# plus one (cheap) inference per grid point, plus the `aap` run's own mapping.
python "$HERE/generate_sweep_samplesheet.py" "$OUTDIR"

nextflow run "$REPO/main.nf" \
    -profile docker \
    -c "$HERE/benchmark.config" \
    --step profile \
    --input "$HERE/sweep_samplesheet.yaml" \
    --outdir "$OUTDIR" \
    --seed 42

# --- Score ---------------------------------------------------------------------------
# Per genus - the one space all four arms share. One row per benchmark dir x method,
# with `arm` and `error_model` columns.
python "$HERE/scripts/score_sweep.py" "$OUTDIR" > "$OUTDIR/silva_sweep_scores.csv"
echo "Scores: $OUTDIR/silva_sweep_scores.csv"
