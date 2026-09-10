#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
OUTDIR="$REPO/results/sr_amplicon_gtdb_sweep"

# --- Phase 1: draw the communities, generate the reads, map them once -------------
# `--step all` trains the error model, draws the negative-binomial communities,
# generates their V4 amplicon reads at every subsample depth, and runs
# superresolution-amplicon once per depth against the pre-built reference set named in
# config.yaml's `database.path:`. Nothing is built: the pipeline resolves
# <path>/<name>_ssu.sr_refs.fasta and hands it to the nested runs. That run publishes
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
# `--step profile` re-profiles those same reads once per grid point. Each row hands back
# the phase-1 mapseq output as `mseq:`, so no grid point re-maps anything against GTDB;
# the run's real cost is one mis-mapping matrix per distinct combination of the matrix
# knobs, plus one (cheap) inference per grid point.
python "$HERE/generate_sweep_samplesheet.py" "$OUTDIR"

nextflow run "$REPO/main.nf" \
    -profile docker \
    -c "$HERE/benchmark.config" \
    --step profile \
    --input "$HERE/sweep_samplesheet.yaml" \
    --outdir "$OUTDIR" \
    --seed 42

# Each grid point's profile lands next to the truth it is scored against:
#   $OUTDIR/S01.amplicon_16s.515-YF-806BR/S01.amplicon_16s.515-YF-806BR.<point>.sr_profile.tsv
# and each matrix under $OUTDIR/mismapping/<reference set>/.
