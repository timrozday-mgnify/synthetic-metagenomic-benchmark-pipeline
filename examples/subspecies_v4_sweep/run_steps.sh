#!/usr/bin/env bash
# Like run.sh, but drives the pipeline through its three phases separately
# (train -> generate -> profile), emitting an intermediate samplesheet per phase.
# Lets each phase be re-run in isolation; all phases share one $OUTDIR so they
# chain through error_models/ and the per-sample benchmark dirs.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
OUTDIR="$REPO/results/subspecies_v4_sweep"

nf() { nextflow run "$REPO/main.nf" -profile docker -c "$HERE/benchmark.config" \
    --outdir "$OUTDIR" --seed 42 "$@"; }

# Base samplesheet.yaml + genomes/*.csv from config.yaml (same as run.sh).
python "$HERE/generate_sweep.py"

# 1. Train the error model once (deduped by train_id). Publishes the model +
#    calibration to $OUTDIR/error_models/<train_id>/.
nf --step train --input "$HERE/samplesheet.yaml"

# 2. Generate samplesheet: samplesheet.yaml + an error_model_dir per sample row
#    (pointing at the just-trained model), so --step generate reuses it instead
#    of retraining. ponytail: ~6-line YAML transform, not worth a 4th generate_*.py.
OUTDIR="$OUTDIR" python - "$HERE/samplesheet.yaml" "$HERE/generate_samplesheet.yaml" <<'PY'
import os, sys, yaml
src, dst = sys.argv[1], sys.argv[2]
outdir = os.environ["OUTDIR"]
doc = yaml.safe_load(open(src).read())
for row in doc["samples"]:
    row["error_model_dir"] = f"{outdir}/error_models/{row['train_id']}"
with open(dst, "w") as fh:
    yaml.safe_dump(doc, fh, sort_keys=False, default_flow_style=False)
print(f"Wrote {dst}: {len(doc['samples'])} samples with error_model_dir")
PY

# 3. Generate reads + ground truth (training skipped, model reused). No DB built.
nf --step generate --input "$HERE/generate_samplesheet.yaml"

# 4. Profile samplesheet: one row per (sample, profiler) + the databases: block.
python "$HERE/generate_profile_samplesheet.py" "$OUTDIR"

# 5. Build the community_v4 DB and profile each sample against it.
nf --step profile --input "$HERE/profile_samplesheet.yaml"
