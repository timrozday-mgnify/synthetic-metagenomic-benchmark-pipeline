#!/usr/bin/env bash
# Pre-fetch the pipeline's container images into the Nextflow Singularity/Apptainer
# cache, once, before launching the pipeline. This avoids the scattered
# GENOME_BLENDER_GENERATE chunk tasks racing to pull the same image into a shared
# cache (which corrupts it -> "unexpected end of JSON input").
#
# Run AFTER the CI rebuild that disables provenance attestations has landed; pulling
# a provenance-decorated :latest fails the same way regardless of the cache.
#
# Usage:
#   export NXF_SINGULARITY_CACHEDIR=/scratch/$USER/nxf_singularity_cache
#   export APPTAINER_TMPDIR=/scratch/$USER/apptainer_tmp   # keep SIF builds off $HOME/tmp
#   ./bin/prefetch-singularity.sh
#
# Env overrides: SMB_SKIVER_TAG, SMB_GENOME_BLENDER_TAG (default: latest).
set -euo pipefail

cache="${NXF_SINGULARITY_CACHEDIR:-${NXF_APPTAINER_CACHEDIR:-$PWD/singularity_cache}}"
mkdir -p "$cache"

# Keep OCI->SIF scratch off a small/quota'd $HOME or /tmp if the caller set a target.
if [[ -n "${APPTAINER_TMPDIR:-}" ]]; then mkdir -p "$APPTAINER_TMPDIR"; fi
if [[ -n "${SINGULARITY_TMPDIR:-}" ]]; then mkdir -p "$SINGULARITY_TMPDIR"; fi

if command -v apptainer >/dev/null 2>&1; then runtime=apptainer
elif command -v singularity >/dev/null 2>&1; then runtime=singularity
else echo "error: neither apptainer nor singularity found on PATH" >&2; exit 1
fi

# ponytail: image list hardcoded because the pipeline uses exactly these two ghcr
# images. If more are added, switch to `nextflow inspect <entry> -profile singularity`.
images=(
  "ghcr.io/timrozday-mgnify/smb-skiver:${SMB_SKIVER_TAG:-latest}"
  "ghcr.io/timrozday-mgnify/smb-genome-blender:${SMB_GENOME_BLENDER_TAG:-latest}"
)

# Nextflow's cache filename: strip protocol, replace ':' and '/' with '-', add .img.
# (Verify against your Nextflow version if it stops reusing the cache — the naming
# scheme is the one fragile assumption here.)
nxf_cache_name() { echo "${1//[:\/]/-}.img"; }

for img in "${images[@]}"; do
  dest="$cache/$(nxf_cache_name "$img")"
  if [[ -f "$dest" ]]; then
    echo "already cached: $dest"
    continue
  fi
  echo "pulling $img -> $dest"
  "$runtime" pull "$dest" "docker://$img"
done

echo "done. point Nextflow at: NXF_SINGULARITY_CACHEDIR=$cache"
