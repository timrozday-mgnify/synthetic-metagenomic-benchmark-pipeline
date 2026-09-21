# smb-genome-blender: genome-blender (`generate-reads`) + skiver-generate on PATH.
# genome-blender shells out to `skiver-generate`, so both must live in one image.
# Also used by GROUND_TRUTH (ships pysam + python).
# Build context is the repo root:
#   docker build -f containers/genome-blender/Dockerfile -t smb-genome-blender .
FROM python:3.11-slim-bookworm

# procps provides `ps`, required by Nextflow to collect task metrics.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4 zlib1g libbz2-1.0 liblzma5 procps \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

# CPU-only torch first so the heavy deps of skiver-generate / genome-blender
# (pyro, etc.) resolve against it instead of pulling the GPU build.
RUN pip install --no-cache-dir torch --index-url https://download.pytorch.org/whl/cpu

# skiver-generate (the generative error-model CLI genome-blender invokes).
COPY vendor/skiver/scripts/ /opt/skiver/scripts/
RUN pip install --no-cache-dir /opt/skiver/scripts
# The platform presets (`--error-model illumina|pacbio|nanopore`, a samplesheet
# `error_model:`). lib/presets.py looks two levels above its own file, which once
# pip-installed is the interpreter's lib dir, not the skiver checkout.
# ponytail: tied to python3.11's layout; ship the presets as package data if the base
# image's python moves.
COPY vendor/skiver/context_error_models/ /usr/local/lib/python3.11/context_error_models/

# genome-blender itself (pure-python; provides `generate-reads`).
# cli.py imports `click` directly, so pin typer + click explicitly.
COPY vendor/genome-blender/ /opt/genome-blender/
RUN pip install --no-cache-dir "typer>=0.9" click \
    && pip install --no-cache-dir /opt/genome-blender

RUN generate-reads --help >/dev/null && skiver-generate --help >/dev/null \
    && normalise-fasta --help >/dev/null \
    && python -c "from lib.presets import available_presets as a, resolve_preset as r; [r(p) for p in a()]" \
    && python -c "import pysam; print('pysam', pysam.__version__)"
