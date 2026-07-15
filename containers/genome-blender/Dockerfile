# smb-genome-blender: genome-blender (`generate-reads`) + skiver-generate on PATH.
# genome-blender shells out to `skiver-generate`, so both must live in one image.
# Also used by GROUND_TRUTH (ships pysam + python).
# Build context is the repo root:
#   docker build -f containers/genome-blender/Dockerfile -t smb-genome-blender .
FROM python:3.11-slim-bookworm

RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4 zlib1g libbz2-1.0 liblzma5 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

# CPU-only torch first so the heavy deps of skiver-generate / genome-blender
# (pyro, etc.) resolve against it instead of pulling the GPU build.
RUN pip install --no-cache-dir torch --index-url https://download.pytorch.org/whl/cpu

# skiver-generate (the generative error-model CLI genome-blender invokes).
COPY vendor/skiver/scripts/ /opt/skiver/scripts/
RUN pip install --no-cache-dir /opt/skiver/scripts

# genome-blender itself (pure-python; provides `generate-reads`).
# cli.py imports `click` directly, so pin typer + click explicitly.
COPY vendor/genome-blender/ /opt/genome-blender/
RUN pip install --no-cache-dir "typer>=0.9" click \
    && pip install --no-cache-dir /opt/genome-blender

RUN generate-reads --help >/dev/null && skiver-generate --help >/dev/null \
    && python -c "import pysam; print('pysam', pysam.__version__)"
