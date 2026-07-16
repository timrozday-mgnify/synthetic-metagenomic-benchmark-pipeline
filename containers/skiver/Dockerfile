# smb-skiver: skiver Rust binary (`skiver dump`) + Python training/calibration env.
# Build context is the repo root so the vendor/skiver submodule is visible:
#   docker build -f containers/skiver/Dockerfile -t smb-skiver .

# ---- stage 1: build the skiver Rust binary --------------------------------
FROM rust:1-bookworm AS builder

# rust-htslib links the system htslib and uses libclang for bindgen.
# libhts-dev pulls the curl dev package it needs (gnutls variant) itself, so we
# don't pin libcurl*-dev here (avoids the gnutls/openssl conflict).
RUN apt-get update && apt-get install -y --no-install-recommends \
        clang libclang-dev pkg-config cmake make \
        libhts-dev zlib1g-dev libbz2-dev liblzma-dev libssl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY vendor/skiver/ ./
RUN cargo build --release --bin skiver && cp target/release/skiver /usr/local/bin/skiver

# ---- stage 2: runtime -----------------------------------------------------
FROM python:3.11-slim-bookworm

ARG QUARTO_VERSION=1.5.57

# htslib shared library for the dynamically linked binary (apt resolves the
# matching libcurl runtime dep); procps provides `ps`, required by Nextflow;
# curl fetches the Quarto .deb (report rendering).
RUN apt-get update && apt-get install -y --no-install-recommends \
        libhts3 procps curl \
    && curl -fsSL -o /tmp/quarto.deb \
        "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb" \
    && dpkg -i /tmp/quarto.deb \
    && rm -f /tmp/quarto.deb \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/bin/skiver /usr/local/bin/skiver

# Training / calibration scripts (run via `python $SKIVER_SCRIPTS/<script>.py`).
ENV SKIVER_SCRIPTS=/opt/skiver/scripts
COPY vendor/skiver/scripts/ ${SKIVER_SCRIPTS}/

# Generic per-platform Quarto error-model report, rendered by SKIVER_REPORT.
ENV SKIVER_NOTEBOOKS=/opt/skiver/notebooks
COPY vendor/skiver/notebooks/context_error_model_*.qmd vendor/skiver/notebooks/_context_error_model_body.qmd ${SKIVER_NOTEBOOKS}/

# CPU-only torch keeps the image lean; pyro + the rest are the training deps
# declared in skiver's requirements.txt (which is gitignored, so pinned here).
# plotly/ipykernel/papermill are needed to render the Quarto error-model report
# (papermill executes the -P/--execute-param substitutions Quarto passes through).
RUN pip install --no-cache-dir torch --index-url https://download.pytorch.org/whl/cpu \
    && pip install --no-cache-dir \
        numpy scipy pandas h5py pyyaml matplotlib rich "pyro-ppl>=1.9" \
        plotly ipykernel jupyter papermill \
    && python -m ipykernel install --name skiver --sys-prefix

# Nextflow's docker executor runs the container as the host UID, not root, so
# $HOME (unset for that UID) must point somewhere writable for Quarto/Deno's
# cache dir (~/.cache/quarto) and Jupyter's runtime dir.
ENV HOME=/tmp
RUN chmod 1777 /tmp

RUN skiver --version && quarto --version

# The version checks above ran as root and left root-owned junk under /tmp
# (quarto's cache dir, jupyter kernel scratch dirs) that a later non-root
# runtime UID can't write into/around; clear it so /tmp starts empty.
RUN rm -rf /tmp/.cache /tmp/tmp*
