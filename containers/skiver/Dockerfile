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

# htslib shared library for the dynamically linked binary (apt resolves the
# matching libcurl runtime dep); procps provides `ps`, required by Nextflow.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libhts3 procps \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/bin/skiver /usr/local/bin/skiver

# Training / calibration scripts (run via `python $SKIVER_SCRIPTS/<script>.py`).
ENV SKIVER_SCRIPTS=/opt/skiver/scripts
COPY vendor/skiver/scripts/ ${SKIVER_SCRIPTS}/

# CPU-only torch keeps the image lean; pyro + the rest are the training deps
# declared in skiver's requirements.txt (which is gitignored, so pinned here).
RUN pip install --no-cache-dir torch --index-url https://download.pytorch.org/whl/cpu \
    && pip install --no-cache-dir \
        numpy scipy pandas h5py pyyaml matplotlib "pyro-ppl>=1.9"

RUN skiver --version
