# AmpliconHunter — in-silico PCR tool for large microbial genome collections
#
# Source: https://github.com/rhowardstone/AmpliconHunter
# Licence: MIT
#
# Dependencies (from setup.py):
#   Python  >= 3.7
#   hyperscan >= 0.7.0   — Python bindings for Intel Hyperscan SIMD regex engine
#   biopython >= 1.78
#   numpy >= 1.19.0
#   pandas >= 1.0.0
#   matplotlib >= 3.3.0
#   seaborn >= 0.11.0
#
# Optional (for --hmm flag):
#   hmmer   (nhmmer + hmmbuild)
#   mafft
#
# Architecture note: Intel Hyperscan uses SSE/AVX SIMD and is x86_64 only.
# On ARM hosts use `docker buildx build --platform linux/amd64 ...` and run
# under emulation, or replace hyperscan with vectorscan.

FROM python:3.11-slim

LABEL org.opencontainers.image.title="AmpliconHunter" \
      org.opencontainers.image.description="High-performance in-silico PCR for large microbial genome collections" \
      org.opencontainers.image.source="https://github.com/rhowardstone/AmpliconHunter" \
      org.opencontainers.image.licenses="MIT"

# ── Environment ────────────────────────────────────────────────────────────────
ENV DEBIAN_FRONTEND=noninteractive \
    # Non-interactive matplotlib backend — no display in a container
    MPLBACKEND=Agg \
    # Writable matplotlib cache for any UID (Nextflow: -u $(id -u):$(id -g))
    MPLCONFIGDIR=/tmp \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# ── System packages ────────────────────────────────────────────────────────────
# hmmer + mafft — optional, only needed for --hmm mode.
# No system hyperscan library is required: the 'hyperscan' PyPI package ships
# a manylinux wheel that bundles Hyperscan/VectorScan statically.  Forcing a
# source build (--no-binary) would require autoreconf + boost headers to
# compile VectorScan from source, adding significant complexity for no gain.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        hmmer \
        mafft \
        git \
        # procps provides 'ps', required by Nextflow for task metrics collection
        procps \
    && rm -rf /var/lib/apt/lists/*

# ── Python dependencies ────────────────────────────────────────────────────────
# Install in a single invocation so pip can resolve the full dependency graph.
# hyperscan is installed last so the heavy packages above are cached as a
# separate layer and a hyperscan version bump doesn't invalidate them.
RUN pip install --no-cache-dir \
        "biopython>=1.78" \
        "numpy>=1.19.0" \
        "pandas>=1.0.0" \
        "matplotlib>=3.3.0" \
        "seaborn>=0.11.0" \
        "hyperscan>=0.7.0"

# ── AmpliconHunter ─────────────────────────────────────────────────────────────
# Install from the GitHub main branch.  Pin to a specific commit SHA for
# fully reproducible builds once a release tag is published:
#   pip install git+https://github.com/rhowardstone/AmpliconHunter.git@<SHA>
RUN pip install --no-cache-dir \
        "git+https://github.com/rhowardstone/AmpliconHunter.git"

# ── Smoke test ─────────────────────────────────────────────────────────────────
RUN ampliconhunter --help

CMD ["ampliconhunter", "--help"]
