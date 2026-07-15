# synthetic-metagenomic-benchmark-pipeline

A [Nextflow](https://www.nextflow.io/) DSL2 pipeline that builds **synthetic
metagenomes with a realistic, data-derived sequencing error profile**, for
benchmarking metagenomic profilers.

For each sample it:

1. **Trains a sequencing error profile** from a *natural* (real, non-synthetic)
   metagenome using [skiver](https://github.com/timrozday-mgnify/skiver) —
   reference-free, from k-mer consensus.
2. **Generates synthetic reads** from a chosen set of reference genomes at chosen
   relative abundances and sequencing depth, applying that error profile via
   [genome-blender](https://github.com/timrozday-mgnify/genome-blender).
3. **Publishes the benchmark ground truth**: the reads, a BAM giving the true
   read→reference mapping, and a ground-truth abundance profile.

> The downstream profiling / evaluation step (e.g. sylph, an amplicon pipeline)
> is intentionally **out of scope** — this pipeline produces the inputs and the
> ground truth a profiler is scored against.

## Quick start

```bash
nextflow run main.nf \
    -profile docker \
    --input samplesheet.csv \
    --outdir results
```

Requires Nextflow (>=25) and a container engine (Docker / Singularity / Apptainer).

## Samplesheet

One row per synthetic sample (`tests/samplesheets/test.csv` is a working example):

```csv
sample,train_id,train_fastq_1,train_fastq_2,platform,genomes_csv,num_reads,mode
S1,natA,natural_R1.fastq.gz,natural_R2.fastq.gz,hq-illumina,genomes_S1.csv,1000000,paired
S2,natA,natural_R1.fastq.gz,natural_R2.fastq.gz,hq-illumina,genomes_S2.csv,1000000,paired
```

| Column | Description |
|--------|-------------|
| `sample` | Unique synthetic-sample ID (names the output dir). |
| `train_id` | Groups error-model training. Rows sharing a `train_id` train the profile **once** and reuse it. |
| `train_fastq_1` | Natural-metagenome reads to learn the error profile from. |
| `train_fastq_2` | Optional mate (paired training reads). Leave blank for single-end. |
| `platform` | `hq-illumina` \| `lq-illumina` \| `ont` \| `pacbio`. |
| `genomes_csv` | A genome-blender input CSV: `genome_id,fasta_path,abundance`. |
| `num_reads` | Sequencing depth as total reads (read pairs × 2 for paired mode). |
| `mode` | `paired` \| `single` \| `long` \| `amplicon`. |

Relative `genomes_csv` / FASTA / FASTQ paths resolve against the pipeline
directory; absolute paths and `scheme://` URLs pass through.

### The genomes CSV

Referenced by each row's `genomes_csv` column — the reference genomes and their
relative abundances (normalised internally):

```csv
genome_id,fasta_path,abundance
genomeA,tests/data/genomeA.fasta,0.7
genomeB,tests/data/genomeB.fasta,0.3
```

## Outputs

Published under `results/<sample>/`:

| File | Description |
|------|-------------|
| `<sample>_R1.fastq.gz`, `<sample>_R2.fastq.gz` | Synthetic reads (single file `<sample>.fastq.gz` for single/long/amplicon). |
| `<sample>.sorted.bam` (+ `.bai`) | Ground-truth read→reference mapping. References are named `<genome_id>:<contig_id>`. |
| `<sample>.truth.tsv` | Ground-truth profile: `genome_id, target_rel_abundance, realized_n_reads, realized_rel_abundance`. |

And under `results/error_models/<train_id>/`: the trained `<train_id>.model.pt`,
`<train_id>.phred_calibration.json`, and `<train_id>.context_model_aic.csv`.

`target_rel_abundance` is the normalised requested abundance; `realized_*` come
from `samtools idxstats` on the true BAM (what was actually generated).

## Parameters

| Param | Default | Description |
|-------|---------|-------------|
| `--input` | – | Samplesheet CSV (required). |
| `--outdir` | `./results` | Output directory. |
| `--seed` | `42` | Global RNG seed (training + generation). |
| `--error_model_candidates` | `AdditiveContext(5),AdditiveContext(7),AdditiveContext(9)` | Candidate skiver component strings; the min-AIC model is kept. |
| `--error_model_components` | `null` | Force a single component string (skips the AIC search). |
| `--smb_skiver_tag` / `--smb_genome_blender_tag` | `latest` | Container image tags. |
| `--max_cpus` / `--max_memory` / `--max_time` | `16` / `128.GB` / `240.h` | Resource caps (override per HPC profile). |

## Containers

Two images are built from the vendored submodules (`vendor/skiver`,
`vendor/genome-blender`) — see `containers/`:

- **`smb-skiver`** — the `skiver` Rust binary (`skiver dump`) plus the Python
  training/calibration environment.
- **`smb-genome-blender`** — `generate-reads` plus `skiver-generate` on `PATH`
  (genome-blender shells out to it), and `pysam` (used by the ground-truth step).

Build locally:

```bash
git submodule update --init --recursive
docker build -f containers/skiver/Dockerfile         -t ghcr.io/timrozday-mgnify/smb-skiver:latest .
docker build -f containers/genome-blender/Dockerfile -t ghcr.io/timrozday-mgnify/smb-genome-blender:latest .
```

The `.github/workflows/build-images.yml` workflow builds and pushes them to GHCR.

### HPC / Singularity

Run with `-profile singularity` (or `apptainer`). Singularity converts the
Docker images automatically. Override resources with a site config, e.g.:

```bash
nextflow run main.nf -profile singularity -c my_hpc.config --input samplesheet.csv
```

## Testing

Three tiers, via [nf-test](https://www.nf-test.com/):

```bash
nf-test test modules/ tests/default.nf.test --tag stub   # fast, host-only, no images
nf-test test tests/default.nf.test --profile docker --tag e2e   # full run (needs images)
pytest tests/bin/test_bin.py                              # bin/ helper unit tests
```

CI (`.github/workflows/ci.yml`) runs pre-commit, the `bin/` unit tests, and the
stub nf-tests on every PR. Image builds and the full e2e run on `main`.

## Development

Set up pre-commit:

```bash
pip install pre-commit && pre-commit install
```

Test fixtures are regenerated with `python tests/data/generate_fixtures.py`.
