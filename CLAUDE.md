# CLAUDE.md

Guidance for Claude Code / developers working in this repo.

## What this is

A Nextflow DSL2 pipeline that trains a skiver sequencing-error profile from a
natural metagenome and applies it (via genome-blender) to reference genomes to
emit synthetic reads + a ground-truth BAM + a ground-truth abundance profile.
The profiling/evaluation stage is deliberately not included.

Built to the conventions of `~/Documents/subspecies-phylogeny` (biocontainer/
GHCR container-per-module, `conf/modules.config` publish layer, three-tier
nf-test, pre-commit + GitHub Actions).

## DAG

```
[samplesheet row] --(dedup per train_id)--> SKIVER_DUMP -> SKIVER_TRAIN
                                                              |  \_ phred_calibration.json
                                                              |  \_ model.pt (min-AIC)
[genomes_csv + fastas + num_reads] --(join by train_id)-----> GENOME_BLENDER_GENERATE
                                                              |  \_ reads.fastq.gz
                                                              |  \_ truth.bam
                                                              v
                                                          GROUND_TRUTH -> sorted bam+bai, truth.tsv
```

## Layout

- `main.nf` — inline samplesheet parse; builds `ch_samples` and a train_id-deduped
  `ch_train`, calls the workflow. Resolves the FASTAs referenced by each
  `genomes_csv` so Nextflow stages them.
- `workflows/synthetic_metagenomic_benchmark.nf` — top workflow; joins trained
  models back to samples by `train_id`.
- `subworkflows/local/train_error_model/` — SKIVER_DUMP → SKIVER_TRAIN.
- `modules/local/` — SKIVER_DUMP, SKIVER_TRAIN, GENOME_BLENDER_GENERATE, GROUND_TRUTH.
- `bin/` — helper scripts staged onto PATH: `build_model_config.py`,
  `pick_best_model.py`, `rewrite_genomes_csv.py`, `build_truth.py`.
- `containers/` — Dockerfiles for the two images (built from `vendor/` submodules).
- `conf/` — `base.config` (resource labels), `modules.config` (publish layer).

## Container ↔ module

| Module | Image | Tool |
|--------|-------|------|
| SKIVER_DUMP, SKIVER_TRAIN | `smb-skiver` | `skiver dump`, `train_context_error_models.py`, `fit_phred_calibration.py` |
| GENOME_BLENDER_GENERATE | `smb-genome-blender` | `generate-reads` (+ `skiver-generate` subprocess) |
| GROUND_TRUTH | `smb-genome-blender` | `pysam` (sort/index/idxstats) — reused, no separate samtools image |

Images are `ghcr.io/timrozday-mgnify/<image>:${params.<image>_tag}`.

## Key implementation notes

- **Error-model training is model-config driven.** The pinned skiver commit has
  no greedy AIC selection (that lives in unmerged skiver WIP), so `SKIVER_TRAIN`
  fits a candidate list (`params.error_model_candidates`) via a generated
  `model_config.json` and keeps the **min-AIC MLE** model (`bin/pick_best_model.py`).
  `params.error_model_components` forces a single model.
- **Training is reference-free** — the natural sample has no reference; `skiver
  dump --base` builds consensus from k-mer coverage. Fixtures therefore need
  *high coverage* (`tests/data/generate_fixtures.py` simulates ~400× over two
  small genomes).
- **Ground truth** is derived from the true BAM: `@SQ` refs are `genome_id:contig_id`,
  so `bin/build_truth.py` aggregates `idxstats` by the part before `:`.
- **Stub tests run on the host** (no `--profile`, so no container engine); stub
  blocks must use only coreutils (no tool calls). Real/e2e tests use
  `--profile docker` + `--tag e2e`; stub selection is `--tag stub`.

## Submodules

`vendor/skiver` and `vendor/genome-blender` are git submodules. The pipeline
pins the skiver commit that carries the `skiver-generate` CLI. After cloning:
`git submodule update --init --recursive`.

## Dev commands

```bash
nf-test test modules/ tests/default.nf.test --tag stub      # host stub tests
nf-test test tests/default.nf.test --profile docker --tag e2e   # full (needs images)
pytest tests/bin/test_bin.py                                 # bin unit tests
pre-commit run --all-files
nextflow run main.nf -preview --input tests/samplesheets/test.csv  # parse/DAG check
```
