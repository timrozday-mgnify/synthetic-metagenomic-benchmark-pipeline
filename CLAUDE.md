# CLAUDE.md

Guidance for Claude Code / developers working in this repo.

## What this is

A Nextflow DSL2 pipeline that trains a skiver sequencing-error profile from a
natural metagenome and applies it (via genome-blender) to reference genomes to
emit synthetic reads + a ground-truth BAM + a ground-truth abundance profile, and
(optionally, `--step all`/`profile`) taxonomically profiles the reads next to that
ground truth with sylph (WGS) or the EBI amplicon-analysis-pipeline (amplicon).
Profiler databases can be built in-pipeline from named sequence collections defined
in a samplesheet `databases:` block, or supplied pre-built.

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

- `main.nf` — inline samplesheet parse (bare list, or a map with `samples:` +
  optional `databases:` block); builds `ch_samples`, a train_id-deduped `ch_train`,
  and `ch_db_specs` (named collections referenced by a sample), calls the workflow.
  Resolves the FASTAs referenced by each `genomes_csv`/collection so Nextflow stages them.
- `workflows/synthetic_metagenomic_benchmark.nf` — top workflow; joins trained
  models back to samples by `train_id`; calls BUILD_DATABASES then PROFILE.
- `subworkflows/local/train_error_model/` — SKIVER_DUMP → SKIVER_TRAIN.
- `subworkflows/local/build_databases/` — build (or resolve pre-built) profiler DBs
  from named sequence collections; emits sylph/mapseq DBs keyed by collection name.
- `subworkflows/local/profile/` — per-sample profiling; selects a DB by the sample's
  `database` name (built collection → params.sylph_databases/aap_config fallback).
- `modules/local/` — SKIVER_DUMP, SKIVER_TRAIN, GENOME_BLENDER_GENERATE, GROUND_TRUTH,
  sylph/build_db (SYLPH_BUILD_DB), amplicon_analysis (RUN_AAP), and the mapseq DB
  builders mapseq/prep (MAPSEQ_PREP), mapseq/build_db (MAPSEQ_CLUSTER), mapseq/otu (MAPSEQ_OTU).
- `bin/` — helper scripts staged onto PATH: `build_model_config.py`,
  `pick_best_model.py`, `rewrite_genomes_csv.py`, `build_truth.py`,
  `normalize_sylph_profile.py`, `build_mapseq_refs.py`, `build_mapseq_otu.py`,
  `write_aap_config.py`.
- `containers/` — Dockerfiles for the two images (built from `vendor/` submodules).
- `conf/` — `base.config` (resource labels), `modules.config` (publish layer).

## Container ↔ module

| Module | Image | Tool |
|--------|-------|------|
| SKIVER_DUMP, SKIVER_TRAIN | `smb-skiver` | `skiver dump`, `train_context_error_models.py`, `fit_phred_calibration.py` |
| GENOME_BLENDER_GENERATE | `smb-genome-blender` | `generate-reads` (+ `skiver-generate` subprocess) |
| GROUND_TRUTH | `smb-genome-blender` | `pysam` (sort/index/idxstats) — reused, no separate samtools image |
| MAPSEQ_PREP, MAPSEQ_OTU | `smb-skiver` | `build_mapseq_refs.py` / `build_mapseq_otu.py` (stdlib python; bin/ on PATH) |

Images above are `ghcr.io/timrozday-mgnify/<image>:${params.<image>_tag}`. Third-party
tools use pinned biocontainers via the singularity/docker ternary: `SYLPH_BUILD_DB`
(sylph 0.9.0) and `MAPSEQ_CLUSTER` (mapseq 2.1.1b). `RUN_AAP` runs on the host
(`executor local`) and launches the nested EBI amplicon-analysis-pipeline itself.

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
- **Profiler DB selection** keys off a sample's `database` name. Names defined in the
  samplesheet `databases:` block are built (or their `path:` dir resolved) by
  BUILD_DATABASES and joined into PROFILE by name; `builtSylphNames`/`builtMapseqNames`
  (computed in `main.nf`) decide built-vs-fallback so `params.sylph_databases` / `self`
  / `params.aap_config` still work for unlisted names. Built DBs publish to
  `${outdir}/databases/<name>/` in the exact layout a `path:` prebuilt entry expects,
  so a run's output is reusable as another run's input. Only the DB type a sample
  actually references (`sylph` vs `aap`) is built. mapseq collections need explicit
  per-sequence `taxonomy` + a pre-extracted `ssu` (no barrnap step in-pipeline).
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
nextflow run main.nf -preview --input tests/samplesheets/test.yaml  # parse/DAG check
```
