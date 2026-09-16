# Benchmark reports

Parameterised [Quarto](https://quarto.org) report that communicates the results of a
synthetic-metagenomic-benchmark-pipeline run: **mis-mapping between ground-truth
genomes**, **ground-truth vs detected abundance**, an optional **swept-pair** view, and
a **performance-across-subsampling** summary. Interactive Plotly plots with tooltips.

The template is generic — it works for amplicon (mapseq) and WGS (competitive mapping)
runs, single or multiple amplified regions, with or without an abundance sweep. When a
run has several **assays** (e.g. multiple amplicon primer sets and/or a WGS arm, each a
`<sweep_point>.<assay>` dir), every visual section is split into one tab per assay so the
same panel can be compared across primers/WGS. Sections and sources that have no data
(e.g. no swept pair, no per-read profiler assignments, no sylph) are omitted automatically.

## Layout

```
reports/
  Taskfile.yml                          # preprocess + render, one task per run
  scripts/benchmark_preprocess.py       # heavy file/BAM/mseq walking -> tidy CSVs
  templates/benchmark_report_template.qmd
  envs/quarto-requirements.txt          # driver env (papermill)
  envs/basicpython-requirements.txt     # execution kernel (pandas, plotly, pysam, ...)
```

Run data and the rendered HTML live **outside** this repo, in the sibling
`../synthetic-metagenomic-benchmark-pipeline_runs/<RUN>/` directory. Nothing run-specific
is written back into the repo.

## Workflow

Preprocessing walks the run tree once and writes four tidy tables into the run dir:

| File | Contents |
|------|----------|
| `abundance.csv` | long: `sample, assay, depth, sweep_x, genome_id, is_sweep_pair, target/realized_rel_abundance, detected_profiler_rel_abundance, detected_reference_rel_abundance, detected_sylph_rel_abundance` |
| `mismapping.csv` | long: `sample, assay, depth, source, truth_genome, assigned_genome, reads, frac_of_truth` |
| `summary.csv` | per `assay × depth × source`: `l1_error_per_sample, pearson_r, mismapping_rate, pair_mismapping_rate, n_reads` |
| `meta.json` | swept pair, assays, depths, samples, sample→sweep_x map, detection sources |

Superresolution's simulated mis-mapping matrices are pipeline outputs rather than
preprocessed tables. They are published once per reference set at
`<pipeline-dir>/mismapping/<reference-set>/mismapping_matrix.{csv,npz}`; the benchmark
report discovers and summarises them from that location.

The report then reads those CSVs (all heavy compute is already done) and builds the plots.

Each sample dir is a `<sweep_point>.<assay>` cell (e.g. `S10_a0.42.amplicon_515YF-806BR_16s.515-YF-806BR`
or `S10_a0.42.wgs`); the sweep point is shared across assays, the assay label (`515-YF-806BR`,
`WGS`, …) is the facet. Up to three detection sources are joined against the ground truth
per cell:

- **profiler** — per-read mapseq assignments (`*.mseq.gz`); `query` encodes the origin
  genome, `dbhit` the assigned genome. Amplicon (amplicon-analysis-pipeline) path.
- **reference** — the ground-truth `*.sorted.bam`; read names encode the origin genome
  and each reference contig is `<genome>:<contig>`. Works for WGS too. In runs where
  reads are mapped to their origin reference this matrix is diagonal by construction; it
  becomes informative under competitive mapping.
- **sylph** — the WGS sylph profile (`*.sylph_profile.tsv`); abundance-only (no per-read
  assignments, so no confusion matrix). Sylph reports reference **accessions**, mapped
  back to community `genome_id`s via the run's `samplesheet.yaml` (`id`↔`genome` fasta
  basename); `--samplesheet` overrides its location.

Every source is optional per cell, so a run missing one (or a whole assay missing a
source) still produces what it can.

## Running it

The `Taskfile.yml` drives everything (needs [go-task](https://taskfile.dev)):

```bash
cd reports
task subspecies-v4          # preprocess + render the bundled example run
```

Other entry points:

```bash
task preprocess RUN_DIR=<run> PIPELINE_DIR=results/<sweep> RUN_LABEL="..."   # CSVs only
task render     RUN_DIR=<run> PIPELINE_DIR=results/<sweep> RUN_LABEL="..."   # re-render only
task subspecies-v4 DIR=/path/to/other/runs                                   # point at a different runs dir
```

Output: `<run>/benchmark_report_template.html` — a single self-contained file
(`embed-resources: true`).

### Adding a run

Copy the `subspecies-v4` task and change its `vars` (`RUN_DIR`, `PIPELINE_DIR`,
`RUN_LABEL`, and optionally `MSEQ_GLOB` or `SWEEP_PAIR` to override auto-detection).

### Preprocessing parameters (`benchmark_preprocess.py`)

| Flag | Default | Meaning |
|------|---------|---------|
| `--run-dir` | — | Run root (contains `--pipeline-dir`). |
| `--pipeline-dir` | `results/subspecies_v4_sweep` | Where the per-sample dirs live. |
| `--mseq-glob` | `profiling/aap/*/taxonomy-summary/*/*.mseq.gz` | Per-cell glob for mapseq files. |
| `--sweep-pair` | auto-detect | `id1,id2` override for the swept genomes. |
| `--run-label` | `""` | Human-readable label stored in `meta.json`. |
| `--output-dir` | `--run-dir` | Where the CSVs are written. |
| `--demo` | — | Self-check the parsing/detection logic and exit. |

The swept pair is auto-detected as the genome_ids whose target abundance varies across
samples — no config needed for the common single-pair sweep.

### Report parameters (`-P` to `quarto render`)

`run_label`, `run_dir`, `pipeline_dir`, `abundance_csv`, `mismapping_csv`, `summary_csv`,
`meta_json`, `output_dir`, `presence_prob_thresholds`, `subspecies`, and
`presence_confusion_threshold`. The CSV params are resolved relative to `run_dir`.

`presence_prob_thresholds` is a YAML list of superresolution posterior cutoffs used for the
presence Jaccard curve. Set `subspecies` to comma-separated genome IDs to add a 2×2
truth-versus-call presence/absence matrix for each one; `presence_confusion_threshold`
sets their calling cutoff (default `0.5`). With the Taskfile, pass the equivalents as
`PRESENCE_THRESHOLDS`, `SUBSPECIES`, and `PRESENCE_CONFUSION_THRESHOLD`.

## Superresolution settings sweep (`sr_settings_report.qmd`)

A second, narrower report for runs that fan a sample out across `sr_settings`: how well
each parameter set resolves a near-identical **sub-species pair**, and what mis-mapping
matrix each one infers against.

```bash
task sr-settings RUN_DIR=<run> PIPELINE_DIR=results/<sweep> RUN_LABEL="..."
```

It reads the run tree directly (per-sample `*.truth.tsv`, `profiling/sr/*.inferred_composition.csv`
and `*.inference_diagnostics.csv`, plus the published `mismapping/*/mismapping_matrix.npz`)
— there is no preprocessing step, because a settings sweep is a few hundred rows.

Settings are the `<sample>.<setting>` suffix on the profile files; the un-suffixed run is
reported as `default`. Matrices are matched to the settings that used them by the
`mean_diagonal` recorded in `inference_diagnostics.csv`, and shown aggregated to genomes
(`genome_mismapping_matrix` in `scripts/mismapping_plots.py`) as well as at their
published per-16S-copy resolution.

Parameters: `run_label`, `run_dir`, `pipeline_dir`, `samplesheet` (for the settings table),
`superresolution_mismapping_dir`, `subspecies` (comma-separated; defaults to whichever
genomes' target abundance varies across the sweep), `presence_threshold`. Via the
Taskfile: `SAMPLESHEET`, `SUBSPECIES`, `PRESENCE_THRESHOLD`.

Output: `<run>/sr_settings_report.html`.

## GTDB panel reinterpretation vs a custom database (`sr_panel_sweep_report.qmd`)

The report for `examples/sr_amplicon_panel_sweep`: the same reads profiled against a
custom database of the panel (`custom.*` grid points) and through GTDB labels
reinterpreted via the panel (`gtdb_panel.*`), all scored against one ground truth.

```bash
task sr-panel-sweep                       # the sr_amp_panel_sweep run, default layout
task sr-panel-sweep RUN_DIR=<run> PIPELINE_DIR=<dir> SAMPLESHEET=<yaml>
```

Like the settings-sweep report it reads the run tree directly: per-dir `*.truth.tsv`,
`profiling/sr/<id>.<setting>.{inferred_composition,inference_diagnostics}.csv` and
`fit_diagnostics.json`, the published `community` matrices under `mismapping/`, and
`pipeline_info/execution_trace_*.txt`. It never opens the reads, BAMs, `*.mseq.gz`,
`posterior_draws.npz` or phase 1's GTDB-scale `*.map.*` outputs, so from a results
tarball extract only those tables. Settings and the panel come from the sweep
samplesheet's row-level `sr_settings:` and `databases:` block.

Scores come from the example's own `scripts/score_sweep.py` (`example_scripts_dir`), so
they match `panel_sweep_scores.csv`. The observed baselines are `custom.mapseq_only` and
`gtdb_panel.home_labels`: the `observed_rel_abundance` of each arm, i.e. no inference. It
also flags grid points whose compositions are identical in every dir, meaning a knob never
reached the nested run.

Parameters: `run_label`, `run_dir`, `pipeline_dir`, `samplesheet`,
`superresolution_mismapping_dir`, `example_scripts_dir`, `strain_pair` (comma-separated;
defaults to the panel id ending `_strain…` and the id it extends), `reference_method` (the
`custom` point the paired TV deltas are taken against; defaults to the best by median TV),
`presence_threshold`. Via the Taskfile: `STRAIN_PAIR`, `REFERENCE_METHOD`,
`PRESENCE_THRESHOLD`.

Output: `<run>/sr_panel_sweep_report.html`. Figures are rendered without plotly.js and the
library is embedded once, which keeps the file around 9 MB rather than ~5 MB per figure.

## Environments

Two conda/mamba envs, matching the `fermentor-run-reports` convention:

```bash
# 1. Quarto CLI itself: https://quarto.org/docs/get-started/
# 2. Driver env — runs `quarto`, needs papermill to apply -P overrides
mamba create -n quarto python=3.12
mamba run -n quarto pip install -r envs/quarto-requirements.txt
# 3. Execution kernel — runs the report's code cells and the preprocess script
mamba create -n basicpython python=3.12
mamba run -n basicpython pip install -r envs/basicpython-requirements.txt
mamba run -n basicpython python -m ipykernel install --user \
  --name basicpython --display-name "Python (basicpython)"
```

The Taskfile references both by absolute path (`QUARTO_ENV`, `BASICPYTHON_ENV`) — adjust
those vars if your envs live elsewhere.

## Gotchas

- **Render via the `quarto` env, not `basicpython`.** `quarto` runs the binary and needs
  `papermill` to apply `-P`; rendering from an env without it fails with *"The papermill
  package is required for processing --execute-params"*.
- **`basicpython` is the execution kernel** (named in the template's `jupyter:`
  frontmatter) and also runs `benchmark_preprocess.py`. It needs `pysam` for the BAM path.
- **Pass absolute paths.** Quarto executes the template from its own directory, so
  `run_dir`/`--output-dir` should be absolute (the Taskfile builds them from `ROOT_DIR`).
- **Glob `-P` values starting with `*`** must be quoted `"'*/...'"` — Quarto parses `-P`
  as YAML and a bare leading `*` is a YAML alias. (Only relevant if you override a glob.)
- `embed-resources: true` yields one standalone HTML; the output keeps the template
  filename: `<run>/benchmark_report_template.html`.
