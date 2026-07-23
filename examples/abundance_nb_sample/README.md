# Negative-binomial abundance sample

Benchmarks profiling accuracy across **random community compositions**. Unlike the
`subspecies_v4_sweep` example (which fixes every species at equal abundance and sweeps
the split of one same-species pair), here **every panel genome is drawn independently
per sample**: a Bernoulli decides whether it's present, and if so a negative-binomial
sets its abundance. Each of the 20 samples is thus a different random mixture — some
genomes dominant, some rare, some absent. Draws are seeded (`sampling.seed`) so a rerun
reproduces the same communities.

- **Reads:** one set per entry in `config.yaml`'s `generation_modes:`. The shipped config
  sweeps **two** modes — `wgs` (shotgun from genomes, 150 bp, profiled by sylph) and
  `amplicon_16s` (V4 515-YF/806BR amplicons extracted from full-length 16S by in-silico PCR,
  2×300 bp, profiled by aap) — so each of the 20 samples yields both a `.wgs` and a
  `.amplicon_16s` sample (40 total). Both modes of a sample share its drawn abundances.
- **Error model:** trained once from the real reads named in `config.yaml` (shared `train_id`).

## The abundance model

Each genome is drawn independently per sample in two steps:

1. **Presence** — `present ~ Bernoulli(presence)`. `sampling.presence` (default 1.0 = always
   present) is the probability the genome shows up at all; lower it to make presence/absence
   itself vary across samples.
2. **Abundance** — if present, `abundance ~ NegativeBinomial(mean, dispersion)`, parameterized
   by mean + dispersion (size *r*): variance = mean + mean²/*r*. If absent, abundance is 0.

- small dispersion (*r* ≈ 0.1–1) → overdispersed/patchy: a few genomes dominate, others near zero;
- large dispersion (*r* ≫ mean) → Poisson-like: abundances cluster near the mean.

Set the globals under `sampling` (`presence`, `negative_binomial: {mean, dispersion}`), and
override `presence` and/or `nb: {mean?, dispersion?}` per genome inline in `panel` (e.g. an
always-present dominant taxon, or a rare one). A sample that ends up all-zero is redrawn.
Integer draws are written straight to the genomes CSV as relative abundances (genome-blender
normalizes internally).

## Fill in before running

Everything is driven by **`config.yaml`** — no paths or metadata are hard-coded in the Python
scripts. Edit:

- `train.fastq_1` / `train.fastq_2` — the real R1/R2 the error model is trained from.
- `sampling` — `n_samples`, `seed`, `presence` (Bernoulli present-probability), and
  `negative_binomial: {mean, dispersion}`.
- `panel` — one entry per genome (`id`, `species`, and the FASTA fields your modes/profilers
  use). A member needs only the fields referenced: `genome` (full assembly — for `source: genome`
  shotgun/primer-PCR generation and the sylph DB), `ssu` (full-length 16S — for `source: ssu`
  generation and the `aap`/mapseq DB), and/or `amplicon` (pre-trimmed FASTA — for `source: amplicon`
  direct generation). `id`s must be unique; species may repeat freely (no sweep pair here).
  `taxonomy` is derived from the `species` slug unless given as a 3-rank `Kingdom;Genus;species`
  string (e.g. `taxonomy: "Bacteria;Bacteroides;fragilis"`); add `kingdom: archaea` for archaea.
  Add `presence:` and/or `nb: {mean?, dispersion?}` to override the global presence/NB for
  that genome.
- `generation_modes` — the read-generation styles to sweep; **every sample is emitted once per
  mode**. Each mode sets `name` (sample/CSV suffix), `source` (`genome`/`ssu`/`amplicon`),
  `mode` (`shotgun`/`amplicon`/`long`), `profiler`, optional `reads:` overrides, and — when
  extracting amplicons from `genome`/`ssu` — a `primers:` list (or TSV path) run as in-silico PCR.
  The four supported combos: genomes+shotgun; genomes+primers+amplicon; 16S+primers+amplicon;
  pre-trimmed amplicon (`source: amplicon`, no primers). `database.profilers` must be the union of
  the modes' profilers. Omit `generation_modes:` and set a single `reads.mode` for the legacy
  single-style behaviour.

The scripts require **PyYAML + numpy** — run them with `python` (the same interpreter `run.sh`
uses).

## Run

```bash
./run.sh          # config.yaml -> samplesheet.yaml, then one `--step all` run
```

`generate_samplesheet.py` reads `config.yaml`, draws the NB abundance matrix (seeded), and
writes a single `samplesheet.yaml` containing a top-level `databases: community_nb` block
**and** the samples — 20 samples × `len(generation_modes)` rows (each with inline `train_*`,
per-mode read-geometry columns, `mode`, `profiler`, and `database: community_nb`), plus
`genomes/sample_NN.<mode>.csv`. `run.sh` then does one `nextflow run --step all`, which —
because training is deduped by `train_id` — trains the error model once, generates every
sample, **builds the `community_nb` profiler DB in-pipeline from the samplesheet's `databases:`
block**, and profiles each sample against it. No out-of-band DB build step is needed.

## Notes

- `reads.num_reads = 1000000` = 500k pairs (pipeline counts paired reads as pairs×2). Halve if
  genome-blender counts pairs directly.
- `train.platform = hq-illumina`; 2×300 MiSeq amplicon tails may fit `lq-illumina` — one line in
  `config.yaml`.
- Read geometry lives in `config.yaml`'s shared `reads:` block (`paired_end`, `read_length_mean`,
  `read_length_variance`, `subsample`), which each `generation_modes` entry may override (e.g. the
  `wgs` mode drops to 150 bp); the merged values are emitted per row — no `ext.args`.
- The three scripts share `scripts/nb_config.py` (config loader + NB sampler + `databases:`-block
  builder); run any script's `--selfcheck` to sanity-test its logic without the pipeline.

## Re-profiling without regenerating

To re-profile already-generated reads (e.g. to compare profilers) without regenerating them:

```bash
python generate_profile_samplesheet.py "$OUTDIR"   # reads config.yaml + samplesheet.yaml
nextflow run ../../main.nf -profile docker -c benchmark.config \
    --step profile --input profile_samplesheet.yaml --outdir "$OUTDIR" --seed 42
```

`generate_profile_samplesheet.py` emits the same `databases:` block plus one row per already-
generated sample, re-profiling each with the profiler its generation mode used (read from
`samplesheet.yaml`), with `benchmark_dir` at `<results_dir>/<sample>`. See the root README's
"Named sequence collections (`databases:`)".

## Building the profiler DB out-of-band (optional)

`scripts/build_profiling_dbs.py` builds the same sylph/mapseq DBs **outside** the pipeline from
the same `config.yaml`, writing `sylph_databases.config` / `aap.config` for a config-based
`database:` run. It runs under Docker by default or Singularity/Apptainer for HPC without Docker:

```bash
python scripts/build_profiling_dbs.py --runtime singularity   # or --runtime docker (default)
```

This is **not** needed for the default in-pipeline path above — use it only if you specifically
want the DB prebuilt. Which DBs are built follows `config.yaml`'s `database.profilers`.
