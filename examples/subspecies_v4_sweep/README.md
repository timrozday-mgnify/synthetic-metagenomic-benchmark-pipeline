# Sub-species V4 abundance sweep

Benchmarks how well a method resolves **two genomes of the same species** as their relative
abundance varies. 20 background species (one genome each, fixed at abundance 1) plus one
additional genome of the same species as one member. Every species carries the same target
abundance; only the split *between the two strains* changes across the 20 samples — swept
from 0/1 to 1/0 with logistic spacing (denser near the extremes; endpoints exactly 0/1), always
summing to 1 (one species' worth, like every other species).

- **Reads:** V4 amplicon (515-YF / 806BR), 2×300 bp paired, 500k pairs/sample.
- **Error model:** trained once from the real reads named in `config.yaml` (shared `train_id`).

## Fill in before running

Everything is driven by **`config.yaml`** — no paths or metadata are hard-coded in the Python
scripts. Edit:

- `train.fastq_1` / `train.fastq_2` — the real R1/R2 the error model is trained from.
- `panel` — one entry per genome (`id`, `species`, `amplicon`, `ssu`). Point `amplicon` at each
  genome's **pre-trimmed V4 amplicon FASTA** (mode=amplicon has no primer logic — it uses each
  record as a fragment directly) and `ssu` at its **pre-extracted full-length 16S FASTA** (needed
  to build the `aap`/mapseq profiler DB — V4 fragments are not valid mapseq input). Exactly one
  `species` must appear twice: that pair is the sweep target. `id`s must be unique. `taxonomy` is
  derived from the `species` slug unless given as a 3-rank `Kingdom;Genus;species` string (e.g.
  `taxonomy: "Bacteria;Bacteroides;fragilis"`); add `kingdom: archaea` for archaea. Add `genome:`
  per entry (and `sylph` to `database.profilers`) to also build a sylph WGS DB.

The scripts require **PyYAML** — run them with `python` (the same interpreter `run.sh` uses).

## Run

```bash
./run.sh          # config.yaml -> samplesheet.yaml, then one `--step all` run
```

`generate_sweep.py` reads `config.yaml` and writes a single `samplesheet.yaml` containing a
top-level `databases: community_v4` block **and** the 20 sweep samples (each with inline
`train_*`, read-geometry columns, `profiler`, and `database: community_v4`). `run.sh` then does
one `nextflow run --step all`, which — because training is deduped by `train_id` — trains the
error model once, generates every sample, **builds the `community_v4` profiler DB in-pipeline
from the samplesheet's `databases:` block**, and profiles each sample against it. No out-of-band
DB build step is needed.

## Notes

- `reads.num_reads = 1000000` = 500k pairs (pipeline counts paired reads as pairs×2). Halve if
  genome-blender counts pairs directly.
- `train.platform = hq-illumina`; 2×300 MiSeq amplicon tails may fit `lq-illumina` — one line in
  `config.yaml`.
- Paired 2×300 read geometry lives in `config.yaml`'s `reads:` block (`mode`, `paired_end`,
  `read_length_mean`, `read_length_variance`), emitted into the samplesheet — no `ext.args`.
- The three scripts share `scripts/sweep_config.py` (config loader + `databases:`-block builder);
  run any script's `--selfcheck` to sanity-test its logic without the pipeline.

## Re-profiling without regenerating

To re-profile already-generated reads (e.g. to compare profilers) without regenerating them:

```bash
python generate_profile_samplesheet.py "$OUTDIR"   # reads config.yaml + samplesheet.yaml
nextflow run ../../main.nf -profile docker -c benchmark.config \
    --step profile --input profile_samplesheet.yaml --outdir "$OUTDIR" --seed 42
```

`generate_profile_samplesheet.py` emits the same `databases:` block plus one row per
(sample, profiler) — every profiler in `database.profilers` — with `benchmark_dir` at
`<results_dir>/<sample>`. See the root README's "Named sequence collections (`databases:`)".

## Building the profiler DB out-of-band (optional)

`scripts/build_profiling_dbs.py` builds the same sylph/mapseq DBs **outside** the pipeline from
the same `config.yaml`, writing `sylph_databases.config` / `aap.config` for a config-based
`database:` run. It runs under Docker by default or Singularity/Apptainer for HPC without Docker:

```bash
python scripts/build_profiling_dbs.py --runtime singularity   # or --runtime docker (default)
```

This is **not** needed for the default in-pipeline path above — use it only if you specifically
want the DB prebuilt. Which DBs are built follows `config.yaml`'s `database.profilers`.
