# Sub-species V4 abundance sweep

Benchmarks how well a method resolves **two genomes of the same species** as their relative
abundance varies. 20 background species (one genome each, fixed at abundance 1) plus one
additional genome of the same species as one member. Every species carries the same target
abundance; only the split *between the two strains* changes across the 20 samples — swept
from 0/1 to 1/0 with logistic spacing (denser near the extremes; endpoints exactly 0/1), always
summing to 1 (one species' worth, like every other species).

- **Reads:** V4 amplicon (515-YF / 806BR), 2×300 bp paired, 500k pairs/sample.
- **Error model:** trained once from the real SC2200627 Illumina reads (shared `train_id`).

## Fill in before running

Edit `generate_sweep.py`:
- `TRAIN_FASTQ_1/2` — the real R1/R2 (defaults point at the SC2200627 lane1C5 pair).
- `PANEL` — 21 rows `(genome_id, species, v4_fasta_path)`. Point `fasta_path` at each genome's
  **pre-trimmed V4 amplicon FASTA** (mode=amplicon has no primer logic — it uses each record as a
  fragment directly). The last row is the additional genome; set its species to match one member
  and give it the second strain's V4 FASTA. `genome_id`s must be unique.

## Run

```bash
./run.sh          # regenerates inputs, then runs train -> generate -> profile in sequence
```

`run.sh` runs the three steps explicitly so the error model is trained once and
reused (and the profiling step is included). `generate_sweep.py` writes
`train_samplesheet.yaml` (consumed by `--step train`) and a `samplesheet.yaml`
whose rows carry `error_model_dir` pointing at the trained model, so
`--step generate` skips retraining. Fill in the reference genomes in
`scripts/build_profiling_dbs.py`'s `GENOMES` before running (see "Profiling this sweep").

## Notes

- `num_reads = 1_000_000` = 500k pairs (pipeline counts paired reads as pairs×2). Halve if
  genome-blender counts pairs directly.
- `platform = hq-illumina`; 2×300 MiSeq amplicon tails may fit `lq-illumina` — one word in the script.
- Paired 2×300 read geometry lives in the samplesheet columns (`mode=amplicon`, `paired_end=true`,
  `read_length_mean=300`, `read_length_variance=0`), written by `generate_sweep.py` — no `ext.args`.

## Profiling this sweep

After a `--step generate` (or `all`) run, each sample is profiled against a custom
`community_v4` database built from this sweep's own reference genomes.

- **Reference genomes** — fill in `GENOMES` at the top of `scripts/build_profiling_dbs.py`:
  a full genome FASTA (for sylph) and a full-length 16S FASTA + taxonomy (for aap/mapseq)
  per `PANEL` genome_id. This dict is the single source of truth for the DB; the profile
  samplesheet reads it directly. V4 amplicon fragments are NOT valid mapseq input — it needs
  full-length 16S. `ssu_fasta=None` is only OK if you drop `aap` from `PROFILERS` (the
  in-pipeline path has no barrnap step, unlike the standalone build script).
- `generate_profile_samplesheet.py [results_dir]` — builds `profile_samplesheet.yaml` for a
  `--step profile` re-run. It emits a top-level `databases: community_v4` block (a sequence
  collection built from `GENOMES`) and one row per sample per profiler (`sylph` + `aap` by
  default, set by `PROFILERS`), each with `database: community_v4` and `benchmark_dir` at
  `<results_dir>/<sample>`. The pipeline builds the sylph `.syldb` and mapseq DB itself during
  the run — no separate build step or config. See the root README's "Named sequence
  collections (`databases:`)" section.
- `scripts/build_profiling_dbs.py` (+ `_singularity.py`) — optional out-of-band builder that
  produces the same DBs plus `sylph_databases.config` / `aap.config` for a config-based run
  (e.g. on HPC where you'd rather prebuild). Not needed for the in-pipeline path above, but it
  still holds the `GENOMES` reference-genome registry both paths share.
