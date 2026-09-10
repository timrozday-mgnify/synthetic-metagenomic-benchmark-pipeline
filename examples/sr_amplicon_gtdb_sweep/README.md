# superresolution-amplicon parameter sweep against GTDB

Asks whether a **large external reference set is fit for purpose**: the synthetic
community is 20 known genomes, the reference set the reads are profiled against is
GTDB (or any other pre-built SSU set), and the superresolution-amplicon parameters are
swept over the same reads so the settings can be compared against one ground truth.

It is two existing examples spliced together, and the scripts import them rather than
copy them:

- the **communities** are `examples/abundance_nb_sample`'s — every panel genome
  presence-gated (`Bernoulli(0.6667)`) then negative-binomial drawn, independently per
  sample, so each of the 20 samples is a different random mixture with about a third of
  the panel genuinely absent;
- the **sweep** is `examples/sr_amplicon_param_sweep`'s — a `sr_sweep.grid:` expanded
  into the samplesheet's `sr_settings:` list, two phases so the reads are mapped once.

What is new is the database: `database.path:` is a **pre-built** reference set, so
nothing is built during the run.

## Why this community rather than a sub-species sweep

`examples/sr_amplicon_param_sweep` sweeps the split of one near-identical strain pair —
the right question for a reference set that holds only the community. Against GTDB the
harder question is the opposite one: every reference that is *not* one of the 20 genomes
is absent from every sample, and a setting is only usable if it says so. The NB
community's own absent genomes (about 7 of 20 per sample, differing per sample) give the
same measurement inside the panel, where the truth is unambiguous.

## The database (`path:`, no build step)

A pre-built database is a directory laid out the way the pipeline publishes its own
(`<outdir>/databases/<name>/`). For `sr_amplicon` that is exactly one file:

```
<database.path>/<database.name>_ssu.sr_refs.fasta      >{genome_id}|{copy}|{original header}
```

`BUILD_DATABASES` resolves it by glob and hands it to the nested runs. Build it once, out
of band, and reuse it: a GTDB-scale FASTA is not something to rebuild per run.

**The panel ids must be the reference set's ids.** `truth.tsv` is written per panel `id`
and a profile is written per reference `genome_id` (the part before the first `|`). If
they disagree the run still succeeds and every genome scores zero. When the reference set
is GTDB, that means the panel ids are GTDB accessions — the species slugs shipped in
`config.yaml` are placeholders. `generate_samplesheet.py` reads the FASTA's headers and
warns about ids it cannot find, whenever the directory exists.

To run the same sweep against a small collection the pipeline builds itself (a cheap
control, on the same code path), delete `path:` and set `sequences_from_panel: true`.

## Fill in before running

Everything is driven by **`config.yaml`** — no paths or metadata are hard-coded in the
Python scripts. Edit:

- `train.fastq_1` / `train.fastq_2` — the real R1/R2 the error model is trained from.
- `database.path` — the pre-built reference set directory, and `database.name` to match
  the `<name>_ssu.sr_refs.fasta` inside it.
- `panel` — one entry per community genome: `id` (**the reference set's id for it**),
  `species`, and `ssu` (a full-length 16S FASTA, which the V4 amplicons are cut from).
- `sampling` — `n_samples`, `seed`, `presence`, `negative_binomial: {mean, dispersion}`.
- `sr_sweep.grid` — the sweep itself, each axis a list of named knob maps, the grid their
  cartesian product. Valid knobs are the pipeline's `sr_settings` keys.
- `primers` (under `generation_modes`) — superresolution-amplicon cuts its own reference
  amplicons with the same pair, so a mismatch means not one read hits a reference.

The scripts need **PyYAML + numpy**. `python scripts/gtdb_sweep.py --selfcheck` exercises
the database block, the depth parsing, the id check and the grid expansion without
touching the pipeline.

## Run

```bash
./run.sh
```

Two runs, back to back:

1. `generate_samplesheet.py` → `samplesheet.yaml` (+ one `genomes/sample_NN.amplicon_16s.csv`
   per community), then `--step all`. Trains the error model, generates the amplicon reads
   at every `reads.subsample` depth, and profiles each depth **once** against the
   pre-built reference set. This is the run that publishes each depth's
   `profiling/sr/<id>.obs.mseq.gz`.
2. `generate_sweep_samplesheet.py` → `sweep_samplesheet.yaml`, then `--step profile`.
   One row per benchmark dir, each carrying that dir's `mseq:`, plus a top-level
   `sr_settings:` list — the expanded grid — that the pipeline fans every row out over.

## What it costs

The shipped config is 20 samples × 2 depths = 40 benchmark dirs × 6 grid points = **240
profiles**, off **3** mis-mapping matrices. Three things keep that affordable, and all
three are about the size of the reference set rather than the size of the grid:

- **The reads are mapped once.** mapseq against a GTDB-scale set is by far the most
  expensive stage and no swept knob changes it. Phase 1 maps; phase 2 hands the
  classification back through `mseq:` and every grid point skips straight to the parts
  that differ. A dir with no phase-1 mseq is not an error — the nested run maps it
  itself, and the generator warns.
- **The matrix is shared where it can be.** Grid points agreeing on every matrix knob
  (`mismapping_method`, `align_backend`, `align_tau`, `matrix_args`) build one matrix
  between them and split only at the inference run — which is why the presence prior is
  swept here rather than in a second run.
- **`mismapping_method: simulate` is left out of the grid.** It simulates and maps
  `sim_n_per_ref` reads for *every reference in the set*: fine for a 20-genome
  collection, not for GTDB. The grid sweeps the alignment backends instead
  (`exact-hash`, and `kmer` at τ=1 and τ=2). `benchmark.config` pins phase 1 to the
  `exact` point for the same reason — phase 1 has no `sr_settings` to fan out over, so
  it would otherwise use the nested pipeline's `simulate` default.

Phase 1's matrix is built separately from phase 2's even when the knobs match: the
reference-set key mixes in the matrix mode only when `sr_settings` is in play, so the two
runs do not share a nested cache. That is one extra `exact` matrix, not one extra
mapping pass.

## Output

```
results/sr_amplicon_gtdb_sweep/
  S01.amplicon_16s.515-YF-806BR/                       <- one per community
    truth.tsv
    S01.amplicon_16s.515-YF-806BR.exact.p01.sr_profile.tsv
    S01.amplicon_16s.515-YF-806BR.kmer1.p01.sr_profile.tsv
    ...
    profiling/sr/S01.amplicon_16s.515-YF-806BR.obs.mseq.gz   <- phase 1; reused by phase 2
    subsample_100000/...                                     <- same, per depth
  S02.amplicon_16s.515-YF-806BR/ ...
  mismapping/
    gtdb_ssu_r220_sr_amplicon_515-YF-806BR-<digest>/   <- one per matrix knob combo
      mismapping_matrix.npz
      mismapping_provenance.json                       <- how that matrix was built
```

No `databases/` directory: nothing was built.

## Reporting

`reports/` renders the settings comparison from that tree:

```bash
cd ../../reports
task sr-settings RUN_DIR=<run> PIPELINE_DIR=results/sr_amplicon_gtdb_sweep \
    RUN_LABEL="GTDB sr_amplicon parameter sweep"
```

Its sub-species sections auto-detect genomes whose target abundance varies across the
sweep; there are none here (every genome varies, independently), so pass `SUBSPECIES=`
a pair of ids to focus them, or read the per-setting scores and the mis-mapping matrices,
which need no sweep pair.

## Notes

- Only the **matrix** knobs split the reference set. Adding an `infer:` point costs one
  extra inference run per benchmark dir and no extra matrix.
- A knob set both in `benchmark.config` (as `sr_amplicon_<knob>`) and in a grid point
  takes the grid point's value — which is why phase 2 ignores the pins phase 1 needs.
- Re-running phase 2 alone is fine and cheap: `python generate_sweep_samplesheet.py
  <outdir>` then `--step profile`. Widen the grid and the nested runs `-resume` the work
  the existing points already did.
- `sr_settings` is emitted top-level here (the default for every row). It also works per
  row, if different samples should be swept differently.
