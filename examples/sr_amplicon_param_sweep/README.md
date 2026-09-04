# superresolution-amplicon parameter sweep

Compares **superresolution-amplicon settings against one another** on a fixed synthetic
community, rather than comparing profilers. Every panel member sits at equal abundance;
what varies is the parameter grid, and every grid point is scored against the same
ground truth.

The two things that make this cheap enough to be a grid rather than a series of separate
runs:

- **The reads are mapped once.** mapseq is by far the most expensive stage of a
  superresolution-amplicon run, and no knob in the grid changes its result. Phase 1 maps
  them; phase 2 hands the classification back with the samplesheet's `mseq:` column and
  every grid point skips straight to the parts that differ.
- **The mis-mapping matrix is shared where it can be.** Grid points agreeing on every
  matrix knob (`mismapping_method`, `align_backend`, `align_tau`, `matrix_args`) belong
  to one reference set and build one matrix between them; they split only at the
  inference run. The shipped 3x2 grid is 6 profiles per benchmark dir off **3** matrices.

## Fill in before running

Everything is driven by **`config.yaml`** — no paths or metadata are hard-coded in the
Python scripts. Edit:

- `train.fastq_1` / `train.fastq_2` — the real R1/R2 the error model is trained from.
- `panel` — one entry per genome, each needing only `ssu` (a full-length 16S FASTA):
  it feeds both the in-silico PCR that generates the reads and the `community_v4`
  reference set superresolution-amplicon resolves against. `taxonomy` is derived from
  the `id` slug unless given; add `kingdom: archaea` for archaea. Two entries of the same
  species (a strain pair) are what make the sub-species question interesting.
- `primers` — the pair the reads are amplified with. superresolution-amplicon cuts its
  own reference amplicons with the same pair, so a mismatch here means not one read hits
  a reference.
- `sr_sweep.grid` — the sweep itself. Each axis is a list of named knob maps and the grid
  is their cartesian product; a point's name (the axis names joined with `.`) becomes
  part of its output filename. Valid knobs are the pipeline's `sr_settings` keys:
  `mismapping_method`, `align_backend`, `align_tau`, `matrix_args`, `infer_presence`,
  `infer_presence_prior`, `infer_presence_temp`, `inference_args`.

The scripts require **PyYAML** — run them with `python` (the same interpreter `run.sh`
uses). `python scripts/sr_sweep.py --selfcheck` sanity-tests the grid expansion without
touching the pipeline.

## Run

```bash
./run.sh
```

Two runs, back to back:

1. `generate_samplesheet.py` → `samplesheet.yaml` (+ `genomes/community.csv`), then
   `--step all`. Trains the error model, generates the amplicon reads at every
   `reads.subsample` depth, builds `community_v4` in-pipeline from the samplesheet's
   `databases:` block, and profiles each depth **once** at the nested pipeline's own
   defaults. This is the run that publishes each depth's
   `profiling/sr/<id>.obs.mseq.gz`.
2. `generate_sweep_samplesheet.py` → `sweep_samplesheet.yaml`, then `--step profile`.
   One row per benchmark dir, each carrying that dir's `mseq:`, plus a top-level
   `sr_settings:` list — the expanded grid — that the pipeline fans every row out over.

## Output

Per grid point, next to the truth it is scored against:

```
results/sr_amplicon_param_sweep/
  community.515-YF-806BR/
    truth.tsv
    community.515-YF-806BR.simulate.p01.sr_profile.tsv
    community.515-YF-806BR.exact.p01.sr_profile.tsv
    ...
    profiling/sr/community.515-YF-806BR.obs.mseq.gz     <- phase 1; reused by phase 2
    subsample_100000/...                                <- same, per depth
  mismapping/
    community_v4_sr_amplicon_515-YF-806BR-<digest>/     <- one per matrix knob combo
      mismapping_matrix.npz
      mismapping_provenance.json                        <- how that matrix was built
```

The `<digest>` in a reference-set directory is of the matrix knobs; `mismapping_provenance.json`
inside it spells them out, so a grid point is attributable without re-deriving it from
the config.

## Notes

- Only the **matrix** knobs split the reference set. Adding a third `infer:` point costs
  one extra inference run per benchmark dir and no extra matrix; adding a third `matrix:`
  point costs a matrix per benchmark-dir-independent reference set plus its inference runs.
- `sr_settings` is emitted top-level here (the default for every row). It also works per
  row, if different samples should be swept differently.
- A knob set both in `benchmark.config` (as `sr_amplicon_<knob>`) and in a grid point
  takes the grid point's value. `benchmark.config` deliberately sets none of them.
- Re-running phase 2 alone is fine and cheap: `python generate_sweep_samplesheet.py
  <outdir>` then `--step profile`. Widen the grid and the nested runs `-resume` the work
  the existing points already did.
- If phase 1's mapseq output is missing for a dir, that row simply omits `mseq:` and the
  nested run maps the reads itself — the script warns rather than failing.
