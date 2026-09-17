# superresolution-amplicon parameter sweep against SILVA

Asks whether a **large external reference set is fit for purpose**: the synthetic
community is 20 known genomes, the reference set the reads are profiled against is
SILVA SSU Ref NR99 (or any other pre-built SSU set with a taxonomy), and the
superresolution-amplicon parameters are swept over the same reads so the settings can be
compared against one ground truth.

It is two existing examples spliced together, and the scripts import them rather than
copy them:

- the **communities** are `examples/abundance_nb_sample`'s — every panel genome
  presence-gated (`Bernoulli(0.6667)`) then negative-binomial drawn, independently per
  sample, so each of the 20 samples is a different random mixture with about a third of
  the panel genuinely absent;
- the **sweep** is `examples/sr_amplicon_param_sweep`'s — a `sr_sweep.grid:` expanded
  into the samplesheet's `sr_settings:` list, two phases so the reads are mapped once.

What is new is the database: `database.path:` is a **pre-built** reference set, so
nothing is built during the run. It is a **label space**, not a set of genomes: a SILVA
accession is one rRNA sequence, so every run infers in V4-group space
(`--infer_space v4_group`) and is scored **per genus**, through the SILVA taxonomy.

## Why this community rather than a sub-species sweep

`examples/sr_amplicon_param_sweep` sweeps the split of one near-identical strain pair —
the right question for a reference set that holds only the community. Against SILVA the
harder question is the opposite one: every genus that is *not* in the community is
absent from every sample, and a setting is only usable if it says so. The NB
community's own absent genomes (about 7 of 20 per sample, differing per sample) give the
same measurement inside the panel, where the truth is unambiguous.

## The database (`path:`, no build step)

A pre-built database is a directory laid out the way the pipeline publishes its own
(`<outdir>/databases/<name>/`). For `sr_amplicon` that is two files:

```
<database.path>/<database.name>_ssu.sr_refs.fasta      >{accession}|0|{accession}
<database.path>/<database.name>_ssu.sr_refs.tax        {header}<TAB>{SILVA lineage}
```

Build them once, out of band, with superresolution-amplicon's
`bin/build_mapseq_database.py --silva-fasta SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz`
(it converts U to T, drops the organism name and pads lineages), then rename them to the
pattern above. `BUILD_DATABASES` resolves both by glob. The `.tax` is optional to the
pipeline but not to this example: it reaches every nested run as `--taxonomy`, which
fills the `lca` column the scores are computed from.

## Scoring (per genus)

`scripts/score_sweep.py` reads each grid point's `profiling/sr/<id>.<point>.inferred_v4_groups.csv`:

- **truth**: `truth.tsv` per genome, rolled up to genus through `panel[].taxonomy` (the
  genome's SILVA lineage down to genus, in `config.yaml`);
- **prediction**: each V4 group's `inferred_mean` goes to the genus of its `lca` (the
  longest common lineage of the group's SILVA members). A group whose `lca` stops above
  genus, because its V4 sequence is shared across genera, goes to `unresolved`.

| Column | Meaning |
|---|---|
| `tv` | total variation over genera plus `unresolved`, both sides renormalised |
| `max_abs_error`, `worst_genus` | the largest single-genus error |
| `mass_absent` | predicted share on genera absent from that community |
| `unresolved` | predicted share SILVA cannot place at genus |

`mapseq_only` is MAPseq's own label share (`observed_rel_abundance`), no inference.
Genus is SILVA's 6th rank, which holds for bacteria and archaea only.

The *B. uniformis* strain question of the other examples cannot be asked here: SILVA
has no species rank, and NR99 collapses near-identical sequences.

To run the same sweep against a small collection the pipeline builds itself (a cheap
control, on the same code path), delete `path:` and set `sequences_from_panel: true`.
The panel's `taxonomy:` gives that collection a `.tax` too, so it scores the same way.

## Fill in before running

Everything is driven by **`config.yaml`** — no paths or metadata are hard-coded in the
Python scripts. Edit:

- `train.fastq_1` / `train.fastq_2` — the real R1/R2 the error model is trained from.
- `database.path` — the pre-built reference set directory, and `database.name` to match
  the `<name>_ssu.sr_refs.{fasta,tax}` inside it.
- `panel` — one entry per community genome: `id`, `species`, `ssu` (a full-length 16S
  FASTA, which the V4 amplicons are cut from) and `taxonomy` (its SILVA lineage down to
  genus). The shipped lineages are SILVA 138.2's.
- `sampling` — `n_samples`, `seed`, `presence`, `negative_binomial: {mean, dispersion}`.
- `sr_sweep.grid` — the sweep itself, each axis a list of named knob maps, the grid their
  cartesian product. Valid knobs are the pipeline's `sr_settings` keys.
- `primers` (under `generation_modes`) — superresolution-amplicon cuts its own reference
  amplicons with the same pair, so a mismatch means not one read hits a reference.

The scripts need **PyYAML + numpy**. `python scripts/silva_sweep.py --selfcheck` exercises
the database block, the genus cut, the depth parsing and the grid expansion, and
`python scripts/score_sweep.py --selfcheck` the scoring, without touching the pipeline.

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
3. `scripts/score_sweep.py` → `results/sr_amplicon_silva_sweep/silva_sweep_scores.csv`.

## What it costs

The shipped config is 20 samples × 2 depths = 40 benchmark dirs × 12 grid points = **480
profiles**, off **3** mis-mapping matrices — `kmer1` and `kmer1_latent` differ only in
whether the distance decay is fitted per sample (`infer_distance_decay`), so they share
the matrix and that comparison costs inference runs only. Three things keep that affordable, and all
three are about the size of the reference set rather than the size of the grid:

- **The reads are mapped once.** mapseq against a SILVA-scale set is by far the most
  expensive stage and no swept knob changes it. Phase 1 maps; phase 2 hands the
  classification back through `mseq:` and every grid point skips straight to the parts
  that differ. A dir with no phase-1 mseq is not an error — the nested run maps it
  itself, and the generator warns.
- **The matrix is shared where it can be.** Grid points agreeing on every matrix knob
  (`mismapping_method`, `align_backend`, `align_tau`, `align_distance_decay`,
  `matrix_args`) build one matrix between them and split only at the inference run —
  which is why the presence prior (including no gate, which superresolution-amplicon
  recommends at database scale in V4-group space), and whether the decay is fitted, are
  swept here rather than in separate runs.
- **`mismapping_method: simulate` is left out of the grid.** It simulates and maps
  `sim_n_per_ref` reads for *every reference in the set*: fine for a 20-genome
  collection, not for SILVA. The grid sweeps the alignment backends instead
  (`exact-hash`, and `kmer` at τ=1 and τ=2). `benchmark.config` pins phase 1 to the
  `exact` point for the same reason — phase 1 has no `sr_settings` to fan out over, so
  it would otherwise use the nested pipeline's `simulate` default.

Phase 1's matrix is built separately from phase 2's even when the knobs match: the
reference-set key mixes in the matrix mode only when `sr_settings` is in play, so the two
runs do not share a nested cache. That is one extra `exact` matrix, not one extra
mapping pass.

## Output

```
results/sr_amplicon_silva_sweep/
  silva_sweep_scores.csv                               <- step 3
  S01.amplicon_16s.515-YF-806BR/                       <- one per community
    truth.tsv
    S01.amplicon_16s.515-YF-806BR.exact.p01.sr_profile.tsv
    S01.amplicon_16s.515-YF-806BR.kmer1.p01.sr_profile.tsv
    ...
    profiling/sr/S01.amplicon_16s.515-YF-806BR.obs.mseq.gz   <- phase 1; reused by phase 2
    profiling/sr/<id>.<point>.inferred_v4_groups.csv         <- what is scored
    subsample_100000/...                                     <- same, per depth
  S02.amplicon_16s.515-YF-806BR/ ...
  mismapping/
    silva_138_2_ssu_nr99_sr_amplicon_515-YF-806BR-<digest>/   <- one per matrix knob combo
      mismapping_matrix.npz
      mismapping_provenance.json                       <- how that matrix was built
```

No `databases/` directory: nothing was built.

## Reporting

`reports/` renders the settings comparison from that tree:

```bash
cd ../../reports
task sr-settings RUN_DIR=<run> PIPELINE_DIR=results/sr_amplicon_silva_sweep \
    RUN_LABEL="SILVA sr_amplicon parameter sweep"
```

That report scores the genome-level `sr_profile.tsv`, which against SILVA names
accessions rather than panel genomes; read `silva_sweep_scores.csv` for accuracy and the
report for the mis-mapping matrices. Its sub-species sections auto-detect genomes whose target abundance varies across the
sweep; there are none here (every genome varies, independently), so pass `SUBSPECIES=`
a pair of ids to focus them, or read the per-setting scores and the mis-mapping matrices,
which need no sweep pair.

## Notes

- Only the **matrix** knobs split the reference set. Adding an `infer:` point costs one
  extra inference run per benchmark dir and no extra matrix.
- A knob set both in `benchmark.config` (as `sr_amplicon_<knob>`) and in a grid point
  takes the grid point's value — which is why phase 2 ignores the pins phase 1 needs.
  `inference_args` (`--infer_space v4_group`) is the exception that is kept, because no
  grid point sets it; a point that does must repeat the flag.
- Re-running phase 2 alone is fine and cheap: `python generate_sweep_samplesheet.py
  <outdir>` then `--step profile`. Widen the grid and the nested runs `-resume` the work
  the existing points already did.
- `sr_settings` is emitted top-level here (the default for every row). It also works per
  row, if different samples should be swept differently.
