# superresolution-amplicon parameter sweep against SILVA

Asks whether a **large external reference set is fit for purpose**: the synthetic
community is 20 known genomes, the reference set the reads are profiled against is
SILVA SSU Ref NR99 (or any other pre-built SSU set with a taxonomy, AAP's included), and
the superresolution-amplicon parameters are swept over the same reads so the settings can
be compared against one ground truth.

Two comparison axes sit on top of that, and everything is scored in one genus space:

- **three profiling arms** over the same reads — `silva_panel` (SILVA's labels
  reinterpreted through the community's genomes; the sweep proper), `custom` (a
  collection built from the community itself) and `aap` (**no superresolution at all**:
  the EBI amplicon-analysis-pipeline's own MAPseq labels). The last is the baseline the
  other two have to beat to be worth running;
- **three error-model arms** — the reads are *generated* three times, from the trained
  model, from a context-free one and from skiver's untrained preset, so any conclusion
  can be checked against how much of it the error model is responsible for.

It is two existing examples spliced together, and the scripts import them rather than
copy them:

- the **communities** are `examples/abundance_nb_sample`'s — every panel genome
  presence-gated (`Bernoulli(0.6667)`) then negative-binomial drawn, independently per
  sample, so each of the 20 samples is a different random mixture with about a third of
  the panel genuinely absent;
- the **sweep** is `examples/sr_amplicon_param_sweep`'s — a `sr_sweep.grid:` expanded
  into the samplesheet's `sr_settings:` list, two phases so the reads are mapped once.

What is new is the two comparison axes, and the database: `database.path:` is a
**pre-built** reference set, so nothing is built during the run. It is a **label space**,
not a set of genomes: a SILVA accession is one rRNA sequence, so every run against it
infers over the community's genomes (a `panel:`), and everything is scored **per genus**,
through the SILVA taxonomy. `examples/sr_amplicon_panel_sweep` asks the neighbouring
question over the 20HM panel, with the kernel's error model swept instead of its mode.

## Why this community rather than a sub-species sweep

`examples/sr_amplicon_param_sweep` sweeps the split of one near-identical strain pair —
the right question for a reference set that holds only the community. Against SILVA the
harder question is the opposite one: every genus that is *not* in the community is
absent from every sample, and a setting is only usable if it says so. The NB
community's own absent genomes (about 7 of 20 per sample, differing per sample) give the
same measurement inside the panel, where the truth is unambiguous.

## The three profiling arms

All three turn the same reads into abundances, and `scripts/score_sweep.py` scores them in
the same genus space — the only space they share, since two report panel genomes and one
reports SILVA lineages.

| Arm | Database | What it is | Scored from |
|---|---|---|---|
| `silva_panel` | `database` (SILVA) | SILVA's labels reinterpreted through the panel collection (an `sr_settings` `panel:`), **fanned over the whole grid** — what knowing the community buys when the reads can only be mapped to SILVA | `inferred_composition.csv`, plus a `background` bucket |
| `custom` | `custom_database` | superresolution against that panel collection itself — the upper bound a reference set that *is* the community gives | `inferred_composition.csv`, panel genomes |
| `aap` | `aap_database` | no superresolution: the amplicon-analysis-pipeline classifies with MAPseq against a pre-built mapseq SILVA DB | its krona table, a count per lineage |

Only `silva_panel` runs the grid. `custom` runs at **one** point, `arm_point:` in
`config.yaml` (shipped as `exact.nogate`), so the arm comparison is not multiplied by the
sweep. `custom` maps its own reads — cheap against 20 references, and phase 1's SILVA
classification is not its classification — while `silva_panel` reuses it.

Both superresolution arms name their settings `<arm>.<point>`, because a profile file is
`<id>.<setting>.sr_profile.tsv` and all three arms publish into the one benchmark dir.
There is no panel-less SILVA arm: without a `panel:` the nested run infers over every
SILVA V4 group, which names sequences rather than genomes and builds a kernel over the
whole database.

## The three error-model arms

These change the **reads**, not the profiling, so each is a separate set of samples with
its own `train_id` (training is deduped and keyed by it) and its own benchmark dirs
(`S01.trained.amplicon_16s.…`, `S01.flat.…`, `S01.naive.…`). The community draw is
shared, so the arms differ only in how the reads were damaged.

| Arm | Samplesheet | Model |
|---|---|---|
| `trained` | — | the pipeline default: `error_model_candidates` fitted from `train.fastq_*`, min-AIC kept |
| `flat` | `error_model_components: BaseContext(1)` | fitted from the same reads, but blind to sequence context — the shortest context skiver accepts |
| `naive` | `error_model: illumina` | no training at all: skiver's bundled `hq-illumina` preset |

Both columns are pipeline features this example uses, not example-local glue:
`error_model` skips training for that row's `train_id` entirely and passes
`--error-model` to genome-blender; `error_model_components` overrides the fitted
component string per `train_id`. Cut the `error_models:` block to one entry to go back to
one set of reads.

## The database (`path:`, no build step)

A pre-built database is a directory holding the reference FASTA, its `.tax` and,
optionally, the FASTA's `.mscluster`:

```
<database.path>/<database.name>_ssu.sr_refs.fasta      >{accession}|0|{accession}
<database.path>/<database.name>_ssu.sr_refs.tax        {header}<TAB>{SILVA lineage}
<database.path>/<database.name>_ssu.sr_refs.fasta.mscluster   optional
```

The FASTA **is** the MAPseq database every nested run maps against: the phase-1 reads
and the `silva_panel` kernels' panel sources alike, so the labels mean the same thing
everywhere. A `.mscluster` beside the FASTA's real path is used; without one the nested
run clusters it once and caches it (`--amplicon_cache`). Nothing is built during the run.

The amplicon-analysis-pipeline's own database directory works as shipped:
`BUILD_DATABASES` falls back to a lone `*.{fasta,fa,fna}` + `*.tax` / `*-tax.txt`, so
`SILVA-SSU.fasta`, `SILVA-SSU-tax.txt` and `SILVA-SSU.fasta.mscluster` need no renaming.
The nested run folds its RNA alphabet to DNA and reads its `sk__…;s__Genus_species`
lineages. Otherwise build one with superresolution-amplicon (`$SRA` its checkout):

```bash
python $SRA/bin/build_mapseq_database.py \
    --silva-fasta SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \
    --output-prefix db/silva_138_2_ssu_nr99_ssu.sr_refs
```

`build_mapseq_database.py` converts U to T, drops the organism name and pads lineages.
Any release works, NR99 or the full Ref. The `.tax` is optional to the pipeline but not to
this example: every nested run gets it as `--taxonomy`.

## Scoring (per genus)

`scripts/score_sweep.py` writes one row per benchmark dir x method, with an `arm` and an
`error_model` column. **Truth** is always `truth.tsv` per genome, rolled up to genus
through `panel[].taxonomy` (the genome's SILVA lineage down to genus, in `config.yaml`).
**Predictions** come from each arm's own output:

- `silva_panel` / `custom`: each panel genome's `inferred_mean` goes to its own lineage's
  genus. `silva_panel`'s `background` — reads the panel cannot explain — goes to
  `unresolved`, where the truth has no mass, so it is charged as error.
- `aap`: the krona table's count per lineage, summed per genus.

| Column | Meaning |
|---|---|
| `arm`, `method` | the profiling arm, and the grid point (or `aap`) within it |
| `error_model` | the arm the **reads** came from: `trained`, `flat` or `naive` |
| `tv` | total variation over genera plus `unresolved`, both sides renormalised |
| `max_abs_error`, `worst_genus` | the largest single-genus error |
| `mass_absent` | predicted share on genera absent from that community |
| `unresolved` | predicted share that arm cannot place at genus |

Rows are only comparable **within one `error_model`**: different arms there are different
reads. Genus is SILVA's 6th rank, which holds for bacteria and archaea only.

The `custom` arm is the control, run alongside: a collection the pipeline builds from the
panel, on the same code path. To point `database:` itself at it instead of at SILVA,
delete `database.path:` and set `sequences_from_panel: true` — the panel's `taxonomy:`
gives that collection a `.tax` too, so it scores the same way.

## Fill in before running

Everything is driven by **`config.yaml`** — no paths or metadata are hard-coded in the
Python scripts. Edit:

- `train.fastq_1` / `train.fastq_2` — the real R1/R2 the error model is trained from.
  The `naive` arm never reads them; the other two do.
- `database.path` — the pre-built reference set directory (AAP's `SILVA-SSU/138.1` works
  as is), and `database.name` to name it (see
  [The database](#the-database-path-no-build-step)).
- `aap_database` — the `aap` arm's **mapseq** database (a different thing from the
  superresolution reference set): a directory holding the four MAPseq files, plus the two
  Rfam paths the pipeline requires of any `aap` collection. MAPseq's own distributed
  SILVA database works as-is.
- `custom_database.name` — nothing to fill in: it is built from `panel[]` and publishes
  to `<outdir>/databases/<name>/`, ready to be another run's `path:`.
- `error_models` — the three read-generation arms. Cut it to one entry for one set of reads.
- `arm_point` — which grid point `custom` and phase 1's mapping run use. Must name a point of the grid.
- `panel` — one entry per community genome: `id`, `species`, `ssu` (a full-length 16S
  FASTA, which the V4 amplicons are cut from) and `taxonomy` (its SILVA lineage down to
  genus). The shipped lineages are SILVA 138.2's.
- `sampling` — `n_samples`, `seed`, `presence`, `negative_binomial: {mean, dispersion}`.
- `sr_sweep.grid` — the sweep itself (the `silva_panel` arm), each axis a list of named
  knob maps, the grid their cartesian product. Valid knobs are the pipeline's
  `sr_settings` keys.
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
   per community), then `--step all`. Trains the error models (one per arm that needs
   one), generates the amplicon reads **once per error-model arm** at every
   `reads.subsample` depth, and profiles each depth **once** against the pre-built
   reference set under a `map` setting (`arm_point`'s knobs, over the panel). This is the
   run that publishes each depth's `profiling/sr/<id>.map.obs.mseq.gz`.
2. `generate_sweep_samplesheet.py` → `sweep_samplesheet.yaml`, then `--step profile`.
   Three rows per benchmark dir, one per profiling arm, each with its own `database:`,
   `profilers:` and `sr_settings:`. The `silva_panel` row carries that dir's `mseq:`;
   `custom` and `aap` map their own reads.
3. `scripts/score_sweep.py` → `results/sr_amplicon_silva_sweep/silva_sweep_scores.csv`.

## What it costs

**The shipped config is large.** 20 communities × 3 error-model arms × 2 depths = **120
benchmark dirs**; each runs 12 `silva_panel` grid points + 1 `custom` + 1 `aap` =
**1680 profiles**. Three knobs cut it, in descending order of effect, and each is one
edit: drop `error_models:` to one entry (÷3, and it also cuts read generation and SILVA
mapping by the same factor — that is the expensive third), shrink `sr_sweep.grid`, or
drop a `reads.subsample` depth.

The error-model arms are the one axis that cannot share work: different reads means a
separate generation *and* a separate SILVA mapping pass per arm. Each also costs its own
`skiver dump` + train, even though `flat` and `trained` dump the same FASTQs; that is one
extra dump, not one extra sweep.

Within one error-model arm, the sweep proper is 40 benchmark dirs × 12 grid points = 480
profiles, off **3** panel kernels — `kmer1` and `kmer1_latent` differ only in whether the
distance decay is fitted per sample (`infer_distance_decay`), so they share the kernel and
that comparison costs inference runs only. What keeps that affordable:

- **The reads are mapped once.** mapseq against a SILVA-scale set is by far the most
  expensive stage and no swept knob changes it. Phase 1 maps; phase 2 hands the
  classification back through `mseq:` and every grid point skips straight to the parts
  that differ. A dir with no phase-1 mseq is not an error — the nested run maps it
  itself, and the generator warns.
- **The kernel is shared where it can be,** and it is small. Grid points agreeing on every
  kernel knob (`mismapping_method`, `align_tau`, `align_distance_decay`, `matrix_args`,
  `panel`) build one kernel between them and split only at the inference run. A kernel
  runs from the panel's ~20 V4 sources to SILVA's labels, not over SILVA.
- **`mismapping_method: simulate` is commented out of the grid.** Over this panel it
  would simulate and map 5,000 reads per source against SILVA per kernel: affordable, but
  the grid sweeps the alignment modes (τ=0, 1, 2) instead.

## Output

```
results/sr_amplicon_silva_sweep/
  silva_sweep_scores.csv                               <- step 3
  error_models/sc2200627.{trained,flat}/               <- one per arm that trains
  databases/community_panel/                           <- the `custom` arm's collection
  S01.trained.amplicon_16s.515-YF-806BR/               <- one per community x error model
    truth.tsv
    <id>.silva_panel.exact.p01.sr_profile.tsv          <- `silva_panel` arm, one per grid point
    <id>.silva_panel.exact.p01.sr_background.tsv       <- ... and its unexplained share
    ...
    <id>.custom.exact.nogate.sr_profile.tsv            <- `custom` arm
    profiling/sr/<id>.map.obs.mseq.gz                  <- phase 1; reused by phase 2
    profiling/sr/<id>.<arm>.<point>.inferred_composition.csv   <- what both arms are scored from
    profiling/aap/<id>/taxonomy-summary/               <- `aap` arm; krona table is scored
    subsample_100000/...                               <- same, per depth
  S01.flat.amplicon_16s.515-YF-806BR/ ...              <- same community, flat model
  S01.naive.amplicon_16s.515-YF-806BR/ ...
  mismapping/
    <reference set>-<digest>/                          <- one per kernel knob combo
      panel_kernel/                                    <- reusable as --panel_kernel
      mismapping_provenance.json                       <- how that kernel was built
```

Only `custom_database` is built; SILVA and the `aap` mapseq DB are resolved from their
`path:`.

## Reporting

`reports/` renders the settings comparison from that tree:

```bash
cd ../../reports
task sr-settings RUN_DIR=<run> PIPELINE_DIR=results/sr_amplicon_silva_sweep \
    RUN_LABEL="SILVA sr_amplicon parameter sweep"
```

That report scores the genome-level `sr_profile.tsv` of every arm; read
`silva_sweep_scores.csv` for the per-genus accuracy the arms are compared on. Its sub-species sections auto-detect genomes whose target abundance varies across the
sweep; there are none here (every genome varies, independently), so pass `SUBSPECIES=`
a pair of ids to focus them, or read the per-setting scores and the mis-mapping matrices,
which need no sweep pair.

## Notes

- Only the **kernel** knobs split the reference set. Adding an `infer:` point costs one
  extra inference run per benchmark dir and no extra kernel.
- `benchmark.config` sets no kernel knob: phase 1's `map` setting and phase 2's grid are
  row-level `sr_settings:`, and a knob pinned there would leak into every entry that
  leaves it out.
- Re-running phase 2 alone is fine and cheap: `python generate_sweep_samplesheet.py
  <outdir>` then `--step profile`. Widen the grid and the nested runs `-resume` the work
  the existing points already did.
- `sr_settings` is emitted **per row** here, because the arms sweep differently. A
  top-level list still works when every row wants the same grid.
- The `custom` arm's collection publishes to `databases/community_panel/` in the layout a
  `path:` entry expects.
- **`naive` is not a straw man.** It is what a benchmark that skips the training step
  does, so the honest reading of any sweep result is the gap between its `trained` and
  `naive` rows, not the `trained` row alone.
