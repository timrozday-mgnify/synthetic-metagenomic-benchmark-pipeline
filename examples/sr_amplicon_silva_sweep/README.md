# superresolution-amplicon parameter sweep against SILVA

Asks whether a **large external reference set is fit for purpose**: the synthetic
community is 20 known genomes, the reference set the reads are profiled against is
SILVA SSU Ref NR99 (or any other pre-built SSU set with a taxonomy), and the
superresolution-amplicon parameters are swept over the same reads so the settings can be
compared against one ground truth.

Two comparison axes sit on top of that, and everything is scored in one genus space:

- **four profiling arms** over the same reads — `silva` (the sweep proper), `custom` (a
  collection built from the community itself), `silva_panel` (SILVA reinterpreted through
  that collection) and `aap` (**no superresolution at all**: the EBI
  amplicon-analysis-pipeline's own MAPseq labels). The last is the baseline the other
  three have to beat to be worth running;
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
**pre-built** reference set, so nothing is built during the run. It is a **label space**, not a set of genomes: a SILVA
accession is one rRNA sequence, so every run infers in V4-group space
(`--infer_space v4_group`) and is scored **per genus**, through the SILVA taxonomy.

## Why this community rather than a sub-species sweep

`examples/sr_amplicon_param_sweep` sweeps the split of one near-identical strain pair —
the right question for a reference set that holds only the community. Against SILVA the
harder question is the opposite one: every genus that is *not* in the community is
absent from every sample, and a setting is only usable if it says so. The NB
community's own absent genomes (about 7 of 20 per sample, differing per sample) give the
same measurement inside the panel, where the truth is unambiguous.

## The four profiling arms

All four turn the same reads into abundances, and `scripts/score_sweep.py` scores them in
the same genus space — the only space they share, since two report panel genomes, one
reports SILVA V4 groups and one reports SILVA lineages.

| Arm | Database | What it is | Scored from |
|---|---|---|---|
| `silva` | `database` (SILVA) | superresolution in V4-group space, **fanned over the whole grid** | `inferred_v4_groups.csv`, each group's `lca` |
| `custom` | `custom_database` | superresolution against a collection the pipeline builds from `panel:` itself — the upper bound a reference set that *is* the community gives | `inferred_composition.csv`, panel genomes |
| `silva_panel` | `database` (SILVA) | the same SILVA labels reinterpreted through that panel collection (an `sr_settings` `panel:`) — what knowing the community buys when the reads can only be mapped to SILVA | `inferred_composition.csv`, plus a `background` bucket |
| `aap` | `aap_database` | no superresolution: the amplicon-analysis-pipeline classifies with MAPseq against a pre-built mapseq SILVA DB | its krona table, a count per lineage |

Only `silva` runs the grid. The other three run at **one** point, `arm_point:` in
`config.yaml` (shipped as `exact.nogate`), so the arm comparison is not multiplied by the
sweep. `custom` maps its own reads — cheap against 20 references, and phase 1's SILVA
classification is not its classification — while `silva` and `silva_panel` reuse it.

`custom` and `silva_panel` name their settings `<arm>.<point>`, because a profile file is
`<id>.<setting>.sr_profile.tsv` and all four arms publish into the one benchmark dir.

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

A pre-built database is a directory laid out the way the pipeline publishes its own
(`<outdir>/databases/<name>/`). For `sr_amplicon` that is two files:

```
<database.path>/<database.name>_ssu.sr_refs.fasta      >{accession}|0|{accession}
<database.path>/<database.name>_ssu.sr_refs.tax        {header}<TAB>{SILVA lineage}
```

Build them once, out of band, with superresolution-amplicon's
`bin/build_mapseq_database.py --silva-fasta SILVA_<ver>_SSURef[_NR99]_tax_silva.fasta.gz`
(it converts U to T, drops the organism name and pads lineages). Any release works, NR99
or the full Ref, and the names above are only a convention: `BUILD_DATABASES` falls back
to a lone `*.{fasta,fa,fna}` + `*.tax` in the directory, so there is nothing to rename.
Only the FASTA and its `.tax` may live there. The `.tax` is optional to the
pipeline but not to this example: it reaches every nested run as `--taxonomy`, which
fills the `lca` column the scores are computed from.

## Scoring (per genus)

`scripts/score_sweep.py` writes one row per benchmark dir x method, with an `arm` and an
`error_model` column. **Truth** is always `truth.tsv` per genome, rolled up to genus
through `panel[].taxonomy` (the genome's SILVA lineage down to genus, in `config.yaml`).
**Predictions** come from each arm's own output:

- `silva`: each V4 group's `inferred_mean` goes to the genus of its `lca` (the longest
  common lineage of the group's SILVA members). A group whose `lca` stops above genus,
  because its V4 sequence is shared across genera, goes to `unresolved`.
- `custom` / `silva_panel`: each panel genome's `inferred_mean` goes to its own lineage's
  genus. `silva_panel`'s `background` — reads the panel cannot explain — goes to
  `unresolved`, where the truth has no mass, so it is charged as error.
- `aap`: the krona table's count per lineage, summed per genus.

| Column | Meaning |
|---|---|
| `arm`, `method` | the profiling arm, and the grid point (or `aap` / `mapseq_only`) within it |
| `error_model` | the arm the **reads** came from: `trained`, `flat` or `naive` |
| `tv` | total variation over genera plus `unresolved`, both sides renormalised |
| `max_abs_error`, `worst_genus` | the largest single-genus error |
| `mass_absent` | predicted share on genera absent from that community |
| `unresolved` | predicted share that arm cannot place at genus |

Rows are only comparable **within one `error_model`**: different arms there are different
reads. `mapseq_only` is MAPseq's own label share against SILVA (`observed_rel_abundance`),
no inference. Genus is SILVA's 6th rank, which holds for bacteria and archaea only.

The *B. uniformis* strain question of the other examples cannot be asked here: SILVA
has no species rank, and NR99 collapses near-identical sequences.

The `custom` arm is that control, run alongside rather than instead: a collection the
pipeline builds from the panel, on the same code path. To run the **whole grid** against
it instead of against SILVA, delete `database.path:` and set `sequences_from_panel: true`
— the panel's `taxonomy:` gives that collection a `.tax` too, so it scores the same way.

## Fill in before running

Everything is driven by **`config.yaml`** — no paths or metadata are hard-coded in the
Python scripts. Edit:

- `train.fastq_1` / `train.fastq_2` — the real R1/R2 the error model is trained from.
  The `naive` arm never reads them; the other two do.
- `database.path` — the pre-built reference set directory, and `database.name` to match
  the `<name>_ssu.sr_refs.{fasta,tax}` inside it.
- `aap_database` — the `aap` arm's **mapseq** database (a different thing from the
  superresolution reference set): a directory holding the four MAPseq files, plus the two
  Rfam paths the pipeline requires of any `aap` collection. MAPseq's own distributed
  SILVA database works as-is.
- `custom_database.name` — nothing to fill in: it is built from `panel[]` and publishes
  to `<outdir>/databases/<name>/`, ready to be another run's `path:`.
- `error_models` — the three read-generation arms. Cut it to one entry for one set of reads.
- `arm_point` — which grid point the non-sweep arms run at. Must name a point of the grid.
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
   per community), then `--step all`. Trains the error models (one per arm that needs
   one), generates the amplicon reads **once per error-model arm** at every
   `reads.subsample` depth, and profiles each depth **once** against the pre-built
   reference set. This is the run that publishes each depth's
   `profiling/sr/<id>.obs.mseq.gz`.
2. `generate_sweep_samplesheet.py` → `sweep_samplesheet.yaml`, then `--step profile`.
   Four rows per benchmark dir, one per profiling arm, each with its own `database:`,
   `profilers:` and `sr_settings:`. The `silva` and `silva_panel` rows carry that dir's
   `mseq:`; `custom` and `aap` map their own reads.
3. `scripts/score_sweep.py` → `results/sr_amplicon_silva_sweep/silva_sweep_scores.csv`.

## What it costs

**The shipped config is large.** 20 communities × 3 error-model arms × 2 depths = **120
benchmark dirs**; each runs 12 `silva` grid points + 1 `custom` + 1 `silva_panel` + 1
`aap` = **1800 profiles**. Three knobs cut it, in descending order of effect, and each is
one edit: drop `error_models:` to one entry (÷3, and it also cuts read generation and
SILVA mapping by the same factor — that is the expensive third), shrink `sr_sweep.grid`
(the `silva` arm only), or drop a `reads.subsample` depth. The arm comparison alone —
one error model, one grid point — is 40 dirs × 4 arms.

The error-model arms are the one axis that cannot share work: different reads means a
separate generation *and* a separate SILVA mapping pass per arm. Each also costs its own
`skiver dump` + train, even though `flat` and `trained` dump the same FASTQs; that is one
extra dump, not one extra sweep.

Within one error-model arm, the sweep proper is 40 benchmark dirs × 12 grid points = 480
profiles, off **3** mis-mapping matrices — `kmer1` and `kmer1_latent` differ only in
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
  error_models/sc2200627.{trained,flat}/               <- one per arm that trains
  databases/community_panel/                           <- the `custom` arm's collection
  S01.trained.amplicon_16s.515-YF-806BR/               <- one per community x error model
    truth.tsv
    <id>.exact.p01.sr_profile.tsv                      <- `silva` arm, one per grid point
    <id>.kmer1.p01.sr_profile.tsv
    ...
    <id>.custom.exact.nogate.sr_profile.tsv            <- `custom` arm
    <id>.silva_panel.exact.nogate.sr_profile.tsv       <- `silva_panel` arm
    profiling/sr/<id>.obs.mseq.gz                      <- phase 1; reused by phase 2
    profiling/sr/<id>.<point>.inferred_v4_groups.csv   <- what the `silva` arm is scored from
    profiling/sr/<id>.<arm>.<point>.inferred_composition.csv   <- ... and the genome arms
    profiling/aap/<id>/taxonomy-summary/               <- `aap` arm; krona table is scored
    subsample_100000/...                               <- same, per depth
  S01.flat.amplicon_16s.515-YF-806BR/ ...              <- same community, flat model
  S01.naive.amplicon_16s.515-YF-806BR/ ...
  mismapping/
    silva_138_2_ssu_nr99_sr_amplicon_515-YF-806BR-<digest>/   <- one per matrix knob combo
      mismapping_matrix.npz
      mismapping_provenance.json                       <- how that matrix was built
```

Only `custom_database` is built; SILVA and the mapseq DB are resolved from their `path:`.

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
- `sr_settings` is emitted **per row** here, because the arms sweep differently. A
  top-level list still works when every row wants the same grid.
- The `custom` arm's collection publishes to `databases/community_panel/` in the exact
  layout a `path:` entry expects, so a second run can consume it pre-built.
- **`naive` is not a straw man.** It is what a benchmark that skips the training step
  does, so the honest reading of any sweep result is the gap between its `trained` and
  `naive` rows, not the `trained` row alone.
