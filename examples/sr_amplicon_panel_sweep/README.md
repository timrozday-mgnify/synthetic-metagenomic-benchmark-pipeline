# SILVA panel reinterpretation vs a custom database

Asks how much is lost when amplicon reads **cannot be re-mapped against the community**.
The reads were classified against SILVA, and the community is known to come from a panel:
the 20HM mock plus a second *B. uniformis* strain (23 genomes). There are two ways to get
panel-genome abundances:

- **`custom`**: map the reads against a database built from the panel itself, and infer
  through a square mis-mapping matrix over those 23 references.
- **`generic_panel`**: keep the SILVA labels and reinterpret them through the panel.
  superresolution-amplicon measures a rectangular kernel (panel V4 amplicons x SILVA
  labels) by simulating reads from the panel and mapping them against SILVA. It then
  infers the 23 genomes plus a `background` bucket for reads the panel cannot explain.

Each arm runs under its **own parameter grid** over the same reads, and everything is
scored against one ground truth.

It is three existing examples spliced together, and the scripts import them rather than
copy them:

- the **communities** come from `examples/abundance_nb_sample`: every panel genome is
  presence-gated (`Bernoulli(0.6667)`) and then negative-binomial drawn, so each of the 20
  samples is a different mixture with about a third of the panel absent;
- the **grids** come from `examples/sr_amplicon_param_sweep`: `grid:` axes expanded into
  `sr_settings:`;
- **map SILVA once, sweep afterwards** comes from `examples/sr_amplicon_gtdb_sweep`:
  phase 1's SILVA classification is handed back as `mseq:`.

## The grids

Both live under `sr_sweep:` in `config.yaml`, one block per arm. Each axis is a list of
named knob maps, and the grid is their cartesian product. A point is named
`<arm>.<axis names joined by '.'>`, and that name is part of its output filename. The
shipped grids:

| Arm | Axis | Points |
|---|---|---|
| `custom` | `matrix` | `simulate` (flat error), `kmer1_latent` (align kmer τ=1, decay fitted per sample) |
| | `prior` | `gate` (presence prior 0.01), `nogate` |
| `generic_panel` | `kernel` | `trained` (simulates with the model the reads were made with), `flat` |
| | `prior` | `gate`, `nogate`, `horseshoe` |
| | `steps` | `s3k` (`s10k` commented out) |

Add an axis or a point by editing `config.yaml`. Knobs are the pipeline's `sr_settings`
keys. Any other superresolution-amplicon param goes through `matrix_args` or
`inference_args` (`--sim_error_model`, `--infer_horseshoe`, `--infer_steps`,
`--infer_alpha`, `--panel_sim_n_per_ref`, ...). When two axes of one point both set one
of those two keys, their flags are joined rather than overwritten, which is what lets
`prior` and `steps` both add inference flags.

Every `generic_panel` point gets `panel: <database.name>` from the generator. There is no
unreinterpreted SILVA arm: its profiles would name SILVA sequences, which cannot be scored
against panel ids. SILVA is a label space here; inference over its own references is
sequence-space inference with no biological reading.

## What the pipeline gives this example

Two samplesheet features exist for it:

- **`panel:`** (an `sr_settings` key) names a `databases:` collection to reinterpret the
  row's own database labels with. It becomes superresolution-amplicon's
  `--panel_references`. The collection is built even if no row profiles against it. A
  panel entry builds no shared matrix: the nested run measures its kernel inside the
  inference run, so it receives the matrix flags instead of `--mismapping_matrix`.
  `infer_distance_decay` is refused with it.
- **`sr_error_model:`** (a row key) is a pre-trained skiver `model.pt`, passed to the
  nested run as `error_model:`. The `trained` kernel simulates with it instead of
  training a model per sample.

## Fill in before running

- `config.yaml` → `train.fastq_1` / `train.fastq_2`: the real reads the error model is
  trained from.
- `config.yaml` → `panel[].ssu`: full-length 16S per genome. The shipped names are those of
  the `sr_amp_param_sweep` run's `references/16S/`.
- `config.yaml` → `generic.path`: a directory holding `<generic.name>_ssu.sr_refs.fasta`
  (headers `{accession}|0|{accession}`) and optionally `<generic.name>_ssu.sr_refs.tax`,
  built once and out of band with superresolution-amplicon's
  `build_mapseq_database.py --silva-fasta` from SILVA Ref NR99. The `.tax` reaches every
  run as `--taxonomy`, which only fills the `lca` column. The sequences must
  still carry the primer sites, because superresolution-amplicon cuts its amplicons by
  in-silico PCR.
- `benchmark.config` → `sr_amplicon_repo`: a superresolution-amplicon checkout (or set
  `sr_revision`) that has panel reinterpretation.

The scripts need **PyYAML + numpy**. `python scripts/panel_sweep.py --selfcheck` and
`python scripts/score_sweep.py --selfcheck` test the grid expansion, samplesheet layout
and scoring without touching the pipeline.

## Run

```bash
./run.sh
```

1. `generate_samplesheet.py` → `samplesheet.yaml` (+ `genomes/sample_NN.amplicon_16s.csv`),
   then `--step all`. This trains the error model, generates the reads at every depth, and
   maps each depth against SILVA **once** under `generic.map_setting` (align exact-hash, cheap
   at SILVA scale). It publishes `profiling/sr/<id>.map.obs.mseq.gz` per dir, and
   `error_models/<train_id>/<train_id>.model.pt`.
2. `generate_sweep_samplesheet.py` → `sweep_samplesheet.yaml`, then `--step profile`. There
   are two rows per benchmark dir. The `custom` row has `database: community_20hm` and the
   custom grid. The `generic_panel` row has `database: silva_138_2_ssu_nr99`, phase 1's `mseq:`,
   `sr_error_model:`, and the panel grid.
3. `scripts/score_sweep.py` → `results/sr_amplicon_panel_sweep/panel_sweep_scores.csv`.

## What it costs

The shipped config is 20 communities x 2 depths = 40 benchmark dirs. Across them run
4 `custom` + 6 `generic_panel` points, which is **400 profiles**. The expensive parts:

- **SILVA read mapping, once**, in phase 1. No phase-2 setting re-maps the reads.
- **2 custom matrices** over 23 references. The two `prior` points share each one.
- **6 panel kernels**, one per `generic_panel` point. The benchmark batches every dir of a
  reference set into one nested run per inference setting, and the kernel is built inside
  that run. So a kernel costs about `panel_sim_n_per_ref` (5,000) x about 30 distinct
  panel amplicons simulated reads mapped against SILVA, once per point, whatever the number of
  dirs. Adding a `generic_panel` point adds a kernel; adding a `custom` inference point does
  not add a matrix.
- The `custom` rows map their own reads (no `mseq:`), which is cheap against 23 references.

## Output

```
results/sr_amplicon_panel_sweep/
  panel_sweep_scores.csv                                   <- step 3
  error_models/sc2200627/sc2200627.model.pt                <- phase 1; sr_error_model
  S01.amplicon_16s.515-YF-806BR/
    S01.amplicon_16s.515-YF-806BR.truth.tsv
    S01.amplicon_16s.515-YF-806BR.custom.simulate.gate.sr_profile.tsv
    S01.amplicon_16s.515-YF-806BR.generic_panel.trained.nogate.s3k.sr_profile.tsv
    ...
    profiling/sr/S01.amplicon_16s.515-YF-806BR.map.obs.mseq.gz   <- phase 1; mseq
    profiling/sr/<id>.<point>.inferred_composition.csv           <- incl. `background`
    subsample_100000/...
  mismapping/community_20hm_sr_amplicon_515-YF-806BR-<digest>/   <- custom matrices only
```

## Scores

`panel_sweep_scores.csv` has one row per benchmark dir x method. The methods are every
grid point, plus `custom.mapseq_only` (MAPseq against the panel, no inference).

| Column | Meaning |
|---|---|
| `tv` | genome total variation, both sides renormalised over panel genomes without `background` |
| `max_abs_error`, `worst_genome` | the largest single-genome error |
| `mass_absent` | predicted share on genomes absent from that community |
| `background` | reinterpretation's out-of-panel share. This community has no out-of-panel genome, so anything here is lost panel mass |
| `strain_split_error` | error in the *B. uniformis* strain-2 fraction within the pair (`score.strain_pair`), blank when the pair is below 1% of the truth |

For reference, the development sweep of the same panel reported median TV 0.009–0.012 for
the custom database and 0.015 for the panel arm (trained kernel, no gate) **against GTDB
r232**, where flat kernels were no better than home labels alone
(superresolution-amplicon `dev/panel_reinterpretation.md`). Rerun against SILVA 138.2 NR99
(`dev/panel_silva_sweep.md`), the best panel arm reached median TV 0.024 against GTDB's
0.013, and a flat kernel with a gate 0.029. Most of the loss is the *B. uniformis* strain
split: the pair shares its home label in SILVA. Both sweeps used equal abundances and a
strain split, so these NB communities are a different test, not a replication.

## Caveats

- **`trained` is an oracle.** It simulates the kernel with the model the reads were
  generated with. `flat` is the no-oracle number.
- **The panel is exactly the community**, so `background` measures lost mass, not
  detected contaminants. To test out-of-panel reads, drop a genome from the collection:
  add a second `database:`-style collection without it, and point a `panel:` at it.
- Phase 1's composition (`<id>.map.sr_profile.tsv`) is an unscored by-product of the
  SILVA mapping. It is inferred over SILVA sequences, so it is not a genome profile.
