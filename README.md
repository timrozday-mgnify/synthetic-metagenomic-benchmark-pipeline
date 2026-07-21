# synthetic-metagenomic-benchmark-pipeline

A [Nextflow](https://www.nextflow.io/) DSL2 pipeline that builds **synthetic
metagenomes with a realistic, data-derived sequencing error profile**, for
benchmarking metagenomic profilers.

For each sample it:

1. **Trains a sequencing error profile** from a *natural* (real, non-synthetic)
   metagenome using [skiver](https://github.com/timrozday-mgnify/skiver) —
   reference-free, from k-mer consensus.
2. **Generates synthetic reads** from a chosen set of reference genomes at chosen
   relative abundances and sequencing depth, applying that error profile via
   [genome-blender](https://github.com/timrozday-mgnify/genome-blender).
3. **Publishes the benchmark ground truth**: the reads, a BAM giving the true
   read→reference mapping, and a ground-truth abundance profile.
4. **Optionally profiles the generated reads** and publishes the predicted
   profile next to the ground truth for easy comparison — [sylph](https://github.com/bluenote-1577/sylph)
   for WGS, the [amplicon-analysis-pipeline](https://github.com/EBI-Metagenomics/amplicon-analysis-pipeline)
   (AAP) for amplicon.

Which phases run is controlled by `--step` (`all` | `generate` | `profile` |
`train`), so you can generate only, profile pre-existing benchmark dirs only, do
both in one run, or just train an error model (`train`) once and reuse it across
later `generate` runs (see [Standalone training](#standalone-training--reusing-a-model)).

## Quick start

```bash
nextflow run main.nf \
    -profile docker \
    --input samplesheet.yaml \
    --outdir results
```

Requires Nextflow (>=25) and a container engine (Docker / Singularity / Apptainer).

On a laptop, cap resources so process requests fit your machine (the defaults
target HPC), e.g. `--max_memory 8.GB --max_cpus 4`.

### HPC (Singularity / Apptainer)

Point the container cache and OCI build scratch at a large shared filesystem (not a
quota'd `$HOME`/`/tmp`), and pre-fetch the images once so the parallel generate tasks
don't race to pull the same image (which corrupts the cache and fails with
`unexpected end of JSON input`):

```bash
export NXF_SINGULARITY_CACHEDIR=/scratch/$USER/nxf_singularity_cache
export APPTAINER_TMPDIR=/scratch/$USER/apptainer_tmp
./bin/prefetch-singularity.sh

nextflow run main.nf -profile singularity --input samplesheet.yaml --outdir results
```

The amplicon path runs the EBI amplicon-analysis-pipeline as a nested `nextflow
run` on the host. It does **not** inherit the outer `-profile`, so tell it which
engine to use via `--aap_configs` (an ordered list of `-c` files). Write an
`aap_singularity.config`:

```groovy
singularity.enabled    = true
singularity.autoMounts = true
```

then add `--aap_configs aap_singularity.config` to the run. `NXF_SINGULARITY_CACHEDIR`
(exported above) is inherited by the nested run. For local development the nested
run has no engine by default — pass `--aap_profile docker` (or an `aap_configs`
file with `docker.enabled = true`).

These two settings can also live in a map-form samplesheet (top-level keys, next to
`databases:`/`samples:`), which overrides the params:

```yaml
aap_configs: [/path/to/aap_singularity.config]   # ordered extra -c files
aap_profile: null                                 # optional -profile (usually unset)
databases: { ... }
samples:   [ ... ]
```

## Samplesheet

A YAML list, one entry per synthetic sample (`tests/samplesheets/test.yaml` is a
working example). It may also be a map with a `samples:` list plus an optional
`databases:` block of named sequence collections the pipeline builds into profiler
databases (see [Named sequence collections](#named-sequence-collections-databases)).

```yaml
- sample: S1
  train_id: natA
  train_fastq_1: natural_R1.fastq.gz
  train_fastq_2: natural_R2.fastq.gz
  train_subsample: 200000   # optional: train on 200k reads instead of the full set
  platform: hq-illumina
  genomes_csv: genomes_S1.csv
  num_reads: 1000000
  mode: shotgun
  paired_end: true
  read_length_mean: 150
  read_length_variance: 10
  profiler: sylph
  database: self
  subsample: [none, 100000, 500000]   # full-depth + two subsampled runs

- sample: S2
  train_id: natA
  train_fastq_1: natural_R1.fastq.gz
  train_fastq_2: natural_R2.fastq.gz
  platform: hq-illumina
  genomes_csv: genomes_S2.csv
  num_reads: 1000000
  mode: shotgun
  paired_end: true
  # no `subsample` -> a single full-depth run
```

| Field | Description |
|--------|-------------|
| `sample` | Unique synthetic-sample ID (names the output dir). |
| `train_id` | Groups error-model training. Rows sharing a `train_id` train the profile **once** and reuse it. |
| `train_fastq_1` | Natural-metagenome reads to learn the error profile from. Not needed when `error_model_dir` is set. |
| `train_fastq_2` | Optional mate (paired training reads). Leave blank for single-end. |
| `error_model_dir` | Optional. Path to an `error_models/<train_id>/` directory produced by an earlier `--step train` (or `generate`/`all`) run. When set, the row's `train_id` is **not** trained — the `*.model.pt` + `*.phred_calibration.json` in that dir are used, and `train_fastq_*`/`train_subsample` are ignored. Set it for **all or none** of a `train_id`'s rows. |
| `train_subsample` | Optional read/pair count to subsample the training reads to before training. `none`/`null`/empty/omitted → train on the full read set. Taken from the first row seen for a given `train_id`, same as `train_fastq_1`/`train_fastq_2`. |
| `platform` | `hq-illumina` \| `lq-illumina` \| `ont` \| `pacbio`. |
| `genomes_csv` | A genome-blender input CSV: `genome_id,fasta_path,abundance`. |
| `num_reads` | Sequencing depth as total reads (read pairs × 2 for paired mode). |
| `mode` | Read structure: `shotgun` \| `amplicon` (also `long`). Default `shotgun`. |
| `paired_end` | Optional. `true` \| `false`. Blank → `params.paired_end` (default `true`). Forced single-end when `mode=long`. |
| `read_length_mean` | Optional. Mean read length. Blank → `params.read_length_mean` (default `150`). |
| `read_length_variance` | Optional. Read-length variance. Blank → `params.read_length_variance` (default `10`). |
| `profiler` | Optional. `sylph` (WGS) or `aap` (amplicon). Blank = generate only, no profiling. |
| `database` | Database to profile against, by name. A name defined in the samplesheet `databases:` block is built (or its prebuilt dir consumed) by the pipeline — works for both `sylph` and `aap`. Otherwise: for `sylph`, a key in `params.sylph_databases`, or `self` to build the DB from this sample's reference genomes; for `aap`, blank (DB comes from `--aap_config`). |
| `subsample` | Optional list of read depths to sweep (absolute read/pair counts). The full draw is generated once, then subsampled to each depth — each gets its own `subsample_<N>/` output dir with its own reads + ground truth + profile. `none`/`null`/empty/omitted → a single full-depth run in `<sample>/`. |
| `chunks` | Optional. Split generation of `num_reads` across N parallel `generate-reads` calls (merged back into one reads-set + BAM before subsampling/ground truth), useful for large `num_reads`. Blank → `params.chunks` (default `1`, no chunking). |

Relative `genomes_csv` / FASTA / FASTQ paths resolve against the pipeline
directory; absolute paths and `scheme://` URLs pass through.

### Named sequence collections (`databases:`)

Instead of pre-building profiler databases out-of-band, define named **sequence
collections** in a top-level `databases:` block and reference them from a sample's
`database` column. The pipeline builds the DB each named collection needs — a sylph
`.syldb` for `profiler: sylph`, a mapseq quartet for `profiler: aap` — once per run,
and profiles the referencing samples against it. Only the DB type actually referenced
by a sample is built (`build_databases` subworkflow → `SYLPH_BUILD_DB` /
`MAPSEQ_PREP → MAPSEQ_CLUSTER → MAPSEQ_OTU`). This replaces the standalone
`examples/subspecies_v4_sweep/scripts/build_profiling_dbs.py` for in-pipeline runs.

```yaml
databases:
  community_v4:                     # build from sequences
    rfam_covariance_model: /dbs/rfam/ribo              # aap only: Rfam .cm dir
    rfam_claninfo: /dbs/rfam/ribo/ribo.clan_info       # aap only: Rfam clan-info file
    sequences:
      - id: bacteroides_fragilis
        genome: references/genomes/CR626927.1.fasta.gz   # for a sylph DB
        ssu: references/16S/CR626927.1_SSU.fasta          # for a mapseq DB
        taxonomy: "Bacteria;Bacteroides;fragilis"          # for a mapseq DB (explicit)
      # ... more sequences ...
  gtdb_r220:                        # OR point at a pre-built DB directory
    path: /dbs/databases/gtdb_r220

samples:
  - sample: S1
    profiler: sylph
    database: community_v4
    # ... the usual sample fields ...
```

Per-entry rules:

- Each named entry is **either** `sequences:` (build) **or** `path:` (a pre-built
  directory), not both.
- `sequences[].genome` is required for a `sylph` collection; `sequences[].ssu` and
  `sequences[].taxonomy` (a `Kingdom;Genus;Species` string) are required for an `aap`
  (mapseq) collection. One collection can serve both if every entry has all three.
  `ssu` must be a pre-extracted full-length 16S FASTA (no barrnap step).
- `rfam_covariance_model` (an Rfam `.cm` directory) and `rfam_claninfo` (the
  `ribo.clan_info` file) are **required for an `aap` collection** — the nested
  amplicon-analysis-pipeline needs them for rRNA detection and aborts with
  `file() ... cannot be empty` if unset. Both are pass-through host paths.
- `path:` points at a directory laid out exactly like this pipeline publishes to
  `<outdir>/databases/<name>/` (`<name>.syldb` and/or `<name>.mapseq.{fasta,tax,otu}`
  + `<name>.mapseq.fasta.mscluster`), so a `databases/<name>/` dir from a prior run is
  directly reusable.

Built DBs are published to `<outdir>/databases/<name>/` for reuse.

### Profile-only samplesheet (`--step profile`)

To profile reads generated by an earlier run (skipping generation), pass
`--step profile` with a samplesheet that points at each benchmark dir:

```yaml
- sample: S1
  profiler: sylph
  benchmark_dir: results/S1
  database: gtdb_r220
- sample: S2
  profiler: aap
  benchmark_dir: results/S2
  database:            # blank -> no configured DB
```

Subsampling is a generate-stage feature; the profile-only step ignores any
`subsample` field.

`benchmark_dir` is a directory containing the reads (`*.fastq.gz`). The predicted
profile is published to `<outdir>/<sample>/`, so point `--outdir` at the benchmark
root to co-locate it with the existing `truth.tsv`. `database=self` is not
available here (no reference genomes) — use a configured database.

### Standalone training / reusing a model (`--step train`)

Training the error model is normally coupled to `generate`, but it can be run on
its own — train once, inspect the report, then point any number of later
`generate` runs at the result instead of retraining each time.

Train-only samplesheet (only the training fields; one model per `train_id`):

```yaml
- train_id: natA
  train_fastq_1: natural_R1.fastq.gz
  train_fastq_2: natural_R2.fastq.gz
  train_subsample: 200000   # optional
  platform: hq-illumina
```

```bash
nextflow run main.nf -profile docker --step train \
    --input train.yaml --outdir results
```

This writes `results/error_models/natA/` (`.model.pt`, `.phred_calibration.json`,
`.context_model_aic.csv`, `.error_model_report.html`, `training_reads/`) and runs
nothing else. Later, a `generate` samplesheet reuses it via `error_model_dir`
(no `train_fastq_*` needed):

```yaml
- sample: S1
  train_id: natA
  error_model_dir: results/error_models/natA
  platform: hq-illumina
  genomes_csv: genomes_S1.csv
  num_reads: 1000000
```

### The genomes CSV

Referenced by each row's `genomes_csv` column — the reference genomes and their
relative abundances (normalised internally):

```csv
genome_id,fasta_path,abundance
genomeA,tests/data/genomeA.fasta,0.7
genomeB,tests/data/genomeB.fasta,0.3
```

## Taxonomic profiling

Set a `profiler` per sample to profile the generated reads and drop the predicted
profile next to `truth.tsv`.

### sylph (WGS)

Uses the nf-core `sylph/profile` module. The database is chosen by the `database`
column:

- **A built collection** — a `database` name defined in the samplesheet
  `databases:` block; the pipeline builds its `.syldb` from the collection's genomes
  once per run (see [Named sequence collections](#named-sequence-collections-databases)).
  A `path:` entry reuses a pre-built DB dir instead.

- **A configured database** — define named databases (mapseq-style) in config and
  reference the key:

  ```groovy
  // nextflow.config / -c my.config
  params.sylph_databases = [
      gtdb_r220: [ syldb: '/dbs/gtdb-r220.syldb', label: 'GTDB-r220' ],
  ]
  ```

  Build a `.syldb` once with `sylph sketch -g genome*.fa -o gtdb-r220`.

- **`self`** — the pipeline builds a `.syldb` from the sample's own reference
  genomes, so sylph's genome-level profile lines up exactly with the ground truth.

Want a fixed "community" database (all of a sample's genomes, but built once and
reused across samples/subsamples, unlike `self`) instead of a production DB?
`examples/subspecies_v4_sweep/scripts/build_profiling_dbs.py` builds one with
`sylph sketch` run via the same docker image as `modules/local/sylph/build_db` —
adapt the genome list for your own reference genomes.

Output `<sample>.sylph_profile.tsv` (`genome_id, predicted_rel_abundance,
predicted_tax_rel_abundance`) is published next to `truth.tsv`; the raw sylph TSV
goes under `<sample>/profiling/sylph/`.

### amplicon (AAP)

`profiler=aap` runs the EBI-Metagenomics amplicon-analysis-pipeline via a nested
`nextflow run` (pinned by `--aap_revision`). The MAPseq DB can come from a built
collection or a config:

- **A built collection** — a `database` name defined in the samplesheet `databases:`
  block; the pipeline builds the mapseq quartet (`fasta`/`tax`/`otu`/`mscluster`) from
  the collection's 16S + explicit taxonomy and writes the `aap.config` for the nested
  run automatically (no `--aap_config` needed).
- **A config** — supply the MAPseq/reference databases through an extra config passed
  with `--aap_config`, e.g.:

```groovy
// aap.config
params.mapseq_databases {
    gtdb_r220 {
        fasta     = '/dbs/gtdb-r220.fasta'
        tax       = '/dbs/gtdb-r220-tax.txt'
        otu       = '/dbs/gtdb-r220.otu'
        mscluster = '/dbs/gtdb-r220.fasta.mscluster'
        label     = 'GTDB-r220'
        run_otu   = true
        run_asv   = false
    }
}
```

Outputs land under `<sample>/profiling/aap/`. The wrapper runs on the host
(`executor local`, no container), so **host `nextflow` and the container engine
must be available** to the task. Samples that share a DB config (same `database`,
`aap_configs`, and `aap_profile`) are batched into a **single** nested AAP run —
AAP's samplesheet is multi-row and namespaces output per sample — so the nested
Nextflow starts once per DB, not once per sample (study-level aggregate outputs
from the batched run are not published).

No production MAPseq database (e.g. SILVA) to hand? Build a small one from your
own reference genomes' full-length 16S rRNA sequences — NOT the amplicon
fragments used to generate the reads, mapseq needs full-length sequences to
classify correctly. The in-pipeline way is an `aap` collection in the `databases:`
block (each entry an `ssu` full-length 16S FASTA + explicit `taxonomy`): the
pipeline concatenates each entry's 16S into one fasta, writes a matching `.tax`,
runs `mapseq` once against itself to build the `.fasta.mscluster` clustering cache,
then majority-votes each cluster's taxonomy into the `.otu` table. The standalone
`examples/subspecies_v4_sweep/scripts/build_profiling_dbs.py` does the same via
docker and additionally runs `barrnap` to predict 16S for any genome missing a
pre-extracted copy.

## Outputs

Published under `results/<sample>/`:

| File | Description |
|------|-------------|
| `<sample>_R1.fastq.gz`, `<sample>_R2.fastq.gz` | Synthetic reads (single file `<sample>.fastq.gz` for single/long/amplicon). |
| `<sample>.sorted.bam` (+ `.bai`) | Ground-truth read→reference mapping. References are named `<genome_id>:<contig_id>`. |
| `<sample>.truth.tsv` | Ground-truth profile: `genome_id, target_rel_abundance, realized_n_reads, realized_rel_abundance`. |
| `<sample>.sylph_profile.tsv` | Predicted profile (sylph), when `profiler=sylph`. |
| `<sample>/profiling/` | Raw profiler outputs (`sylph/`, `aap/`). |

Built profiler databases (from a samplesheet `databases:` block) are published under
`results/databases/<name>/` — `<name>.syldb` (sylph) and/or
`<name>.mapseq.{fasta,tax,otu}` + `<name>.mapseq.fasta.mscluster` (mapseq). This dir
is directly reusable as a `path:` pre-built database entry.

And under `results/error_models/<train_id>/`:

| File | Description |
|------|-------------|
| `<train_id>.error_model_report.html` | Human-readable summary — **open this first**. Model-selection table (AIC/BIC/AICc vs context length), winning model's parameters with uncertainty, per-base error-rate/Weibull survival curve. |
| `<train_id>.model.pt` | The serialized model actually consumed by `genome-blender`/`skiver-generate` to apply the error profile; not human-readable. |
| `<train_id>.context_model_aic.csv` | Full candidate comparison. Check the `aic` column (`maximum_likelihood` rows) across candidates — lower is better; compare `train_log_likelihood` vs `test_log_likelihood` for the same model as a quick overfitting check. |
| `<train_id>.phred_calibration.json` | Empirical `P(Q \| error_type)` counts/probs. Compare the empirical error rate per reported Q against the Phred-implied rate (`10^(-Q/10)`) to sanity-check the sequencer's quality scores. |

`target_rel_abundance` is the normalised requested abundance; `realized_*` come
from `samtools idxstats` on the true BAM (what was actually generated).

## Parameters

| Param | Default | Description |
|-------|---------|-------------|
| `--input` | – | Samplesheet CSV (required). |
| `--outdir` | `./results` | Output directory. |
| `--step` | `all` | `all` (generate + profile) \| `generate` \| `profile` \| `train` (train the error model only). |
| `--seed` | `42` | Global RNG seed (training + generation). |
| `--chunks` | `1` | Default number of chunks to split generation into (per-sample override via samplesheet `chunks`). |
| `--sylph_databases` | `[:]` | Named sylph databases (see Profiling). |
| `--aap_revision` | `main` | Revision of the amplicon-analysis-pipeline for the nested run. |
| `--aap_profile` | `null` | Optional `-profile` for the nested AAP run; omitted when null (set the engine via `--aap_configs` instead). |
| `--aap_config` | `null` | Single `-c` DB-passthrough config for the nested AAP run (MAPseq databases etc.). |
| `--aap_configs` | `[]` | Ordered extra `-c` config files forwarded to the nested AAP run (engine/site config); comma-separated on the CLI. |
| `--error_model_candidates` | `AdditiveContext(5),AdditiveContext(7),AdditiveContext(9)` | Candidate skiver component strings; the min-AIC model is kept. |
| `--error_model_components` | `null` | Force a single component string (skips the AIC search). |
| `--smb_skiver_tag` / `--smb_genome_blender_tag` | `latest` | Container image tags. |
| `--max_cpus` / `--max_memory` / `--max_time` | `16` / `128.GB` / `240.h` | Resource caps (override per HPC profile). |

## Containers

Two images are built from the vendored submodules (`vendor/skiver`,
`vendor/genome-blender`) — see `containers/`:

- **`smb-skiver`** — the `skiver` Rust binary (`skiver dump`) plus the Python
  training/calibration environment.
- **`smb-genome-blender`** — `generate-reads` plus `skiver-generate` on `PATH`
  (genome-blender shells out to it), and `pysam` (used by the ground-truth step).

Build locally:

```bash
git submodule update --init --recursive
docker build -f containers/skiver/Dockerfile         -t ghcr.io/timrozday-mgnify/smb-skiver:latest .
docker build -f containers/genome-blender/Dockerfile -t ghcr.io/timrozday-mgnify/smb-genome-blender:latest .
```

The `.github/workflows/build-images.yml` workflow builds and pushes them to GHCR
on `main`/tags. GHCR packages are created **private** on first push — set each
package's visibility to **public** once in its GitHub *Package settings* so HPC
pulls need no authentication (or keep them private and pull with a token).

### HPC / Singularity

Run with `-profile singularity` (or `apptainer`). Singularity converts the
Docker images automatically. Override resources with a site config, e.g.:

```bash
nextflow run main.nf -profile singularity -c my_hpc.config --input samplesheet.yaml
```

## Testing

Three tiers, via [nf-test](https://www.nf-test.com/):

```bash
nf-test test modules/ tests/default.nf.test --tag stub   # fast, host-only, no images
nf-test test tests/default.nf.test --profile docker --tag e2e   # full run (needs images)
pytest tests/bin/test_bin.py                              # bin/ helper unit tests
```

CI (`.github/workflows/ci.yml`) runs pre-commit, the `bin/` unit tests, and the
stub nf-tests on every PR. Image builds and the full e2e run on `main`.

## Development

Set up pre-commit:

```bash
pip install pre-commit && pre-commit install
```

Test fixtures are regenerated with `python tests/data/generate_fixtures.py`.
