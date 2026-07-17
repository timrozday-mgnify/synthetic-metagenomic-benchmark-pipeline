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

## Samplesheet

A YAML list, one entry per synthetic sample (`tests/samplesheets/test.yaml` is a
working example):

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
| `database` | Sylph only: a key in `params.sylph_databases`, or `self` to build the DB from this sample's reference genomes. |
| `subsample` | Optional list of read depths to sweep (absolute read/pair counts). The full draw is generated once, then subsampled to each depth — each gets its own `subsample_<N>/` output dir with its own reads + ground truth + profile. `none`/`null`/empty/omitted → a single full-depth run in `<sample>/`. |
| `chunks` | Optional. Split generation of `num_reads` across N parallel `generate-reads` calls (merged back into one reads-set + BAM before subsampling/ground truth), useful for large `num_reads`. Blank → `params.chunks` (default `1`, no chunking). |

Relative `genomes_csv` / FASTA / FASTQ paths resolve against the pipeline
directory; absolute paths and `scheme://` URLs pass through.

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
`nextflow run` (pinned by `--aap_revision`). Its MAPseq/reference databases are
supplied through an extra config passed with `--aap_config`, e.g.:

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
must be available** to the task.

No production MAPseq database (e.g. SILVA) to hand? Build a small one from your
own reference genomes' full-length 16S rRNA sequences — NOT the amplicon
fragments used to generate the reads, mapseq needs full-length sequences to
classify correctly. `examples/subspecies_v4_sweep/scripts/build_profiling_dbs.py`
does this end-to-end via docker (`mapseq`'s own image, `barrnap` to predict 16S
for any genome missing a pre-extracted copy): concatenates each genome's 16S
into one fasta, writes a matching `.tax` file, runs `mapseq` once against itself
to build the `.fasta.mscluster` clustering cache, then derives the `.otu` table
from that clustering by majority-voting each cluster's genome-level taxonomy.

## Outputs

Published under `results/<sample>/`:

| File | Description |
|------|-------------|
| `<sample>_R1.fastq.gz`, `<sample>_R2.fastq.gz` | Synthetic reads (single file `<sample>.fastq.gz` for single/long/amplicon). |
| `<sample>.sorted.bam` (+ `.bai`) | Ground-truth read→reference mapping. References are named `<genome_id>:<contig_id>`. |
| `<sample>.truth.tsv` | Ground-truth profile: `genome_id, target_rel_abundance, realized_n_reads, realized_rel_abundance`. |
| `<sample>.sylph_profile.tsv` | Predicted profile (sylph), when `profiler=sylph`. |
| `<sample>/profiling/` | Raw profiler outputs (`sylph/`, `aap/`). |

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
| `--aap_profile` | `docker` | `-profile` for the nested AAP run. |
| `--aap_config` | `null` | Extra `-c` config for the nested AAP run (MAPseq databases etc.). |
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
