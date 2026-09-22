# sr_amplicon on AAP output: the Phase 4 truth sample

Phase 4 of `docs/aap_silva_compat_plan.md`. One synthetic V4 community is profiled by the
amplicon-analysis-pipeline (AAP) and by superresolution-amplicon (SR). Both map against
AAP's own SILVA-SSU 138.1 directory, and the results are scored against the known truth.

## The sample (4.1)

- **Genomes:** 16 of the 22 20HM genomes at abundances from 400 down to 1 (`genomes.csv`).
  The other 6 are absent, but they stay in the panel, so false positives show up as
  `mass_absent`.
- **Reads:** 200k pairs, 2x310, 515F/806R, from full-length 16S by in-silico PCR. Primers
  stay on the reads, which AAP needs.
- **Error model:** skiver, trained on a raw fermentor run (SC2200627 lane C1, 2x310, Finn).
  No raw Nov2025 FASTQs were at hand, so this run stands in for them.
- **No chimeras.** genome-blender doesn't make them, so the chimera risk in the plan is
  not sized here.

## The arms (4.2, 4.3)

| arm | name | reads | kernel |
|---|---|---|---|
| 1 | `aap` | AAP's own SILVA labels | none |
| 2 | `sim_merged` | AAP merged reads + `.mseq` | simulated merged reads, flat at 0.2's rates |
| 3 | `sim_pairs` | AAP merged reads + `.mseq` | simulated mate pairs through AAP's fastp, generating model (an oracle) |
| 4 | `align` | AAP merged reads + `.mseq` | align, tau 1, sub and indel decays `auto` from 0.2's rates |
| 5 | `raw` | raw pairs, SR's own merge and trim | simulated, flat at SC2200627's calibrated rate |

Arms 2–4 use the `aap_reads` sr_settings knob. Every SR arm reinterprets through the 22-genome
panel. Each arm builds its own kernel.

## Run

    ./run.sh

It trains the model (`--step train`), then generates the sample and profiles it
(`--step all`), then writes `results/sr_amplicon_aap/scores.csv`. The configs cap memory
at 18 GB for a 24 GB machine: `benchmark.config` for the outer tasks, `nested.config` for
the nested runs. Paths in `samplesheet.yaml` are this machine's. Edit them elsewhere.

## Score

`score.py` writes one row per arm:
- `tv`: species TV over the panel, renormalised over the panel on both sides;
- `background`: the share outside the panel;
- `mass_absent`: the share on absent genomes.

Each genome is its own species here, so genome and species TV coincide.

`species.tsv` maps each genome to its SILVA 138.1 species label. Six genomes have none
(e.g. *R. gnavus*, *E. rectale*): AAP's labels cannot reach them, and arm 1 scores that.

**Acceptance:**
- arm 2 beats arm 1 on species TV;
- arm 3 is not better than arm 2 by more than 0.01, otherwise `pairs` becomes the default.
