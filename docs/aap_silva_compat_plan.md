# superresolution-amplicon on amplicon-analysis-pipeline output and its SILVA SSU database

Status: plan, 2026-09-17. Nothing implemented. Spans `timrozday-mgnify/superresolution-amplicon`
(upstream, most of the work) and this repo (validation in Phase 4, contract changes in
Phase P). Order: Phases 0 and P are independent; P lands before Phase 1.

## Goal

Run superresolution-amplicon (SR) on what the EBI amplicon-analysis-pipeline (AAP, v6.1.5)
produces, labelling reads with the MAPseq database AAP uses (`silva-ssu-138.1`, db label
`SILVA-SSU`). Three things must hold:

1. SR's observed labels are AAP's own OTU-path labels. An AAP `.mseq` can then go straight
   into SR as `mseq:`, and a sample without one is re-mapped exactly as AAP would map it.
2. The panel kernel `K` (panel sources × SILVA labels) describes reads that went through **AAP's**
   read preparation (fastp merge, no primer trim, cmsearch clip), not SR's own. Both
   `simulate` and `align` modes have to account for the merge.
3. SR has one inference method: the panel. A sample without a panel uses the whole
   reference database as its panel. The square mis-mapping matrix and everything that
   exists only for it are removed (Phase P).

Out of scope: AAP's ASV path (DADA2), except as a comparison baseline. Also out of scope:
LSU/ITS, runs with more than one amplified region, and runs whose primers differ within a
batch.

## What AAP does to a read before MAPseq

From `~/Documents/Pipelines/amplicon-analysis-pipeline` (commit `4727842`):

| step | where | what |
|---|---|---|
| fastp | `workflows/amplicon_pipeline.nf:161` (`READS_QC_MERGE`), args `conf/modules.config:747` | `-m --merged_out --detect_adapter_for_pe`, `--cut_front/--cut_tail 3`, `--cut_right` window 4 at mean Q15, `-l 100`. fastp's merge defaults apply: overlap ≥ 30, ≤ 5 mismatches, ≤ 20%. |
| drop unmerged | `subworkflows/ebi-metagenomics/reads_qc/main.nf` (`ch_reads_se_and_merged`) | For PE runs only merged reads go on. Unmerged pairs never reach MAPseq. |
| SSU extraction | `DETECT_RNA`: `cmsearch --hmmonly --cut_ga`, then `EXTRACTCOORDS` | Each read is cut to its SSU hit. **Primers are not trimmed on the OTU path**: cutadapt only feeds DADA2. |
| MAPseq | `MAPSEQ_OTU_KRONA`, args `conf/modules.config:803` | mapseq 2.1.1b (`--h3ab3c3b_0`), `-seed 12 -tophits 80 -topotus 40 -outfmt simple`. Target is the full-length SILVA-SSU FASTA, with the shipped `.mscluster`. |

### Measured on the Nov2025 example run

Sample `SC2189280-SC3-1-26s000344`, from
`~/Documents/mimicc/fermentor-run-reports_rendered/Nov2025/amplicon-pipeline-results`:

- Reads are 2×310 cycles: 533,552 pairs, of which 517,471 merged (97.0%). fastp corrected
  13,105 bases in 7,857 reads, so correction is on. That is ~9×10⁻⁵ per merged base.
- Merged reads are 291–302 bp (mode 299). Each starts with the 515F oligo: the two most
  common starts, `GTGCCAGCAGCC` and `GTGCCAGCCGCC`, are the two options of its `M` code, at
  about 65:35. Each ends with rc(806R). The insert (≈ 292 bp) is shorter than a read, so
  the mates overlap along the whole amplicon and read into the adapter.
- The cmsearch hit removes 0–5 bases at the 5′ end and 0–2 at the 3′ end. A typical hit
  spans read coordinates 6–298 on a 299-bp read.
- Residual error, from merged reads aligned to this sample's DADA2 ASVs (edlib, 15k reads,
  rough, a lower bound):
  - Edits pile up in the first and last 30 bp: primer degeneracy plus a 0–3 bp overhang
    at the 3′ end.
  - The interior rate is about **3×10⁻⁴ per base**. That is 17× below `--flat_sub_rate`
    (5×10⁻³) and 10× below the 3.1×10⁻³ calibrated on SC2200627
    (`dev/error_rate_calibration.md`).
- DADA2 flagged 36% of reads as chimeric. The OTU path keeps chimeras.
- The primers AAP identified (`515-YF GTGYCAGCMGCCGCGGTAA`, `806BR GGACTACNVGGGTWTCTAAT`)
  are SR's defaults.
- The rendered results have no `.mseq`. `PUBLISH_OTU_RESULTS` publishes one, and the
  benchmark's `RUN_AAP` keeps it (`modules/local/amplicon_analysis/main.nf:99`), so both
  cases occur.

### The database

`~/Downloads/silva-ssu-138.1.tar.gz` unpacks to `silva-ssu-138.1/138.1/`:

- `SILVA-SSU.fasta`: 2,098,433 sequences, 2.9 GB, one line per sequence. This is SILVA
  138.1 **SSU Ref, not NR99**, all three domains. It uses the **RNA alphabet (`U`)**, and
  headers are bare `acc.start.stop`.
- `SILVA-SSU-tax.txt`: a MAPseq header with 8 levels and per-level `#cutoff`. Lineages look
  like `sk__Bacteria;k__;p__…;g__Collinsella;s__Genus_species`.
- `SILVA-SSU.otu` and `SILVA-SSU.fasta.mscluster` are shipped.

In-silico PCR on a 0.1% random sample (2,071 sequences):

- As shipped, `subspecies_infer.extract_v4` amplifies **none** of them. `U` is not in
  `_IUPAC`, so every `T` in a primer counts as a mismatch.
- After U→T it amplifies 98.6% of Bacteria (1,842/1,868), 100% of Archaea (46/46) and 54%
  of Eukaryota (84/157).

## How the merge changes M

`M[a, b]` = P(MAPseq label `b` | read from source `a`). Every step between template and
MAPseq is part of `M`. AAP's chain is: PCR → sequence each mate → quality-trim → fastp
merge (correct, or drop the pair) → cmsearch clip → MAPseq. The merge changes `M` in five
ways.

> **Phase 0.1 (2026-09-17, upstream `dev/aap_merge_effects.md`) corrects effects 1, 2 and 4.**
> fastp checks its ≤ 5 mismatch limit only on the first 50 bases of the overlap (the start
> of the fragment here). Past those, pairs merge with any number of substitutions or an
> indel. The merged read is R1: R1 substitutions and indels pass unchanged, and an R2 indel
> becomes a 1-base length error at the 3′ end. Only indels in fragment bases ~8–40 drop the
> pair. Correction needs R1 ≤ Q14 against R2 ≥ Q30, i.e. only Q12 against Q38 in the binned
> Nov2025 reads. The text below is the pre-measurement expectation.

1. **The error rate falls, but the merged read is not a consensus.**
   - Across the overlap, fastp's merged sequence takes R1's bases; rc(R2) contributes only
     beyond R1's end.
   - A mismatch is corrected only when one base is high quality and the other low (fastp's
     `BaseCorrector`: ≥ Q30 against ≤ Q14; verify in Phase 0.1). An R1 error at moderate
     quality survives even when R2 is right. That is why only 9×10⁻⁵/base was corrected.
   - The bigger effect is the **drop filter**. A pair with more than 5 overlap mismatches
     is not merged. Here the overlap is the whole insert, so R1 and R2 errors together are
     capped at about 5. That cuts off the tail of the error-count distribution.
   - Result: less leak onto 1–2-edit neighbours and a heavier diagonal.
2. **Sequencing indels are almost entirely removed.**
   - fastp's overlap is ungapped. After an indel in one mate, each base mismatches with
     probability 3/4, so the pair fails the ≤ 5 limit unless the indel is within ~6 bp of
     the overlap's end.
   - The indels that survive are in both mates, which means PCR indels.
   - This matters because indels send simulated reads to distant labels: the flat indel
     rate was 10× too high on GTDB (`dev/error_rate_calibration.md`). For merged reads the
     right indel rate is the PCR rate, not the sequencing rate.
3. **Errors shared by both mates pass untouched.**
   - Polymerase substitutions and chimeras agree across the overlap, so they merge. They
     set a floor under the leak.
   - A chimera depends on the second parent's abundance, not on the source alone. No `M`
     can hold chimeras.
4. **Where errors fall matters, not just how many.**
   - What survives is R1's moderate-quality errors, which lean toward R1's 3′ end, plus
     rc(R2) beyond R1.
   - The sites that tell V4 neighbours apart are not spread evenly along the amplicon. So
     the leak from `a` to `b` depends on where they differ, not just on their edit
     distance, and a uniform per-base rate puts the leak in the wrong places.
5. **Yield depends on the source.**
   - A pair merges only when its overlap after trimming is at least 30 bp, i.e. when
     insert ≤ 2L − trim − 30.
   - At 2×310 over 291–302 bp this never binds: 97% merge, and the failures are driven by
     errors.
   - At V3–V4 on 2×250 (≈ 460-bp insert, ≈ 40-bp overlap) it does bind. Taxa with longer
     amplicons lose reads before MAPseq.
   - That loss is a per-source efficiency: the row sums to less than 1 over observed
     labels. A row-stochastic `M` cannot show it, and it biases abundance directly.

Two more effects come from AAP's read preparation rather than from the merge:

6. **Primers stay on, and primer bases come from the oligo, not the template.** Degenerate
   positions carry the primer mix (`M` is A:C ≈ 65:35 here), whatever the source. MAPseq
   scores these bases against full-length references. `simulate_amplicon_reads.flank`
   resolves every degenerate code to its first option, which is wrong for this.
7. **cmsearch clipping** (0–5 bp) probably does not change the top hit. Phase 0.4 checks.

What each mode needs as a result:

- **simulate**: each simulated read has to pass through the same steps as a real one. That
  means either an error model of the *merged* read, or simulated mates run through the real
  fastp. Either way the reads are mapped against the same database with the same MAPseq
  build and arguments.
- **align**: `c**d` stands for "d independent errors, each at the right base". Under the
  merge, `c` becomes the merged-read residual rate (effect 1). Indel edits need a much
  smaller weight than substitutions (effect 2), and position can matter (effect 4). All of
  these can come from one error profile measured on the reads. Yield (effect 5) cannot come
  from distances at all.

## Decisions

1. **Map against AAP's database as shipped, not an extracted amplicon database.**
   - Observed, simulated and panel reads all use the same FASTA, `.tax`, `.mscluster`,
     MAPseq version and arguments.
   - Labels are then identical to AAP's, so its `.mseq` can be reused.
   - SR avoids clustering a 2.1M-sequence amplicon DB, since AAP ships the `.mscluster`.
   - Amplicon extraction still runs. It provides V4 groups, sources, the translation
     table and the `align` input, but it is no longer the MAPseq target.
2. **Use AAP's merged reads and do not merge again.** SR's `reads_to_fasta.merge_pair`
   (R1 wins, overlap ≥ 20, ≤ 25% mismatches, no qualities) accepts pairs fastp rejects,
   indel-bearing pairs in particular. Re-merging raw reads would give a different read set
   with a different `M`.
3. **The default simulation is a merged-read error model.** skiver trains on merged reads
   reference-free, as SR already does for single-end `reads:`, with position covariates.
   Simulated mates run through real fastp (`pairs`) are the validation reference, and the
   choice when yield binds.
4. **Inference is over a panel.** SILVA is a label space (`silva_taxon_panel_plan.md`
   §1.1). Only the panel's sources are simulated: its genome V4 amplicons, plus the SILVA
   V4 groups under each taxon entry (≤ `--panel_taxon_max_sources`). Those reads are mapped
   against all 2.1M references. The result is the rectangular kernel `K` (sources × SILVA
   labels). Reads on labels no source reaches go to `background`. No database reference is
   simulated unless it is a taxon-entry source.
5. **The panel is the only method** (Phase P).
   - With no `panel_references`/`panel_taxa`, the panel is the whole reference database.
     Its sources are the database's distinct amplicons, and its entries are the
     database's genomes.
   - The kernel builders, inference and fit check are the panel ones for every sample.
     The square `M`, its backends, `--mismapping_matrix` and `--infer_space v4_group` go.
   - **The whole-database panel is never tested or swept against full SILVA-SSU or GTDB.**
     It is tested only on small reference sets: the bundled fixture, the 81-reference
     *B. uniformis* set and the 20HM database. Against a generic database it would simulate
     every distinct V4 group and fit a genome space V4 cannot identify. The AAP phases
     below always name a panel.

## Phases

### Phase 0 — measure before coding (`dev/aap_merge_effects.{py,md}` upstream)

Run on the 20 Nov2025 samples. Raw reads come from the sequencing facility or ENA.

0.1 **fastp behaviour.** *Done: `dev/aap_merge_effects.md`, see the note under "How the merge
    changes M".* Simulate pairs from one known amplicon with placed errors:
    substitutions at Q20 and Q35 in either mate, and single indels at positions along the
    overlap. Run fastp 1.0.1 with AAP's arguments. Record which base the merged read keeps,
    the correction thresholds, and merge success against indel position and mismatch
    count. This confirms or corrects effects 1–2.
0.2 **Residual error profile.** *Done: `dev/aap_merge_effects.md`. Aligned to the 20HM genomes' V4
    amplicons, not ASVs. Interior substitutions 3.4×10⁻⁴/base (median column, stable across
    runs), sequencing indels ~4×10⁻⁶ (~100× lower), 20-bp window max/min ≤ 2.1, so 3.3 is
    not triggered and 3.2 is justified.* Align every merged read to its nearest ASV (edlib) and
    count substitution, insertion and deletion rates per amplicon position, excluding
    primer columns. Do this per sample and pooled. Output: `position_error_profile.tsv`
    (`pos sub ins del`). It feeds Phases 2.2 and 3.
0.3 **Primer-site mix.** *Done: `dev/aap_merge_effects.md`. The oligos are EMP 515F
    `GTGCCAGCMGCCGCGGTAA` / 806R `GGACTACHVGGGTWTCTAAT` (Y reads C, N never G), with
    heterogeneity spacers of 0/4–7 bp before 515F and 0–3 bp after 806R. The mix does not
    depend on the source genome. The reverse mix depends on overhang length and its sites
    co-vary, so 2.1 should sample the joint `ends.tsv`.* Base frequencies at each degenerate primer position, from read
    starts and ends. Output: `primer_mix.tsv`.
0.4 **cmsearch clip.** *Done: `dev/aap_merge_effects.md`. Not negligible. Clipping changes the
    MAPseq top hit for 19.2% of reads (replicate: 0%). 98.8% of those changes stay inside one
    V4 group, so V4-group agreement is 99.77–99.81%, below the 99.9% threshold. Emulate the
    clip, and compare labels as V4 groups. Also: the Nov2025 run mapped against the older
    `SSU.fasta`/`slv_ssu_filtered2.txt`, not `silva-ssu-138.1`.* Map one sample's merged reads with AAP's arguments against the full
    DB, both unclipped and clipped to the deoverlapped coordinates. If top hits agree on
    ≥ 99.9% of reads, clipping is not emulated anywhere.
0.5 **MAPseq build and format.** *Done: `dev/aap_merge_effects.md`, all 517k clipped reads. With AAP's
    arguments the two builds give identical output. With SR's default `mapseq_args` ('')
    the top hit differs for 20.2% of reads and the V4 group for 2.8%, so 1.3 must pass AAP's
    `-seed 12 -tophits 80 -topotus 40`. Columns 1, 2 and 4 are the same in the simple and
    confidences formats.* Map the same 10k reads with `2.1.1b--h3ab3c3b_0` (AAP) and
    `2.1.1b--hc47f52e_1` (SR); top hits must be identical. Check that `-outfmt simple`
    columns 1, 2 and 4 are what `subspecies_infer.iter_mseq` reads.
0.6 **Label coverage.** *Done: `dev/aap_merge_effects.md`, all 20 runs. 0.29% pooled (at most
    0.56% per run) of reads hit a reference without an amplicon, below the 1% trigger, so 2.4's
    relaxed extraction is not needed. 99.96% of those reads hit references with a
    one-base indel in the primer site, not truncated ones; three references hold 98%.* Per sample, the fraction of reads whose top hit is a
    non-amplifiable reference (after U→T). This decides Phase 2.4.

**Acceptance:** a number for each item, and Decision 3 plus the Phase 2.4 and 3.3 triggers
confirmed or revised in this plan.

### Phase P — panel-only pipeline (upstream, plus the benchmark contract)

Today a sample runs one of two paths:
- **square:** `MATRIX_KEY` → `SIMULATE_READS`/`MAPSEQ_SIM`/`BUILD_MISMAPPING`, or
  `ALIGN_MISMAPPING`/`GROUPED_MISMAPPING`/`MINIMAP2_*` → `PUBLISH_MISMAPPING` →
  `infer_composition.run`;
- **panel:** `PANEL_PREPARE` → `MAPSEQ_PANEL_HOME` + `SIMULATE_PANEL_READS`/
  `MAPSEQ_PANEL_SIM` → `PANEL_KERNEL` → `infer_composition.run_panel`.

`main.nf:79–87` limits the panel path to `simulate` and genome space, and
`build_panel_kernel.py align` exists but is not wired. Phase P keeps only the panel path.

P.0 **Snapshot before deleting.** *Done 2026-09-17: upstream
    `dev/panel_only_parity/` holds the fixture and 81-amplicon B. uniformis square-path
    outputs, the existing 20HM genome-panel profiles, and a stored-MAPseq SILVA species
    panel re-score.  See its README for source paths, settings and checksums.* The snapshot
    comprises four runs:
    - the bundled fixture with a square matrix (`simulate`, `align exact-hash`,
      `align kmer tau 1`);
    - the *B. uniformis* 81-reference set (`dev/alignment_mismapping.md`), same three modes;
    - the 20HM database, with its existing genome panel;
    - one species taxon panel from `dev/panel_silva_sweep.md`, re-scored from its stored
      `.mseq`. Do not re-run it against SILVA.

    These are the P.7 references. Nothing against full SILVA or GTDB is re-run.
P.1 **The default panel.** *Done 2026-09-17: `prepare --whole-database` reuses the
    extracted amplicons and translation table; unconfigured samples now carry the
    `database` panel marker and pass the amplicon directory to `PANEL_PREPARE`. Focused
    pytest and stub workflow checks cover both the new default and existing explicit panels.*
    - `build_panel_kernel.py prepare --whole-database <amplicon_dir>` reads
      `EXTRACT_AMPLICONS`' `amplicons.fasta` and `translation_table.tsv`. It writes
      `sources.fasta` (one `v4g_<sha16>` per distinct amplicon), `panel_translation.tsv`
      (genome → source, copy weight as in `T`) and `sources.tsv` (`in_db` all true).
    - It does not repeat the in-silico PCR, so sources are byte-identical to the database
      groups. Today `prepare` cuts panel genomes at `--max-mismatch 2` while
      `EXTRACT_AMPLICONS` uses 3; a second PCR would drift from the labels.
    - It logs a warning with the source count above 10,000, so an accidental full-database
      run is visible early. It does not refuse to run.
    - `main.nf`: a row with no `panel_references`/`panel_taxa` gets `meta.panel = 'database'`.
      The `main.nf:79–87` restrictions go, except that `panel_taxa` still needs
      `--taxonomy`. `PANEL_PREPARE` passes `--whole-database` for that value.
P.2 **Workflow.** *Done 2026-09-17: removed the square-matrix modules and workflow
    branch; panel kernels now group and publish by database/panel/settings identity. The
    new `PANEL_ALIGN` branch invokes `build_panel_kernel.py align`, while the simulation
    branch remains panel-only. Stub coverage passes 13/13, including default, explicit,
    taxon, shared-panel, and alignment cases.*
    - Delete the square branch of `workflows/superresolution_amplicon.nf`:
      `MAPSEQ_CLUSTER_MATRIX`, `SIMULATE_READS`, `MAPSEQ_SIM`, `BUILD_MISMAPPING`,
      `ALIGN_MISMAPPING`, `GROUPED_MISMAPPING`, `MINIMAP2_INDEX`, `MINIMAP2_ALLVSALL`, the
      square `MATRIX_KEY` grouping, and `PUBLISH_MISMAPPING`'s `groups.tsv`. Delete their
      module directories.
    - `MAPSEQ_CLUSTER` stays: observed reads still map against the database.
    - Add `PANEL_ALIGN`, which wires the existing `build_panel_kernel.py align`. Under
      `--mismapping_method align` it replaces `SIMULATE_PANEL_READS` and
      `MAPSEQ_PANEL_SIM`; the home mapping still runs.
    - The kernel group key is (database amplicons sha256, panel or `database`, panel taxa,
      model identity, method and settings). Kernels publish to
      `mismapping/panel_<key>/`, keeping the existing layout.
P.3 **Align on the whole-database panel.** *Done 2026-09-17: whole-database sources
    now use a literal hash join at tau 0 and source-probed pigeonhole candidates with
    bounded IUPAC verification at tau >= 1. `auto` decay again measures either the flat
    rates or a staged pre-trained model; backend/minimap settings are removed. The
    measured-mismapping tests pass 14/14 and the stub workflow suite 13/13.*
    - Removed: `--align_backend`, `--minimap2_args` and `--minimap2_index_args` (one builder
      remains).
    - Kept: `--align_tau`, `--align_distance_decay` (including `auto`),
      `--align_decay_model`, `--align_ambiguity_weight`, `--max_ambiguous_bases`,
      `--max_postings`.
P.4 **A supplied kernel replaces `--mismapping_matrix`.**
    - `--panel_kernel <dir>` takes a published `mismapping/panel_<key>/` bundle
      (`mismapping_matrix.npz`, `panel_translation.tsv`, `sources.tsv`, `provenance.json`)
      and skips kernel building.
    - It is refused unless the bundle's `db_amplicons_sha256` (already written by
      `sm.write_kernel`) matches the sample's extracted amplicons and its panel identity
      matches the sample's.
    - This keeps the benchmark's build-once / infer-many split and cheap `mseq:` sweeps.
P.5 **Inference and fit check.**
    - `infer_composition.run` becomes `run_panel`. Delete the square code:
      `_mismapping_matrix`, `--build-mismapping`, `--build-grouped`, `--active-amplicons`,
      `--infer-space`, `--infer-prune`.
    - Delete the square branch of `check_composition_fit.py`.
    - Carried over, confirm each with a test: presence gate, horseshoe, latent distance
      decay (strata from `align --tau >= 1`), `--no-mismapping`, `low_depth` and
      `no_reference_hits`, and the `resolution` column.
    - Prune only if P.7 shows the whole-database panel is too slow on the 20HM database.
      Pruning would drop sources whose row reaches no observed label, with the pruned mass
      going to the existing zero-count sink.
    - `--sim_n_per_ref` and `--panel_sim_n_per_ref` merge into one `--sim_n_per_ref`,
      default 5,000.
P.6 **Delete what only the square path used.**
    - Delete: `bin/select_active_amplicons.py`;
      `subspecies_infer.build_mismapping`/`build_mismapping_grouped`; the grouped and CSR
      readers in `sparse_matrix.py` (keep `read_kernel`/`write_kernel`); legacy CSV
      matrices; `tests/test_grouped_mismapping.py`; the square parts of
      `tests/test_measured_mismapping.py` and `tests/test_composition_fit.py`.
    - In `build_mismapping_align.py`, keep only what `build_panel_kernel.py` imports:
      `bounded_iupac_distance`, `ambiguity_weights`, `pigeonhole_candidates`,
      `measure_error_rate`. Move them into a `kernel_align.py` if the file is left mostly
      empty.
    - `dev/` scripts that import deleted functions stay as a record and are marked
      unmaintained in their `.md`.
    - `tally_kernel` stays (`build` uses it).
P.7 **Tests.**
    - Stub nf-test: a row with no panel, a genome panel, and a taxon panel, each under
      `simulate` and `align`, plus `--panel_kernel` reuse.
    - Real e2e on the bundled fixture: whole-database panel under `simulate`.
    - pytest: `prepare --whole-database` sources equal the database groups; the
      pigeonhole `align` equals the brute-force `align` on the *B. uniformis* set; a
      `--panel_kernel` sha mismatch is refused.
    - Parity against P.0, whole-database panel against the old square run:
      - genome TV ≤ 0.02 on the fixture and on *B. uniformis* in each mode;
      - the 20HM panel and the taxon-panel re-score match their snapshots exactly, since
        no numerical code on that path changes.

      Any larger difference is explained in `dev/panel_only_parity.md` before the square
      code is deleted. The expected source is the `background` row and home labels,
      neither of which the square path has.
    - No test, sweep or example runs the whole-database panel against full SILVA-SSU or
      GTDB.
P.8 **Benchmark contract** (this repo).
    - `BUILD_SUPERRESOLUTION_MISMAPPING` lifts the `mismapping/panel_<key>/` directory
      instead of `mismapping_matrix.npz`. `RUN_SUPERRESOLUTION` passes `--panel_kernel`.
      The "exactly one bundle per run" rule stays.
    - Panel `sr_settings` stop being a special case: every amplicon run builds a kernel
      first. `srSetKey` keeps the panel in the key, as now.
    - `sr_amplicon_matrix_args` loses the backend and minimap2 flags.
    - `bin/normalize_sr_profile.py` drops the `background` row and renormalises. It also
      writes the background share to a sidecar, so the profile is still over genomes and
      the unexplained fraction is still reported.
    - `examples/sr_amplicon_silva_sweep` runs `--infer_space v4_group` over full SILVA,
      which is gone and is the untested case. Retire it and point to
      `sr_amplicon_panel_sweep` (species taxon panels).
    - Update `examples/sr_amplicon_param_sweep`, `subspecies_v4_sweep`, README and
      `CLAUDE.md` for the removed flags. The shotgun profiler is untouched.

**Acceptance:**
- P.7 passes.
- `git grep -E "build_mismapping_grouped|GROUPED_MISMAPPING|infer_space|mismapping_matrix\b"`
  upstream finds only `dev/` and the kernel file name.
- The benchmark's stub tests pass, and its e2e passes on the fixture panel sweep.

### Phase 1 — read AAP's database and outputs

1.1 **U→T in `subspecies_infer.read_fasta`.** One line. Everything that reads references
    goes through it: extraction, `align`, panel prepare. Add a `U` record to
    `stage_amplicons`'s demo.
1.2 **`build_mapseq_database.py --aap-tax SILVA-SSU-tax.txt -o <prefix>.tax`.** Writes the
    SR-form lineage used by `--taxonomy`, `lca` and panel taxa:
    - strip the `xx__` prefixes;
    - write the species as `Genus epithet` (`s__Bacteroides_fragilis` →
      `Bacteroides fragilis`), matching the SILVA form's species rank;
    - pad empty ranks with `unclassified`, as the SILVA form does.

    MAPseq still gets AAP's own tax file, unchanged.
1.3 **`--mapseq_target amplicons|references`** (default `amplicons`, so nothing changes
    for existing runs).
    - With `references`, MAPSEQ_OBS and MAPSEQ_PANEL_{HOME,SIM} map against the
      `--references` FASTA itself. Its tax file (`<stem>-tax.txt` or `<stem>.tax`) and
      `<fasta>.mscluster` must sit next to it; if either is missing the run errors.
      Clustering 2.1M full-length sequences in the pipeline is not offered, and
      MAPSEQ_CLUSTER is skipped.
    - The panel kernel id and `provenance.json` gain the target plus the sha256 of the FASTA
      and `.mscluster`.
    - `references` implies AAP's read form: observed and simulated reads keep their
      primers. `--trim_primers true` together with `references` is an error in `main.nf`,
      because an AAP `.mseq` of untrimmed reads would otherwise be paired with a matrix of
      trimmed ones.
1.4 **Bare-accession headers.** With no `|`, `genome_of_header` returns the whole
    accession, so each sequence is its own genome. That is what the SILVA form means
    already, so no header rewrite is needed. Add a test that pins it.
1.5 **`bin/aap_samplesheet.py <aap_outdir> [--db-label SILVA-SSU]`.** Writes SR samplesheet
    rows `{id, reads: <id>/qc/<id>.merged.fastq.gz, merged: true, merge_rate,
    mseq: <id>/taxonomy-summary/SILVA-SSU/<id>_SILVA-SSU.mseq(.gz)}`, with `mseq` only when
    the file exists. It:
    - skips runs missing from `qc_passed_runs.csv`;
    - refuses a run whose `amplified-region-inference/<id>.tsv` is not a single 16S V4;
    - refuses a run whose `primer-identification/*.fasta` differs from
      `--fwd_primer/--rev_primer` (primers are global in SR);
    - takes `merge_rate` from `fastp.json` (`merged_and_filtered.total_reads` over input
      pairs) and warns below 80%, the same threshold `reads_to_fasta.py` uses.
1.6 **`merged: true` in the samplesheet.** READS_TO_FASTA never merges such reads and never
    trims their primers under `references`. `obs_max_reads` still applies.

Tests:
- a stub nf-test on a miniature AAP outdir fixture (two runs, one with `.mseq`);
- pytest for `aap_samplesheet.py` and `--aap-tax`;
- the U→T demo.

**Acceptance:** on one Nov2025 sample, an SR run with `--mapseq_target references` and a
small genome panel (20HM) completes. Its `obs.mseq.gz` top hits must equal those of a
direct AAP-arguments MAPseq run on the same reads.

### Phase 2 — simulate mode for merged reads

2.1 **Primer flank from the observed mix.** `simulate_amplicon_reads.py --primer-mix
    primer_mix.tsv`: each degenerate primer position is drawn from the Phase 0.3
    frequencies. Without the file it draws uniformly over the code's options, never the
    first option. Under `references` the primers are not trimmed after errors are applied.
2.2 **Merged-read error model** (default for `merged: true`).
    - `--sim_error_model trained` trains on the merged reads. Add
      `AdditiveContext(7)+Position(2)` and `AdditiveContext(9)+Position(2)` to the
      candidates (skiver's `parse_model_components` accepts `Position(N)`), so the
      3′-weighted residual of effect 4 can be learned.
    - Default to `--trained_error_model_scope pooled` for AAP batches.
    - Under `flat`, the README gives Phase 0.2's pooled sub/ins/del rates as the AAP
      starting point, marked as per-run numbers.
2.3 **`--sim_read_structure merged|pairs`**, with `pairs` for validation and for runs where
    yield binds.
    - New `bin/simulate_amplicon_pairs.py`:
      - draws the primer-flanked fragment;
      - reads `--sim_mate_len` cycles of R1 forward and of R2 on the reverse strand
        (`apply_batch(..., is_forward=False, emit_quality=True)`), continuing into the
        TruSeq adapter past the fragment end so `--detect_adapter_for_pe` has something to
        remove;
      - writes paired FASTQ.
    - New module `FASTP_MERGE` (fastp 1.0.1 biocontainer). Its `ext.args` are AAP's
      READS_QC_MERGE arguments, verbatim, plus `-m --merged_out --detect_adapter_for_pe`.
    - The merged reads go to MAPSEQ_PANEL_SIM. Check that the source-carrying read name survives
      fastp's ` merged_x_y` suffix as `iter_mseq`'s first token.
    - Per-source yield (merged over simulated) is written to `yield.tsv` in the bundle.
    - The model is a per-mate skiver model trained on raw R1/R2
      (`AdditiveContext(7)+Position(2)+Strand`), from `error_model:` or from a new
      `raw_reads:` samplesheet key.
    - SR's own `--paired` merge stays as it is for non-AAP inputs. Its divergence from
      fastp is documented, not fixed.
2.4 **Labels outside the amplifiable set.**
    - Panel reads mapped against the full DB can land on references that have no
      amplicon. Each such reference becomes its own label `nonamp_<acc>` in
      `build_panel_kernel.py`'s label space. It is a column only; no source row is needed,
      because only panel sources have rows.
    - Trigger for more: if Phase 0.6 puts > 1% of observed reads on such labels, add a
      relaxed extraction pass that cuts at one primer and the expected length.
2.5 **What gets simulated.** Only the panel sources: `PANEL_PREPARE`'s distinct V4
    amplicons of the genome entries, plus each taxon entry's SILVA groups (capped by
    `--panel_taxon_max_sources`). They are simulated at `--panel_sim_n_per_ref` (5,000)
    and mapped against all 2.1M references. The existing panel workflow already does this
    (`MAPSEQ_PANEL_HOME`, `SIMULATE_PANEL_READS`, `MAPSEQ_PANEL_SIM`); under
    `--mapseq_target references` those mappings just use the AAP database.
    - Cost: sources × 5,000 reads. A 20-genome panel with a few species entries is
      ~100–500 sources, i.e. 0.5–2.5M reads, the MAPseq work of 1–5 AAP samples.
    - One kernel serves every sample on the same panel, database and model, so nothing
      is batch-specific.
    - A sample with no panel against the AAP database would make the whole-database panel
      all 2.1M references. Per Decision 5 that is never run or tested. `aap_samplesheet.py`
      needs `--panel` or `--panel-taxa` and writes them into every row.

**Acceptance** (Nov2025 batch plus the Phase 4 truth sample):
- `merged` against `pairs` on the same panel: median row L1 of `K` ≤ 0.05. If it fails,
  `pairs` becomes the default and 2.2 is documented as insufficient.
- Predicted mass on labels that drew no observed reads ≤ 0.3% (GTDB calibration reached
  0.24%).
- The panel fit check is `ok` on ≥ 18 of the 20 Nov2025 samples. Where it is not, check
  `unexplained_fraction` and whether the misfit sits on labels the panel reaches (chimeras
  between panel members, see Risks) before tuning anything else.

### Phase 3 — align mode for merged reads

3.1 **Decay from the merged reads.** `--align_distance_decay auto --align_decay_model
    <merged-read model>` already measures `c` by applying the model
    (`build_mismapping_align.measure_error_rate`), so no code is needed here. Add a test
    that `auto` on the 2.2 model lands within 2× of Phase 0.2's interior substitution rate,
    and a README note that for AAP inputs the model must be the merged-read one, not a
    per-mate model.
3.0 **Where align applies.** After Phase P, align means `PANEL_ALIGN`
    (`build_panel_kernel.py align`, pigeonhole candidates). It writes
    `K[s, l] ∝ w(l)·c^d(s, l)` over the SILVA groups within tau of each panel source,
    plus one MAPseq home mapping of the sources against the full DB. Only panel sources
    have rows. 3.1–3.5 change that builder and the shared distance helpers it imports.
3.2 **Separate indel decay.** `--align_indel_decay`, defaulting to `--align_distance_decay`
    so existing kernels do not change.
    - The IUPAC-aware bounded edit distance shared with `build_mismapping_align.py` returns
      `(subs, indels)`, and the weight becomes `c_sub**s * c_indel**i`.
    - `auto` measures `c_indel` from the model's indel rate.
    - Effect 2 is why this is needed: a one-indel neighbour of a merged read leaks at the
      PCR indel rate, not the substitution rate.
3.3 **Position-aware decay** (conditional).
    - Trigger: Phase 0.2's interior profile varies more than 3× (max/min over 20-bp
      windows). Otherwise a scalar `c` is enough; skip it.
    - Mechanism: `--align_position_profile position_error_profile.tsv`. A pair's weight is
      the product of `r(p)` over its differing columns, normalised so a flat profile
      reproduces `c**d`.
    - It needs the differing columns, so verified pairs within tau get an edlib traceback.
      Profile coordinates include the forward primer; offset by its length.
3.4 **Refuse `align` when yield binds.** If any sample's `merge_rate` (1.5) is below 0.8,
    `align` errors with the same style of message as the `--sim_read_len` refusal: use
    `simulate --sim_read_structure pairs`.
3.5 **`--infer_distance_decay`** fits `c_sub` only; `c_indel` stays at its built value.
    Document it.

**Acceptance:** on the Phase 2 panel, align (tau 1, merged-read `auto` decay, indel
decay) against simulate `merged`:
- mean row L1 of `K` ≤ 0.08, the level reached on the *B. uniformis* set;
- panel-entry TV between the two inferences ≤ 0.02 on every sample.

### Phase 4 — validation against truth (this repo)

4.1 **A synthetic sample.** Paired-end amplicon sample (`paired_end: true`,
    `read_length_mean: 310`, V4), error model trained on a raw Nov2025 run, `profilers:
    [aap, sr_amplicon]`. The sample's genomes come from SILVA-represented taxa present in
    the fermentor runs.
4.2 **sr_amplicon reads RUN_AAP's output** for the same sample.
    - PROFILE passes `profiling/aap/<id>/qc/<id>.merged.fastq.gz` as `reads` with
      `merged: true`, and the `.mseq.gz` as `mseq:`.
    - RUN_SUPERRESOLUTION already writes `mseq:` rows.
    - The DB is a `databases:` `path:` entry pointing at the unpacked AAP SILVA directory.
    - Worked example: `examples/sr_amplicon_aap/`.
4.3 **Score** panel-entry (genome and species) TV to truth, one arm each:
    1. raw AAP labels;
    2. simulate `merged`;
    3. simulate `pairs`;
    4. align (Phase 3);
    5. SR from raw reads with its own merge and amplicon DB, which measures what
       Decisions 1–2 buy.

**Acceptance:**
- Arm 2 beats arm 1 on species TV.
- Arm 3 is not better than arm 2 by more than 0.01; otherwise `pairs` becomes the default.
- The Nov2025 comparison against exact-matched DADA2 ASV profiles is reported as
  descriptive only. ASVs are not a gold standard, and 36% of reads were removed as chimeras.

### Phase 5 — documentation

- Upstream README: an "Inputs from amplicon-analysis-pipeline" section covering the
  samplesheet helper, `--mapseq_target references`, the merged-read error model, and the
  align decays.
- Upstream `docs/alignment_mismapping_plan.md` §7: merged reads from a fastp merge satisfy
  the spanning assumption, subject to 3.4.
- This repo's `CLAUDE.md`: one line under Key implementation notes.

## Risks and open questions

- **Chimeras.** The OTU path keeps them: 36% of reads by DADA2's count in the example
  sample, possibly inflated by closely related strains. A chimera on a label no panel
  source reaches goes to `background`, which is harmless. A chimera between two panel
  members lands on a label the panel does reach and gets absorbed into a member. This is
  the known ceiling in the panel plan's "Out-of-panel reads", and Phase 4's truth sample
  should include chimeras to size it.
- **SILVA Ref is not NR99.** Its duplicate groups are larger. Species taxon entries will
  have more sources than NR99's 1–24 groups per species, so re-check
  `--panel_taxon_max_sources`.
- **Extraction cost.** In-silico PCR over 2.1M sequences is pure Python: hours, once,
  cached by path. Publish the amplicon directory and reuse it.
- **fastp and MAPseq versions** must follow AAP's pins (fastp 1.0.1 per `fastp.json`;
  mapseq per 0.5). An AAP upgrade changes `K`, so both versions go in the panel kernel id.
- **Non-V4 or 2×250 runs.** Nothing here has been measured on them. Yield (effect 5) is
  expected to dominate there, and the `pairs` path exists for that case.

## Files touched

Upstream:
- `bin/subspecies_infer.py` (U→T)
- `bin/build_mapseq_database.py` (`--aap-tax`)
- `bin/aap_samplesheet.py` (new)
- `bin/simulate_amplicon_reads.py` (primer mix)
- `bin/simulate_amplicon_pairs.py` (new)
- `bin/build_panel_kernel.py` (`--whole-database`; pigeonhole `align`; non-amplifiable
  labels; indel and position decay)
- `bin/build_mismapping_align.py` (cut to the shared helpers; edit distance returns subs
  and indels)
- `bin/infer_composition.py`, `bin/check_composition_fit.py`, `bin/sparse_matrix.py`,
  `bin/subspecies_infer.py` (square code removed)
- deleted: `bin/select_active_amplicons.py`; `modules/local/{build_mismapping,
  align_mismapping,grouped_mismapping,minimap2}` (`simulate_reads` stays as
  `SIMULATE_PANEL_READS`; `matrix_key` stays because the panel groups key off it); `tests/test_grouped_mismapping.py`
- `modules/local/panel_kernel/align/` (new), `modules/local/fastp/merge/` (new)
- `modules/local/mapseq/{map,cluster}/`
- `workflows/superresolution_amplicon.nf`, `main.nf`, `conf/modules.config`,
  `nextflow.config`
- `README.md`, `dev/aap_merge_effects.{py,md}`, `dev/panel_only_parity/`, tests

This repo:
- `subworkflows/local/profile/main.nf`, `modules/local/superresolution/run/` (kernel bundle,
  `--panel_kernel`)
- `bin/normalize_sr_profile.py` (drop `background`)
- `examples/sr_amplicon_aap/` (new); `examples/sr_amplicon_silva_sweep/` (retired);
  `examples/sr_amplicon_param_sweep/`, `examples/subspecies_v4_sweep/`
- `README.md`, `CLAUDE.md`
