# superresolution-amplicon on amplicon-analysis-pipeline output and its SILVA SSU database

Status: plan, revised 2026-09-21. It was first written 2026-09-17 and now absorbs
`sr_full_length_mapseq_plan.md`. The work spans `timrozday-mgnify/superresolution-amplicon`
(upstream, most of it) and this repo (the benchmark contract and validation).

| Phase | State |
|---|---|
| 0 (measure) | done |
| P.0–P.3 (panel-only) | done, on upstream `main` only |
| 1.1 (U→T) | done, on both upstream lines |
| R (reconcile) | done (#18) |
| F (full-length MAPseq database) | done (#19) |
| P.4 (`--panel_kernel`) | done (#20) |
| P.5 (panel-only inference) | done (#21) |
| P.6 (delete square-only code) | done (#21) |
| P.7 (tests, parity) | done (#21); test at >= 10k reads |
| P.8 (benchmark contract) | done on `panel-only-contract` (#38); upstream #21 merged (`e29bf69`) |
| 1.5–1.6 (AAP samplesheet, `merged:`) | done (upstream #22, `462603b`); acceptance run passed |
| 2.1 (primer mix) | done, upstream #23, not merged |
| 2.2 (merged-read model) | done without `Position`, upstream #24 (stacked on #23), not merged |
| 2.3 (`pairs`) | done without `raw_reads:`, upstream #25 (stacked on #24), not merged |
| everything else | not done |

Order: R → F → P.4–P.8 → 1 → 2 → 3 → 4 → 5. See Decision 7.

## Goal

Run superresolution-amplicon (SR) on what the EBI amplicon-analysis-pipeline (AAP, v6.1.5)
produces, labelling reads with the MAPseq database AAP uses (`silva-ssu-138.1`, db label
`SILVA-SSU`). Four things must hold:

1. **One MAPseq database: the reference FASTA itself.** Every MAPseq mapping goes against
   the full-length reference FASTA, its tax file and its `.mscluster`. That covers
   observed, simulated, panel and home reads. There is no amplicon MAPseq database,
   prebuilt or built, for any reference set.
2. **SR's observed labels are AAP's own OTU-path labels.** An AAP `.mseq` can then go
   straight into SR as `mseq:`. A sample without one is re-mapped exactly as AAP would map
   it.
3. **The panel kernel `K` (panel sources × SILVA labels) describes reads that went through
   AAP's read preparation** (fastp merge, no primer trim, cmsearch clip), not SR's own.
   Both `simulate` and `align` modes have to account for the merge.
4. **SR has one inference method: the panel.** A sample without a panel uses the whole
   reference database as its panel. The square mis-mapping matrix, and everything that
   exists only for it, is removed.

Out of scope:
- AAP's ASV path (DADA2), except as a comparison baseline;
- LSU/ITS;
- runs with more than one amplified region;
- runs whose primers differ within a batch.

## Upstream has two lines

Both branched from #12 (`4733662`):

- **`main`**: #13 (panel-only, P.0–P.3) and #14 (U→T, `5aeffa3`).
- **`silva-fixes-pinned`**: #15–#17, not the panel-only migration. The benchmark pins this
  line (`examples/sr_amplicon_silva_sweep/benchmark.config`, `sr_revision bd9228f`)
  because the SILVA sweep's `silva` arm runs the square `--infer_space v4_group` path that
  panel-only removes. Its commits:
  - `105f484`: U→T, a duplicate of #14;
  - `d448dfa`: vectorised, multi-process amplicon extraction;
  - `6c9390a`: `--amplicon_cache`;
  - `ac9328e`: one-record-query clustering, plus `dev/mapseq_free_labels`;
  - `4d8f665`: AAP's MAPseq flags by default;
  - `c1f2031`: `--align_home_probes`, which only feeds the square grouped matrix;
  - `eca1763`: the prebuilt amplicon MAPseq database, which is what makes AAP's own
    directory unusable today.

Phase R ports what is independent of the line onto `main` and drops the rest.

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

- **Merge.** Reads are 2×310 cycles: 533,552 pairs, of which 517,471 merged (97.0%). fastp
  corrected 13,105 bases in 7,857 reads, so correction is on. That is ~9×10⁻⁵ per merged
  base.
- **Read shape.** Merged reads are 291–302 bp (mode 299).
  - Each starts with the 515F oligo. The two most common starts, `GTGCCAGCAGCC` and
    `GTGCCAGCCGCC`, are the two options of its `M` code, at about 65:35.
  - Each ends with rc(806R).
  - The insert (≈ 292 bp) is shorter than a read, so the mates overlap along the whole
    amplicon and read into the adapter.
- **cmsearch clip.** The hit removes 0–5 bases at the 5′ end and 0–2 at the 3′ end. A
  typical hit spans read coordinates 6–298 on a 299-bp read.
- **Residual error**, from merged reads aligned to this sample's DADA2 ASVs (edlib, 15k
  reads, rough, a lower bound):
  - Edits pile up in the first and last 30 bp: primer degeneracy plus a 0–3 bp overhang
    at the 3′ end.
  - The interior rate is about **3×10⁻⁴ per base**. That is 17× below `--flat_sub_rate`
    (5×10⁻³) and 10× below the 3.1×10⁻³ calibrated on SC2200627
    (`dev/error_rate_calibration.md`).
- **Chimeras.** DADA2 flagged 36% of reads as chimeric. The OTU path keeps chimeras.
- **Primers.** The primers AAP identified (`515-YF GTGYCAGCMGCCGCGGTAA`,
  `806BR GGACTACNVGGGTWTCTAAT`) are SR's defaults.
- **`.mseq` availability.** The rendered results have no `.mseq`. `PUBLISH_OTU_RESULTS`
  publishes one, and the benchmark's `RUN_AAP` keeps it
  (`modules/local/amplicon_analysis/main.nf:99`), so both cases occur.

### The database

`~/Downloads/silva-ssu-138.1.tar.gz` unpacks to `silva-ssu-138.1/138.1/`. The same layout is
at `/hps/nobackup/rdf/metagenomics/service-team/ref-dbs/amplicon_pipeline/SILVA-SSU/138.1`.

- `SILVA-SSU.fasta`: 2,098,433 sequences, 2.9 GB, one line per sequence.
  - This is SILVA 138.1 **SSU Ref, not NR99**, all three domains.
  - It uses the **RNA alphabet (`U`)**.
  - Headers are bare `acc.start.stop`.
- `SILVA-SSU-tax.txt`: a MAPseq header with 8 levels and per-level `#cutoff`. Lineages look
  like `sk__Bacteria;k__;p__…;g__Collinsella;s__Genus_species`.
- `SILVA-SSU.otu` and `SILVA-SSU.fasta.mscluster` are shipped.

In-silico PCR on a 0.1% random sample (2,071 sequences):

- As shipped, `subspecies_infer.extract_v4` amplifies **none** of them. `U` is not in
  `_IUPAC`, so every `T` in a primer counts as a mismatch. Fixed by 1.1.
- After U→T it amplifies 98.6% of Bacteria (1,842/1,868), 100% of Archaea (46/46) and 54%
  of Eukaryota (84/157).

## How the merge changes M

`M[a, b]` = P(MAPseq label `b` | read from source `a`). Every step between template and
MAPseq is part of `M`. AAP's chain is: PCR → sequence each mate → quality-trim → fastp
merge (correct, or drop the pair) → cmsearch clip → MAPseq. The merge changes `M` in five
ways.

> **Phase 0.1 (2026-09-17, upstream `dev/aap_merge_effects.md`) corrects effects 1, 2 and 4.**
> - fastp checks its ≤ 5 mismatch limit only on the first 50 bases of the overlap (the
>   start of the fragment here). Past those, pairs merge with any number of substitutions
>   or an indel.
> - The merged read is R1. R1 substitutions and indels pass unchanged, and an R2 indel
>   becomes a 1-base length error at the 3′ end.
> - Only indels in fragment bases ~8–40 drop the pair.
> - Correction needs R1 ≤ Q14 against R2 ≥ Q30, i.e. only Q12 against Q38 in the binned
>   Nov2025 reads.
>
> The text below is the pre-measurement expectation.

1. **The error rate falls, but the merged read is not a consensus.**
   - Across the overlap, fastp's merged sequence takes R1's bases; rc(R2) contributes only
     beyond R1's end.
   - A mismatch is corrected only when one base is high quality and the other low (fastp's
     `BaseCorrector`: ≥ Q30 against ≤ Q14). An R1 error at moderate quality survives even
     when R2 is right. That is why only 9×10⁻⁵/base was corrected.
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
7. **cmsearch clipping** changes the top hit often, but mostly within one V4 group. See
   0.4.

What each mode needs as a result:

- **simulate**: each simulated read has to pass through the same steps as a real one.
  - That means either an error model of the *merged* read, or simulated mates run
    through the real fastp.
  - Either way the reads are mapped against the same database with the same MAPseq build
    and arguments.
- **align**: `c**d` stands for "d independent errors, each at the right base".
  - Under the merge, `c` becomes the merged-read residual rate (effect 1).
  - Indel edits need a much smaller weight than substitutions (effect 2), and position
    can matter (effect 4). All of these can come from one error profile measured on the
    reads.
  - Yield (effect 5) cannot come from distances at all.

## Decisions

1. **Map against the reference FASTA as shipped. That is the only MAPseq target.**
   - Observed, simulated, panel and home reads all use the same FASTA, tax file,
     `.mscluster`, MAPseq version and arguments.
   - Against AAP's directory the labels are then identical to AAP's, so its `.mseq` can be
     reused.
   - There is no `--mapseq_target` switch and no amplicon MAPseq database. The amplicon
     headers *are* the FASTA's headers, so a full-length hit has the id an amplicon hit
     would have. The id-keyed tallies (`observed_refseq_counts`, `tally_kernel`,
     `_panel_context`) already drop hits outside the amplifiable set.
   - Phase 0.6 put those hits at 0.29% of reads.
   - A missing `.mscluster` is built once and cached (F.4). AAP ships its own, so SR never
     clusters 2.1M sequences.
   - Amplicon extraction still runs. It provides V4 groups, sources, the translation table
     and the `align` input, but it is not the MAPseq target.
2. **V4 groups (`db_groups`) stay as the label space.** They are not scale machinery.
   - Against a full-length database, an error-free V4 read ties across every reference
     that shares its V4, and MAPseq breaks the tie arbitrarily.
   - 0.4 found 98.8% of top-hit changes stay inside one V4 group. Grouping is what makes
     those labels comparable.
   - The groups come from `amplicons.fasta` as sequences, never from the MAPseq FASTA.
3. **Use AAP's merged reads and do not merge again.** SR's `reads_to_fasta.merge_pair`
   (R1 wins, overlap ≥ 20, ≤ 25% mismatches, no qualities) accepts pairs fastp rejects,
   indel-bearing pairs in particular. Re-merging raw reads would give a different read set
   with a different `M`.
4. **The default simulation is a merged-read error model.**
   - skiver trains on merged reads reference-free, as SR already does for single-end
     `reads:`, with position covariates.
   - Simulated mates run through real fastp (`pairs`) are the validation reference, and the
     choice when yield binds.
5. **Inference is over a panel, and the panel is the only method** (Phase P).
   - SILVA is a label space (`silva_taxon_panel_plan.md` §1.1).
   - Only the panel's sources are simulated: its genome V4 amplicons, plus the SILVA V4
     groups under each taxon entry (≤ `--panel_taxon_max_sources`). Those reads are mapped
     against all 2.1M references.
   - The result is the rectangular kernel `K` (sources × SILVA labels). Reads on labels no
     source reaches go to `background`. No database reference is simulated unless it is a
     taxon-entry source.
   - With no `panel_references`/`panel_taxa`, the panel is the whole reference database:
     its sources are the database's distinct amplicons, and its entries are the database's
     genomes.
   - The square `M`, its backends, `--mismapping_matrix` and `--infer_space v4_group` go.
   - **The whole-database panel is never tested or swept against full SILVA-SSU or GTDB.**
     It is tested only on small reference sets: the bundled fixture, the 81-reference
     *B. uniformis* set and the 20HM database. Against a generic database it would simulate
     every distinct V4 group and fit a genome space V4 cannot identify. The AAP phases
     below always name a panel.
6. **Primer trimming follows the read source, not the MAPseq target.**
   - SR's own reads are trimmed or not per `--trim_primers`. Observed and simulated reads
     are prepared the same way, so either works against the full-length FASTA.
   - `merged: true` rows (AAP inputs, 1.6) are never trimmed, because an AAP `.mseq` is of
     untrimmed reads. `--trim_primers true` with a `merged: true` row that supplies `mseq:`
     is an error.
7. **Order: reconcile the upstream lines first (R), then full-length mapping (F), then
   finish P.**
   - This writes F once, on `main`, instead of on `silva-fixes-pinned` followed by a
     conflicting port.
   - Cost: the SILVA sweep's square `silva` arm stays broken until P.8 retires it. If that
     arm is needed sooner, do F on `silva-fixes-pinned` first, and treat R as porting F too.

## Phases

### Phase 0 — measure before coding (`dev/aap_merge_effects.{py,md}` upstream) — done

Run on the 20 Nov2025 samples.

- **0.1 fastp behaviour.** See the note under "How the merge changes M".
- **0.2 Residual error profile.** Aligned to the 20HM genomes' V4 amplicons.
  - Interior substitutions are 3.4×10⁻⁴/base (median column, stable across runs).
  - Sequencing indels are ~4×10⁻⁶ (~100× lower).
  - The 20-bp window max/min is ≤ 2.1, so 3.3 is not triggered and 3.2 is justified.
  - Output: `position_error_profile.tsv` (`pos sub ins del`).
- **0.3 Primer-site mix.**
  - The oligos are EMP 515F `GTGCCAGCMGCCGCGGTAA` / 806R `GGACTACHVGGGTWTCTAAT` (Y reads
    C, N never G).
  - Heterogeneity spacers are 0/4–7 bp before 515F and 0–3 bp after 806R.
  - The mix does not depend on the source genome. The reverse mix depends on overhang
    length and its sites co-vary, so 2.1 samples the joint `ends.tsv`.
  - Output: `primer_mix.tsv`.
- **0.4 cmsearch clip.**
  - Clipping changes the MAPseq top hit for 19.2% of reads (replicate: 0%).
  - 98.8% of those changes stay inside one V4 group, so V4-group agreement is
    99.77–99.81%, below the 99.9% threshold.
  - So: emulate the clip, and compare labels as V4 groups.
  - The Nov2025 run mapped against the older `SSU.fasta`/`slv_ssu_filtered2.txt`, not
    `silva-ssu-138.1`.
- **0.5 MAPseq build and format.** Measured on all 517k clipped reads.
  - With AAP's arguments, AAP's and SR's mapseq builds give identical output.
  - With SR's old default `mapseq_args` ('') the top hit differs for 20.2% of reads and the
    V4 group for 2.8%. AAP's `-seed 12 -tophits 80 -topotus 40` is therefore the default
    (`4d8f665`, ported in R.1).
  - Columns 1, 2 and 4 are the same in the simple and confidences formats.
- **0.6 Label coverage.**
  - 0.29% of reads pooled (at most 0.56% per run) hit a reference without an amplicon.
    That is below the 1% trigger, so 2.4's relaxed extraction is not needed.
  - 99.96% of those reads hit references with a one-base indel in the primer site, not
    truncated ones. Three references hold 98%.

### Phase P — panel-only pipeline

Before P, a sample ran one of two paths:
- **square:** `MATRIX_KEY` → `SIMULATE_READS`/`MAPSEQ_SIM`/`BUILD_MISMAPPING`, or
  `ALIGN_MISMAPPING`/`GROUPED_MISMAPPING`/`MINIMAP2_*` → `PUBLISH_MISMAPPING` →
  `infer_composition.run`;
- **panel:** `PANEL_PREPARE` → `MAPSEQ_PANEL_HOME` + `SIMULATE_PANEL_READS`/
  `MAPSEQ_PANEL_SIM` → `PANEL_KERNEL` → `infer_composition.run_panel`.

Phase P keeps only the panel path.

- **P.0 Snapshot before deleting.** *Done 2026-09-17 (upstream `dev/panel_only_parity/`).*
  - The bundled fixture and the *B. uniformis* 81-reference set, each under `simulate`,
    `align exact-hash` and `align kmer tau 1`.
  - The 20HM genome panel.
  - One species taxon panel from `dev/panel_silva_sweep.md`, re-scored from its stored
    `.mseq` and not re-run against SILVA.
  - These are the P.7 references.
- **P.1 The default panel.** *Done 2026-09-17.*
  - `build_panel_kernel.py prepare --whole-database <amplicon_dir>` reuses the extracted
    amplicons and translation table. It writes one `v4g_<sha16>` source per distinct
    amplicon.
  - A row with no panel carries `meta.panel = 'database'`.
- **P.2 Workflow.** *Done 2026-09-17.*
  - The square modules and branch are gone, and `PANEL_ALIGN` wires
    `build_panel_kernel.py align`.
  - `MAPSEQ_CLUSTER` now clusters the reference FASTA when it has no `.mscluster` beside
    it (F.4), not the amplicons.
  - The kernel group key is (database amplicons sha256, panel or `database`, panel taxa,
    model identity, method and settings), plus the MAPseq target identity (F.5). Kernels
    publish to `mismapping/panel_<key>/`.
- **P.3 Align on the whole-database panel.** *Done 2026-09-17.*
  - Uses a literal hash join at tau 0 and source-probed pigeonhole candidates with bounded
    IUPAC verification at tau ≥ 1.
  - Removed: `--align_backend`, `--minimap2_args`, `--minimap2_index_args`.
  - Kept: `--align_tau`, `--align_distance_decay` (including `auto`),
    `--align_decay_model`, `--align_ambiguity_weight`, `--max_ambiguous_bases`,
    `--max_postings`.

### Phase R — reconcile `silva-fixes-pinned` onto `main` (upstream)

- **R.1 Port onto `main`,** as one PR of cherry-picks:
  - `d448dfa`: vectorised extraction. It matters more now: SSU Ref is 2.1M sequences.
  - `6c9390a`: `--amplicon_cache`. It keys extraction, and the clustering of built sets,
    by reference set + primers + extractor code.
  - `ac9328e`'s `modules/local/mapseq/cluster` change: the one-record query. Its
    `dev/mapseq_free_labels.*` go too, as a record.
  - `4d8f665`: AAP's MAPseq flags by default.
- **R.2 Drop:**
  - `105f484`, since #14 already has U→T;
  - `c1f2031` (`--align_home_probes`), whose grouped square matrix is gone. `PANEL_ALIGN`
    already maps its sources home against the database;
  - `eca1763` (prebuilt amplicon database, `--build_mapseq_db`), which F supersedes.
- **R.3** Delete the `silva-fixes-pinned` branch once the benchmark no longer pins it
  (P.8).

**Acceptance:** `main`'s stub suite and the fixture e2e pass. `--amplicon_cache` reuses an
extraction across two runs.

### Phase F — the reference FASTA is the MAPseq database (upstream `main`)

- **F.1 `ch_db` comes from the reference FASTA** (`workflows/superresolution_amplicon.nf`). *Done, with F.2 and F.4–F.8 (upstream `c323bc1`).*
  - It becomes `[ id, refs_fasta, refs_tax, refs_fasta.mscluster ]`.
  - MAPSEQ_OBS, MAPSEQ_PANEL_HOME and MAPSEQ_PANEL_SIM use it.
  - The amplicon directory still feeds simulation sources, `align`, PANEL_PREPARE,
    `db_groups` and the translation table.
- **F.2 Fix the two places that use the MAPseq FASTA as sequences.**
  - PANEL_KERNEL and PANEL_ALIGN are handed `ch_panel_db`'s `fasta` as the database
    amplicons (the `.join(ch_panel_db.map { … [ meta.id, fasta ] })` lines).
  - Under F.1 that would be full-length, and `db_groups` would group whole SSU sequences.
  - Pass the set's `amplicons.fasta` instead. This is the one correctness trap.
- **F.3 MAPseq always gets a generated tax file.** *Done (upstream `c323bc1`).*
  - EXTRACT_AMPLICONS writes `references.tax` over *all* headers, and every MAPseq task
    uses it. No tax file is looked up beside the FASTA.
  - Measured: 400/400 SILVA queries gave identical top hits and identities under AAP's
    `SILVA-SSU-tax.txt` and under a generated one. The tax file changes only MAPseq's
    taxonomy columns, which nothing reads.
  - Lineages come from `--taxonomy` alone (F.6).
- **F.4 The `.mscluster` is optional and cached.**
  - If `<fasta>.mscluster` exists it is staged beside the FASTA.
  - Otherwise MAPSEQ_CLUSTER clusters the full-length FASTA once, under `--amplicon_cache`.
    That takes seconds for built sets.
  - Log a warning above 100k records.
  - ponytail: an unclustered 2.1M-sequence FASTA clusters silently for hours. Refuse it
    only if that bites.
- **F.5 Keys and provenance.**
  - The panel kernel key and `provenance.json` add the MAPseq FASTA's and `.mscluster`'s
    identity (path, size, mtime, like the set id).
  - `db_amplicons_sha256` stays: it identifies the label space, which is unchanged.
- **F.6 Read AAP-form lineages in `--taxonomy`.** This replaces the old 1.2
  `--aap-tax` converter.
  - `build_panel_kernel._ranks` and the `lca` reader strip `xx__` prefixes, write
    `s__Genus_species` as `Genus species`, and treat an empty rank as `unclassified`.
  - AAP's `SILVA-SSU-tax.txt` then serves as both the MAPseq tax and `--taxonomy`, and
    panel `taxon:` entries keep bare names.
- **F.7 Bare-accession headers** (the old 1.4). With no `|`, `genome_of_header` returns the
  whole accession, so each sequence is its own genome. Add one assert that pins it.
- **F.8 Docs and tests.**
  - README: "The MAPseq database is the reference FASTA", with the directory layout: FASTA,
    tax, optional `.mscluster`.
  - Stub tests: a set with a `.mscluster` beside it, a set without one, and a panel set
    for F.2.
  - pytest for F.6 and F.7.

**Acceptance:**
- The fixture e2e passes.
- On one Nov2025 sample, `obs.mseq.gz` top hits equal a direct AAP-flag MAPseq run against
  SILVA-SSU/138.1 on the same reads. This was the old 1.3 acceptance.
- On the fixture and 20HM, panel-entry TV between the pre-F and post-F runs is ≤ 0.02. Any
  gap is explained by dropped non-amplifiable hits.

### Phase P, continued

- **P.4 A supplied kernel replaces `--mismapping_matrix`.** *Done (upstream `9177cb7`).*
  - `--panel_kernel <dir>` takes a published `mismapping/panel_<key>/` bundle
    (`mismapping_matrix.npz`, `panel_translation.tsv`, `sources.tsv`, `provenance.json`)
    and skips kernel building.
  - It is refused unless the bundle's `db_amplicons_sha256` and MAPseq target identity
    (F.5) match the sample's, and its panel identity matches the sample's.
  - This keeps the benchmark's build-once / infer-many split and cheap `mseq:` sweeps.
  - As built: PANEL_KERNEL/PANEL_ALIGN write the group's provenance into the bundle as
    `provenance.json`. A supplied bundle is checked on `reference_sha256`, `mapseq_db`,
    `panel` and `panel_taxa`, all compared as strings. Its build settings are not
    checked. MATRIX_KEY still runs, and nothing trains.
  - The panel and MAPseq identities are paths, as in the key. The benchmark must hand the
    build run and the inference runs the same absolute paths (P.8).
- **P.5 Inference and fit check.** *Done (upstream `panel-only-inference`).*
  - `infer_composition.run` becomes `run_panel`.
  - Delete the square code: `_mismapping_matrix`, `--build-mismapping`,
    `--build-grouped`, `--active-amplicons`, `--infer-space`, `--infer-prune`. Delete the
    square branch of `check_composition_fit.py`.
  - Carried over, and confirmed each with a test:
    - presence gate, horseshoe, latent distance decay (strata from `align --tau >= 1`);
    - `--no-mismapping`;
    - `low_depth` and `no_reference_hits`;
    - the `resolution` column.
  - Prune only if P.7 shows the whole-database panel is too slow on 20HM. Pruning would
    drop sources whose row reaches no observed label, with the pruned mass going to the
    zero-count sink.
  - `--sim_n_per_ref` and `--panel_sim_n_per_ref` merge into one `--sim_n_per_ref`,
    default 5,000.
  - As built:
    - `infer_composition.py` takes `--amplicon-dir` (the prepare dir),
      `--mismapping-matrix` (the kernel) and `--obs-mseq`, all required. `--sim-mseq`,
      `--no-prune` and `--demo` went with the square flags; `--infer_space` and
      `--infer_prune` left the pipeline, and so did `inferred_v4_groups.csv`.
    - The fit check has one statistic, over the fitted labels, under `label`
      (`gate: "label"`); `raw_reference`/`unique_v4_group` are gone. The sweep report
      follows `gate`, so it reads either form.
    - Carried-over tests are in `tests/test_panel_inference.py` (presence gate,
      horseshoe, `--no-mismapping`, low depth, no hits); resolution and the latent decay
      were already covered in `test_measured_mismapping.py`. `test_zero_hits.py` folded in.
    - Test code that only exercised deleted inference paths went now, not in P.6:
      the square tests in `test_composition_fit.py`, the three `infer_composition.run`
      tests in `test_grouped_mismapping.py`, and `--build-grouped` in
      `test_measured_mismapping.py`.
    - Not removed: `--taxonomy` / `lca_composition.csv`. A panel run passes no lineages,
      so inference never writes that table now. Drop it in P.6 or wire panel lineages.
    - The benchmark's `tests/sr_fast.config` `sim_n_per_ref = 50` now reaches the panel
      simulation, which `panel_sim_n_per_ref` used to hide.
- **P.6 Delete what only the square path used.** *Done (upstream `f18b1d7`).*
  - `bin/select_active_amplicons.py`;
  - `subspecies_infer.build_mismapping`/`build_mismapping_grouped`;
  - the grouped and CSR readers in `sparse_matrix.py` (keep `read_kernel`/`write_kernel`);
  - legacy CSV matrices;
  - `tests/test_grouped_mismapping.py`, and the square parts of
    `tests/test_measured_mismapping.py` and `tests/test_composition_fit.py`;
  - in `build_mismapping_align.py`, everything except what `build_panel_kernel.py`
    imports. Move what is kept into `kernel_align.py` if the file is left mostly empty.

  `tally_kernel` stays. `dev/` scripts that import deleted functions stay as a record,
  marked unmaintained.

  As built:
  - `kernel_align.py` holds what `build_panel_kernel.py align` used from
    `build_mismapping_align.py` (pigeonhole candidates, bounded IUPAC distance,
    `measure_error_rate`, `ambiguity_weights`); the rest of that file is gone.
  - Also deleted: `write_mismapping_bundle.py` (no caller) and
    `subspecies_infer.observed_refseq_counts`. The README lost the align-backend and
    minimap2 rows P.3 had left behind.
  - The MATRIX_KEY `kernel_version` guard moved to `test_measured_mismapping.py`.
  - Left for later: the torch model (`subspecies_infer._apply_mismapping`,
    `DecayKernel`) still accepts dense and grouped `M`, and `subspecies_infer.py --demo`
    exercises those forms. Only the rectangular form is reached from the pipeline.
- **P.7 Tests.** *Done (upstream `fdc5e06`–`b0f6d1e`).*
  - Stub nf-test: a row with no panel, a genome panel, and a taxon panel, each under
    `simulate` and `align`, plus `--panel_kernel` reuse.
  - Real e2e on the bundled fixture: the whole-database panel under `simulate`.
  - pytest:
    - `prepare --whole-database` sources equal the database groups;
    - pigeonhole `align` equals brute force on *B. uniformis*;
    - a `--panel_kernel` sha mismatch is refused.
  - Parity against P.0:
    - genome TV ≤ 0.02 on the fixture and *B. uniformis* in each mode;
    - the 20HM panel and the taxon-panel re-score match their snapshots up to F's dropped
      non-amplifiable hits.

    Any larger difference is explained in `dev/panel_only_parity.md`. The expected sources
    are the `background` row and home labels.
  - Nothing runs the whole-database panel against full SILVA-SSU or GTDB.
  - As built (upstream `dev/panel_only_parity.md`):
    - Stub: align-mode tests for a genome panel and taxon panels were added; the rest
      existed. The `--panel_kernel` sha refusal is a workflow check, covered by its stub
      test rather than pytest.
    - e2e: `tests/default.nf.test` gained its first `e2e` test, the fixture on the
      whole-database panel under `simulate`. It passes.
    - pytest: pigeonhole `align` equals brute force on *B. uniformis* at tau 1 and 2.
      `prepare --whole-database` was already covered.
    - Parity: fixture TV 0.013 in all three modes. The 20HM (S01–S20) and SILVA species
      (S05) re-scores are bit-identical, fit `ok`.
    - **Open: *B. uniformis* misses the bar at every depth.**
      - At 1,175 reads over 11 V4-group labels, the default presence gate collapses:
        `conc_frac` drops to 0.014 and nothing is called present. TV to P.0 is 0.24.
      - A depth sweep (`dev/panel_only_parity.py depth`) refits the same observation
        at 1x to 128x its reads with the P.0 square code and the panel code. The gate
        recovers by 8x (10,000 reads).
      - The P.0-to-panel gap then plateaus at TV 0.06–0.07 (V4-inseparable pairs
        merged), or 0.11 unmerged, up to 160,000 reads.
      - The gap is P.0's error. The observation's truth is known, and at 128x the panel
        is 0.10 from it against P.0's 0.16. P.0 halves *D. formicigenerans*; the rest
        is the prior-decided *B. uniformis* strain split.
      - **Decided 2026-09-21:** `--infer_presence` stays on by default. Test inference
        at >= 10,000 reads from now on. The 1,175-read result measures the low-depth
        failure, not the code. At 10k+ reads the remaining gap is P.0's error, so P.7 is
        closed.
- **P.8 Benchmark contract (this repo).** This lands with the `sr_revision` bump to `main`. *Done on branch `panel-only-contract`.*
  - **Delete the prebuilt-MAPseq-database machinery** (b57f775):
    - `srPrebuiltMapseqDb` and its check (`subworkflows/local/profile/main.nf`);
    - `srBuildDb` / `meta.sr_build_db`;
    - `srBuildDbArgs` and its three call sites (`modules/local/superresolution/run/main.nf`);
    - the `prebuilt` 5th element of `ch_sr_dbs` (the BUILD_DATABASES emit and PROFILE's
      tuple);
    - `tests/samplesheets/sr_prebuilt_no_mapseq_db.yaml` and its test;
    - the `silvaNoMapseq` and `silvaLoose/*_amplicons/` fixtures.
  - **Accept AAP's directory as is.** Add `'*-tax.txt'` to `taxPats`
    (`subworkflows/local/build_databases/main.nf`). The FASTA pattern already matches
    `SILVA-SSU.fasta`. The nested run finds the `.mscluster` beside `realpath refs` itself.
  - **Keep `srCacheArgs`.** Extraction is the costly per-set step now.
  - **Kernel bundle.** `BUILD_SUPERRESOLUTION_MISMAPPING` lifts the `mismapping/panel_<key>/`
    directory instead of `mismapping_matrix.npz`, and `RUN_SUPERRESOLUTION` passes
    `--panel_kernel`. The "exactly one bundle per run" rule stays.
  - **Panel `sr_settings` stop being a special case:** every amplicon run builds a kernel
    first. `srSetKey` keeps the panel in the key.
  - **`sr_amplicon_matrix_args`** loses the backend and minimap2 flags.
  - **`bin/normalize_sr_profile.py`** drops the `background` row and renormalises. It
    writes the background share to a sidecar, so the profile is still over genomes and the
    unexplained fraction is still reported.
  - **`examples/sr_amplicon_silva_sweep`.** Its `silva` arm (the square `v4_group` grid over
    full SILVA) goes. Keep its `custom`, `silva_panel` and `aap` arms, or fold them into
    `sr_amplicon_panel_sweep`, whichever is less duplication. The `database.path:` docs
    become "FASTA + tax (+ `.mscluster`)", and AAP's SILVA-SSU/138.1 works unmodified.
  - **Other docs:** update `examples/sr_amplicon_param_sweep`, `subspecies_v4_sweep`,
    README and `CLAUDE.md` (the prebuilt-MAPseq bullet goes, and the amplicon-cache
    bullet shrinks) for the removed flags. The shotgun profiler is untouched.
  - As built:
    - Every `sr_amplicon` set, panel settings included, builds its kernel in
      BUILD_SUPERRESOLUTION_MISMAPPING, which is handed the kernel flags plus the panel
      flags (`--panel_references`, `--panel_taxa`, `--taxonomy`). The bundle is lifted
      whole as `<id>.panel_kernel/` and published to
      `mismapping/<set>/panel_kernel/`.
    - Inference gets `--panel_kernel` plus the same panel flags: the upstream check
      compares them as absolute paths. RUN_SUPERRESOLUTION no longer publishes kernels.
    - `align_backend` left `sr_settings` and the params, and so did the
      "panel kernel records no distances" refusal (untrue since P.3).
    - `silvaLoose` became `silvaAap`: AAP's file names, no `_amplicons/`. A stub run
      against the real `~/Downloads/silva-ssu-138.1/138.1` with a genome panel passes;
      that stands in for the `-preview` acceptance.
    - Sidecar: `<id>.sr_background.tsv` (`background_fraction`). The panel sweep's scorer
      reads it; the SILVA sweep's reads `inferred_composition.csv` and needed nothing.
    - The panel sweep's mapping run and the SILVA sweep's phase 1 now carry `panel:`,
      because a panel-less run against SILVA would build a whole-SILVA kernel
      (Decision 5).
    - The SILVA sweep kept `custom`, `silva_panel` and `aap`, and its grid now fans over
      `silva_panel`. Folding into the panel sweep would have meant porting the `aap` and
      error-model arms.
    - The three SR examples pin upstream `e29bf693e5952c8584d26497f7c202227cb78d52`,
      #21's merge commit on `main`.
    - Not run: the benchmark e2e, which pulls the nested pipeline from GitHub.

**Acceptance (P):**
- P.7 passes.
- `git grep -E "build_mismapping_grouped|GROUPED_MISMAPPING|infer_space|mismapping_matrix\b|build_mapseq_db|_amplicons/"`
  upstream finds only `dev/` and the kernel file name.
- The benchmark's stub tests pass, and its e2e passes on the fixture panel sweep.
- `nextflow run main.nf -preview` on a SILVA panel samplesheet whose `path:` is AAP's
  directory gets past PROFILE.

### Phase 1 — read AAP's outputs

- **1.1 U→T in `subspecies_infer.read_fasta`.** *Done (#14).*
- **1.2 and 1.4** moved to F.6 and F.7.
- **1.3** replaced by Phase F. There is no `--mapseq_target`. Its primer rule is now
  Decision 6 and 1.6.
- **1.5 `bin/aap_samplesheet.py <aap_outdir> [--db-label SILVA-SSU]`.**
  - It writes SR samplesheet rows `{id, reads: <id>/qc/<id>.merged.fastq.gz, merged: true,
    merge_rate, mseq: <id>/taxonomy-summary/SILVA-SSU/<id>_SILVA-SSU.mseq(.gz)}`, with
    `mseq` only when the file exists.
  - It skips runs missing from `qc_passed_runs.csv`.
  - It refuses a run whose `amplified-region-inference/<id>.tsv` is not a single 16S V4.
  - It refuses a run whose `primer-identification/*.fasta` differs from
    `--fwd_primer/--rev_primer` (primers are global in SR).
  - It takes `merge_rate` from `fastp.json` (`merged_and_filtered.total_reads` over input
    pairs) and warns below 80%, the threshold `reads_to_fasta.py` uses.
  - It needs `--panel` or `--panel-taxa` and writes them into every row (Decision 5).
- **1.6 `merged: true` in the samplesheet.**
  - READS_TO_FASTA never merges such reads and never trims their primers.
  - `obs_max_reads` still applies.
  - `--trim_primers true` with a `merged: true` row carrying `mseq:` is an error in
    `main.nf` (Decision 6).

Tests: a stub nf-test on a miniature AAP outdir fixture (two runs, one with `.mseq`), and
pytest for `aap_samplesheet.py`.

As built (upstream #22, `34effd7`), three changes from the above:
- **Primers are checked for compatibility, not equality.** The 20 Nov2025 runs report four
  primer pairs: fwd `GTGYCAGCMGCCGCGGTAA` or `CCAGCAGCCGCGGTAATACG` (the same site, shifted
  3 bases), rev with `N` or `H` at one position. Equality would refuse 18 of 20. A run is
  refused only if its primer cannot overlap SR's by ≥ 15 IUPAC-compatible positions within a
  4-base shift. A differing but compatible primer is logged. The ends of those reads differ
  by a few bases, which is Phase 2.1's primer-flank work.
- **The `.mseq` is found by glob.** AAP writes `taxonomy-summary/<label>/<id>.mseq`, not
  `<id>_<label>.mseq`. The script takes the single `*.mseq*` under the label directory.
- **Every `merged: true` row needs `--trim_primers false`, not only rows with `mseq:`.**
  `--trim_primers` also trims the simulated reads, run-wide. An untrimmed merged row
  against trimmed simulated reads is the same mismatch with or without a supplied `.mseq`.
  With that enforced, READS_TO_FASTA needs no per-row switch. `merge_rate` is written to
  the row but not read yet (3.4 will).

All 20 Nov2025 runs pass the script. Stub nf-tests (21) and pytest (43) pass.

Acceptance run (2026-09-21, upstream `main` at `462603b`, in upstream
`work/phase1_acceptance/`). It passed.
- **Input.** Sample `SC2189280-SC3-1-26s000344`. The rendered AAP results have no `.mseq`, so
  a one-run AAP outdir was staged with Phase 0.4's `work/aap_clip/clipped.mseq` as
  `taxonomy-summary/SILVA-SSU/<id>.mseq`. That file is AAP's MAPseq build and arguments on
  the cmsearch-clipped reads against SILVA-SSU/138.1, i.e. what AAP would write. The panel was
  the 20HM `community_20hm_ssu.sr_refs.fasta` (95 SSU copies, 22 genomes).
- **Run.** `-profile docker`, `--references` AAP's `SILVA-SSU.fasta` unmodified,
  `--trim_primers false`, flat model at 0.2's rates (sub 3.4×10⁻⁴, ins/del 4×10⁻⁶),
  `--amplicon_cache`.
  - It needs `process.resourceLimits = [memory: 18.GB]` on a 24 GB machine, because
    EXTRACT_AMPLICONS asks for 36 GB.
  - EXTRACT_AMPLICONS took 1 min 21 s on the 2.1M sequences, not hours. Its peak RSS was
    24.6 GB, summed over its worker processes.
  - The whole run took 5 min 24 s.
- **Result.** READS_TO_FASTA and MAPSEQ_OBS did not run, so AAP's `.mseq` was reused without
  re-mapping. Inference used 514,626 reads with fit `ok` and 16 of 22 genomes present.
  - `obs_max_reads` does not cap a supplied `.mseq`: all of its reads were used.
  - `background` is 0.79 of the profile. 225 of 295 observed labels are reached by no
    panel source. This sample is mostly reads outside the 20HM panel (Phase 0.4's top
    label changes are within Clostridiaceae). It is a property of the sample, not of the
    Phase 1 plumbing, and Phase 2's fit check is the place to judge it.

**Acceptance:** on one Nov2025 sample, an SR run from `aap_samplesheet.py` output with a
20HM genome panel completes, and it reuses AAP's `.mseq` without re-mapping.

### Phase 2 — simulate mode for merged reads

- **2.1 Primer flank from the observed mix.**
  - `simulate_amplicon_reads.py --primer-mix primer_mix.tsv` draws each degenerate primer
    position from the 0.3 frequencies, sampling the joint `ends.tsv` for the reverse end.
  - Without the file it draws uniformly over the code's options, never the first option.
  - For `merged: true` samples the primers are not trimmed after errors are applied.
  - As built (upstream `50b86e6`):
    - Before this, untrimmed simulated reads were bare amplicons: no primers at all, while
      the observed merged reads keep both.
    - The table is `--primer_mix` (`fwd rev reads`, concrete oligos, rev on its own
      strand), not `ends.tsv`. `assets/primer_mix_emp_v4.tsv` is `ends.tsv`'s pooled
      rows summed over spacer and overhang, which keeps the joint reverse triple. The
      0.15% of reads with an off-code base were dropped: the error model supplies those.
    - Each oligo must be an instance of `--fwd_primer`/`--rev_primer`, else refused.
    - Spacers and the overhang are not simulated. The clip removes most of both (0.4).
    - Without a table, untrimmed reads draw uniformly over each code, with a warning.
      Trimmed runs are unchanged unless a table is given.
    - `primer_mix` joins the kernel key as a path.
- **2.2 Merged-read error model** (default for `merged: true`).
  - `--sim_error_model trained` trains on the merged reads.
  - Add `AdditiveContext(7)+Position(2)` and `AdditiveContext(9)+Position(2)` to the
    candidates, so the 3′-weighted residual of effect 4 can be learned.
  - Default to `--trained_error_model_scope pooled` for AAP batches.
  - Under `flat`, the README gives 0.2's pooled sub/ins/del rates as the AAP starting
    point, marked as per-run numbers.
  - As built (upstream #24):
    - Training on merged reads needed no change: `merged: true` rows reach
      TRAIN_ERROR_MODEL like any single-end reads.
    - `--trained_error_model_scope` defaults to null. `main.nf` resolves it to `pooled`
      when any row is `merged: true`, else `per-sample`.
    - `--sim_error_model` stays `flat` by default. The acceptance comparison below decides
      whether `trained` should be the default for merged rows.
    - **The `Position(2)` candidates were not added.** skiver's `Position(N)` is not a
      per-base term: training sums `read_pos`/`dist_to_end` per context and fits
      `mean_pos(context) @ position_weights`. The pinned `apply_batch` also refuses any
      model carrying `position_weights`, so an AIC win would crash SIMULATE_PANEL_READS.
      0.2 put interior variation at ≤ 2.1× (under 3.3's 3× trigger), and the end pile-ups
      were primer degeneracy (2.1) and overhang (clipped). Add a per-base generative
      position term to skiver only if the fit check or 2.3's `merged`-vs-`pairs`
      comparison shows a 3′ bias.
- **2.3 `--sim_read_structure merged|pairs`,** with `pairs` for validation and for runs
  where yield binds.
  - A new `bin/simulate_amplicon_pairs.py`:
    - draws the primer-flanked fragment;
    - reads `--sim_mate_len` cycles of R1 forward and of R2 on the reverse strand
      (`apply_batch(..., is_forward=False, emit_quality=True)`), continuing into the
      TruSeq adapter;
    - writes paired FASTQ.
  - A new `FASTP_MERGE` module (fastp 1.0.1 biocontainer), whose `ext.args` are AAP's
    READS_QC_MERGE arguments verbatim. The merged reads go to MAPSEQ_PANEL_SIM.
  - Check that the source-carrying read name survives fastp's ` merged_x_y` suffix as
    `iter_mseq`'s first token.
  - Per-source yield (merged over simulated) goes to `yield.tsv` in the bundle.
  - The model is a per-mate skiver model trained on raw R1/R2
    (`AdditiveContext(7)+Position(2)+Strand`), from `error_model:` or a new `raw_reads:`
    key.
  - SR's own `--paired` merge stays for non-AAP inputs. Its divergence from fastp is
    documented, not fixed.
  - As built (upstream #25):
    - `simulate_amplicon_reads.py --mate-len` instead of a new script: it already holds the
      flank, the primer mix and both samplers. R2 is read off rc(fragment) with
      `is_forward=False`, and each mate runs into its TruSeq adapter, as in
      `dev/aap_merge_effects.py`. Flat mates are Q38 throughout.
    - FASTP_MERGE writes the merged reads as FASTA named by their first token, which drops
      the ` merged_x_y` suffix, so `iter_mseq` sees the source name. It also writes
      `yield.tsv` (`source simulated merged yield`), which PANEL_KERNEL copies into the
      bundle. Nothing reads it yet (3.4).
    - `pairs` needs `--trim_primers false`. Under `trained`, every row needs `error_model`,
      because a model trained on merged reads is not a mate model. **`raw_reads:` was not
      added**, and the plan's mate model can't have `Position` (2.2).
    - No cmsearch clip: the mates carry no spacer or overhang (2.1), so a merged read
      already spans primer to primer.
    - Without skiver's Phred calibration, trained qualities come from the context error
      rate, so fastp's Q ≤ 14 against Q ≥ 30 correction rarely fires. It was 9×10⁻⁵ per
      base in the real run.
    - The read structure and mate length join the kernel key and `provenance.json`.
    - A docker run on the upstream AAP fixture completes. At the default flat rates,
      90–94% of simulated pairs merge.
- **2.4 Labels outside the amplifiable set.** Not needed: 0.6 is below its 1% trigger.
  - Under F such hits are dropped as off-target.
  - If a later batch exceeds 1%, each such reference becomes its own `nonamp_<acc>` label
    column in `build_panel_kernel.py`.
- **2.5 What gets simulated.**
  - Only the panel sources, at `--sim_n_per_ref` (5,000), mapped against all 2.1M
    references.
  - Cost: sources × 5,000 reads. A 20-genome panel with a few species entries is ~100–500
    sources, i.e. the MAPseq work of 1–5 AAP samples.
  - One kernel serves every sample on the same panel, database and model.

**Acceptance** (Nov2025 batch plus the Phase 4 truth sample):
- `merged` against `pairs` on the same panel: median row L1 of `K` ≤ 0.05. If not,
  `pairs` becomes the default and 2.2 is documented as insufficient.
- Predicted mass on labels that drew no observed reads ≤ 0.3% (GTDB calibration reached
  0.24%).
- The panel fit check is `ok` on ≥ 18 of the 20 Nov2025 samples. Where it is not, check
  `unexplained_fraction` and whether the misfit sits on labels the panel reaches (chimeras
  between panel members, see Risks) before tuning anything else.

### Phase 3 — align mode for merged reads

`align` means `PANEL_ALIGN` (`build_panel_kernel.py align`, pigeonhole candidates).
- It writes `K[s, l] ∝ w(l)·c^d(s, l)` over the SILVA groups within tau of each panel
  source, plus one MAPseq home mapping of the sources against the full-length database.
- Only panel sources have rows.

- **3.1 Decay from the merged reads.**
  - `--align_distance_decay auto --align_decay_model <merged-read model>` already measures
    `c`.
  - Add a test that `auto` on the 2.2 model lands within 2× of 0.2's interior substitution
    rate.
  - Add a README note that for AAP inputs the model must be the merged-read one.
- **3.2 Separate indel decay.** `--align_indel_decay`, defaulting to
  `--align_distance_decay`.
  - The IUPAC-aware bounded edit distance returns `(subs, indels)`, and the weight becomes
    `c_sub**s * c_indel**i`.
  - `auto` measures `c_indel` from the model's indel rate.
  - The reason is effect 2: one-indel neighbours of merged reads leak at the PCR indel rate.
- **3.3 Position-aware decay.** Not triggered: 0.2's max/min is ≤ 2.1, against a 3×
  trigger. Skip it.
- **3.4 Refuse `align` when yield binds.** If any sample's `merge_rate` (1.5) is below
  0.8, `align` errors: use `simulate --sim_read_structure pairs`.
- **3.5 `--infer_distance_decay`** fits `c_sub` only; `c_indel` stays at its built value.
  Document it.

**Acceptance:** on the Phase 2 panel, align (tau 1, merged-read `auto` decay, indel decay)
against simulate `merged`:
- mean row L1 of `K` ≤ 0.08;
- panel-entry TV ≤ 0.02 on every sample.

### Phase 4 — validation against truth (this repo)

- **4.1 A synthetic sample.** Paired-end amplicon (`paired_end: true`,
  `read_length_mean: 310`, V4), error model trained on a raw Nov2025 run,
  `profilers: [aap, sr_amplicon]`. Its genomes come from SILVA-represented taxa present in
  the fermentor runs.
- **4.2 sr_amplicon reads RUN_AAP's output.**
  - PROFILE passes `profiling/aap/<id>/qc/<id>.merged.fastq.gz` as `reads` with
    `merged: true`, and the `.mseq.gz` as `mseq:`.
  - The DB is a `databases:` `path:` entry pointing at the unpacked AAP SILVA directory,
    unmodified.
  - Worked example: `examples/sr_amplicon_aap/`.
- **4.3 Score** panel-entry (genome and species) TV to truth, one arm each:
  1. raw AAP labels;
  2. simulate `merged`;
  3. simulate `pairs`;
  4. align (Phase 3);
  5. SR from raw reads with its own merge and read prep, against the same full-length
     database. This measures what Decision 3 buys.

**Acceptance:**
- Arm 2 beats arm 1 on species TV.
- Arm 3 is not better than arm 2 by more than 0.01; otherwise `pairs` becomes the default.
- The Nov2025 comparison against exact-matched DADA2 ASV profiles is descriptive only.

### Phase 5 — documentation

- Upstream README: an "Inputs from amplicon-analysis-pipeline" section covering the
  samplesheet helper, the full-length MAPseq database, the merged-read error model and the
  align decays.
- Upstream `docs/alignment_mismapping_plan.md` §7: merged reads satisfy the spanning
  assumption, subject to 3.4.
- This repo's `CLAUDE.md`: one line under Key implementation notes.

## Risks and open questions

- **Chimeras.** The OTU path keeps them: 36% of reads by DADA2's count in the example
  sample. A chimera on a label no panel source reaches goes to `background`, which is
  harmless. A chimera between two panel members is absorbed into a member. Phase 4's truth
  sample should include chimeras to size it.
- **SILVA Ref is not NR99.**
  - It has 2.1M sequences, 4× NR99, and its duplicate groups are larger.
  - Re-check `--panel_taxon_max_sources`.
  - Extraction is hours the first time, then cached (R.1). Publish the amplicon directory
    and reuse it.
- **Supplied `mseq:` across the F boundary.** A `.mseq` made against a pre-F amplicon
  database has the same ids but different tie-breaks.
  - The kernel key (F.5) keeps `-resume` and `--panel_kernel` from mixing them.
  - A user-supplied `mseq:` is not checked; document it.
- **fastp and MAPseq versions** must follow AAP's pins (fastp 1.0.1; mapseq per 0.5). An
  AAP upgrade changes `K`, so both versions go in the panel kernel id.
- **Non-V4 or 2×250 runs.** Nothing here has been measured on them. Yield (effect 5) is
  expected to dominate there, and the `pairs` path exists for that case.
- **The SILVA sweep's square arm is broken from now until P.8 retires it** (Decision 7).

## Files touched

Upstream:
- R: cherry-picks of `d448dfa`, `6c9390a`, `ac9328e` (cluster module), `4d8f665`.
- F: `workflows/superresolution_amplicon.nf`, `modules/local/{extract_amplicons,
  mapseq/cluster,panel_kernel/build}`, `bin/subspecies_infer.py` (tax writer),
  `bin/build_panel_kernel.py` and `bin/infer_composition.py` (lineage reader),
  `conf/modules.config`, `nextflow.config`, README, tests.
- P.4–P.7: `bin/infer_composition.py`, `bin/check_composition_fit.py`,
  `bin/sparse_matrix.py`, `bin/subspecies_infer.py`, `bin/build_mismapping_align.py`
  (square code removed); `bin/select_active_amplicons.py` and
  `tests/test_grouped_mismapping.py` deleted.
- 1–3: `bin/aap_samplesheet.py` (new), `bin/simulate_amplicon_reads.py` (primer mix),
  `bin/simulate_amplicon_pairs.py` (new), `bin/build_panel_kernel.py` (indel decay),
  `modules/local/fastp/merge/` (new), `main.nf`, `dev/aap_merge_effects.*`.

This repo:
- P.8: `subworkflows/local/{profile,build_databases}/main.nf`,
  `modules/local/superresolution/run/`, `bin/normalize_sr_profile.py`,
  `tests/default.nf.test`, the `tests/samplesheets` and `tests/data/prebuilt_db` fixtures,
  `examples/sr_amplicon_silva_sweep/` (square arm removed),
  `examples/sr_amplicon_param_sweep/`, `examples/subspecies_v4_sweep/`, `README.md`,
  `CLAUDE.md`.
- 4: `examples/sr_amplicon_aap/` (new).
