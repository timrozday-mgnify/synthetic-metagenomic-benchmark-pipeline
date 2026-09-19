# SILVA SSU as the generic database, and taxon panels

Status: 2026-09-19. **Complete.** Phases 1, 2 and V upstream (`silva-ssu-database`,
superresolution-amplicon#11); phase 3 in benchmark PR #26 and phase 4 in PR #27, both
merged. Deviations: the built `.tax` mirrors upstream's genome-collection form (padded
lineage plus a `Genome` level); the SILVA sweep grid gained a `nogate` point; and
`<id>.inferred_panel_members.csv` needed no publish wiring — upstream writes it into
`composition/<id>/`, which `RUN_SUPERRESOLUTION`'s existing `sr_out/composition/**`
output and publish pattern already carry. Phase V (`dev/panel_silva_sweep.md` upstream)
**limited Phase 4 to species**: taxon panels fail at genus scale (8,587 members) and work
at species scale (88 members, genome TV 0.034 against 0.030 for the genome panel on the
same kernel). The simulated reads now get the observed reads' primer trim
(`simulate_amplicon_reads.py --trim-primers`), which brings a SILVA genome panel to TV
0.024 (emulated on the existing trained reads; not yet resimulated with the oracle model).

Left open, both upstream: refusing `uncultured` / `Incertae Sedis` taxa, and validating
the presence gate with taxa (a 50-group entry is 50 gated parameters) before documenting
gate + taxa as supported.

Two changes, spanning this repo and `timrozday-mgnify/superresolution-amplicon` (upstream):

1. **SILVA SSU replaces GTDB** as the generic large MAPseq reference set in the examples, in
   the documentation of the mis-mapping matrix `M` / panel kernel `K`, and in what that
   database is *for*.
2. **Taxon panels.** A panel entry can be a taxon instead of a genome, or a panel can mix
   both. A taxon entry's sequences come from the generic database's own taxonomy. **Tested
   and supported at species level only** (Phase V).

Upstream lands first; the benchmark pins it through `sr_revision` / `sr_amplicon_repo`.

## Background: what exists today

- The generic database reaches the nested run as a pre-built `path:` dir holding
  `<name>_ssu.sr_refs.fasta` (headers `{genome_id}|{copy}|{orig}`). There is no taxonomy
  sidecar. Upstream writes a synthetic `Bacteria;<genome>;<header>` `.tax` for MAPseq
  (`subspecies_infer.write_mapseq_tax`). The real `--taxonomy` is used only for the `lca`
  column of `v4_group` output.
- `panel:` in `sr_settings` names a `databases:` collection. Its `<name>_ssu.sr_refs.fasta`
  becomes `--panel_references`, via `srPanelArgs` in `subworkflows/local/profile/main.nf`.
- Upstream `build_panel_kernel.py prepare` cuts each panel genome's V4 copies into sources
  (`v4g_<sha16>`) and writes `panel_translation.tsv` (genome → source, copy weight).
  `build` measures `K` (sources × database V4 groups), and `infer_composition.run_panel`
  fits θ over the panel genomes plus `background`.
- GTDB shaped two assumptions:
  - Database ids are genome accessions, so `examples/sr_amplicon_gtdb_sweep` scores
    profiles against panel ids that *must be GTDB accessions*.
  - Genome space over the database means something.

  SILVA has neither. A SILVA "genome" is one rRNA sequence, and its taxonomy stops at genus
  (the organism name is free text).

## Part 1: SILVA SSU replaces GTDB

### 1.1 What the generic database is for (the design statement to write down)

With GTDB gone, the generic database is a **label space only**:

- MAPseq labels reads with database V4 groups. The groups' lineages come from the SILVA
  taxonomy.
- Abundances are never reported *over the database's own references*. They are reported
  over one of two things:
  - a panel, via the rectangular `K` (panel sources × SILVA V4 groups), or
  - taxa, via `v4_group` inference plus `lca` (the `sr_amplicon_gtdb_sweep` use case,
    Part 1.4).
- Square `M` over the generic database is still used, but only inside `v4_group` inference
  (one source per SILVA V4 group). Genome-space inference against SILVA is sequence-space
  inference and has no biological reading.

Write this into:

- upstream `README.md`: the "GTDB-scale safety" section becomes "Generic-database scale";
  also the `--infer_space`, `--taxonomy` and panel reinterpretation sections;
- upstream `nextflow.config` comments (lines ~87–93 and ~193–212);
- upstream `docs/panel_reinterpretation_plan.md` "Goal" (it already says "GTDB, SILVA");
- the benchmark `README.md` panel paragraph (~line 534–545) and `CLAUDE.md` "Key
  implementation notes" (one sentence on the generic DB being a label space).

Dev results measured on GTDB r232 (`dev/panel_reinterpretation.md`,
`dev/gtdb_*`) stay as they are and are cited as GTDB numbers. They are not replaced by
SILVA claims until Phase V below reruns them.

### 1.2 Upstream: build the SILVA reference set

`bin/build_mapseq_database.py` already has a GTDB form (`--ssu-fasta`). Add `--silva-fasta`
beside it, not a generic header-parser abstraction:

- Input: `SILVA_<ver>_SSURef_NR99_tax_silva.fasta(.gz)`. Header:
  `>AB000001.1.1500 Bacteria;Bacillota;...;Genus;organism name`.
- Convert **U→T**. SILVA ships RNA, and `extract_v4` and MAPseq expect DNA; this is the
  easiest thing to miss.
- Header written: `{accession}|0|{accession}` (one "copy" per sequence; 0-based like the
  other two forms).
- `.tax` line: `{header}\t{lineage without the trailing organism-name field}`.
- Default to Ref **NR99**. Full Ref is roughly 4× larger for near-duplicate labels that
  MAPseq would split arbitrarily anyway.
- Partial sequences that lack a primer site drop out at in-silico PCR as now (logged
  count). Nothing new to handle.

Status: done on upstream branch `silva-ssu-database` (uncommitted). The staged-dir taxonomy
sidecar was **dropped**: `--taxonomy` is already a global input to INFER_COMPOSITION, and
Phase 2 can hand the same `params.taxonomy` to PANEL_PREPARE. MAPseq keeps its synthetic
`.tax`; only field 2 is used.

Tests: `--demo` covers the header, U→T, organism-name drop, padding, duplicate accessions
and a lineage-less header. `tests/test_build_mapseq_database.py` gains a Docker-gated
MAPseq acceptance test for the SILVA form.

### 1.3 Benchmark: a pre-built dir carries its taxonomy

- The pre-built layout gains an **optional** `<name>_ssu.sr_refs.tax`, in MAPseq `.tax`
  format as written by 1.2. In `subworkflows/local/build_databases/main.nf`, resolve it
  next to the `sr_refs.fasta` glob (optional, `[]` when absent) and emit it with `sr_dbs`
  as `[ "<name>:ssu", refs_fasta, tax|[] ]`.
- Built collections: `SR_BUILD_REFS` (`bin/build_sr_refs.py`) writes the same `.tax` when
  every sequence has `taxonomy:`, and none otherwise. The `taxonomy:` field already exists
  for AAP collections, so there is no new samplesheet key.
- `PROFILE` passes `--taxonomy <tax>` to the nested run whenever the row's database has
  one. That is harmless for existing runs: it only fills the `lca` column.
- Add a unit test in `tests/bin/test_bin.py` for the `build_sr_refs.py` `.tax` output.

### 1.4 Examples

| Where | Change |
|---|---|
| `examples/sr_amplicon_panel_sweep` | `gtdb:` block → `generic:` (values: `silva_138_2_ssu_nr99`, a `path:`). Arm `gtdb_panel` → `generic_panel`. Update `scripts/panel_sweep.py` (`ARMS`, block validation), `score_sweep.py`, `generate_*`, `run.sh` and README prose. The "development sweep reported 0.015" paragraph stays but is labelled as a GTDB r232 result. |
| `examples/sr_amplicon_gtdb_sweep` | Rename it `sr_amplicon_silva_sweep`. Its scoring **cannot survive** as written, because panel ids can no longer be database accessions. Switch it to `--infer_space v4_group` with `--taxonomy`, and score at **genus**. Truth genomes carry a SILVA lineage in `config.yaml` (`panel[].taxonomy`). The prediction is v4-group mass summed by `lca` truncated to genus, and mass whose `lca` stops above genus is reported as `unresolved`. Drop the "panel ids must be the reference set's ids" check and its selfcheck. Its question ("is the large set fit for purpose") is unchanged. |
| Benchmark `README.md` | Lines 202/241 (`gtdb_r220 path:` example for a pre-built DB) and 394–399 (mapseq/AAP config): use a SILVA name. **Sylph examples stay GTDB** (lines 359–363, `nextflow.config:27`), because sylph needs genomes and that is the reason for the switch. |
| `tests/samplesheets/sr_panel.yaml` | Comment: "stands in for the generic database (SILVA SSU)". |

Leave `examples/*/scripts/build_profiling_dbs.py` alone: its "production GTDB db" comment
refers to sylph.

## Part 2: taxon panels

### 2.1 Semantics

A panel is a list of **entries**. Each entry is one of two kinds:

- **genome**: `{id, ssu}`. As today: its V4 copies are sources with fixed copy weights.
- **taxon**: `{id, taxon}`. Its sources are every distinct V4 group of the generic
  database whose members' lineage falls under `taxon`. **Each group is a free source**, and
  the entry's abundance is the sum over its groups.

Rules:

1. **Matching.** `taxon` is a lineage prefix matched at `;` boundaries
   (`Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides`). A bare
   name (`Bacteroides`) is accepted only if it resolves to exactly one prefix. Zero
   matches, or an ambiguous name, is an error that lists the candidates. SILVA has
   cross-kingdom homonyms, so silent first-match is not acceptable.
2. **Group membership.** A V4 group belongs to a taxon if **any** member sequence's
   lineage is under it. A group shared across taxa, such as a V4 sequence identical
   between two genera, goes to each such entry as a shared source. The existing
   shared-source path (the *B. uniformis* pair) already handles that, and the entries are
   then reported `not_identifiable` as today.
3. **Precedence.** A group that is a genome entry's source belongs to that genome, not to
   any taxon. When taxa nest, a group belongs to the most specific matching entry. So
   `g__Bacteroides` next to a `B. fragilis` genome means "other *Bacteroides*".
4. **Species only.** Benchmark samplesheets, tests and examples use species taxa. Higher
   ranks are not refused upstream, but they are untested: inference collapses above a few
   hundred members, and a gut genus has hundreds to thousands of SILVA V4 groups.
5. **Cost guard.** `panel_taxon_max_sources` (default 200). An entry resolving to more V4
   groups fails with its count and a hint to use a lower rank. Every source costs
   `panel_sim_n_per_ref` simulated reads mapped against SILVA.
   `# ponytail: hard cap; sub-sample or cluster groups if broad taxa turn out to be needed.`

Rejected alternative: a taxon as one "genome" whose copy weights are its groups'
representation in SILVA. It needs no inference change, but a sample holds one or a few
strains of a genus, not SILVA's average. A fixed mixture puts the wrong label distribution
on the row, and the misfit leaks into `background` and neighbours.

### 2.2 Upstream (superresolution-amplicon)

Inputs:

- New param **`panel_taxa`**, plus a per-row samplesheet column overriding it (as
  `panel_references` does). It is a TSV `id<TAB>taxon`.
- `panel_references` becomes optional when `panel_taxa` is set. At least one of the two is
  required for a panel sample.
- `--taxonomy` is required whenever `panel_taxa` is set. Fail in `main.nf` validation.

`bin/build_panel_kernel.py prepare`:

- New `--panel-taxa` and `--db-taxonomy`. `db_groups()` already returns headers, group
  index and group sequences. Add a header → lineage map and resolve each taxon to groups
  under rules 1–3 and 5.
- Taxon sources are the database group sequences themselves. They are already insert-only,
  so no `extract_v4`, and `in_db` is always true.
- `panel_translation.tsv`: one row per taxon member, `genome_id = <entry>::<v4g>`, weight
  1. Genome entries are unchanged.
- New `panel_entries.tsv` (`member`, `entry`, `kind`). Genome members map to themselves.
- `sources.tsv` gains `entries`.
- `build` and `align` are unchanged: they already work per source.

`bin/infer_composition.py` `run_panel`:

- Fit θ over members exactly as now; the model is unchanged.
- Aggregate **posterior draws** by `panel_entries.tsv` before summarising: entry mean,
  `inferred_lo`/`inferred_hi` and `presence_prob` come from summed draws, not from summed
  intervals.
- `inferred_composition.csv` is reported at entry level. The member table is published
  as `<id>.inferred_panel_members.csv`.
- Identifiability (`_panel_identifiability`) is computed on entry rows (T·K summed per
  entry). Members of one taxon colliding with each other is expected and harmless.
- `check_composition_fit.py` reads the member-level table (it needs θ over kernel rows).
  Point it at the members file when it exists.

Status: done (uncommitted, same upstream branch). Deviations from the text above:

- No `panel_entries.tsv` and no `entries` column in `sources.tsv`. The entry is the part of
  a member id before `::` (`infer_composition._panel_entries`), so the kernel dir needs no
  new file. Prepare rejects ids containing `::` or equal to `background`.
- The fit check is not re-pointed at the members file. It keeps reading
  `inferred_composition.csv` and compares it to member draws summed per entry, so the
  workflow wiring is unchanged. `<id>.inferred_panel_members.csv` is still published.
- Entry `presence_prob` is `1 − Π(1 − p_member)`. That is exact under the mean-field
  presence guide, which has no per-draw gate samples to sum.
- Rule 3 is applied per database sequence, not per group: each sequence counts for the most
  specific entry above it. A group whose sequences sit under `Bacteroides` and under other
  `Bacteroidales` goes to both entries, not only to `Bacteroides`.
- An entry left with no group (none amplified, or all taken by genomes or narrower taxa)
  is dropped with a warning, like an unamplifiable genome.
- Not done: refusing `uncultured` / `Incertae Sedis` taxa (still an open question below).

Tests: `tests/test_panel_taxa.py` (7) and two stubs in `tests/default.nf.test` (taxa-only
plus mixed, and missing `--taxonomy`).

Priors: the presence gate and horseshoe act per member, so a 50-group genus is 50 gated
parameters. The panel recommendation is already no gate. Validate that before documenting
gate + taxa as supported (Phase V).

Tests:

- `tests/test_panel_taxa.py`:
  - a toy DB + `.tax`: prefix match, bare-name ambiguity error, nested-taxon precedence,
    genome-over-taxon exclusion, cap error;
  - draw aggregation (entry interval from summed draws, not the sum of intervals);
  - a regression case: a taxon whose only group equals a genome's single source gives the
    same θ as that genome entry.
- `tests/default.nf.test`: stubs for a taxa-only and a mixed panel.

### 2.3 Benchmark

Samplesheet:

```yaml
databases:
  silva_nr99:
    path: /dbs/silva_138_2_ssu_nr99        # sr_refs.fasta + sr_refs.tax
  panel_mixed:
    sequences:
      - {id: bacteroides_uniformis,         ssu: references/16S/FNPN01_SSU.fasta}
      - {id: bacteroides_uniformis_strain2, ssu: references/16S/BU_JCM13286_NT5170.1_SSU.fasta}
      - {id: bacteroides_fragilis,     taxon: "Bacteroides fragilis"}
      - {id: streptococcus_salivarius, taxon: "Bacteria;Bacillota;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus;Streptococcus salivarius"}
```

- The database `.tax` must carry a species rank. `--silva-fasta` drops the organism name
  today; it gains the rule from Phase V (`dev/panel_silva_sweep.species_taxonomy` upstream):
  a genus-level Bacteria/Archaea lineage gets `<genus> <epithet>` when the organism name is
  a binomial whose epithet is not `sp.`/`bacterium`/…. This is an upstream change and lands
  first.
- `main.nf`:
  - a sequence has `taxon:` **or** `genome`/`ssu`, not both;
  - the existing "`sr_amplicon` requires `ssu`" check skips taxon entries;
  - a collection with any taxon entry may only be named by `panel:`. It is an error if a
    row's `database` names it, or if any non-`sr_amplicon` profiler would build it.
  - Extend the seqs map with `taxon`.
- `BUILD_DATABASES` / `SR_BUILD_REFS`:
  - a panel collection publishes `<name>_ssu.sr_refs.fasta` from its genome entries (absent
    when taxa-only) plus `<name>.panel_taxa.tsv` from its taxon entries;
  - a `path:` dir reads both optionally, so the published dir stays reusable;
  - key stays `"<name>:ssu"`, with value `[ refs|[], taxa|[], tax|[] ]` (the tax slot from
    1.3).
- `PROFILE` (`ch_sr_panel`, `srPanelArgs`):
  - carry `panel_refs` and `panel_taxa` in meta;
  - emit `--panel_references` and/or `--panel_taxa`;
  - when the panel has taxa, require that the row's own database provides a tax. Fail
    with the sample, setting and database named. `database: self` has none, so a taxon
    panel over `self` is an error.
  - `srSetKey` already includes the panel name. No batching change.
- `bin/normalize_sr_profile.py`: entry ids pass through unchanged. Add a test row with a
  taxon id and `background`.
- `RUN_SUPERRESOLUTION` publish: add `<id>.inferred_panel_members.csv` next to
  `inferred_composition.csv`.

**Truth for taxon entries.** `truth.tsv` stays per genome; the pipeline does not roll it
up. The example's scorer does, using the same rules as 2.1:

- a truth genome that is a genome entry maps to itself;
- otherwise it maps to the most specific taxon entry its lineage falls under;
- otherwise it goes to `background`.

The genome's lineage is an explicit SILVA lineage in the example config.
`# ponytail: pipeline-level truth rollup deferred until a second consumer needs it.`

Example: add a third arm to `examples/sr_amplicon_panel_sweep`, `generic_taxa`.

- It maps against SILVA (reuses phase 1's `mseq:`).
- Its panel keeps the *B. uniformis* pair as genomes and makes each other genome a taxon
  entry for its SILVA species (the 20 names in `dev/panel_silva_sweep.SPECIES` upstream).
- Its grid is `generic_panel`'s kernel × `{nogate, horseshoe}`.
- Scoring uses genome TV (each entry stands for one truth genome) plus the BU split.
- Report the resolved source count per entry in the generator; Phase V saw 1–24.

Tests:

- `tests/samplesheets/sr_panel.yaml`: add a built generic collection with `taxonomy:`
  (standing in for SILVA), a taxa-only panel setting and a mixed one.
- Validation errors each get a `-preview`-level nf-test: taxon collection as `database`,
  taxon panel over `self`, and `taxon` + `ssu` on one entry.
- Test fixtures use species-rank taxa.

## Phases

| # | Repo | Work | Depends |
|---|---|---|---|
| 1 | upstream | 1.2 SILVA builder + taxonomy through staging | – |
| 2 | upstream | 2.2 taxon panels (prepare, aggregation, fit-check wiring, tests) | 1 |
| V | upstream `dev/` | Rerun `dev/panel_reinterpretation_sweep.py` against SILVA NR99 V4: genome panel (compare to GTDB 0.015), taxon panel at genus, gate vs no gate with taxa | 1, 2 (done: genome panel 0.025 with a length-matched kernel against GTDB 0.013; genus panel entry TV 0.21–0.35, inference collapses above a few hundred members; species panel genome TV 0.034 horseshoe / 0.039 no gate, 0/20 misfits) |
| 3 | benchmark | 1.3 tax sidecar + `--taxonomy`, 1.1/1.4 docs and example rename | 1 (done) |
| 4a | upstream | species rank in the `--silva-fasta` `.tax` (superresolution-amplicon#11) | V |
| 4 | benchmark | 2.3 species taxon entries, PROFILE wiring, `generic_taxa` arm, tests | 2, 3, 4a (done: branch `sr-taxon-panels`) |

Phases 3 and 4 can land as one PR if upstream is ready. Write the design statement
(1.1) after Phase V, so its numbers are SILVA numbers.

## Risks and open questions

- **SILVA lineages are uneven**: `uncultured`, `Incertae Sedis`, and genus-less paths. A
  taxon at `...;uncultured` is legal but meaningless. Proposal: refuse entries whose last
  field is `uncultured` or `Incertae Sedis`.
- **No species rank in SILVA.** It is derived from organism names (4a); only 119k of 510k
  NR99 sequences get one, and a species entry owns only groups with at least one named
  sequence. Below species, use genome entries.
- **NR99 clustering** removes near-identical references, so there are fewer, coarser V4
  labels than GTDB. Kernels get smaller; inference over the generic DB is less
  discriminating. Phase V measures it.
- **Truth lineage source.** The plan uses an explicit lineage per community genome. The
  alternative is deriving it from the DB by MAPseq-ing each truth genome's amplicon (its
  home label's `lca`). That is closer to "taxonomy from the reference database" but ties
  truth to MAPseq behaviour. **Decided: explicit** (`panel[].silva_taxon` in the example
  config). With species entries named after their genomes the scorer needs no rollup.
- **Broad taxa cost.** The 200-source cap never binds at species (max 24 in Phase V). Keep
  it as the guard against genus-scale entries.
