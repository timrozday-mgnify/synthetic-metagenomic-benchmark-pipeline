#!/usr/bin/env python
"""Pre-process a synthetic-metagenomic-benchmark-pipeline run into tidy CSVs.

Walks one completed run tree once and emits the flat tables the benchmark Quarto
report reads via ``-P`` paths, so all the heavy file-walking / BAM / mseq parsing
happens here (in the ``basicpython`` env) rather than inside the report.

For every ``(sample, depth)`` cell it joins the ground-truth abundance table
against two independent detection sources:

* **profiler** — per-read mapseq assignments (``*.mseq.gz``); the ``query`` field
  encodes the true origin genome and ``dbhit`` the assigned genome, so this gives
  both per-genome detected abundance and a truth x assigned confusion matrix. This
  is the amplicon (amplicon-analysis-pipeline) path.
* **reference** — the ground-truth read->reference BAM (``*.sorted.bam``). Read
  names encode the origin genome and each reference contig is ``<genome>:<contig>``,
  so primary alignments give a general (amplicon *and* WGS) mis-mapping signal.

Outputs (written into ``--output-dir``): ``abundance.csv``, ``mismapping.csv``,
``summary.csv``, ``meta.json``. Both detection sources are optional per cell, so a
WGS/sylph run with no mseq (or a run with no BAM) still produces the tables it can.

Run:
    benchmark_preprocess.py --run-dir <run> --output-dir <run>
    benchmark_preprocess.py --demo        # self-check, no run needed
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# Globs are relative to a single (sample, depth) cell directory.
TRUTH_GLOB = "*.truth.tsv"
MSEQ_GLOB_DEFAULT = "profiling/aap/*/taxonomy-summary/*/*.mseq.gz"
BAM_GLOB = "*.sorted.bam"
SYLPH_GLOB = "*.sylph_profile.tsv"

# `S10_a0.42` -> the swept fraction encoded in the sample-dir suffix.
SWEEP_X_RE = re.compile(r"_a([0-9]*\.?[0-9]+)$")

# A pipeline sample dir is `<sweep_point>.<assay>`, e.g.
# `S10_a0.42.amplicon_515YF-806BR_16s.515-YF-806BR` or `S10_a0.42.wgs`. The sweep
# point (`S10_a0.42`) is shared across assays; the assay suffix must start with a
# letter so the sweep fraction's own decimals aren't mistaken for it.
SAMPLE_ASSAY_RE = re.compile(r"^(?P<sample>.*_a[0-9]*\.?[0-9]+)\.(?P<assay>[A-Za-z].*)$")

# Abundance is compared for up to three independent detection sources.
SOURCE_ORDER = ["profiler", "reference", "sylph"]
ABUND_COL = {
    "profiler": "detected_profiler_rel_abundance",
    "reference": "detected_reference_rel_abundance",
    "sylph": "detected_sylph_rel_abundance",
}


def assay_label(raw: str) -> str:
    """Assay dir suffix -> short label. `wgs` -> `WGS`; amplicon suffixes end in the
    primer token (e.g. `...16s.515-YF-806BR` -> `515-YF-806BR`)."""
    if raw == "wgs":
        return "WGS"
    return raw.rsplit(".", 1)[-1]


def split_sample_assay(dirname: str) -> tuple[str, str]:
    """Sample dir name -> (sweep_point, assay_label). Runs without an assay suffix
    (older single-assay layout) fall back to one `reads` assay."""
    m = SAMPLE_ASSAY_RE.match(dirname)
    if m:
        return m.group("sample"), assay_label(m.group("assay"))
    return dirname, "reads"


def truth_genome_from_query(query: str) -> str:
    """mapseq `#query` / BAM read name -> origin genome id (2nd `:`-field)."""
    # e.g. "S10_a0.42.chunk0:bacteroides_uniformis:ENA|...".
    parts = query.split(":")
    return parts[1] if len(parts) > 1 else parts[0]


def assigned_genome_from_dbhit(dbhit: str) -> str:
    """mapseq `dbhit` -> assigned genome id (before first `|`)."""
    return dbhit.split("|", 1)[0]


def ref_genome_from_name(ref_name: str) -> str:
    """BAM reference contig `<genome>:<contig>` -> genome id (before first `:`)."""
    return ref_name.split(":", 1)[0]


def read_truth(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    need = {"genome_id", "target_rel_abundance", "realized_n_reads", "realized_rel_abundance"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"{path}: truth table missing columns {sorted(missing)}")
    return df


def _parse_mseq_lines(handle) -> tuple[Counter, Counter]:
    """mseq lines -> (assigned-genome counts, (truth,assigned) counts). Handle-based
    so the streaming file parser and the in-memory demo share one body."""
    detected: Counter = Counter()
    confusion: Counter = Counter()
    for line in handle:
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 2 or not fields[1]:
            continue
        truth = truth_genome_from_query(fields[0])
        assigned = assigned_genome_from_dbhit(fields[1])
        detected[assigned] += 1
        confusion[(truth, assigned)] += 1
    return detected, confusion


def parse_mseq(path: Path) -> tuple[Counter, Counter]:
    """Stream a gzipped mseq file -> (assigned-genome counts, (truth,assigned) counts)."""
    with gzip.open(path, "rt") as fh:
        return _parse_mseq_lines(fh)


def parse_bam(path: Path) -> tuple[Counter, Counter]:
    """Primary alignments of a BAM -> (reference-genome counts, (truth,ref) counts)."""
    import pysam  # local import: only the reference path needs it

    detected: Counter = Counter()
    confusion: Counter = Counter()
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for aln in bam:
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            truth = truth_genome_from_query(aln.query_name)
            ref = ref_genome_from_name(bam.get_reference_name(aln.reference_id))
            detected[ref] += 1
            confusion[(truth, ref)] += 1
    return detected, confusion


def load_accession_map(samplesheet: Path | None) -> dict[str, str]:
    """Pipeline samplesheet `id`/`genome` pairs -> {ref_accession: genome_id}. The
    accession is the genome fasta basename sans extension — exactly what sylph reports
    in its `genome_id` column — so it maps sylph rows back to community genome ids."""
    if not samplesheet or not samplesheet.exists():
        return {}
    acc2genome: dict[str, str] = {}
    last_id = None
    for line in samplesheet.read_text().splitlines():
        s = line.strip()
        m = re.match(r"-?\s*id:\s*(\S+)", s)
        if m:
            last_id = m.group(1)
            continue
        m = re.match(r"genome:\s*(\S+)", s)
        if m and last_id:
            acc2genome[Path(m.group(1)).stem] = last_id
    return acc2genome


def parse_sylph(path: Path, acc2genome: dict[str, str]) -> Counter:
    """A sylph profile TSV -> per-genome detected abundance (sequence abundance),
    accessions folded to community genome ids. Abundance-only: sylph gives no per-read
    assignments, so there is no confusion matrix."""
    df = pd.read_csv(path, sep="\t")
    counts: Counter = Counter()
    for acc, ab in zip(df["genome_id"].astype(str), df["predicted_rel_abundance"].astype(float)):
        counts[acc2genome.get(acc, acc)] += ab
    return counts


def find_cells(pipeline_root: Path) -> list[dict]:
    """Discover (sample, assay, depth) cells: each sample dir + its subsample_<N> subdirs.
    The dir name splits into a sweep point (shared across assays) and an assay label."""
    cells = []
    for sample_dir in sorted(p for p in pipeline_root.iterdir() if p.is_dir()):
        if not list(sample_dir.glob(TRUTH_GLOB)):
            continue  # skip databases/, error_models/, pipeline_info/ etc.
        sample, assay = split_sample_assay(sample_dir.name)
        m = SWEEP_X_RE.search(sample)
        sweep_x = float(m.group(1)) if m else None
        base = {"sample": sample, "assay": assay, "sweep_x": sweep_x}
        cells.append({**base, "depth": "full", "dir": sample_dir})
        for sub in sorted(sample_dir.glob("subsample_*")):
            if sub.is_dir() and list(sub.glob(TRUTH_GLOB)):
                cells.append({**base, "depth": sub.name.replace("subsample_", "sub"), "dir": sub})
    return cells


def rel(counts: Counter) -> dict[str, float]:
    total = sum(counts.values())
    return {g: c / total for g, c in counts.items()} if total else {}


def detect_sweep_pair(abundance: pd.DataFrame, override: str | None) -> list[str]:
    """Genome_ids whose target abundance varies across sweep points (the swept pair).
    Variation is checked *within* each assay, then unioned: different assays can
    normalise target abundance differently (e.g. amplicon vs WGS truth tables), so
    pooling assays would make every genome look like it varies."""
    if override:
        return [g.strip() for g in override.split(",") if g.strip()]
    full = abundance[abundance["depth"] == "full"].copy()
    if "assay" not in full.columns:
        full["assay"] = "reads"
    varying = {
        g for (_assay, g), sub in full.groupby(["assay", "genome_id"])
        if sub["target_rel_abundance"].dropna().nunique() > 1
    }
    return sorted(varying)


def build_tables(run_dir: Path, pipeline_dir: str, mseq_glob: str,
                 sweep_pair_override: str | None, run_label: str,
                 acc2genome: dict[str, str] | None = None) -> dict:
    pipeline_root = run_dir / pipeline_dir
    if not pipeline_root.is_dir():
        raise ValueError(f"pipeline dir not found: {pipeline_root} (set --pipeline-dir)")
    cells = find_cells(pipeline_root)
    if not cells:
        raise ValueError(f"no cells with a {TRUTH_GLOB} found under {pipeline_root}")
    acc2genome = acc2genome or {}

    abundance_rows: list[dict] = []
    mismap_rows: list[dict] = []

    for cell in cells:
        cdir, sample, assay = cell["dir"], cell["sample"], cell["assay"]
        depth, sweep_x = cell["depth"], cell["sweep_x"]
        truth = read_truth(next(iter(cdir.glob(TRUTH_GLOB))))

        prof_detected, prof_conf = Counter(), Counter()
        for mseq in cdir.glob(mseq_glob):
            d, c = parse_mseq(mseq)
            prof_detected += d
            prof_conf += c
        ref_detected, ref_conf = Counter(), Counter()
        for bam in cdir.glob(BAM_GLOB):
            d, c = parse_bam(bam)
            ref_detected += d
            ref_conf += c
        syl_detected: Counter = Counter()
        for syl in cdir.glob(SYLPH_GLOB):
            syl_detected += parse_sylph(syl, acc2genome)

        prof_rel, ref_rel, syl_rel = rel(prof_detected), rel(ref_detected), rel(syl_detected)
        genomes = set(truth["genome_id"]) | set(prof_rel) | set(ref_rel) | set(syl_rel)
        tmap = truth.set_index("genome_id")
        for g in sorted(genomes):
            row = {"sample": sample, "assay": assay, "depth": depth, "sweep_x": sweep_x,
                   "genome_id": g}
            if g in tmap.index:
                row["target_rel_abundance"] = float(tmap.at[g, "target_rel_abundance"])
                row["realized_rel_abundance"] = float(tmap.at[g, "realized_rel_abundance"])
            else:
                row["target_rel_abundance"] = 0.0
                row["realized_rel_abundance"] = 0.0
            row["detected_profiler_rel_abundance"] = prof_rel.get(g) if prof_detected else None
            row["detected_reference_rel_abundance"] = ref_rel.get(g) if ref_detected else None
            row["detected_sylph_rel_abundance"] = syl_rel.get(g) if syl_detected else None
            abundance_rows.append(row)

        for source, conf in (("profiler", prof_conf), ("reference", ref_conf)):
            per_truth = defaultdict(int)
            for (t, _a), n in conf.items():
                per_truth[t] += n
            for (t, a), n in conf.items():
                mismap_rows.append({
                    "sample": sample, "assay": assay, "depth": depth, "source": source,
                    "truth_genome": t, "assigned_genome": a, "reads": n,
                    "frac_of_truth": n / per_truth[t] if per_truth[t] else 0.0,
                })

    abundance = pd.DataFrame(abundance_rows)
    mismapping = pd.DataFrame(mismap_rows)
    sweep_pair = detect_sweep_pair(abundance, sweep_pair_override)
    abundance["is_sweep_pair"] = abundance["genome_id"].isin(sweep_pair)

    summary = build_summary(abundance, mismapping, sweep_pair)
    sources = [s for s in SOURCE_ORDER
               if abundance[ABUND_COL[s]].notna().any()
               or (not mismapping.empty and s in set(mismapping["source"]))]
    assays = sorted(abundance["assay"].unique().tolist())
    meta = {
        "run_label": run_label,
        "run_dir": str(run_dir),
        "pipeline_dir": pipeline_dir,
        "sweep_pair": sweep_pair,
        "assays": assays,
        "depths": sorted(abundance["depth"].unique().tolist()),
        "samples": sorted(abundance["sample"].unique().tolist()),
        "sample_sweep_x": (
            abundance.dropna(subset=["sweep_x"]).drop_duplicates("sample")
            .set_index("sample")["sweep_x"].to_dict()
        ),
        "sources": sources,
        "n_genomes": int(abundance["genome_id"].nunique()),
    }
    return {"abundance": abundance, "mismapping": mismapping, "summary": summary, "meta": meta}


def build_summary(abundance: pd.DataFrame, mismapping: pd.DataFrame,
                  sweep_pair: list[str]) -> pd.DataFrame:
    """Per (assay, depth, source): abundance accuracy, overall mis-mapping rate, and the
    within-pair mis-mapping rate (swept genomes assigned to their partner) across the sweep.
    Sources without a confusion matrix (e.g. sylph) get abundance metrics only."""
    pair = set(sweep_pair)
    rows = []
    for assay in sorted(abundance["assay"].unique()):
        aass = abundance[abundance["assay"] == assay]
        for depth in sorted(aass["depth"].unique()):
            adep = aass[aass["depth"] == depth]
            for source in SOURCE_ORDER:
                col = ABUND_COL[source]
                pts = adep.dropna(subset=[col])
                if pts.empty:
                    continue
                detected = pts[col].to_numpy(dtype=float)
                realized = pts["realized_rel_abundance"].to_numpy(dtype=float)
                l1 = 0.5 * abs(detected - realized).sum() / max(pts["sample"].nunique(), 1)
                pearson = pts[[col, "realized_rel_abundance"]].corr().iloc[0, 1]
                mm = mismapping[(mismapping["assay"] == assay) & (mismapping["depth"] == depth)
                                & (mismapping["source"] == source)] if not mismapping.empty \
                    else mismapping
                total = mm["reads"].sum() if "reads" in mm.columns else 0
                offdiag = mm[mm["truth_genome"] != mm["assigned_genome"]]["reads"].sum()
                # Within-pair mis-mapping: swept-pair reads assigned to the *other* member.
                pair_reads = mm[mm["truth_genome"].isin(pair)]["reads"].sum() if pair and total else 0
                pair_cross = mm[mm["truth_genome"].isin(pair) & mm["assigned_genome"].isin(pair)
                                & (mm["truth_genome"] != mm["assigned_genome"])]["reads"].sum() \
                    if total else 0
                rows.append({
                    "assay": assay, "depth": depth, "source": source,
                    "l1_error_per_sample": l1,
                    "pearson_r": pearson,
                    "mismapping_rate": (offdiag / total) if total else np.nan,
                    "pair_mismapping_rate": (pair_cross / pair_reads) if pair_reads else np.nan,
                    "n_reads": int(total),
                })
    return pd.DataFrame(rows)


def write_outputs(tables: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    tables["abundance"].to_csv(output_dir / "abundance.csv", index=False)
    tables["mismapping"].to_csv(output_dir / "mismapping.csv", index=False)
    tables["summary"].to_csv(output_dir / "summary.csv", index=False)
    (output_dir / "meta.json").write_text(json.dumps(tables["meta"], indent=2))


def print_summary(tables: dict) -> None:
    meta, summary = tables["meta"], tables["summary"]
    print(f"run: {meta['run_label'] or meta['run_dir']}")
    print(f"cells: {len(meta['samples'])} samples x {len(meta['assays'])} assays "
          f"x {len(meta['depths'])} depths ({meta['depths']}), {meta['n_genomes']} genomes")
    print(f"assays: {meta['assays']}")
    print(f"sweep pair: {meta['sweep_pair'] or 'none detected'}")
    print(f"detection sources: {meta['sources']}")
    if not summary.empty:
        print("summary (per depth x source):")
        print(summary.to_string(index=False))


def run_demo() -> None:
    """Self-check the mseq/BAM parsing, cross-tab and sweep-pair detection inline."""
    assert truth_genome_from_query("S10_a0.42.chunk0:bacteroides_uniformis:ENA|X") == "bacteroides_uniformis"
    assert assigned_genome_from_dbhit("bacteroides_uniformis_strain2|2|ENA|X") == "bacteroides_uniformis_strain2"
    assert ref_genome_from_name("bacteroides_uniformis:ENA|CR626927|CR626927.1_3") == "bacteroides_uniformis"

    # Tiny mseq: 2 reads from A stay, 1 read from A mis-assigned to B.
    import io
    text = (
        "#query\tdbhit\n"
        "s.c0:A:x\tA|1|z\n"
        "s.c0:A:x\tA|1|z\n"
        "s.c0:A:x\tB|1|z\n"
        "s.c0:B:x\tB|1|z\n"
    )
    detected, conf = _parse_mseq_lines(io.StringIO(text))
    assert detected == Counter({"A": 2, "B": 2}), detected
    assert conf[("A", "B")] == 1 and conf[("A", "A")] == 2, conf

    # Sweep-pair detection: only genome P varies across samples.
    ab = pd.DataFrame([
        {"depth": "full", "sample": "S1", "genome_id": "P", "target_rel_abundance": 0.1},
        {"depth": "full", "sample": "S2", "genome_id": "P", "target_rel_abundance": 0.9},
        {"depth": "full", "sample": "S1", "genome_id": "Q", "target_rel_abundance": 0.5},
        {"depth": "full", "sample": "S2", "genome_id": "Q", "target_rel_abundance": 0.5},
    ])
    assert detect_sweep_pair(ab, None) == ["P"], detect_sweep_pair(ab, None)
    assert detect_sweep_pair(ab, "Q,R") == ["Q", "R"]

    # Sample-dir -> (sweep point, assay) split, and the older no-suffix layout.
    assert split_sample_assay("S10_a0.42.amplicon_515YF-806BR_16s.515-YF-806BR") == \
        ("S10_a0.42", "515-YF-806BR"), split_sample_assay("S10_a0.42.amplicon_515YF-806BR_16s.515-YF-806BR")
    assert split_sample_assay("S10_a0.42.wgs") == ("S10_a0.42", "WGS")
    assert split_sample_assay("S10_a0.42") == ("S10_a0.42", "reads")

    # sylph accession -> genome_id via a samplesheet, incl. a `.fa` (not `.fasta`) genome.
    import tempfile
    ss = Path(tempfile.mkdtemp()) / "samplesheet.yaml"
    ss.write_text(
        "    - id: bacteroides_uniformis\n"
        "      genome: /refs/genomes/FNPN01.fasta\n"
        "    - id: bacteroides_uniformis_strain2\n"
        "      genome: /refs/genomes/BU_JCM13286_NT5170.1.fa\n"
    )
    amap = load_accession_map(ss)
    assert amap == {"FNPN01": "bacteroides_uniformis",
                    "BU_JCM13286_NT5170.1": "bacteroides_uniformis_strain2"}, amap
    print("demo: OK")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", type=Path, help="Run root (contains --pipeline-dir).")
    ap.add_argument("--pipeline-dir", default="results/subspecies_v4_sweep",
                    help="Sample dirs live here, relative to --run-dir.")
    ap.add_argument("--mseq-glob", default=MSEQ_GLOB_DEFAULT,
                    help="Per-cell glob for gzipped mapseq files.")
    ap.add_argument("--sweep-pair", default=None,
                    help="Override auto-detected swept pair, e.g. 'id1,id2'.")
    ap.add_argument("--samplesheet", type=Path, default=None,
                    help="Pipeline samplesheet.yaml for sylph accession->genome_id mapping "
                         "(default: <run-dir>/samplesheet.yaml if present).")
    ap.add_argument("--run-label", default="", help="Human-readable run label for meta.json.")
    ap.add_argument("--output-dir", type=Path, help="Where to write the CSVs (default: --run-dir).")
    ap.add_argument("--demo", action="store_true", help="Run self-check and exit.")
    args = ap.parse_args()

    if args.demo:
        run_demo()
        return
    if not args.run_dir:
        ap.error("--run-dir is required (or use --demo)")

    samplesheet = args.samplesheet or (args.run_dir / "samplesheet.yaml")
    acc2genome = load_accession_map(samplesheet)
    tables = build_tables(args.run_dir, args.pipeline_dir, args.mseq_glob,
                          args.sweep_pair, args.run_label, acc2genome)
    write_outputs(tables, args.output_dir or args.run_dir)
    print_summary(tables)


if __name__ == "__main__":
    main()
