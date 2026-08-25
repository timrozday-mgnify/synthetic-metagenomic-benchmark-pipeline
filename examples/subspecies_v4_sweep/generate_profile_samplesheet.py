#!/usr/bin/env python3
"""Generate a profile-only samplesheet (`--step profile`) to re-profile reads
already written by an earlier run, without regenerating them.

Reuses the same config.yaml as generate_sweep.py: it emits the same `databases:`
block (so the pipeline builds the `community_v4` DB) and one `samples:` row per
already-generated *benchmark directory* x profiler. A generate run fans each
samplesheet row out over its primer pairs and its subsampling depths, so the dirs
holding the reads are `<sample>.<pair>/` and `<sample>.<pair>/subsample_<N>/`, not
`<sample>/` - a row per pre-fan-out sample points at a directory that does not
exist. Each emitted row therefore carries:

  benchmark_dir : the dir the reads are actually in (pair suffix + depth subdir)
  subsample     : that depth, so the profile is named and published like the
                  generate step named it (`<sample>.subN.*` under `subsample_N/`)
  primers       : the ONE pair those reads were amplified with - sr_amplicon cuts
                  its reference amplicons with its own PCR, defaulting to V4, and
                  silently matches nothing against reads from another region

    python generate_profile_samplesheet.py [results_dir] [config.yaml]
    nextflow run ../../main.nf -profile docker --step profile \\
        --input profile_samplesheet.yaml --outdir <results_dir>
"""
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))
import sweep_config as sc

HERE = Path(__file__).resolve().parent
DEFAULT_RESULTS_DIR = HERE.parent.parent / "results" / "subspecies_v4_sweep"


def primer_pairs(primers):
    """A generate row's `primers:` (inline entries or a TSV path) -> [pair, ...], each
    kept in the samplesheet's own shape so it can be re-emitted verbatim. [None] when
    the row has no primers (no in-silico PCR, so no pair suffix on its dir)."""
    if not primers:
        return [None]
    if isinstance(primers, str):
        lines = [l for l in Path(primers).read_text().splitlines() if l.strip()]
        head = [c.strip() for c in lines[0].split("\t")]
        return [dict(zip(head, [c.strip() for c in l.split("\t")])) for l in lines[1:]]
    return list(primers)


def pair_id(pair):
    return pair["pair_id"] if isinstance(pair, dict) else pair[0]


def generated_dirs(samplesheet, results_dir):
    """One (sample_id, benchmark_dir, depth, pair, profiler) per generated benchmark
    directory x profiler, mirroring the generate step's fan-out: a row with `primers:`
    becomes one sample per pair (`<sample>.<pair>`), and each `subsample:` depth
    publishes into `<sample>/subsample_<N>/`."""
    if not samplesheet.exists():
        sys.exit(f"{samplesheet} not found - run generate_sweep.py first.")
    doc = yaml.safe_load(samplesheet.read_text())
    out = []
    for s in doc["samples"]:
        for pair in primer_pairs(s.get("primers")):
            sample = f"{s['sample']}.{pair_id(pair)}" if pair else s["sample"]
            for depth in s.get("subsample") or ["none"]:
                depth = None if depth in (None, "none", "") else depth
                d = results_dir / sample
                if depth:
                    d = d / f"subsample_{depth}"
                for profiler in s["profilers"]:
                    out.append((sample, d, depth, pair, profiler))
    return out


def main():
    results_dir = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else DEFAULT_RESULTS_DIR
    cfg_path = sys.argv[2] if len(sys.argv) > 2 else HERE / "config.yaml"
    cfg = sc.load_config(cfg_path)

    entries = generated_dirs(HERE / "samplesheet.yaml", results_dir)
    db_name = cfg["database"]["name"]

    rows = []
    for sample, bench_dir, depth, pair, profiler in entries:
        row = {"sample": sample, "profilers": [profiler],
               "benchmark_dir": str(bench_dir), "database": db_name}
        if depth:
            row["subsample"] = depth
        if pair:
            row["primers"] = [pair]
        rows.append(row)
    doc = {"databases": sc.database_block(cfg), **sc.aap_settings(cfg), "samples": rows}
    with open(HERE / "profile_samplesheet.yaml", "w") as fh:
        sc.dump_yaml(doc, fh)

    missing = sorted({r["benchmark_dir"] for r in rows
                      if not list(Path(r["benchmark_dir"]).glob("*.fastq.gz"))})
    print(f"Wrote profile_samplesheet.yaml: {len(rows)} rows (one per generated "
          f"benchmark dir x profiler) against '{db_name}', "
          f"benchmark_dir root {results_dir}")
    if missing:
        print(f"WARNING: {len(missing)} benchmark dir(s) hold no *.fastq.gz, e.g. "
              f"{missing[0]} - the profile run will fail on them.", file=sys.stderr)


if __name__ == "__main__":
    main()
