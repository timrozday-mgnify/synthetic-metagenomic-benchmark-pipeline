#!/usr/bin/env python3
"""Phase 2: the parameter sweep, as a profile-only samplesheet.

Emits one row per benchmark directory phase 1 generated (primer pair x subsample depth),
each carrying:

  benchmark_dir : where the reads actually are (`<sample>.<pair>/[subsample_<N>/]`)
  subsample     : that depth, so the profile is named and published the way the
                  generate step named it
  primers       : the ONE pair those reads were amplified with (sr_amplicon cuts its
                  reference amplicons with the same PCR)
  mseq          : that dir's mapseq classification from the phase-1 run. Mapping reads
                  against a GTDB-sized reference set is by far the most expensive stage
                  and no swept knob changes it, so every grid point reuses this file.

plus a top-level `sr_settings:` list - the expanded grid - which the pipeline fans each
row out over. So the run is `rows x grid points` inference runs but only
`distinct matrix knob combinations` mis-mapping matrices.

    python generate_sweep_samplesheet.py [results_dir] [config.yaml]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))
import gtdb_sweep as gs

HERE = Path(__file__).resolve().parent
DEFAULT_RESULTS_DIR = HERE.parent.parent / "results" / "sr_amplicon_gtdb_sweep"


def benchmark_dirs(cfg, results_dir):
    """(sample, dir, depth, pair) per generated benchmark dir, mirroring the generate
    step's fan-out: one sample per community x mode x primer pair, one subdir per depth."""
    for gm in gs.generation_modes(cfg):
        for pair in gm.get("primers") or [None]:
            for i in range(1, cfg["sampling"]["n_samples"] + 1):
                sample = f"S{i:02d}.{gm['name']}"
                if pair:
                    sample = f"{sample}.{pair['pair_id']}"
                for depth in gs.depths(cfg):
                    directory = results_dir / sample
                    if depth:
                        directory = directory / f"subsample_{depth}"
                    yield sample, directory, depth, pair


def main():
    results_dir = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else DEFAULT_RESULTS_DIR
    cfg = gs.load_config(sys.argv[2] if len(sys.argv) > 2 else HERE / "config.yaml")
    settings = gs.settings(cfg)

    rows, missing = [], []
    for sample, directory, depth, pair in benchmark_dirs(cfg, results_dir):
        # Phase 1 published the classification under the generate step's own run id.
        run_id = f"{sample}.sub{depth}" if depth else sample
        mseq = directory / "profiling" / "sr" / f"{run_id}.obs.mseq.gz"
        row = {"sample": sample, "profilers": ["sr_amplicon"],
               "benchmark_dir": str(directory), "database": cfg["database"]["name"]}
        if depth:
            row["subsample"] = depth
        if pair:
            row["primers"] = [dict(pair)]
        if mseq.exists():
            row["mseq"] = str(mseq)
        else:
            missing.append(str(mseq))
        rows.append(row)

    doc = {"databases": gs.database_block(cfg), "sr_settings": settings, "samples": rows}
    with open(HERE / "sweep_samplesheet.yaml", "w") as fh:
        gs.dump_yaml(doc, fh)

    print(f"Wrote sweep_samplesheet.yaml: {len(rows)} benchmark dir(s) x "
          f"{len(settings)} grid point(s) = {len(rows) * len(settings)} profiles, "
          f"from {gs.n_matrices(settings)} mis-mapping matrix/matrices over "
          f"'{cfg['database']['name']}'")
    if missing:
        print(f"NOTE: {len(missing)} benchmark dir(s) have no phase-1 mapseq output "
              f"(e.g. {missing[0]}); those rows omit `mseq:` and the nested run will map "
              "their reads itself - which against a GTDB-sized reference set is the "
              "expensive thing this example exists to do once.", file=sys.stderr)


if __name__ == "__main__":
    main()
