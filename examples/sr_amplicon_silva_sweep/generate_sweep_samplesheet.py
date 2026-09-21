#!/usr/bin/env python3
"""Phase 2: the parameter sweep and the arm comparison, as a profile-only samplesheet.

Emits one row per benchmark directory phase 1 generated (error-model arm x primer pair x
subsample depth) PER PROFILING ARM. The arms are `silva_panel` (SILVA reinterpreted
through the panel: the sweep proper, fanned over the whole grid), `custom` (the same reads
against a collection built from the panel) and `aap` (no superresolution at all). They
differ only in `database:`, `profilers:` and `sr_settings:`, so all three land in the
benchmark dir the reads are in and are scored against one truth.tsv.

Each row carries:

  benchmark_dir : where the reads actually are (`<sample>.<pair>/[subsample_<N>/]`)
  subsample     : that depth, so the profile is named and published the way the
                  generate step named it
  primers       : the ONE pair those reads were amplified with (sr_amplicon cuts its
                  reference amplicons with the same PCR)
  mseq          : that dir's mapseq classification from the phase-1 run. Mapping reads
                  against a SILVA-sized reference set is by far the most expensive stage
                  and no swept knob changes it, so every grid point reuses this file.

Every superresolution row carries its arm's own `sr_settings:` list, which the pipeline
fans that row out over. So the run is `rows x that arm's points` inference runs but only
`distinct kernel knob combinations` kernels.

    python generate_sweep_samplesheet.py [results_dir] [config.yaml]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))
import silva_sweep as gs

HERE = Path(__file__).resolve().parent
DEFAULT_RESULTS_DIR = HERE.parent.parent / "results" / "sr_amplicon_silva_sweep"


def main():
    results_dir = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else DEFAULT_RESULTS_DIR
    cfg = gs.load_config(sys.argv[2] if len(sys.argv) > 2 else HERE / "config.yaml")
    settings = gs.settings(cfg)

    rows, missing, n_dirs = [], [], 0
    for sample, directory, depth, pair, _em in gs.benchmark_dirs(cfg, results_dir):
        n_dirs += 1
        # Phase 1 published the classification under the generate step's own run id.
        run_id = f"{sample}.sub{depth}" if depth else sample
        mseq = directory / "profiling" / "sr" / f"{run_id}.map.obs.mseq.gz"
        if not mseq.exists():
            missing.append(str(mseq))
        for arm in gs.ARMS:
            row = {"sample": sample,
                   "profilers": ["aap"] if arm == "aap" else ["sr_amplicon"],
                   "benchmark_dir": str(directory),
                   "database": gs.arm_database(cfg, arm)}
            if arm != "aap":
                row["sr_settings"] = gs.arm_settings(cfg, arm)
            if depth:
                row["subsample"] = depth
            if pair:
                row["primers"] = [dict(pair)]
            # `custom` maps against its own 20-reference collection, so phase 1's SILVA
            # classification is not its classification. `aap` does its own mapping too -
            # that is the baseline being measured.
            if arm == "silva_panel" and mseq.exists():
                row["mseq"] = str(mseq)
            rows.append(row)

    doc = {"databases": gs.database_block(cfg), "samples": rows}
    with open(HERE / "sweep_samplesheet.yaml", "w") as fh:
        gs.dump_yaml(doc, fh)

    profiles = sum(len(r.get("sr_settings", [None])) for r in rows)
    print(f"Wrote sweep_samplesheet.yaml: {n_dirs} benchmark dir(s) x "
          f"{len(gs.ARMS)} arm(s) = {profiles} profiles "
          f"({len(settings)} grid point(s) on the 'silva_panel' arm, the rest at "
          f"'{gs.arm_point(cfg)['name']}'), from "
          f"{gs.n_matrices(gs.arm_settings(cfg, 'silva_panel'))} panel kernel(s) over "
          f"'{cfg['database']['name']}'")
    if missing:
        print(f"NOTE: {len(missing)} benchmark dir(s) have no phase-1 mapseq output "
              f"(e.g. {missing[0]}); those rows omit `mseq:` and the nested run will map "
              "their reads itself - which against a SILVA-sized reference set is the "
              "expensive thing this example exists to do once.", file=sys.stderr)


if __name__ == "__main__":
    main()
