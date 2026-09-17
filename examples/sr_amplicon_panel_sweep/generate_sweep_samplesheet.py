#!/usr/bin/env python3
"""Phase 2: both arms of the sweep, as one profile-only samplesheet.

Two rows per benchmark dir phase 1 generated, each with its own row-level `sr_settings:`:

  custom      database: the collection built from `panel:`. No `mseq:`: phase 1 mapped
              against SILVA, the wrong reference set for this row, and mapping against 23
              references is cheap.
  generic_panel  database: SILVA. `mseq:` is phase 1's classification of this dir against it,
              and `sr_error_model:` is the skiver model phase 1 trained. Every setting
              carries `panel:` (the custom collection), so inference runs over the panel
              genomes plus `background`, not over SILVA's sequences.

    python generate_sweep_samplesheet.py [results_dir] [config.yaml]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))
import panel_sweep as ps

HERE = Path(__file__).resolve().parent
DEFAULT_RESULTS_DIR = HERE.parent.parent / "results" / "sr_amplicon_panel_sweep"


def main():
    results_dir = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else DEFAULT_RESULTS_DIR
    cfg = ps.load_config(sys.argv[2] if len(sys.argv) > 2 else HERE / "config.yaml")
    custom, panel = ps.settings(cfg, "custom"), ps.settings(cfg, "generic_panel")
    map_name = cfg["generic"]["map_setting"]["name"]
    model = ps.trained_model(cfg, results_dir)

    rows, missing = [], []
    for sample, run_id, directory, depth, pair in ps.benchmark_dirs(cfg, results_dir):
        base = {"sample": sample, "profilers": ["sr_amplicon"], "benchmark_dir": str(directory),
                **({"subsample": depth} if depth else {}),
                **({"primers": [dict(pair)]} if pair else {})}
        rows.append({**base, "database": cfg["database"]["name"], "sr_settings": custom})
        generic = {**base, "database": cfg["generic"]["name"], "sr_settings": panel}
        # Phase 1 fanned out over its one mapping setting, so its id carries that name.
        mseq = directory / "profiling" / "sr" / f"{run_id}.{map_name}.obs.mseq.gz"
        if mseq.exists():
            generic["mseq"] = str(mseq)
        else:
            missing.append(str(mseq))
        if model:
            generic["sr_error_model"] = model
        rows.append(generic)

    doc = {"databases": ps.databases_block(cfg), "samples": rows}
    with open(HERE / "sweep_samplesheet.yaml", "w") as fh:
        ps.dump_yaml(doc, fh)

    n_dirs = len(rows) // 2
    print(f"Wrote sweep_samplesheet.yaml: {n_dirs} benchmark dir(s) x ({len(custom)} custom "
          f"+ {len(panel)} generic_panel) setting(s) = {n_dirs * (len(custom) + len(panel))} "
          f"profiles, from {ps.n_matrices(custom)} custom matrix/matrices and "
          f"{len(panel)} panel kernel(s)")
    if missing:
        print(f"NOTE: {len(missing)} benchmark dir(s) have no phase-1 SILVA mapping "
              f"(e.g. {missing[0]}); their generic_panel rows omit `mseq:` and every panel "
              "setting will map those reads against SILVA again.", file=sys.stderr)
    if not model:
        print(f"WARNING: no single *.model.pt under {results_dir}/error_models/"
              f"{cfg['train']['id']}/; `--sim_error_model trained` panel settings will "
              "train their own model per sample instead of using the one the reads were "
              "generated with.", file=sys.stderr)


if __name__ == "__main__":
    main()
