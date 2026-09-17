#!/usr/bin/env python3
"""Shared config glue for the SILVA panel-reinterpretation vs custom-database sweep.

Two sibling examples already do most of this, so they are imported rather than copied:

* ``abundance_nb_sample/scripts/nb_config.py`` — config schema, path resolution, the
  seeded negative-binomial communities and the custom collection's ``databases:`` entry.
* ``sr_amplicon_param_sweep/scripts/sr_sweep.py`` — expansion of a ``grid:`` into
  ``sr_settings:`` entries, and the count of mis-mapping matrices a grid costs.

What is new is that the same reads are profiled two ways, each under a grid of its own:
``custom`` maps them against a collection built from the panel, and ``generic_panel`` maps
them against SILVA and reinterprets the labels through that same collection (the
``panel`` sr_settings knob).

    python scripts/panel_sweep.py --selfcheck
"""
import sys
from pathlib import Path

_EXAMPLES = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(_EXAMPLES / "abundance_nb_sample" / "scripts"))
sys.path.insert(0, str(_EXAMPLES / "sr_amplicon_param_sweep" / "scripts"))
import nb_config as nc          # noqa: E402
import sr_sweep as sw           # noqa: E402

# Re-exported so the generator scripts import one module, not three.
dump_yaml = nc.dump_yaml
sample_abundances = nc.sample_abundances
generation_modes = nc.generation_modes
mode_reads = nc.mode_reads
depths = sw.depths
n_matrices = sw.n_matrices

ARMS = ("custom", "generic_panel")


def load_config(path):
    """`nb_config.load_config`, plus this example's `generic:`, `score:` and two-arm `sr_sweep:`."""
    cfg = nc.load_config(path)
    generic = cfg.get("generic") or {}
    for key in ("name", "path", "map_setting"):
        if not generic.get(key):
            sys.exit(f"config.yaml: generic needs '{key}:'")
    if generic["name"] == cfg["database"]["name"]:
        sys.exit("config.yaml: generic.name and database.name must differ")
    generic["path"] = str(Path(str(generic["path"])).expanduser())
    if sorted(cfg.get("sr_sweep") or {}) != sorted(ARMS):
        sys.exit(f"config.yaml: sr_sweep needs exactly the arms {list(ARMS)}")
    ids = {m["id"] for m in cfg["panel"]}
    if not set((cfg.get("score") or {}).get("strain_pair") or []) <= ids:
        sys.exit("config.yaml: score.strain_pair must name two panel ids")
    for arm in ARMS:
        settings(cfg, arm)              # expands (and validates) both grids up front
    return cfg


def settings(cfg, arm):
    """One arm's grid as `sr_settings:` entries named `<arm>.<point>`, so both arms'
    profiles can sit in one benchmark dir. Every `generic_panel` point reinterprets SILVA's
    labels through the custom database's collection."""
    out = sw.settings({"sr_sweep": cfg["sr_sweep"][arm]})
    for s in out:
        s["name"] = f"{arm}.{s['name']}"
        if arm == "generic_panel":
            s["panel"] = cfg["database"]["name"]
    return out


def databases_block(cfg):
    """Both reference sets: the collection built from `panel:` (the custom database, and
    the reinterpretation panel) and SILVA, pre-built. The pipeline builds only the ones a
    row or a `panel:` knob references, so phase 1 builds neither."""
    return {**nc.database_block(cfg), cfg["generic"]["name"]: {"path": cfg["generic"]["path"]}}


def benchmark_dirs(cfg, results_dir):
    """(sample, run id, dir, depth, primer pair) per benchmark dir the generate step
    publishes: community x mode x primer pair, one subdir per depth."""
    for gm in generation_modes(cfg):
        for pair in gm.get("primers") or [None]:
            for i in range(1, cfg["sampling"]["n_samples"] + 1):
                sample = f"S{i:02d}.{gm['name']}" + (f".{pair['pair_id']}" if pair else "")
                for depth in depths(cfg):
                    directory = results_dir / sample / (f"subsample_{depth}" if depth else "")
                    run_id = f"{sample}.sub{depth}" if depth else sample
                    yield sample, run_id, directory, depth, pair


def trained_model(cfg, results_dir):
    """The skiver model phase 1 generated the reads with, or None. The `trained` panel
    kernel simulates with it, which makes that arm oracle-optimistic (see README)."""
    found = sorted((results_dir / "error_models" / cfg["train"]["id"]).glob("*.model.pt"))
    return str(found[0]) if len(found) == 1 else None


def _selfcheck():
    cfg = {"database": {"name": "custom", "profilers": ["sr_amplicon"]},
           "generic": {"name": "silva", "path": "/x/silva"},
           "panel": [{"id": "a", "ssu": "/x/a.fa"}],
           "sampling": {"n_samples": 2},
           "reads": {"subsample": ["none", 100]},
           "generation_modes": [{"name": "amp", "primers": [{"pair_id": "V4"}]}],
           "sr_sweep": {
               "custom": {"grid": {"matrix": [{"name": "sim", "mismapping_method": "simulate"}]}},
               "generic_panel": {"grid": {
                   "kernel": [{"name": "trained", "mismapping_method": "simulate",
                               "matrix_args": "--sim_error_model trained"}],
                   "prior": [{"name": "nogate", "infer_presence": False},
                             {"name": "hs", "infer_presence": False,
                              "inference_args": "--infer_horseshoe true"}],
                   "steps": [{"name": "s10k", "inference_args": "--infer_steps 10000"}]}}}}

    custom = settings(cfg, "custom")
    assert [s["name"] for s in custom] == ["custom.sim"] and "panel" not in custom[0], custom
    panel = settings(cfg, "generic_panel")
    assert [s["name"] for s in panel] == ["generic_panel.trained.nogate.s10k",
                                          "generic_panel.trained.hs.s10k"], panel
    assert all(s["panel"] == "custom" for s in panel), panel
    # Flags from two axes add up; the steps axis must not drop the horseshoe flag.
    assert panel[0]["inference_args"] == "--infer_steps 10000", panel[0]
    assert panel[1]["inference_args"] == "--infer_horseshoe true --infer_steps 10000", panel[1]

    block = databases_block(cfg)
    assert block["silva"] == {"path": "/x/silva"}, block
    assert block["custom"] == {"sequences": [{"id": "a", "ssu": "/x/a.fa"}]}, block

    dirs = list(benchmark_dirs(cfg, Path("/r")))
    assert len(dirs) == 4, dirs
    assert dirs[0][:4] == ("S01.amp.V4", "S01.amp.V4", Path("/r/S01.amp.V4"), None), dirs[0]
    assert dirs[1][:4] == ("S01.amp.V4", "S01.amp.V4.sub100",
                           Path("/r/S01.amp.V4/subsample_100"), 100), dirs[1]
    assert trained_model({"train": {"id": "t"}}, Path("/nonexistent")) is None
    print("panel_sweep selfcheck OK")


if __name__ == "__main__":
    if "--selfcheck" in sys.argv:
        _selfcheck()
    else:
        print(__doc__)
