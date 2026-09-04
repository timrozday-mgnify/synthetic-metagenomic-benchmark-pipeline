#!/usr/bin/env python3
"""Shared config loader for the superresolution-amplicon parameter sweep.

One `config.yaml` is the single source of truth for both generator scripts:
`generate_samplesheet.py` (phase 1, generate the reads) and
`generate_sweep_samplesheet.py` (phase 2, re-profile them under every parameter
setting). Relative paths in the config resolve against the config file's directory.

Run `python scripts/sr_sweep.py --selfcheck` to exercise the grid expansion and the
database block without touching the pipeline.
"""
import itertools
import sys
from pathlib import Path

import yaml

# Knobs an sr_settings entry may carry (mirrors srSettingKeys() in ../../main.nf).
SETTING_KEYS = {
    "mismapping_method", "align_backend", "align_tau", "matrix_args",
    "infer_presence", "infer_presence_prior", "infer_presence_temp", "inference_args",
}
# The subset that changes the mis-mapping matrix. Settings agreeing on all of these
# share one matrix and only re-run the (cheap) inference — which is what makes a grid
# with several inference points affordable.
MATRIX_KEYS = {"mismapping_method", "align_backend", "align_tau", "matrix_args"}


def dump_yaml(doc, fh):
    yaml.safe_dump(doc, fh, sort_keys=False, default_flow_style=False, width=10**6)


def load_config(path):
    path = Path(path).resolve()
    cfg = yaml.safe_load(path.read_text())
    cfg["_dir"] = path.parent
    _validate(cfg)
    return cfg


def resolve(cfg, p):
    """Config-relative path -> absolute string (the pipeline resolves relative paths
    against its own projectDir, not the samplesheet's, so everything is emitted absolute)."""
    p = Path(str(p))
    return str(p if p.is_absolute() else (cfg["_dir"] / p).resolve())


def _validate(cfg):
    for key in ("train", "reads", "panel", "database", "sr_sweep", "primers"):
        if key not in cfg:
            sys.exit(f"config.yaml: missing required '{key}:' block")
    if not cfg["panel"]:
        sys.exit("config.yaml: 'panel:' is empty")
    ids = [m["id"] for m in cfg["panel"]]
    if len(set(ids)) != len(ids):
        sys.exit("config.yaml: panel 'id's must be unique")
    for m in cfg["panel"]:
        if not m.get("ssu"):
            sys.exit(f"config.yaml: panel member {m['id']} needs 'ssu' "
                     "(16S drives both amplicon generation and the sr_amplicon references)")
    settings(cfg)  # expands (and validates) the grid up front


def taxonomy(member):
    """3-rank 'Kingdom;Genus;species', derived from the id unless given."""
    if member.get("taxonomy"):
        return member["taxonomy"]
    kingdom = "Archaea" if member.get("kingdom") == "archaea" else "Bacteria"
    genus, _, species = member["id"].partition("_")
    return f"{kingdom};{genus.capitalize()};{species or 'sp'}"


def database_block(cfg):
    """The samplesheet `databases:` block: sr_amplicon needs only `ssu` per sequence."""
    db = cfg["database"]
    return {db["name"]: {"sequences": [
        {"id": m["id"], "ssu": resolve(cfg, m["ssu"]), "taxonomy": taxonomy(m)}
        for m in cfg["panel"]
    ]}}


def settings(cfg):
    """Expand `sr_sweep.grid` (a map of axis -> list of named knob maps) into the flat
    `sr_settings:` list the pipeline fans out over. One entry per combination across
    axes; its name is the axis entry names joined with '.', which becomes part of the
    output filename (`<id>.<name>.sr_profile.tsv`)."""
    grid = cfg["sr_sweep"].get("grid") or {}
    if not grid:
        sys.exit("config.yaml: 'sr_sweep.grid:' is empty — nothing to sweep")
    axes = []
    for axis, entries in grid.items():
        if not entries:
            sys.exit(f"config.yaml: sr_sweep axis '{axis}' is empty")
        for e in entries:
            if not e.get("name"):
                sys.exit(f"config.yaml: sr_sweep axis '{axis}' entry {e} needs a 'name'")
            unknown = set(e) - SETTING_KEYS - {"name"}
            if unknown:
                sys.exit(f"config.yaml: sr_sweep '{e['name']}': unknown knob(s) "
                         f"{sorted(unknown)}; expected {sorted(SETTING_KEYS)}")
        axes.append(entries)
    out = []
    for combo in itertools.product(*axes):
        merged = {"name": ".".join(e["name"] for e in combo)}
        for e in combo:
            merged.update({k: v for k, v in e.items() if k != "name"})
        out.append(merged)
    names = [s["name"] for s in out]
    if len(set(names)) != len(names):
        sys.exit(f"config.yaml: sr_sweep produced duplicate setting names: {names}")
    return out


def n_matrices(sweep_settings):
    """How many mis-mapping matrices the grid actually costs: settings agreeing on every
    matrix knob share one reference set, and so one matrix."""
    return len({tuple(sorted((k, str(v)) for k, v in s.items() if k in MATRIX_KEYS))
                for s in sweep_settings})


def depths(cfg):
    """Subsample depths, `None` for full depth. One benchmark dir per depth."""
    v = cfg["reads"].get("subsample") or [None]
    return [None if d in (None, "none", "") else int(d) for d in v]


def _selfcheck():
    cfg = {"sr_sweep": {"grid": {
        "mismapping": [{"name": "sim", "mismapping_method": "simulate"},
                       {"name": "exact", "mismapping_method": "align",
                        "align_backend": "exact-hash", "align_tau": 0}],
        "inference": [{"name": "p01", "infer_presence_prior": 0.01},
                      {"name": "p001", "infer_presence_prior": 0.001}],
    }}, "reads": {"subsample": ["none", 100000]}}
    s = settings(cfg)
    assert [e["name"] for e in s] == ["sim.p01", "sim.p001", "exact.p01", "exact.p001"], s
    assert s[3] == {"name": "exact.p001", "mismapping_method": "align",
                    "align_backend": "exact-hash", "align_tau": 0,
                    "infer_presence_prior": 0.001}, s[3]
    # Two matrix modes x two inference points = 4 runs but only 2 matrices.
    assert n_matrices(s) == 2, n_matrices(s)
    assert depths(cfg) == [None, 100000], depths(cfg)
    assert taxonomy({"id": "bacteroides_fragilis"}) == "Bacteria;Bacteroides;fragilis"
    assert taxonomy({"id": "methanobrevibacter_smithii", "kingdom": "archaea"}) \
        == "Archaea;Methanobrevibacter;smithii"
    print("sr_sweep selfcheck OK")


if __name__ == "__main__":
    if "--selfcheck" in sys.argv:
        _selfcheck()
    else:
        print(__doc__)
