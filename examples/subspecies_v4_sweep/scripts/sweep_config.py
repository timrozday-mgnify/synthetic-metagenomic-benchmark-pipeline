#!/usr/bin/env python3
"""Shared config loader for the subspecies_v4_sweep example scripts.

A single `config.yaml` (train + reads + sweep + database + panel) drives
generate_sweep.py, generate_profile_samplesheet.py and (optionally)
scripts/build_profiling_dbs.py, so no paths or metadata are hard-coded in Python.
Requires PyYAML - run the scripts with `python` (not necessarily `python3`).

Schema (see ../config.yaml for a filled-in template):

    train:    {id, fastq_1, fastq_2, platform}
    reads:    {num_reads, mode, paired_end, read_length_mean, read_length_variance}
    sweep:    {n_samples, steepness}
    database: {name, profilers: [aap|sylph, ...],
               rfam_covariance_model?, rfam_claninfo?}  # required if 'aap' in profilers
    aap:      {configs?: [path, ...], profile?}         # optional; nested-AAP engine/-c files
    panel:    [ {id, species, amplicon, ssu?, genome?, taxonomy?, kingdom?}, ... ]

Exactly one `species` must appear twice: that pair is the abundance-sweep target.
"""
import math
import sys
from collections import Counter
from pathlib import Path

import yaml


def load_config(path):
    """Parse config.yaml, resolve panel/train paths (relative -> against the
    config file's dir), validate, and return the dict (with `_dir` added)."""
    path = Path(path).resolve()
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    base = path.parent
    cfg["_dir"] = base

    def resolve(p):
        p = Path(str(p)).expanduser()
        return str(p if p.is_absolute() else base / p)

    for key in ("fastq_1", "fastq_2"):
        if cfg.get("train", {}).get(key):
            cfg["train"][key] = resolve(cfg["train"][key])
    for member in cfg.get("panel") or []:
        for key in ("amplicon", "ssu", "genome"):
            if member.get(key):
                member[key] = resolve(member[key])
    for key in ("rfam_covariance_model", "rfam_claninfo"):
        if cfg.get("database", {}).get(key):
            cfg["database"][key] = resolve(cfg["database"][key])
    aap_configs = (cfg.get("aap") or {}).get("configs")
    if aap_configs:
        cfg["aap"]["configs"] = [resolve(p) for p in aap_configs]

    _validate(cfg)
    return cfg


def _validate(cfg):
    panel = cfg.get("panel") or []
    if not panel:
        sys.exit("config: `panel:` is empty")

    ids = [m["id"] for m in panel]
    dup_ids = [i for i, n in Counter(ids).items() if n > 1]
    if dup_ids:
        sys.exit(f"config: duplicate panel id(s): {dup_ids}")

    species_counts = Counter(m["species"] for m in panel)
    doubled = [sp for sp, n in species_counts.items() if n == 2]
    over = [sp for sp, n in species_counts.items() if n > 2]
    if over or len(doubled) != 1:
        sys.exit(f"config: expected exactly one species appearing twice (the sweep "
                 f"pair); doubled={doubled}, >2={over}")

    profilers = cfg["database"]["profilers"]
    if "aap" in profilers:
        for key in ("rfam_covariance_model", "rfam_claninfo"):
            if not cfg["database"].get(key):
                sys.exit(f"config: database needs `{key}:` (aap rRNA detection)")
    for m in panel:
        if not m.get("amplicon"):
            sys.exit(f"config: panel member '{m['id']}' needs `amplicon:`")
        if "aap" in profilers and not m.get("ssu"):
            sys.exit(f"config: panel member '{m['id']}' needs `ssu:` (aap profiling)")
        if "sylph" in profilers and not m.get("genome"):
            sys.exit(f"config: panel member '{m['id']}' needs `genome:` (sylph profiling)")


def sweep_pair(cfg):
    """(doubled_species, major_id, minor_id) - major is the first-listed of the two."""
    panel = cfg["panel"]
    counts = Counter(m["species"] for m in panel)
    doubled = next(sp for sp, n in counts.items() if n == 2)
    major, minor = (m["id"] for m in panel if m["species"] == doubled)
    return doubled, major, minor


def taxonomy(member):
    """`Kingdom;Genus;species` from member.taxonomy, else derived from the species
    slug (`genus_species`) and optional `kingdom: archaea`."""
    if member.get("taxonomy"):
        return member["taxonomy"]
    genus, species = member["species"].split("_", 1)
    kingdom = "Archaea" if str(member.get("kingdom", "")).lower() == "archaea" else "Bacteria"
    return f"{kingdom};{genus.capitalize()};{species}"


def database_block(cfg):
    """`{name: {sequences: [...]}}` for the config's database, emitting only the
    fields each requested profiler needs (genome for sylph; ssu+taxonomy for aap)."""
    db = cfg["database"]
    profilers = db["profilers"]
    sequences = []
    for m in cfg["panel"]:
        seq = {"id": m["id"]}
        if "sylph" in profilers:
            seq["genome"] = m["genome"]
        if "aap" in profilers:
            seq["ssu"] = m["ssu"]
            seq["taxonomy"] = taxonomy(m)
        sequences.append(seq)
    entry = {"sequences": sequences}
    # aap needs the Rfam rRNA-detection DBs (main.nf requires both for an aap DB).
    if "aap" in profilers:
        entry["rfam_covariance_model"] = db["rfam_covariance_model"]
        entry["rfam_claninfo"] = db["rfam_claninfo"]
    return {db["name"]: entry}


def aap_settings(cfg):
    """Top-level samplesheet keys for the nested amplicon-analysis-pipeline run
    (container engine / extra `-c` configs). Emitted only when `aap:` is set in
    config.yaml; the pipeline falls back to its params otherwise."""
    aap = cfg.get("aap") or {}
    out = {}
    if aap.get("configs"):
        out["aap_configs"] = aap["configs"]
    if aap.get("profile"):
        out["aap_profile"] = aap["profile"]
    return out


def logistic_fracs(n, k):
    """n fractions in [0,1], symmetric, denser near the 0/1 extremes (logistic
    spacing). Endpoints are rescaled to land exactly on 0 and 1."""
    s = [1 / (1 + math.exp(-k * (2 * i / (n - 1) - 1))) for i in range(n)]
    lo, hi = s[0], s[-1]
    return [(v - lo) / (hi - lo) for v in s]


def _selfcheck():
    """Round-trip a fixture config; check path resolution, database_block, sweep
    pair, taxonomy derivation, and the aap/sylph validation errors."""
    import tempfile
    import textwrap

    cfg_text = textwrap.dedent("""
        train: {id: t, fastq_1: reads/r1.fq, fastq_2: /abs/r2.fq, platform: hq-illumina}
        reads: {num_reads: 100, mode: amplicon, paired_end: true, read_length_mean: 300, read_length_variance: 0}
        sweep: {n_samples: 4, steepness: 6.0}
        database: {name: community_v4, profilers: [aap], rfam_covariance_model: /abs/ribo, rfam_claninfo: /abs/ribo.clan}
        aap: {configs: [engine.config, /abs/site.config], profile: singularity}
        panel:
          - {id: a, species: genus_a, amplicon: refs/a.fa, ssu: refs/a.16s.fa}
          - {id: b, species: genus_b, amplicon: refs/b.fa, ssu: refs/b.16s.fa}
          - {id: b2, species: genus_b, amplicon: refs/b2.fa, ssu: refs/b2.16s.fa}
    """)
    with tempfile.TemporaryDirectory() as d:
        base = Path(d).resolve()  # macOS /var -> /private/var symlink
        p = base / "config.yaml"
        p.write_text(cfg_text)
        cfg = load_config(p)

        assert cfg["panel"][0]["ssu"] == str(base / "refs/a.16s.fa"), "relative -> config dir"
        assert cfg["train"]["fastq_2"] == "/abs/r2.fq", "absolute path untouched"

        seqs = database_block(cfg)["community_v4"]["sequences"]
        assert [s["id"] for s in seqs] == ["a", "b", "b2"], seqs
        assert seqs[0]["taxonomy"] == "Bacteria;Genus;a" and "genome" not in seqs[0], seqs

        assert sweep_pair(cfg) == ("genus_b", "b", "b2"), sweep_pair(cfg)

        aap = aap_settings(cfg)
        assert aap["aap_configs"] == [str(base / "engine.config"), "/abs/site.config"], aap
        assert aap["aap_profile"] == "singularity", aap

        fr = logistic_fracs(4, 6.0)
        assert abs(fr[0]) < 1e-9 and abs(fr[-1] - 1) < 1e-9, fr

        # sylph requested but no `genome:` -> validation exits.
        p.write_text(cfg_text.replace("profilers: [aap]", "profilers: [sylph]"))
        try:
            load_config(p)
        except SystemExit as e:
            assert "genome" in str(e), e
        else:
            raise AssertionError("expected SystemExit for missing genome")
    print("selfcheck OK")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--selfcheck":
        _selfcheck()
    else:
        sys.exit("sweep_config is a library; run a generate_*/build_* script instead "
                 "(or `--selfcheck`).")
