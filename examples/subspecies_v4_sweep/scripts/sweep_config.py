#!/usr/bin/env python3
"""Shared config loader for the subspecies_v4_sweep example scripts.

A single `config.yaml` (train + reads + sweep + database + panel) drives
generate_sweep.py, generate_profile_samplesheet.py and (optionally)
scripts/build_profiling_dbs.py, so no paths or metadata are hard-coded in Python.
Requires PyYAML - run the scripts with `python` (not necessarily `python3`).

Schema (see ../config.yaml for a filled-in template):

    train:    {id, fastq_1, fastq_2, platform}
    reads:    {num_reads, subsample?, paired_end, read_length_mean, read_length_variance}
              # shared read-generation defaults; a generation mode may override any field.
    sweep:    {n_samples, steepness}
    database: {name, profilers: [aap|sylph, ...],
               rfam_covariance_model?, rfam_claninfo?}  # required if 'aap' in profilers
    aap:      {configs?: [path, ...], profile?}         # optional; nested-AAP engine/-c files
    generation_modes:                                   # each sweep sample is emitted once per mode
              [ {name, source: genome|ssu|amplicon, mode: shotgun|amplicon|long,
                 profiler?, primers?: [...]|path-to-TSV, reads?: {...}}, ... ]
    panel:    [ {id, species, amplicon?, ssu?, genome?, taxonomy?, kingdom?}, ... ]

Exactly one `species` must appear twice: that pair is the abundance-sweep target.

A generation mode's `source` picks which panel FASTA feeds read generation
(`genome`/`ssu`/`amplicon`). `mode: amplicon` from `genome`/`ssu` requires
`primers:` (in-silico PCR, each pair run as its own benchmark); `source: amplicon`
uses pre-trimmed FASTAs directly (no primers). `generation_modes:` is optional: a
legacy config with `reads.mode` (+ optional top-level `primers:`) synthesizes one.
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
    # `primers:` (top-level, legacy) or per-mode may be an inline list of pairs
    # (passed through) or a path to a TSV.
    if isinstance(cfg.get("primers"), str):
        cfg["primers"] = resolve(cfg["primers"])
    for gm in cfg.get("generation_modes") or []:
        if isinstance(gm.get("primers"), str):
            gm["primers"] = resolve(gm["primers"])

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

    # Each generation mode picks a panel FASTA source; validate its shape and that
    # every member carries the source field it needs.
    for gm in generation_modes(cfg):
        name = gm.get("name", "?")
        src = gm.get("source")
        if src not in ("genome", "ssu", "amplicon"):
            sys.exit(f"config: generation mode '{name}' needs `source:` one of genome|ssu|amplicon")
        rmode = gm.get("mode", "shotgun")
        has_primers = bool(gm.get("primers"))
        if rmode == "amplicon" and src in ("genome", "ssu") and not has_primers:
            sys.exit(f"config: generation mode '{name}' (amplicon from {src}) needs `primers:`")
        if src == "amplicon" and has_primers:
            sys.exit(f"config: generation mode '{name}' uses pre-trimmed amplicons; remove `primers:`")
        if gm.get("profiler") and gm["profiler"] not in profilers:
            sys.exit(f"config: generation mode '{name}' profiler '{gm['profiler']}' "
                     f"not in database.profilers {profilers}")
        for m in panel:
            if not m.get(src):
                sys.exit(f"config: panel member '{m['id']}' needs `{src}:` (generation mode '{name}')")

    for m in panel:
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


def generation_modes(cfg):
    """Normalized list of read-generation modes (each sweep sample is emitted once
    per mode). Falls back to a single mode synthesized from the legacy `reads.mode`
    (+ optional top-level `primers:`) when `generation_modes:` is absent."""
    modes = cfg.get("generation_modes")
    if modes:
        return modes
    primers = cfg.get("primers")
    rmode = cfg["reads"].get("mode", "shotgun")
    return [{
        "name": rmode,
        "source": "genome" if primers else ("genome" if rmode == "shotgun" else "amplicon"),
        "mode": rmode,
        "profiler": cfg["database"]["profilers"][0],
        **({"primers": primers} if primers else {}),
    }]


def mode_source_field(m):
    """Panel FASTA field feeding read generation for this mode (genome|ssu|amplicon)."""
    return m["source"]


def mode_reads(cfg, m):
    """Global `reads:` defaults overlaid with this mode's optional `reads:` override."""
    reads = dict(cfg["reads"])
    reads.update(m.get("reads") or {})
    return reads


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

        # `primers:` set: genome required (extraction source), primers path resolved.
        primed = cfg_text.replace(
            "panel:",
            "primers: [{pair_id: v4, forward: GTGYCAG, reverse: GGACTAC}]\npanel:",
        ).replace("amplicon: refs/a.fa", "genome: refs/a.fna").replace(
            "amplicon: refs/b.fa", "genome: refs/b.fna"
        ).replace("amplicon: refs/b2.fa", "genome: refs/b2.fna")
        p.write_text(primed)
        cfg2 = load_config(p)
        assert cfg2["primers"][0]["pair_id"] == "v4", cfg2["primers"]
        assert cfg2["panel"][0]["genome"] == str(base / "refs/a.fna"), cfg2["panel"][0]
        # primers set but a member lacks genome -> exits.
        p.write_text(primed.replace("genome: refs/b2.fna", "amplicon: refs/b2.fa"))
        try:
            load_config(p)
        except SystemExit as e:
            assert "genome" in str(e), e
        else:
            raise AssertionError("expected SystemExit for missing genome under primers")

        # generation_modes: fan-out source selection + per-mode reads merge.
        gm_text = textwrap.dedent("""
            train: {id: t, fastq_1: /abs/r1.fq, fastq_2: /abs/r2.fq, platform: hq-illumina}
            reads: {num_reads: 100, subsample: none, paired_end: true, read_length_mean: 300, read_length_variance: 0}
            sweep: {n_samples: 4, steepness: 6.0}
            database: {name: db, profilers: [sylph, aap], rfam_covariance_model: /abs/ribo, rfam_claninfo: /abs/ribo.clan}
            generation_modes:
              - {name: wgs, source: genome, mode: shotgun, profiler: sylph, reads: {read_length_mean: 150}}
              - {name: amp16s, source: ssu, mode: amplicon, profiler: aap, primers: [{pair_id: v4, forward: GTGYCAG, reverse: GGACTAC}]}
            panel:
              - {id: a, species: genus_a, genome: refs/a.fna, ssu: refs/a.16s.fa, amplicon: refs/a.amp.fa}
              - {id: b, species: genus_b, genome: refs/b.fna, ssu: refs/b.16s.fa, amplicon: refs/b.amp.fa}
              - {id: b2, species: genus_b, genome: refs/b2.fna, ssu: refs/b2.16s.fa, amplicon: refs/b2.amp.fa}
        """)
        p.write_text(gm_text)
        cfg3 = load_config(p)
        modes = generation_modes(cfg3)
        assert [gm["name"] for gm in modes] == ["wgs", "amp16s"], modes
        assert mode_source_field(modes[0]) == "genome" and mode_source_field(modes[1]) == "ssu"
        assert mode_reads(cfg3, modes[0])["read_length_mean"] == 150, "per-mode override"
        assert mode_reads(cfg3, modes[1])["read_length_mean"] == 300, "falls back to global"

        # amplicon-from-ssu without primers -> exits.
        p.write_text(gm_text.replace(
            ", primers: [{pair_id: v4, forward: GTGYCAG, reverse: GGACTAC}]", ""))
        try:
            load_config(p)
        except SystemExit as e:
            assert "primers" in str(e), e
        else:
            raise AssertionError("expected SystemExit for amplicon-from-ssu without primers")

        # legacy fallback: no generation_modes -> single mode from reads.mode.
        assert generation_modes(cfg)[0]["source"] == "amplicon", "legacy amplicon fallback"
    print("selfcheck OK")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--selfcheck":
        _selfcheck()
    else:
        sys.exit("sweep_config is a library; run a generate_*/build_* script instead "
                 "(or `--selfcheck`).")
