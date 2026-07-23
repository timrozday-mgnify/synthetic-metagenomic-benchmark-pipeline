#!/usr/bin/env python3
"""Shared config loader for the abundance_nb_sample example scripts.

A single `config.yaml` (train + reads + sampling + database + panel) drives
generate_samplesheet.py, generate_profile_samplesheet.py and (optionally)
scripts/build_profiling_dbs.py, so no paths or metadata are hard-coded in Python.
Requires PyYAML + numpy - run the scripts with `python` (not necessarily `python3`).

Unlike the subspecies_v4_sweep example (which sweeps the split of one same-species
pair), here EVERY genome is drawn independently per sample: a Bernoulli(`presence`)
decides whether it is present, and if so a negative-binomial(mean, dispersion) sets
its abundance. Each sample is thus a different random community composition. The
draws are seeded (`sampling.seed`) for reproducibility.

Schema (see ../config.yaml for a filled-in template):

    train:    {id, fastq_1, fastq_2, platform}
    reads:    {num_reads, subsample?, paired_end, read_length_mean, read_length_variance}
              # shared read-generation defaults; a generation mode may override any field.
    sampling: {n_samples, seed, presence?, negative_binomial: {mean, dispersion}}
              # Per genome per sample: Bernoulli(presence) decides present/absent
              #   (default 1.0), then NB(mean, dispersion) sets abundance if present.
              #   NB variance = mean + mean^2/r. A panel member may override
              #   `presence:` and/or `nb: {mean?, dispersion?}` inline.
    database: {name, profilers: [aap|sylph, ...],
               rfam_covariance_model?, rfam_claninfo?}  # required if 'aap' in profilers
    aap:      {configs?: [path, ...], profile?}         # optional; nested-AAP engine/-c files
    generation_modes:                                   # each sample is emitted once per mode
              [ {name, source: genome|ssu|amplicon, mode: shotgun|amplicon|long,
                 profiler?, primers?: [...]|path-to-TSV, reads?: {...}}, ... ]
    panel:    [ {id, species, amplicon?, ssu?, genome?, taxonomy?, kingdom?,
                 presence?, nb?: {mean?, dispersion?}}, ... ]

A generation mode's `source` picks which panel FASTA feeds read generation
(`genome`/`ssu`/`amplicon`). `mode: amplicon` from `genome`/`ssu` requires
`primers:` (in-silico PCR, each pair run as its own benchmark); `source: amplicon`
uses pre-trimmed FASTAs directly (no primers). `generation_modes:` is optional: a
legacy config with `reads.mode` (+ optional top-level `primers:`) synthesizes one.
"""
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import yaml


class _NoAliasDumper(yaml.SafeDumper):
    # PyYAML emits *anchors for reused dict objects (e.g. the shared `subsample`
    # block across all samples); Nextflow's SnakeYAML caps non-scalar aliases at
    # 50 and dies with >50. Inlining every occurrence sidesteps the limit.
    def ignore_aliases(self, data):
        return True


def dump_yaml(doc, fh):
    yaml.dump(doc, fh, Dumper=_NoAliasDumper, sort_keys=False,
              default_flow_style=False)


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

    samp = cfg.get("sampling") or {}
    if not samp.get("n_samples"):
        sys.exit("config: `sampling.n_samples:` is required")
    nb = samp.get("negative_binomial") or {}
    if nb.get("mean") is None or nb.get("dispersion") is None:
        sys.exit("config: `sampling.negative_binomial:` needs `mean:` and `dispersion:`")
    for m in panel:
        mean, disp = _nb_params(cfg, m)
        if mean < 0:
            sys.exit(f"config: panel member '{m['id']}' nb.mean must be >= 0")
        if disp <= 0:
            sys.exit(f"config: panel member '{m['id']}' nb.dispersion must be > 0")
        if not 0 <= _presence_param(cfg, m) <= 1:
            sys.exit(f"config: panel member '{m['id']}' presence must be in [0, 1]")
    if max(_presence_param(cfg, m) for m in panel) == 0:
        sys.exit("config: all genomes have presence 0 - no community is possible")

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


def _nb_params(cfg, member):
    """(mean, dispersion) for a genome: the global `sampling.negative_binomial`
    defaults, overlaid with the member's optional inline `nb:` override."""
    nb = cfg["sampling"]["negative_binomial"]
    over = member.get("nb") or {}
    return (float(over.get("mean", nb["mean"])),
            float(over.get("dispersion", nb["dispersion"])))


def _presence_param(cfg, member):
    """Bernoulli presence probability for a genome: global `sampling.presence`
    (default 1.0 = always present), overlaid with the member's optional
    `presence:` override."""
    default = cfg["sampling"].get("presence", 1.0)
    return float(member.get("presence", default))


def sample_abundances(cfg):
    """(n_samples, n_genomes) integer matrix of abundances, one row per sample.
    Two independent draws per genome per sample: a Bernoulli(`presence`) decides
    whether the genome is present at all, and if so a NegBinom(mean, dispersion)
    sets its abundance (absent genomes are 0). Seeded by `sampling.seed`; a row
    that ends up all-zero (empty community) is redrawn."""
    panel = cfg["panel"]
    samp = cfg["sampling"]
    means = np.array([_nb_params(cfg, m)[0] for m in panel])
    disps = np.array([_nb_params(cfg, m)[1] for m in panel])
    pres = np.array([_presence_param(cfg, m) for m in panel])
    # numpy's NB is parameterized by (n=size r, p); mean = n*(1-p)/p, so with
    # dispersion r as n: p = r / (r + mean).
    p = disps / (disps + means)
    rng = np.random.default_rng(samp.get("seed", 0))
    rows = []
    for _ in range(samp["n_samples"]):
        while True:
            present = rng.random(len(panel)) < pres  # Bernoulli presence gate
            draw = present * rng.negative_binomial(disps, p)
            if draw.sum() > 0:
                rows.append(draw)
                break
    return np.array(rows, dtype=int)


def generation_modes(cfg):
    """Normalized list of read-generation modes (each sample is emitted once per
    mode). Falls back to a single mode synthesized from the legacy `reads.mode`
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


def _selfcheck():
    """Round-trip a fixture config; check path resolution, database_block, NB
    sampling (seeded/reproducible, no all-zero rows, per-genome override), taxonomy
    derivation, and the aap/sylph validation errors."""
    import tempfile
    import textwrap

    cfg_text = textwrap.dedent("""
        train: {id: t, fastq_1: reads/r1.fq, fastq_2: /abs/r2.fq, platform: hq-illumina}
        reads: {num_reads: 100, mode: amplicon, paired_end: true, read_length_mean: 300, read_length_variance: 0}
        sampling: {n_samples: 8, seed: 7, negative_binomial: {mean: 20, dispersion: 1.0}}
        database: {name: community, profilers: [aap], rfam_covariance_model: /abs/ribo, rfam_claninfo: /abs/ribo.clan}
        aap: {configs: [engine.config, /abs/site.config], profile: singularity}
        panel:
          - {id: a, species: genus_a, amplicon: refs/a.fa, ssu: refs/a.16s.fa, nb: {mean: 200}}
          - {id: b, species: genus_b, amplicon: refs/b.fa, ssu: refs/b.16s.fa}
          - {id: c, species: genus_c, amplicon: refs/c.fa, ssu: refs/c.16s.fa}
    """)
    with tempfile.TemporaryDirectory() as d:
        base = Path(d).resolve()  # macOS /var -> /private/var symlink
        p = base / "config.yaml"
        p.write_text(cfg_text)
        cfg = load_config(p)

        assert cfg["panel"][0]["ssu"] == str(base / "refs/a.16s.fa"), "relative -> config dir"
        assert cfg["train"]["fastq_2"] == "/abs/r2.fq", "absolute path untouched"

        seqs = database_block(cfg)["community"]["sequences"]
        assert [s["id"] for s in seqs] == ["a", "b", "c"], seqs
        assert seqs[0]["taxonomy"] == "Bacteria;Genus;a" and "genome" not in seqs[0], seqs

        m = sample_abundances(cfg)
        assert m.shape == (8, 3), m.shape
        assert (m.sum(axis=1) > 0).all(), "no all-zero community rows"
        assert np.array_equal(m, sample_abundances(cfg)), "seeded draw is reproducible"
        # per-genome override: 'a' has 10x the mean, so its column should dominate.
        assert m[:, 0].mean() > m[:, 1].mean(), "nb override lifts genome a"

        # Bernoulli presence gate: a per-genome `presence: 0` is always absent; a
        # global `presence: 0.5` still yields non-empty communities.
        p.write_text(cfg_text.replace(
            "sampling: {n_samples: 8, seed: 7,",
            "sampling: {n_samples: 30, seed: 7, presence: 0.5,").replace(
            "ssu: refs/c.16s.fa}", "ssu: refs/c.16s.fa, presence: 0}"))
        cfg_pres = load_config(p)
        mp = sample_abundances(cfg_pres)
        assert (mp[:, 2] == 0).all(), "presence 0 -> genome c always absent"
        assert (mp[:, :2] == 0).any(), "presence 0.5 -> some genomes gated absent"
        assert (mp.sum(axis=1) > 0).all(), "no all-zero community rows"

        aap = aap_settings(cfg)
        assert aap["aap_configs"] == [str(base / "engine.config"), "/abs/site.config"], aap
        assert aap["aap_profile"] == "singularity", aap

        # sylph requested but no `genome:` -> validation exits.
        p.write_text(cfg_text.replace("profilers: [aap]", "profilers: [sylph]"))
        try:
            load_config(p)
        except SystemExit as e:
            assert "genome" in str(e), e
        else:
            raise AssertionError("expected SystemExit for missing genome")

        # bad NB dispersion -> exits.
        p.write_text(cfg_text.replace("dispersion: 1.0", "dispersion: 0"))
        try:
            load_config(p)
        except SystemExit as e:
            assert "dispersion" in str(e), e
        else:
            raise AssertionError("expected SystemExit for dispersion <= 0")

        # generation_modes: fan-out source selection + per-mode reads merge.
        gm_text = textwrap.dedent("""
            train: {id: t, fastq_1: /abs/r1.fq, fastq_2: /abs/r2.fq, platform: hq-illumina}
            reads: {num_reads: 100, subsample: none, paired_end: true, read_length_mean: 300, read_length_variance: 0}
            sampling: {n_samples: 4, seed: 1, negative_binomial: {mean: 20, dispersion: 1.0}}
            database: {name: db, profilers: [sylph, aap], rfam_covariance_model: /abs/ribo, rfam_claninfo: /abs/ribo.clan}
            generation_modes:
              - {name: wgs, source: genome, mode: shotgun, profiler: sylph, reads: {read_length_mean: 150}}
              - {name: amp16s, source: ssu, mode: amplicon, profiler: aap, primers: [{pair_id: v4, forward: GTGYCAG, reverse: GGACTAC}]}
            panel:
              - {id: a, species: genus_a, genome: refs/a.fna, ssu: refs/a.16s.fa, amplicon: refs/a.amp.fa}
              - {id: b, species: genus_b, genome: refs/b.fna, ssu: refs/b.16s.fa, amplicon: refs/b.amp.fa}
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
    print("selfcheck OK")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--selfcheck":
        _selfcheck()
    else:
        sys.exit("nb_config is a library; run a generate_*/build_* script instead "
                 "(or `--selfcheck`).")
