#!/usr/bin/env python3
"""Shared config glue for the SILVA superresolution-amplicon parameter sweep.

Two sibling examples already do most of this, so they are imported rather than copied:

* ``abundance_nb_sample/scripts/nb_config.py`` — config schema, path resolution and the
  seeded negative-binomial community sampler. The communities here are that example's,
  unchanged: every panel genome presence-gated then NB-drawn, independently per sample.
* ``sr_amplicon_param_sweep/scripts/sr_sweep.py`` — expansion of ``sr_sweep.grid`` into
  the flat ``sr_settings:`` list the pipeline fans each row out over, and the count of
  distinct mis-mapping matrices that grid actually costs.

What is new is the database: this example profiles against a **pre-built** generic
reference set (a `path:` entry, SILVA SSU) instead of a collection the pipeline builds
from the samplesheet. SILVA is a label space, not a set of genomes, so inference runs in
V4-group space and is scored per genus (`scripts/score_sweep.py`): every panel genome
carries its SILVA lineage in `config.yaml`.

Two comparison axes sit on top of that one reference set, and both are scored in the same
genus space:

* **profiling arms** (`arms`) — the SILVA sweep itself, the same reads against a
  collection built from the panel (`custom_database`), SILVA reinterpreted through that
  panel (an `sr_settings` `panel:`), and no superresolution at all (the EBI
  amplicon-analysis-pipeline against `aap_database`).
* **error-model arms** (`error_models`) — the reads are *generated* three times, once per
  error model: the AIC-selected trained model, a context-free one, and skiver's bundled
  preset with no training. Each arm needs its own `train_id`, which is what the pipeline
  keys training by.

    python scripts/silva_sweep.py --selfcheck
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
mode_profilers = nc.mode_profilers
mode_source_field = nc.mode_source_field
mode_reads = nc.mode_reads
settings = sw.settings
n_matrices = sw.n_matrices

# SILVA puts a prokaryote's genus at the 6th rank: domain;phylum;class;order;family;genus.
# ponytail: prokaryotes only; eukaryote lineages are deeper and would need a rank map.
GENUS_RANK = 6


def genus(lineage):
    """The genus-level lineage prefix, or None when `lineage` stops above genus."""
    ranks = lineage.split(";")[:GENUS_RANK]
    if len(ranks) < GENUS_RANK or ranks[-1] in ("", "unclassified"):
        return None
    return ";".join(ranks)


# The four ways this example turns one set of reads into a genus profile. `silva` is the
# sweep proper (fanned over the whole grid); the other three run at one grid point each,
# so the arm comparison is not multiplied by the sweep.
#   custom      : sr_amplicon against a collection built from the panel - the upper bound
#                 a reference set that *is* the community gives.
#   silva_panel : sr_amplicon against SILVA, its labels reinterpreted through that same
#                 panel collection (an sr_settings `panel:`).
#   aap         : no superresolution - the EBI amplicon-analysis-pipeline's own MAPseq
#                 labels against a pre-built mapseq SILVA database.
SR_ARMS = ("silva", "custom", "silva_panel")
ARMS = SR_ARMS + ("aap",)


def error_models(cfg):
    """The error-model arms, `[{name, train_id, components?, preset?}]`.

    Absent `error_models:` = one unnamed arm, which is the pre-amendment behaviour: the
    train_id is `train.id` itself and the samplesheet carries no new columns.
    """
    raw = cfg.get("error_models")
    base = cfg["train"]["id"]
    if not raw:
        return [{"name": None, "train_id": base}]
    out = []
    for entry in raw:
        name = str(entry["name"])
        if entry.get("preset") and entry.get("components"):
            sys.exit(f"config.yaml: error_models '{name}' sets both preset and components")
        # One train_id per arm: training is deduped by it, so two arms sharing one id
        # would silently share one model.
        out.append({k: v for k, v in
                    [("name", name), ("train_id", f"{base}.{name}"),
                     ("components", entry.get("components")),
                     ("preset", entry.get("preset"))] if v is not None})
    names = [e["name"] for e in out]
    if len(set(names)) != len(names):
        sys.exit(f"config.yaml: duplicate error_models name(s) in {names}")
    return out


def arm_point(cfg):
    """The grid point the non-sweep arms run at, as a settings entry."""
    name = str(cfg.get("arm_point") or settings(cfg)[0]["name"])
    for entry in settings(cfg):
        if entry["name"] == name:
            return entry
    sys.exit(f"config.yaml: arm_point '{name}' is not a point of sr_sweep.grid "
             f"({', '.join(e['name'] for e in settings(cfg))})")


def arm_settings(cfg, arm):
    """The `sr_settings:` list a row of `arm` carries.

    `silva` takes the whole grid. `custom` and `silva_panel` take `arm_point` alone,
    renamed so their profiles do not collide with the sweep's in the shared benchmark
    dir - `<id>.<setting>.sr_profile.tsv` is keyed by the setting name.
    """
    if arm == "silva":
        return settings(cfg)
    point = dict(arm_point(cfg))
    point["name"] = f"{arm}.{point['name']}"
    if arm == "silva_panel":
        # A panel kernel is measured inside the inference run and records no distances.
        point.pop("infer_distance_decay", None)
        point["panel"] = cfg["custom_database"]["name"]
    return [point]


def arm_database(cfg, arm):
    """The `database:` name a row of `arm` profiles against."""
    return {"silva": cfg["database"]["name"],
            "silva_panel": cfg["database"]["name"],
            "custom": cfg["custom_database"]["name"],
            "aap": cfg["aap_database"]["name"]}[arm]


def load_config(path):
    """`nb_config.load_config`, plus this example's own required blocks."""
    cfg = nc.load_config(path)
    if "sr_sweep" not in cfg:
        sys.exit("config.yaml: missing required 'sr_sweep:' block")
    settings(cfg)                       # expands (and validates) the grid up front
    db = cfg["database"]
    if bool(db.get("path")) == bool(db.get("sequences_from_panel")):
        sys.exit("config.yaml: database needs either 'path:' (pre-built, e.g. SILVA) or "
                 "'sequences_from_panel: true' (built in-pipeline from the panel)")
    if db.get("path"):
        cfg["database"]["path"] = str(Path(str(db["path"])).expanduser())
    for m in cfg["panel"]:
        if not genus(str(m.get("taxonomy") or "")):
            sys.exit(f"config.yaml: panel '{m['id']}' needs a SILVA 'taxonomy:' lineage "
                     "down to genus")
    for block in ("custom_database", "aap_database"):
        if block not in cfg:
            sys.exit(f"config.yaml: missing required '{block}:' block")
    aap = cfg["aap_database"]
    for key in ("name", "path", "rfam_covariance_model", "rfam_claninfo"):
        if not aap.get(key):
            sys.exit(f"config.yaml: aap_database needs '{key}:' "
                     "(the pipeline requires both Rfam files for an 'aap' database)")
    aap["path"] = str(Path(str(aap["path"])).expanduser())
    error_models(cfg)                   # validates the arms up front
    arm_point(cfg)                      # ... and that the arms' grid point exists
    return cfg


def database_block(cfg):
    """The samplesheet `databases:` block.

    A `path:` database is handed to the pipeline as a pre-built directory: nothing is
    built, and BUILD_DATABASES resolves `<name>_ssu.sr_refs.{fasta,tax}` inside it.
    `sequences_from_panel: true` falls back to the sibling examples' in-pipeline build (a
    cheap control on the same code path), with each genome's SILVA lineage as its
    `taxonomy:`, so the control gets a `.tax` and scores the same way.
    """
    db = cfg["database"]
    lineage = {m["id"]: m["taxonomy"] for m in cfg["panel"]}
    if db.get("path"):
        block = {db["name"]: {"path": db["path"]}}
    else:
        block = nc.database_block(cfg)
        for seq in block[db["name"]]["sequences"]:
            seq["taxonomy"] = lineage[seq["id"]]
    # The panel collection: profiled against directly by the `custom` arm, and named as
    # the `panel:` of the `silva_panel` arm. One collection serves both.
    block[cfg["custom_database"]["name"]] = {
        "profilers": ["sr_amplicon"],
        "sequences": [{"id": m["id"], "ssu": m["ssu"], "taxonomy": m["taxonomy"]}
                      for m in cfg["panel"]],
    }
    # The mapseq database the amplicon-analysis-pipeline classifies against. Pre-built
    # out of band, like the superresolution reference set, plus the two Rfam files
    # main.nf requires of any 'aap' collection.
    aap = cfg["aap_database"]
    block[aap["name"]] = {"profilers": ["aap"], "path": aap["path"],
                          "rfam_covariance_model": aap["rfam_covariance_model"],
                          "rfam_claninfo": aap["rfam_claninfo"]}
    return block


def depths(cfg):
    """Subsample depths, `None` for full depth. One benchmark dir per depth."""
    value = cfg["reads"].get("subsample") or [None]
    if not isinstance(value, list):
        value = [value]
    return [None if d in (None, "none", "") else int(d) for d in value]


def sample_id(cfg, em, gm, i):
    """The phase-1 sample name. The error-model arm is part of it: each arm generates
    its own reads, so it is a sample of its own, not a profiling setting."""
    parts = [f"S{i:02d}"] + ([em["name"]] if em["name"] else []) + [gm["name"]]
    return ".".join(parts)


def benchmark_dirs(cfg, results_dir):
    """(sample, dir, depth, pair, error_model) per generated benchmark dir, mirroring
    the generate step's fan-out: one sample per error-model arm x community x mode x
    primer pair, one subdir per depth."""
    for em in error_models(cfg):
        for gm in generation_modes(cfg):
            for pair in gm.get("primers") or [None]:
                for i in range(1, cfg["sampling"]["n_samples"] + 1):
                    sample = sample_id(cfg, em, gm, i)
                    if pair:
                        sample = f"{sample}.{pair['pair_id']}"
                    for depth in depths(cfg):
                        directory = results_dir / sample
                        if depth:
                            directory = directory / f"subsample_{depth}"
                        yield sample, directory, depth, pair, em


def _selfcheck():
    bacteroides = "Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides"
    assert genus(bacteroides + ";unclassified;unclassified") == bacteroides
    assert genus("Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae") is None
    assert genus("Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;"
                 "unclassified") is None
    assert genus("unclassified_v4_group") is None

    cfg = {"database": {"name": "silva", "path": "/nonexistent/silva"},
           "custom_database": {"name": "ctl"},
           "aap_database": {"name": "silva_ms", "path": "/nonexistent/ms",
                            "rfam_covariance_model": "/rfam.cm",
                            "rfam_claninfo": "/rfam.claninfo"},
           "train": {"id": "t"},
           "reads": {"subsample": ["none", 100000]},
           "panel": [{"id": "bf", "taxonomy": bacteroides, "ssu": "bf.fasta"}]}
    block = database_block(cfg)
    panel_seqs = [{"id": "bf", "ssu": "bf.fasta", "taxonomy": bacteroides}]
    assert block == {"silva": {"path": "/nonexistent/silva"},
                     "ctl": {"profilers": ["sr_amplicon"], "sequences": panel_seqs},
                     "silva_ms": {"profilers": ["aap"], "path": "/nonexistent/ms",
                                  "rfam_covariance_model": "/rfam.cm",
                                  "rfam_claninfo": "/rfam.claninfo"}}, block
    assert depths(cfg) == [None, 100000], depths(cfg)
    assert depths({"reads": {}}) == [None]

    # Error-model arms: absent block = one unnamed arm on the bare train_id, which is
    # what keeps a samplesheet without `error_models:` byte-identical.
    assert error_models(cfg) == [{"name": None, "train_id": "t"}], error_models(cfg)
    cfg["error_models"] = [{"name": "trained"},
                           {"name": "flat", "components": "BaseContext(1)"},
                           {"name": "naive", "preset": "illumina"}]
    arms = error_models(cfg)
    assert [a["train_id"] for a in arms] == ["t.trained", "t.flat", "t.naive"], arms
    assert arms[1]["components"] == "BaseContext(1)" and arms[2]["preset"] == "illumina"
    gm = {"name": "amplicon_16s"}
    assert sample_id(cfg, arms[0], gm, 3) == "S03.trained.amplicon_16s"
    assert sample_id(cfg, {"name": None}, gm, 3) == "S03.amplicon_16s"

    # The grid expansion and matrix accounting are sr_sweep's, exercised here so a
    # change in the sibling example fails this example's selfcheck too.
    grid = {"sr_sweep": {"grid": {
        "matrix": [{"name": "exact", "mismapping_method": "align",
                    "align_backend": "exact-hash", "align_tau": 0},
                   {"name": "kmer1", "mismapping_method": "align",
                    "align_backend": "kmer", "align_tau": 1}],
        "infer": [{"name": "p01", "infer_presence_prior": 0.01},
                  {"name": "p001", "infer_presence_prior": 0.001}]}}}
    expanded = settings(grid)
    assert [e["name"] for e in expanded] == ["exact.p01", "exact.p001",
                                             "kmer1.p01", "kmer1.p001"], expanded
    assert n_matrices(expanded) == 2, n_matrices(expanded)

    # Profiling arms. Only `silva` takes the whole grid; the others take one point,
    # renamed so their profiles do not collide with the sweep's in the same dir.
    grid["arm_point"] = "kmer1.p001"
    grid["custom_database"] = {"name": "ctl"}
    grid["database"] = {"name": "silva"}
    grid["aap_database"] = {"name": "silva_ms"}
    assert [e["name"] for e in arm_settings(grid, "silva")] == [e["name"] for e in expanded]
    assert [e["name"] for e in arm_settings(grid, "custom")] == ["custom.kmer1.p001"]
    panel_pt = arm_settings(grid, "silva_panel")[0]
    assert panel_pt["name"] == "silva_panel.kmer1.p001" and panel_pt["panel"] == "ctl", panel_pt
    assert arm_database(grid, "aap") == "silva_ms"
    assert arm_database(grid, "silva_panel") == "silva"
    print("silva_sweep selfcheck OK")


if __name__ == "__main__":
    if "--selfcheck" in sys.argv:
        _selfcheck()
    else:
        print(__doc__)
