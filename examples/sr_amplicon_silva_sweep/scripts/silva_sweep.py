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
    if db.get("path"):
        return {db["name"]: {"path": db["path"]}}
    block = nc.database_block(cfg)
    lineage = {m["id"]: m["taxonomy"] for m in cfg["panel"]}
    for seq in block[db["name"]]["sequences"]:
        seq["taxonomy"] = lineage[seq["id"]]
    return block


def depths(cfg):
    """Subsample depths, `None` for full depth. One benchmark dir per depth."""
    value = cfg["reads"].get("subsample") or [None]
    if not isinstance(value, list):
        value = [value]
    return [None if d in (None, "none", "") else int(d) for d in value]


def benchmark_dirs(cfg, results_dir):
    """(sample, dir, depth, pair) per generated benchmark dir, mirroring the generate
    step's fan-out: one sample per community x mode x primer pair, one subdir per depth."""
    for gm in generation_modes(cfg):
        for pair in gm.get("primers") or [None]:
            for i in range(1, cfg["sampling"]["n_samples"] + 1):
                sample = f"S{i:02d}.{gm['name']}"
                if pair:
                    sample = f"{sample}.{pair['pair_id']}"
                for depth in depths(cfg):
                    directory = results_dir / sample
                    if depth:
                        directory = directory / f"subsample_{depth}"
                    yield sample, directory, depth, pair


def _selfcheck():
    bacteroides = "Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides"
    assert genus(bacteroides + ";unclassified;unclassified") == bacteroides
    assert genus("Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae") is None
    assert genus("Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;"
                 "unclassified") is None
    assert genus("unclassified_v4_group") is None

    cfg = {"database": {"name": "silva", "path": "/nonexistent/silva"},
           "reads": {"subsample": ["none", 100000]},
           "panel": [{"id": "bf", "taxonomy": bacteroides, "ssu": "bf.fasta"}]}
    assert database_block(cfg) == {"silva": {"path": "/nonexistent/silva"}}
    cfg["database"] = {"name": "ctl", "profilers": ["sr_amplicon"], "sequences_from_panel": True}
    assert database_block(cfg) == {"ctl": {"sequences": [
        {"id": "bf", "ssu": "bf.fasta", "taxonomy": bacteroides}]}}, database_block(cfg)
    assert depths(cfg) == [None, 100000], depths(cfg)
    assert depths({"reads": {}}) == [None]

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
    print("silva_sweep selfcheck OK")


if __name__ == "__main__":
    if "--selfcheck" in sys.argv:
        _selfcheck()
    else:
        print(__doc__)
