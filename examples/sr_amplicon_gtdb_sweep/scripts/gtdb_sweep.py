#!/usr/bin/env python3
"""Shared config glue for the GTDB superresolution-amplicon parameter sweep.

Two sibling examples already do most of this, so they are imported rather than copied:

* ``abundance_nb_sample/scripts/nb_config.py`` — config schema, path resolution and the
  seeded negative-binomial community sampler. The communities here are that example's,
  unchanged: every panel genome presence-gated then NB-drawn, independently per sample.
* ``sr_amplicon_param_sweep/scripts/sr_sweep.py`` — expansion of ``sr_sweep.grid`` into
  the flat ``sr_settings:`` list the pipeline fans each row out over, and the count of
  distinct mis-mapping matrices that grid actually costs.

What is genuinely new is the database: this example profiles against a **pre-built**
reference set (a `path:` entry, i.e. GTDB) instead of a collection the pipeline builds
from the samplesheet, and the panel ids therefore have to be the ids that reference set
uses — see `check_reference_ids`.

    python scripts/gtdb_sweep.py --selfcheck
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

# superresolution-amplicon's reference FASTA inside a pre-built database directory.
# BUILD_DATABASES resolves it by this glob (`*_<source>.sr_refs.fasta`), so the file
# must be named after the collection: `<name>_ssu.sr_refs.fasta`.
SR_REFS_GLOB = "*_ssu.sr_refs.fasta"


def load_config(path):
    """`nb_config.load_config`, plus this example's own required blocks."""
    cfg = nc.load_config(path)
    if "sr_sweep" not in cfg:
        sys.exit("config.yaml: missing required 'sr_sweep:' block")
    settings(cfg)                       # expands (and validates) the grid up front
    db = cfg["database"]
    if bool(db.get("path")) == bool(db.get("sequences_from_panel")):
        sys.exit("config.yaml: database needs either 'path:' (pre-built, e.g. GTDB) or "
                 "'sequences_from_panel: true' (built in-pipeline from the panel)")
    if db.get("path"):
        cfg["database"]["path"] = str(Path(str(db["path"])).expanduser())
    return cfg


def database_block(cfg):
    """The samplesheet `databases:` block.

    A `path:` database is handed to the pipeline as a pre-built directory: nothing is
    built, and BUILD_DATABASES just resolves `<name>_ssu.sr_refs.fasta` inside it. That
    is the whole point of this example — a reference set far too large to build per run.
    `sequences_from_panel: true` falls back to the sibling examples' in-pipeline build,
    which is the way to get a small control run against the same code path.
    """
    db = cfg["database"]
    if db.get("path"):
        return {db["name"]: {"path": db["path"]}}
    return nc.database_block(cfg)


def depths(cfg):
    """Subsample depths, `None` for full depth. One benchmark dir per depth."""
    value = cfg["reads"].get("subsample") or [None]
    if not isinstance(value, list):
        value = [value]
    return [None if d in (None, "none", "") else int(d) for d in value]


def check_reference_ids(cfg):
    """Warn when a panel id is not a reference id in the pre-built FASTA.

    Truth is written per panel `id`, and a superresolution profile is written per
    reference id (the `{genome_id}|{copy}|{orig}` prefix in the refs FASTA). If the two
    disagree, the run still succeeds and every single genome scores zero — the panel has
    to be named in the reference set's own ids (GTDB accessions, for GTDB).

    Returns the ids that are missing, or None when the check could not run (no `path:`,
    directory absent — the shipped dummy config).
    """
    db = cfg["database"]
    if not db.get("path"):
        return None
    matches = sorted(Path(db["path"]).glob(SR_REFS_GLOB)) if Path(db["path"]).is_dir() else []
    if len(matches) != 1:
        return None
    reference_ids = set()
    with open(matches[0]) as fh:
        for line in fh:
            if line.startswith(">"):
                reference_ids.add(line[1:].split("|", 1)[0].strip())
    return [m["id"] for m in cfg["panel"] if m["id"] not in reference_ids]


def _selfcheck():
    cfg = {"database": {"name": "gtdb_ssu", "path": "/nonexistent/gtdb"},
           "reads": {"subsample": ["none", 100000]},
           "panel": [{"id": "GCF_1"}, {"id": "GCF_2"}]}
    assert database_block(cfg) == {"gtdb_ssu": {"path": "/nonexistent/gtdb"}}
    assert depths(cfg) == [None, 100000], depths(cfg)
    assert depths({"reads": {}}) == [None]
    # No reference FASTA to read => "could not check", which is not "all ids missing".
    assert check_reference_ids(cfg) is None

    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        (Path(tmp) / "gtdb_ssu_ssu.sr_refs.fasta").write_text(
            ">GCF_1|0|orig\nACGT\n>GCF_3|0|orig\nACGT\n")
        cfg["database"]["path"] = tmp
        assert check_reference_ids(cfg) == ["GCF_2"], check_reference_ids(cfg)

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
    print("gtdb_sweep selfcheck OK")


if __name__ == "__main__":
    if "--selfcheck" in sys.argv:
        _selfcheck()
    else:
        print(__doc__)
