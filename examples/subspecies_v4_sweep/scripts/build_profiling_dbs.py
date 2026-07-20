#!/usr/bin/env python3
"""OPTIONAL / out-of-band builder for a custom sylph (WGS) and/or mapseq (AAP
amplicon) profiler database, driven by ../config.yaml.

The DEFAULT path does NOT need this script: generate_sweep.py emits a `databases:`
block into samplesheet.yaml and the pipeline (BUILD_DATABASES) builds the DB itself
during the run. Use this only if you want the DB prebuilt out-of-band; it also
writes sylph_databases.config / aap.config for a config-based `database:` run.

Containers run under Docker (default) or Singularity/Apptainer (`--runtime
singularity`, for HPC without Docker). Under Singularity the `docker://` images are
pulled to a SIF cache once (many HPC installs won't resolve docker:// at exec
time); for an air-gapped cluster, pre-pull the SIF and point the *_IMAGE constants
at the local .sif paths.

Which DBs are built follows config.yaml's `database.profilers`:
- sylph DB (if 'sylph'): sketched from each panel member's `genome`.
- mapseq DB (if 'aap'): built from each member's `ssu` (full-length 16S; NOT the V4
  amplicon fragments). If a member has no `ssu` file on disk, its 16S is predicted
  from `genome` with barrnap.

    python scripts/build_profiling_dbs.py [config.yaml] [--runtime docker|singularity]
Writes genome_references/, mapseq_references/, mapseq_db/, sylph_db/, and (for the
built DBs) sylph_databases.config / aap.config next to this script's parent.
"""
import argparse
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sweep_config as sc

HERE = Path(__file__).resolve().parent
RUN_DIR = HERE.parent
CONFIG = RUN_DIR / "config.yaml"
GENOME_REFS_DIR = RUN_DIR / "genome_references"
SSU_REFS_DIR = RUN_DIR / "mapseq_references"
MAPSEQ_DB_DIR = RUN_DIR / "mapseq_db"
SYLPH_DB_DIR = RUN_DIR / "sylph_db"

# docker:// biocontainers (match modules/ebi-metagenomics/mapseq and
# modules/local/sylph/build_db). Under --runtime docker the prefix is stripped;
# under singularity it's pulled to a SIF. Swap for local .sif paths if air-gapped.
MAPSEQ_IMAGE = "docker://quay.io/biocontainers/mapseq:2.1.1b--h3ab3c3b_0"
SYLPH_IMAGE = "docker://quay.io/biocontainers/sylph:0.9.0--ha6fb395_0"
BARRNAP_IMAGE = "docker://quay.io/biocontainers/barrnap:0.9--hdfd78af_4"  # not an nf-core module here

# Container runtime, set from --runtime in main(). `singularity` may be `apptainer`.
RUNTIME = "docker"
SINGULARITY = "singularity"
# Pulled .sif cache (singularity only). Reuses $SINGULARITY_CACHEDIR if set, else
# sif_cache/ in the run dir. Point at a shared path on HPC to avoid re-pulling.
SIF_CACHE_DIR = Path(os.environ.get("SINGULARITY_CACHEDIR", RUN_DIR / "sif_cache"))

# genome_id -> config panel member; populated in main(), used by tax_string().
BY_ID = {}


def sh(cmd, **kw):
    print("+", " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True, **kw)


def _sif(image):
    """Local path Singularity can exec: local .sif passes through; docker:// URIs
    are pulled to SIF_CACHE_DIR once."""
    if not image.startswith("docker://"):
        return image
    SIF_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    sif = SIF_CACHE_DIR / (image.removeprefix("docker://").replace("/", "_").replace(":", "_") + ".sif")
    if not sif.exists():
        sh([SINGULARITY, "pull", str(sif), image])
    return str(sif)


def container_run(image, args, workdir):
    """Run `args` in `image` with `workdir` bound to /data (cwd /data), under the
    selected runtime. Docker and Singularity get the same /data view so callers'
    /data paths are identical."""
    if RUNTIME == "singularity":
        sh([SINGULARITY, "exec", "--bind", f"{workdir}:/data", "--pwd", "/data", _sif(image), *args])
    else:
        sh(["docker", "run", "--rm", "-v", f"{workdir}:/data", "-w", "/data",
            image.removeprefix("docker://"), *args])


def parse_fasta(path):
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks)


def stage_genome(member):
    """Decompress/copy this genome's assembly into genome_references/{id}.fasta."""
    out = GENOME_REFS_DIR / f"{member['id']}.fasta"
    src = Path(member["genome"])
    if src.suffix == ".gz":
        with open(out, "wb") as fh:
            sh(["gzip", "-dc", str(src)], stdout=fh)
    else:
        shutil.copy(src, out)
    return out


def stage_ssu(member):
    """Write mapseq_references/{id}.ssu.fasta: full-length 16S copies, headers
    prefixed with genome_id so cluster/taxonomy assignment is traceable per genome.
    Falls back to barrnap on the genome if no pre-extracted 16S file is present."""
    gid = member["id"]
    out = SSU_REFS_DIR / f"{gid}.ssu.fasta"
    ssu = member.get("ssu")

    if ssu and Path(ssu).exists():
        records = list(parse_fasta(ssu))
    else:
        if not member.get("genome"):
            raise RuntimeError(f"{gid}: no ssu file at {ssu!r} and no `genome:` to run barrnap")
        print(f"  no pre-extracted SSU for {gid} ({ssu}); running barrnap on the genome")
        genome_fasta = stage_genome(member)
        rrna_fa = GENOME_REFS_DIR / f"{gid}.barrnap.fasta"
        kingdom = "arc" if str(member.get("kingdom", "")).lower() == "archaea" else "bac"
        container_run(BARRNAP_IMAGE, [
            "barrnap", "--kingdom", kingdom, "--quiet",
            "--outseq", f"/data/{rrna_fa.name}", f"/data/{genome_fasta.name}",
        ], workdir=str(GENOME_REFS_DIR))
        records = [(h, s) for h, s in parse_fasta(rrna_fa) if h.startswith("16S_rRNA")]
        if not records:
            raise RuntimeError(f"barrnap found no 16S_rRNA hits for {gid}")
        rrna_fa.unlink()
        (genome_fasta.with_suffix(".fasta.fai")).unlink(missing_ok=True)

    with open(out, "w") as fh:
        for i, (header, seq) in enumerate(records):
            fh.write(f">{gid}|{i}|{header}\n{seq}\n")
    return out


def build_sylph_db(genome_fastas):
    SYLPH_DB_DIR.mkdir(exist_ok=True)
    prefix = "community"
    for f in genome_fastas:
        shutil.copy(f, SYLPH_DB_DIR / f.name)
    container_run(SYLPH_IMAGE, [
        "sylph", "sketch", "-g", *[f.name for f in genome_fastas], "-o", prefix, "-t", "4",
    ], workdir=str(SYLPH_DB_DIR))
    for f in genome_fastas:
        (SYLPH_DB_DIR / f.name).unlink()
    return SYLPH_DB_DIR / f"{prefix}.syldb"


TAX_LEVELS = ["Kingdom", "Genus", "Species"]
TAX_CUTOFFS = "0.00:0.08 0.85:0.65 0.95:0.85"  # loose defaults; only 3 ranks used here


def tax_string(genome_id):
    return sc.taxonomy(BY_ID[genome_id])


def build_mapseq_tax(headers_by_genome):
    """headers_by_genome: genome_id -> [fasta header, ...] (already genome_id-prefixed)."""
    tax_path = MAPSEQ_DB_DIR / "mapseq_db.tax"
    with open(tax_path, "w") as fh:
        fh.write(f"#cutoff: {TAX_CUTOFFS}\n")
        fh.write("#name: custom\n")
        fh.write(f"#levels: {' '.join(TAX_LEVELS)}\n")
        for genome_id, headers in headers_by_genome.items():
            for h in headers:
                fh.write(f"{h}\t{tax_string(genome_id)}\n")
    return tax_path


def parse_mscluster(path):
    """Yield (cluster_id, [member seq_index, ...]) from a mapseq .mscluster file
    (one line per cluster: `<cluster_id> <member seq_index> <member seq_index> ...`)."""
    with open(path) as fh:
        for line in fh:
            cluster_i, *member_is = (int(x) for x in line.split())
            yield cluster_i, member_is


def build_mapseq_otu(mscluster_path, fasta_order, genome_of_header):
    """Majority-vote each mapseq-assigned cluster to one genome's taxonomy string."""
    members = {
        cluster_i: [genome_of_header[fasta_order[seq_i]] for seq_i in member_is]
        for cluster_i, member_is in parse_mscluster(mscluster_path)
    }
    otu_path = MAPSEQ_DB_DIR / "mapseq_db.otu"
    with open(otu_path, "w") as fh:
        for cluster_i in sorted(members):
            genome_ids = members[cluster_i]
            winner = max(set(genome_ids), key=genome_ids.count)
            # ponytail: real SILVA .otu files carry an NCBI taxid in the 3rd column;
            # we have no such ID for a custom db, so this is just a stable placeholder.
            fh.write(f"{cluster_i}\t{tax_string(winner)}\t{cluster_i}\n")
    return otu_path


def build_mapseq_db(ssu_fastas):
    MAPSEQ_DB_DIR.mkdir(exist_ok=True)
    db_fasta = MAPSEQ_DB_DIR / "mapseq_db.fasta"
    fasta_order, genome_of_header, headers_by_genome = [], {}, defaultdict(list)
    with open(db_fasta, "w") as out:
        for f in ssu_fastas:
            genome_id = f.name.removesuffix(".ssu.fasta")
            for header, seq in parse_fasta(f):
                out.write(f">{header}\n{seq}\n")
                fasta_order.append(header)
                genome_of_header[header] = genome_id
                headers_by_genome[genome_id].append(header)

    build_mapseq_tax(headers_by_genome)

    # Self-search to force mapseq to build & cache `mapseq_db.fasta.mscluster`
    # (no clustering file -> mapseq clusters the db on this first run).
    container_run(MAPSEQ_IMAGE, [
        "mapseq", db_fasta.name, db_fasta.name, "mapseq_db.tax", "-nthreads", "4",
    ], workdir=str(MAPSEQ_DB_DIR))
    mscluster_path = MAPSEQ_DB_DIR / "mapseq_db.fasta.mscluster"

    build_mapseq_otu(mscluster_path, fasta_order, genome_of_header)
    return db_fasta, MAPSEQ_DB_DIR / "mapseq_db.tax", MAPSEQ_DB_DIR / "mapseq_db.otu", mscluster_path


def write_configs(db_name, syldb_path, mapseq_paths):
    if syldb_path:
        (RUN_DIR / "sylph_databases.config").write_text(f"""\
// Config for `--step profile` (or `all`) with a per-sample `database` column set
// to '{db_name}' instead of 'self' or a production GTDB db.
params {{
    sylph_databases = [
        {db_name}: [ syldb: '{syldb_path}', label: 'subspecies_v4_sweep community' ],
    ]
}}
""")
    if mapseq_paths:
        mapseq_fasta, mapseq_tax, mapseq_otu, mapseq_mscluster = mapseq_paths
        (RUN_DIR / "aap.config").write_text(f"""\
// --aap_config for `profiler=aap`: a custom mapseq DB built from this example's own
// reference genomes' full-length 16S, instead of production SILVA.
params {{
    mapseq_databases {{
        {db_name} {{
            fasta = '{mapseq_fasta}'
            tax = '{mapseq_tax}'
            otu = '{mapseq_otu}'
            mscluster = '{mapseq_mscluster}'
            label = 'subspecies_v4_sweep community'
            run_otu = true
            run_asv = false
        }}
    }}
}}
""")


def _selfcheck():
    """Runnable check for the .mscluster-parsing / majority-vote OTU logic (no container)."""
    import tempfile

    global MAPSEQ_DB_DIR, BY_ID
    BY_ID = {gid: {"id": gid, "species": gid} for gid in
             ("bacteroides_fragilis", "clostridium_bolteae", "clostridium_ramosum")}
    fasta_order = ["h0", "h1", "h2", "h3"]
    genome_of_header = {"h0": "bacteroides_fragilis", "h1": "bacteroides_fragilis",
                        "h2": "clostridium_bolteae", "h3": "clostridium_ramosum"}

    with tempfile.TemporaryDirectory() as d:
        mscluster_path = Path(d) / "test.mscluster"
        mscluster_path.write_text("0 0 1 2\n1 3\n")

        clusters = dict(parse_mscluster(mscluster_path))
        assert clusters == {0: [0, 1, 2], 1: [3]}, clusters

        saved, MAPSEQ_DB_DIR = MAPSEQ_DB_DIR, Path(d)
        otu_path = build_mapseq_otu(mscluster_path, fasta_order, genome_of_header)
        MAPSEQ_DB_DIR = saved

        rows = [line.split("\t") for line in otu_path.read_text().splitlines()]
    # cluster 0 = [fragilis, fragilis, bolteae] -> majority fragilis; cluster 1 = [ramosum].
    assert rows[0][0] == "0" and rows[0][1] == "Bacteria;Bacteroides;fragilis", rows
    assert rows[1][0] == "1" and rows[1][1] == "Bacteria;Clostridium;ramosum", rows
    print("selfcheck OK")


def main():
    ap = argparse.ArgumentParser(description="Build custom sylph/mapseq profiler DBs from config.yaml.")
    ap.add_argument("config", nargs="?", default=str(CONFIG), help="config.yaml (default: ../config.yaml)")
    ap.add_argument("--runtime", choices=["docker", "singularity"], default="docker",
                    help="container runtime (default: docker)")
    ap.add_argument("--selfcheck", action="store_true", help="run the OTU-logic self-check and exit")
    args = ap.parse_args()

    if args.selfcheck:
        _selfcheck()
        return

    global RUNTIME, BY_ID
    RUNTIME = args.runtime
    cfg = sc.load_config(args.config)
    BY_ID = {m["id"]: m for m in cfg["panel"]}
    profilers = cfg["database"]["profilers"]
    db_name = cfg["database"]["name"]
    for d in (GENOME_REFS_DIR, SSU_REFS_DIR):
        d.mkdir(exist_ok=True)

    genome_fastas = []
    if "sylph" in profilers:
        for m in cfg["panel"]:
            print(f"staging genome {m['id']}")
            genome_fastas.append(stage_genome(m))
    ssu_fastas = []
    if "aap" in profilers:
        for m in cfg["panel"]:
            print(f"staging 16S {m['id']}")
            ssu_fastas.append(stage_ssu(m))

    syldb_path = build_sylph_db(genome_fastas) if genome_fastas else None
    mapseq_paths = build_mapseq_db(ssu_fastas) if ssu_fastas else None
    write_configs(db_name, syldb_path, mapseq_paths)

    if syldb_path:
        print(f"\nWrote {syldb_path} (+ sylph_databases.config)")
    if mapseq_paths:
        print(f"Wrote {mapseq_paths[0]} (+ .tax/.otu/.mscluster, aap.config)")


if __name__ == "__main__":
    main()
