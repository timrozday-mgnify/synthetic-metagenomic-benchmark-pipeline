#!/usr/bin/env python3
"""Singularity/Apptainer version of build_profiling_dbs.py, for HPC where Docker
isn't available. Identical behaviour - builds a custom sylph (WGS) DB and a
custom mapseq (AAP amplicon) DB from this example's reference genomes - but runs
mapseq/sylph/barrnap via `singularity exec` instead of `docker run`.

Containers are pulled straight from the same quay.io biocontainers via the
`docker://` URI (Singularity/Apptainer converts them to SIF on first use). On an
air-gapped cluster, pre-build the SIF once on a login node
(`singularity pull mapseq.sif docker://quay.io/...`) and set the *_IMAGE
constants below to the local .sif paths instead.

- sylph DB: sketched directly from the full genome assemblies.
- mapseq DB: built from each genome's full-length 16S rRNA sequence(s) - NOT the
  V4 amplicon fragments in ../references/*.amplicons.fasta, mapseq needs
  full-length 16S to classify correctly. If you don't have a pre-extracted
  full-length 16S for a genome, leave its GENOMES entry's ssu_fasta as None and
  this script predicts one from the genome with barrnap.

Fill in GENOMES below, then run once:
    python scripts/build_profiling_dbs_singularity.py
Writes genome_references/, mapseq_references/, mapseq_db/, sylph_db/, and two
example config snippets (sylph_databases.config, aap.config) next to this
script's parent directory.
"""
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUN_DIR = HERE.parent
GENOME_REFS_DIR = RUN_DIR / "genome_references"
SSU_REFS_DIR = RUN_DIR / "mapseq_references"
MAPSEQ_DB_DIR = RUN_DIR / "mapseq_db"
SYLPH_DB_DIR = RUN_DIR / "sylph_db"

# docker:// URIs (converted to SIF on first use); swap for local .sif paths on
# an air-gapped cluster. Same biocontainers as the Docker version.
MAPSEQ_IMAGE = "docker://quay.io/biocontainers/mapseq:2.1.1b--h3ab3c3b_0"     # matches modules/ebi-metagenomics/mapseq
SYLPH_IMAGE = "docker://quay.io/biocontainers/sylph:0.9.0--ha6fb395_0"        # matches modules/local/sylph/build_db
BARRNAP_IMAGE = "docker://quay.io/biocontainers/barrnap:0.9--hdfd78af_4"      # fallback 16S extraction, not an nf-core module here

# `singularity` or `apptainer` - both accept the same exec/bind/pwd flags.
SINGULARITY = "singularity"

# --- Fill these in ---------------------------------------------------------
# One entry per genome used in ../generate_sweep.py's PANEL. genome_id is
# "genus_species" (matches PANEL genome_ids); override GENUS_SPECIES below for
# any genome_id that doesn't split cleanly (e.g. an extra strain like
# "species_strain2"). ssu_fasta=None -> predict full-length 16S with barrnap.
GENOMES = {
    "bacteroides_fragilis": dict(genome_fasta="/PATH/TO/bacteroides_fragilis.fna", ssu_fasta=None),
    "bacteroides_thetaiotaomicron": dict(genome_fasta="/PATH/TO/bacteroides_thetaiotaomicron.fna", ssu_fasta=None),
    "bacteroides_uniformis": dict(genome_fasta="/PATH/TO/bacteroides_uniformis.fna", ssu_fasta=None),
    # ... one entry per PANEL genome_id ...
    "bacteroides_uniformis_strain2": dict(
        genome_fasta="/PATH/TO/bacteroides_uniformis_strain2.fna", ssu_fasta=None,
    ),
}
GENUS_SPECIES_OVERRIDES = {"bacteroides_uniformis_strain2": ("bacteroides", "uniformis")}
KINGDOM_OVERRIDES = {}  # genome_id -> "arc" for Archaea; rest default to "bac" (barrnap --kingdom)
# ---------------------------------------------------------------------------


def sh(cmd, **kw):
    print("+", " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True, **kw)


def container_run(image, args, workdir):
    # Bind workdir to /data and set pwd there so the callers' /data paths work
    # unchanged, exactly as the Docker `-v {workdir}:/data -w /data` version did.
    sh([SINGULARITY, "exec", "--bind", f"{workdir}:/data", "--pwd", "/data", image, *args])


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


def genus_species(genome_id):
    if genome_id in GENUS_SPECIES_OVERRIDES:
        return GENUS_SPECIES_OVERRIDES[genome_id]
    genus, species = genome_id.split("_", 1)
    return genus, species


def stage_genome(genome_id):
    """Decompress/copy this genome's assembly into genome_references/{genome_id}.fasta."""
    out = GENOME_REFS_DIR / f"{genome_id}.fasta"
    src = Path(str(GENOMES[genome_id]["genome_fasta"])).expanduser()
    if src.suffix == ".gz":
        with open(out, "wb") as fh:
            sh(["gzip", "-dc", str(src)], stdout=fh)
    else:
        shutil.copy(src, out)
    return out


def stage_ssu(genome_id, genome_fasta):
    """Write mapseq_references/{genome_id}.ssu.fasta: full-length 16S copies, headers
    prefixed with genome_id so cluster/taxonomy assignment is traceable per genome."""
    out = SSU_REFS_DIR / f"{genome_id}.ssu.fasta"
    ssu_fasta = GENOMES[genome_id]["ssu_fasta"]

    if ssu_fasta and Path(ssu_fasta).expanduser().exists():
        records = list(parse_fasta(Path(ssu_fasta).expanduser()))
    else:
        # No pre-extracted full-length 16S for this genome - predict it directly.
        print(f"  no ssu_fasta for {genome_id}; running barrnap on the genome")
        gff = GENOME_REFS_DIR / f"{genome_id}.barrnap.gff3"
        rrna_fa = GENOME_REFS_DIR / f"{genome_id}.barrnap.fasta"
        container_run(BARRNAP_IMAGE, [
            "barrnap", "--kingdom", KINGDOM_OVERRIDES.get(genome_id, "bac"), "--quiet",
            "--outseq", f"/data/{rrna_fa.name}", f"/data/{genome_fasta.name}",
        ], workdir=str(GENOME_REFS_DIR))
        records = [(h, s) for h, s in parse_fasta(rrna_fa) if h.startswith("16S_rRNA")]
        if not records:
            raise RuntimeError(f"barrnap found no 16S_rRNA hits for {genome_id}")
        gff.unlink(missing_ok=True)
        rrna_fa.unlink()
        (genome_fasta.with_suffix(".fasta.fai")).unlink(missing_ok=True)

    with open(out, "w") as fh:
        for i, (header, seq) in enumerate(records):
            fh.write(f">{genome_id}|{i}|{header}\n{seq}\n")
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
    genus, species = genus_species(genome_id)
    kingdom = "Archaea" if KINGDOM_OVERRIDES.get(genome_id) == "arc" else "Bacteria"
    return f"{kingdom};{genus.capitalize()};{species}"


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


def write_configs(syldb_path, mapseq_paths):
    mapseq_fasta, mapseq_tax, mapseq_otu, mapseq_mscluster = mapseq_paths

    (RUN_DIR / "sylph_databases.config").write_text(f"""\
// Example config for `--step profile` (or `all`) with a per-sample `database`
// column set to 'community_v4' instead of 'self' or a production GTDB db.
params {{
    sylph_databases = [
        community_v4: [ syldb: '{syldb_path}', label: 'subspecies_v4_sweep community' ],
    ]
}}
""")

    (RUN_DIR / "aap.config").write_text(f"""\
// Example --aap_config for `profiler=aap`: a custom mapseq DB built from this
// example's own reference genomes' full-length 16S, instead of production SILVA.
params {{
    mapseq_databases {{
        community_v4 {{
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

    global MAPSEQ_DB_DIR
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
    if len(sys.argv) > 1 and sys.argv[1] == "--selfcheck":
        _selfcheck()
        return
    for d in (GENOME_REFS_DIR, SSU_REFS_DIR):
        d.mkdir(exist_ok=True)

    genome_fastas, ssu_fastas = [], []
    for genome_id in GENOMES:
        print(f"staging {genome_id}")
        genome_fasta = stage_genome(genome_id)
        genome_fastas.append(genome_fasta)
        ssu_fastas.append(stage_ssu(genome_id, genome_fasta))

    print("building sylph db")
    syldb_path = build_sylph_db(genome_fastas)

    print("building mapseq db")
    mapseq_paths = build_mapseq_db(ssu_fastas)

    write_configs(syldb_path, mapseq_paths)
    print(f"\nWrote {syldb_path}")
    print(f"Wrote {mapseq_paths[0]} (+ .tax/.otu/.mscluster)")
    print("Wrote sylph_databases.config, aap.config")


if __name__ == "__main__":
    main()
