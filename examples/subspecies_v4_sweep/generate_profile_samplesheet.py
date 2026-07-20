#!/usr/bin/env python3
"""Generate a profile-only samplesheet (`--step profile`) for reads already
written by an earlier `--step generate` run, profiling every sample against a
custom `community_v4` database built from this example's own reference genomes.

The database ships *inside* the samplesheet as a top-level `databases:` block of
named sequence collections; the pipeline builds the sylph `.syldb` and the mapseq
DB itself during the profile run (no separate build_profiling_dbs.py / docker /
singularity step). Each sample is profiled with every profiler in PROFILERS.

Per-genome genome/ssu/taxonomy come from scripts/build_profiling_dbs.py's GENOMES
(the single place this example lists its reference genomes) - fill that in first.
Then:
    python generate_profile_samplesheet.py [results_dir]
    nextflow run ../../main.nf -profile docker --step profile \\
        --input profile_samplesheet.yaml --outdir <results_dir>
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "scripts"))
import build_profiling_dbs as refs  # GENOMES + tax_string(); no side effects on import

DEFAULT_RESULTS_DIR = HERE.parent.parent / "results" / "subspecies_v4_sweep"

# --- Fill these in ---------------------------------------------------------
DATABASE = "community_v4"        # name of the in-samplesheet collection
PROFILERS = ["sylph", "aap"]     # sylph uses `genome`; aap uses `ssu` + `taxonomy`
# ---------------------------------------------------------------------------


def sample_names(samplesheet):
    """Sample ids from a generate samplesheet.yaml (its `- sample:` lines)."""
    return [
        line.split(":", 1)[1].strip()
        for line in samplesheet.read_text().splitlines()
        if line.startswith("- sample:")
    ]


def database_block(genomes, tax_of, database, need_ssu):
    """`databases: <database>:` sequence collection built from `genomes`
    ({genome_id: {genome_fasta, ssu_fasta}}). `ssu`/`taxonomy` are only emitted
    (and required) when an aap/mapseq DB is needed."""
    lines = ["databases:", f"  {database}:", "    sequences:"]
    for gid, entry in genomes.items():
        genome = Path(str(entry["genome_fasta"])).expanduser()
        lines += [f"      - id: {gid}", f"        genome: {genome}"]
        if need_ssu:
            ssu = entry["ssu_fasta"]
            if not ssu:
                sys.exit(f"'aap' profiling needs a pre-extracted full-length 16S per "
                         f"genome, but '{gid}' has ssu_fasta=None in "
                         f"scripts/build_profiling_dbs.py GENOMES. Provide one, or drop "
                         f"'aap' from PROFILERS.")
            lines += [f"        ssu: {Path(str(ssu)).expanduser()}",
                      f'        taxonomy: "{tax_of(gid)}"']
    return lines


def sample_rows(samples, profilers, results_dir, database):
    """One `samples:` row per (sample, profiler); reads live in results_dir/sample."""
    lines = ["samples:"]
    for sample in samples:
        for profiler in profilers:
            lines += [f"  - sample: {sample}",
                      f"    profiler: {profiler}",
                      f"    benchmark_dir: {Path(results_dir) / sample}",
                      f"    database: {database}",
                      ""]
    return lines


def main():
    results_dir = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else DEFAULT_RESULTS_DIR
    samples = sample_names(HERE / "samplesheet.yaml")

    lines = database_block(refs.GENOMES, refs.tax_string, DATABASE, need_ssu="aap" in PROFILERS)
    lines.append("")
    lines += sample_rows(samples, PROFILERS, results_dir, DATABASE)

    (HERE / "profile_samplesheet.yaml").write_text("\n".join(lines) + "\n")
    print(f"Wrote profile_samplesheet.yaml: {len(samples)} samples x {PROFILERS} against "
          f"'{DATABASE}' (built in-pipeline from {len(refs.GENOMES)} genomes), "
          f"benchmark_dir root {results_dir}")


def _selfcheck():
    """Runnable check for the databases-block / sample-row emission (no pipeline)."""
    genomes = {"g1": dict(genome_fasta="/g1.fna", ssu_fasta="/g1.16s.fa"),
               "g2": dict(genome_fasta="/g2.fna", ssu_fasta="/g2.16s.fa")}
    tax_of = lambda gid: f"Bacteria;Genus;{gid}"

    block = database_block(genomes, tax_of, "community_v4", need_ssu=True)
    assert block[:3] == ["databases:", "  community_v4:", "    sequences:"], block
    assert "      - id: g1" in block and "        genome: /g1.fna" in block, block
    assert '        taxonomy: "Bacteria;Genus;g1"' in block, block
    # sylph-only run omits ssu/taxonomy.
    assert not any("ssu" in l for l in database_block(genomes, tax_of, "d", need_ssu=False))

    rows = sample_rows(["S01", "S02"], ["sylph", "aap"], "/res", "community_v4")
    assert sum(l.startswith("  - sample:") for l in rows) == 4, rows  # 2 samples x 2 profilers
    assert "    profiler: sylph" in rows and "    profiler: aap" in rows, rows
    assert "    benchmark_dir: /res/S01" in rows, rows
    print("selfcheck OK")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--selfcheck":
        _selfcheck()
    else:
        main()
