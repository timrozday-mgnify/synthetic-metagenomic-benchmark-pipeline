#!/usr/bin/env python3
"""Split a batched AmpliconHunter amplicons.fa back into per-genome FASTAs.

AmpliconHunter runs once per primer pair over ALL of a sample's genomes and emits
one amplicons.fa whose headers carry `source=<fasta_basename>`. The benchmark needs
per-genome fragments to preserve each genome's abundance and ground truth, so this
groups amplicons by their source genome, rewrites headers to unique
`>{genome_id}_amplicon{N}` names, and emits a new genomes CSV pointing each row at
that genome's amplicon FASTA (carrying the original abundance forward).

    build_amplicon_genomes_csv.py AMPLICONS_FA GENOMES_CSV OUT_CSV OUT_DIR
"""

import csv
import os
import re
import sys

_SOURCE_RE = re.compile(r"source=([^.\s]+(?:\.[^.\s]+)*?)\.coordinates=")


def _strip_gz(name):
    return name[:-3] if name.endswith(".gz") else name


def _source_key(header):
    """Extract the source genome basename from an AmpliconHunter FASTA header.

    Header form: `>ENA|...| ...source=AAXE02.fasta.coordinates=2445-2737...`.
    Falls back to the first `source=` token split on `.coordinates` if the tight
    regex misses (defensive against header-format drift).
    """
    m = _SOURCE_RE.search(header)
    if m:
        return _strip_gz(m.group(1))
    if "source=" in header:
        tail = header.split("source=", 1)[1]
        return _strip_gz(tail.split(".coordinates", 1)[0].split()[0])
    return None


def _iter_fasta(lines):
    """Yield (header_without_gt, sequence) from FASTA text lines."""
    header, seq = None, []
    for line in lines:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq)
            header, seq = line[1:], []
        elif header is not None:
            seq.append(line)
    if header is not None:
        yield header, "".join(seq)


def build(amplicons_lines, genomes_csv_rows):
    """Return (out_header, out_rows, per_genome_fastas).

    genomes_csv_rows: list of dict rows with keys genome_id, fasta_path, abundance.
    per_genome_fastas: {genome_id: fasta_text} for genomes with >=1 amplicon.
    """
    # basename(fasta_path) (gz-stripped) -> (genome_id, abundance)
    by_source = {}
    for row in genomes_csv_rows:
        base = _strip_gz(os.path.basename(row["fasta_path"].strip()))
        by_source[base] = (row["genome_id"].strip(), row["abundance"].strip())

    seqs_by_genome = {}  # genome_id -> [seq, ...] (in file order)
    for header, seq in _iter_fasta(amplicons_lines):
        source = _source_key(header)
        if source is None or source not in by_source:
            continue  # off-target/unmapped source: skip rather than guess
        genome_id = by_source[source][0]
        seqs_by_genome.setdefault(genome_id, []).append(seq)

    if not seqs_by_genome:
        raise SystemExit(
            "build_amplicon_genomes_csv: no amplicons mapped to any input genome "
            "(empty amplicons.fa, or source= headers don't match the genomes CSV)"
        )

    fastas = {}
    for genome_id, seqs in seqs_by_genome.items():
        fastas[genome_id] = "".join(
            f">{genome_id}_amplicon{i}\n{seq}\n" for i, seq in enumerate(seqs)
        )

    out_header = ["genome_id", "fasta_path", "abundance"]
    out_rows = []
    for row in genomes_csv_rows:  # preserve input ordering
        genome_id = row["genome_id"].strip()
        if genome_id not in fastas:
            continue
        abundance = row["abundance"].strip()
        out_rows.append([genome_id, f"amplicons/{genome_id}.fa", abundance])
    return out_header, out_rows, fastas


def main(argv):
    amplicons_fa, genomes_csv, out_csv, out_dir = argv[1:5]
    with open(amplicons_fa) as fh:
        amplicons_lines = fh.readlines()
    with open(genomes_csv, newline="") as fh:
        genomes_rows = list(csv.DictReader(fh))

    out_header, out_rows, fastas = build(amplicons_lines, genomes_rows)

    os.makedirs(out_dir, exist_ok=True)
    for genome_id, text in fastas.items():
        with open(os.path.join(out_dir, f"{genome_id}.fa"), "w") as fh:
            fh.write(text)
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(out_header)
        w.writerows(out_rows)


def _selfcheck():
    amplicons = [
        ">ENA|X| bla source=genomeA.fasta.coordinates=10-20.orientation=FR.Tm=60\n",
        "ACGT\n",
        ">ENA|Y| bla source=genomeA.fasta.coordinates=30-40.orientation=RF.Tm=61\n",
        "TTGG\n",
        ">ENA|Z| bla source=genomeB.fasta.coordinates=5-9.orientation=FR.Tm=59\n",
        "CCAA\n",
    ]
    rows = [
        {"genome_id": "genomeA", "fasta_path": "d/genomeA.fasta.gz", "abundance": "0.7"},
        {"genome_id": "genomeB", "fasta_path": "d/genomeB.fasta", "abundance": "0.3"},
        {"genome_id": "genomeC", "fasta_path": "d/genomeC.fasta", "abundance": "0.1"},
    ]
    header, out_rows, fastas = build(amplicons, rows)
    assert header == ["genome_id", "fasta_path", "abundance"]
    # genomeC produced no amplicons -> dropped; A before B (input order preserved).
    assert out_rows == [
        ["genomeA", "amplicons/genomeA.fa", "0.7"],
        ["genomeB", "amplicons/genomeB.fa", "0.3"],
    ], out_rows
    # genomeA got both its amplicons with unique, traceable headers.
    assert fastas["genomeA"] == ">genomeA_amplicon0\nACGT\n>genomeA_amplicon1\nTTGG\n"
    assert fastas["genomeB"] == ">genomeB_amplicon0\nCCAA\n"
    assert "genomeC" not in fastas
    print("selfcheck ok")


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--selfcheck":
        _selfcheck()
    else:
        main(sys.argv)
