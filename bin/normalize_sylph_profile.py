#!/usr/bin/env python3
"""Normalise a `sylph profile` TSV into a genome_id-keyed predicted profile that
lines up with the ground-truth `truth.tsv`.

sylph reports one row per detected genome with a `Genome_file` (the FASTA it was
built from) plus `Sequence_abundance` / `Taxonomic_abundance` (percentages). We
map `Genome_file` -> genome_id and renormalise the abundances to sum to 1 over
detected genomes, so `predicted_rel_abundance` is directly comparable to
`realized_rel_abundance` in truth.tsv.

genome_id mapping: with --genomes-csv (genome_id,fasta_path,...) we map by FASTA
basename (the `database = self` case, where the DB is the sample's own genomes).
Without it we fall back to the FASTA basename minus its extension.
"""

import argparse
import csv
import os


def _basename_stem(path: str) -> str:
    base = os.path.basename(path.strip())
    for ext in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"):
        if base.endswith(ext):
            return base[: -len(ext)]
    return os.path.splitext(base)[0]


def load_mapping(genomes_csv: str | None) -> dict[str, str]:
    """basename(fasta_path) -> genome_id from the genomes CSV (empty if none)."""
    if not genomes_csv:
        return {}
    mapping: dict[str, str] = {}
    with open(genomes_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            mapping[os.path.basename(row["fasta_path"].strip())] = row["genome_id"].strip()
    return mapping


def normalize(rows: list[dict[str, str]], mapping: dict[str, str]) -> list[tuple[str, float, float]]:
    """Return [(genome_id, seq_rel, tax_rel)] renormalised to sum to 1 each."""
    seq: dict[str, float] = {}
    tax: dict[str, float] = {}
    for row in rows:
        gfile = row["Genome_file"].strip()
        gid = mapping.get(os.path.basename(gfile)) or _basename_stem(gfile)
        seq[gid] = seq.get(gid, 0.0) + float(row.get("Sequence_abundance") or 0.0)
        tax[gid] = tax.get(gid, 0.0) + float(row.get("Taxonomic_abundance") or 0.0)
    seq_total = sum(seq.values()) or 1.0
    tax_total = sum(tax.values()) or 1.0
    return [(gid, seq[gid] / seq_total, tax[gid] / tax_total) for gid in sorted(seq)]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sylph-tsv", required=True, help="Output of `sylph profile`.")
    p.add_argument("--genomes-csv", help="genome_id,fasta_path[,abundance] mapping (self DB).")
    p.add_argument("--output", required=True, help="Output normalised profile TSV.")
    args = p.parse_args()

    with open(args.sylph_tsv, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    result = normalize(rows, load_mapping(args.genomes_csv))
    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["genome_id", "predicted_rel_abundance", "predicted_tax_rel_abundance"])
        for gid, seq_rel, tax_rel in result:
            writer.writerow([gid, f"{seq_rel:.6f}", f"{tax_rel:.6f}"])
    return 0


def _selfcheck() -> None:
    rows = [
        {"Genome_file": "/db/g1.fasta", "Sequence_abundance": "60", "Taxonomic_abundance": "70"},
        {"Genome_file": "/db/g2.fa", "Sequence_abundance": "20", "Taxonomic_abundance": "30"},
    ]
    mapping = {"g1.fasta": "genomeA", "g2.fa": "genomeB"}
    out = dict((gid, seq) for gid, seq, _ in normalize(rows, mapping))
    assert set(out) == {"genomeA", "genomeB"}, out
    assert abs(sum(out.values()) - 1.0) < 1e-9, out
    assert abs(out["genomeA"] - 0.75) < 1e-9, out  # 60 / (60+20)
    # No mapping -> basename stem is the genome_id.
    nostem = dict((gid, seq) for gid, seq, _ in normalize(rows, {}))
    assert set(nostem) == {"g1", "g2"}, nostem
    print("normalize_sylph_profile self-check ok")


if __name__ == "__main__":
    import sys

    if "--selfcheck" in sys.argv:
        _selfcheck()
    else:
        raise SystemExit(main())
