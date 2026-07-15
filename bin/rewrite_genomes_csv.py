#!/usr/bin/env python3
"""Rewrite a genome-blender input CSV so fasta_path points at the file's
basename (the FASTA is staged into the work dir by Nextflow).

Usage: rewrite_genomes_csv.py in.csv out.csv
Columns: genome_id,fasta_path,abundance
"""

import csv
import os
import sys


def main() -> int:
    src, dst = sys.argv[1], sys.argv[2]
    with open(src, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        sys.exit(f"No genome rows in {src}")
    fields = ["genome_id", "fasta_path", "abundance"]
    with open(dst, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "genome_id": row["genome_id"],
                    "fasta_path": os.path.basename(row["fasta_path"]),
                    "abundance": row["abundance"],
                }
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
