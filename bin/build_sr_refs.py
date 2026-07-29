#!/usr/bin/env python3
"""Build the combined reference FASTA the superresolution pipelines expect.

Both superresolution-amplicon and superresolution-shotgun take one FASTA over the
whole community whose headers encode genome membership as `{genome}|...` — the
text before the first '|' is the genome id, and entries sharing a genome id are
that genome's contigs (shotgun) or 16S copies (amplicon). Emitting
`{genome_id}|{n}|{original header}` satisfies both contracts at once.

Input is a genomes CSV (`genome_id,fasta_path[,...]`, the same file the rest of
the pipeline passes around). Each row's FASTA is resolved by basename in the
working directory, since Nextflow stages inputs there under their basenames.
"""

import argparse
import csv
import gzip
import os


def _open_text(path: str):
    """Open a FASTA that may or may not be gzipped."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_rows(genomes_csv: str) -> list[tuple[str, str]]:
    """[(genome_id, staged_fasta_basename)] from the genomes CSV."""
    with open(genomes_csv, newline="") as fh:
        return [
            (row["genome_id"].strip(), os.path.basename(row["fasta_path"].strip()))
            for row in csv.DictReader(fh)
            if row.get("genome_id")
        ]


def build(rows: list[tuple[str, str]], out_path: str, resolve=lambda p: p) -> int:
    """Write the combined FASTA; return the number of records emitted."""
    written = 0
    with open(out_path, "w") as out:
        for genome_id, fasta in rows:
            path = resolve(fasta)
            if not os.path.exists(path):
                raise SystemExit(f"build_sr_refs: {path} (genome '{genome_id}') not found")
            index = 0
            with _open_text(path) as fh:
                for line in fh:
                    if line.startswith(">"):
                        orig = line[1:].strip()
                        out.write(f">{genome_id}|{index}|{orig}\n")
                        index += 1
                        written += 1
                    else:
                        out.write(line if line.endswith("\n") else line + "\n")
            if index == 0:
                raise SystemExit(f"build_sr_refs: {path} (genome '{genome_id}') has no records")
    if not written:
        raise SystemExit(f"build_sr_refs: no records written from {len(rows)} genome(s)")
    return written


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genomes-csv", required=True, help="genome_id,fasta_path[,abundance] CSV.")
    p.add_argument("--output", required=True, help="Combined reference FASTA to write.")
    args = p.parse_args()

    n = build(read_rows(args.genomes_csv), args.output)
    print(f"build_sr_refs: wrote {n} records to {args.output}")
    return 0


def _selfcheck() -> None:
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        a = os.path.join(d, "a.fasta")
        b = os.path.join(d, "b.fasta.gz")
        with open(a, "w") as fh:
            fh.write(">ctg1 some desc\nACGT\n>ctg2\nTTTT\n")
        with gzip.open(b, "wt") as fh:
            fh.write(">ctg1\nGGGG\n")
        csv_path = os.path.join(d, "genomes.csv")
        with open(csv_path, "w") as fh:
            fh.write("genome_id,fasta_path,abundance\n")
            fh.write(f"genomeA,/elsewhere/{os.path.basename(a)},0.5\n")
            fh.write(f"genomeB,/elsewhere/{os.path.basename(b)},0.5\n")

        out = os.path.join(d, "refs.fasta")
        n = build(read_rows(csv_path), out, resolve=lambda f: os.path.join(d, f))
        assert n == 3, n
        heads = [ln.strip() for ln in open(out) if ln.startswith(">")]
        assert heads == [">genomeA|0|ctg1 some desc", ">genomeA|1|ctg2", ">genomeB|0|ctg1"], heads
        # The genome id is recoverable exactly as both SR pipelines recover it.
        assert [h[1:].split("|", 1)[0] for h in heads] == ["genomeA", "genomeA", "genomeB"]
        seqs = [ln.strip() for ln in open(out) if not ln.startswith(">")]
        assert seqs == ["ACGT", "TTTT", "GGGG"], seqs
    print("build_sr_refs self-check ok")


if __name__ == "__main__":
    import sys

    if "--selfcheck" in sys.argv:
        _selfcheck()
    else:
        raise SystemExit(main())
