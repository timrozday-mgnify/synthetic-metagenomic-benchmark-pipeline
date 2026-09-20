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

When every row also has a `taxonomy` lineage, a MAPseq `.tax` sidecar is written
beside the FASTA in the form superresolution-amplicon's `build_mapseq_database.py`
writes for a genome collection: lineages padded with `unclassified` to the deepest
one, the genome id as a last `Genome` level, keyed by each record's first header
token (the id the nested pipeline reads). It feeds the nested `--taxonomy`.

A row with a `taxon` (and no `fasta_path`) is a taxon panel entry: it is written to
`--taxa-output` as `id<TAB>taxon`, superresolution-amplicon's `--panel_taxa`, and takes no
part in the FASTA or `.tax`. A taxa-only collection writes no FASTA.
"""

import argparse
import csv
import gzip
import os

# ponytail: SILVA (and any other rRNA reference) ships RNA, so primer matching downstream
# finds nothing unless U is folded to T. Everything else passes through untouched.
_RNA_TO_DNA = str.maketrans("Uu", "Tt")


def _open_text(path: str):
    """Open a FASTA that may or may not be gzipped."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_rows(genomes_csv: str) -> list[tuple[str, str, str, str]]:
    """[(genome_id, staged_fasta_basename, taxonomy or '', taxon or '')] from the genomes CSV."""
    with open(genomes_csv, newline="") as fh:
        return [
            (row["genome_id"].strip(), os.path.basename((row.get("fasta_path") or "").strip()),
             (row.get("taxonomy") or "").strip(), (row.get("taxon") or "").strip())
            for row in csv.DictReader(fh)
            if row.get("genome_id")
        ]


def write_tax(rows: list[tuple[str, str, str, str]], headers: list[tuple[str, str]], tax_path: str) -> None:
    """MAPseq `.tax` over `headers` [(genome_id, record id)], lineages from `rows`."""
    lineage = {g: [f.strip() for f in t.split(";")] for g, _, t, _ in rows}
    if any("" in fields for fields in lineage.values()):
        raise SystemExit(f"build_sr_refs: empty rank in a taxonomy lineage: {lineage}")
    depth = max(len(fields) for fields in lineage.values())
    with open(tax_path, "w") as out:
        out.write("#cutoff: " + " ".join(["0.00:0.00"] * (depth + 1)) + "\n")
        out.write("#name: refdb\n")
        out.write("#levels: " + " ".join([f"Taxonomy_{i}" for i in range(1, depth + 1)] + ["Genome"]) + "\n")
        for genome_id, record in headers:
            padded = lineage[genome_id] + ["unclassified"] * (depth - len(lineage[genome_id]))
            out.write(f"{record}\t{';'.join([*padded, genome_id])}\n")


def build(rows: list[tuple[str, str, str, str]], out_path: str, resolve=lambda p: p,
          tax_path: str | None = None, taxa_path: str | None = None) -> int:
    """Write the combined FASTA over genome rows (and `.tax` when each has a lineage), and
    taxon rows to `taxa_path`; return FASTA records."""
    taxa = [(g, taxon) for g, fasta, _, taxon in rows if taxon]
    if any(fasta for g, fasta, _, taxon in rows if taxon):
        raise SystemExit("build_sr_refs: a row has both a taxon and a fasta_path")
    rows = [r for r in rows if not r[3]]
    if taxa:
        if not taxa_path:
            raise SystemExit("build_sr_refs: taxon rows need --taxa-output")
        with open(taxa_path, "w") as out:
            out.write("id\ttaxon\n")
            out.writelines(f"{g}\t{taxon}\n" for g, taxon in taxa)
        if not rows:
            return 0
    written = 0
    headers = []
    with open(out_path, "w") as out:
        for genome_id, fasta, _, _ in rows:
            path = resolve(fasta)
            if not os.path.exists(path):
                raise SystemExit(f"build_sr_refs: {path} (genome '{genome_id}') not found")
            index = 0
            with _open_text(path) as fh:
                for line in fh:
                    if line.startswith(">"):
                        record = f"{genome_id}|{index}|{line[1:].strip()}"
                        out.write(f">{record}\n")
                        headers.append((genome_id, record.split()[0]))
                        index += 1
                        written += 1
                    else:
                        line = line.translate(_RNA_TO_DNA)
                        out.write(line if line.endswith("\n") else line + "\n")
            if index == 0:
                raise SystemExit(f"build_sr_refs: {path} (genome '{genome_id}') has no records")
    if not written:
        raise SystemExit(f"build_sr_refs: no records written from {len(rows)} genome(s)")
    if tax_path and all(t for _, _, t, _ in rows):
        write_tax(rows, headers, tax_path)
    return written


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genomes-csv", required=True, help="genome_id,fasta_path[,abundance] CSV.")
    p.add_argument("--output", required=True, help="Combined reference FASTA to write.")
    p.add_argument("--tax-output", help="MAPseq .tax to write when every row has a taxonomy.")
    p.add_argument("--taxa-output", help="id<TAB>taxon TSV of the rows with a taxon.")
    args = p.parse_args()

    n = build(read_rows(args.genomes_csv), args.output, tax_path=args.tax_output,
              taxa_path=args.taxa_output)
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

        # RNA references (SILVA) are folded to DNA so primers can match.
        rna = os.path.join(d, "rna.fasta")
        with open(rna, "w") as fh:
            fh.write(">r1\nACGUuGU\n")
        rna_csv = os.path.join(d, "rna.csv")
        with open(rna_csv, "w") as fh:
            fh.write("genome_id,fasta_path\nrna,rna.fasta\n")
        rna_out = os.path.join(d, "rna_refs.fasta")
        build(read_rows(rna_csv), rna_out, resolve=lambda f: os.path.join(d, f))
        assert [ln.strip() for ln in open(rna_out) if not ln.startswith(">")] == ["ACGTtGT"]
        assert not os.path.exists(os.path.join(d, "none.tax"))
        build(read_rows(csv_path), out, resolve=lambda f: os.path.join(d, f),
              tax_path=os.path.join(d, "none.tax"))
        assert not os.path.exists(os.path.join(d, "none.tax")), "no taxonomy column, no .tax"

        with open(csv_path, "w") as fh:
            fh.write("genome_id,fasta_path,taxonomy\n")
            fh.write(f'genomeA,{os.path.basename(a)},"Bacteria; Bacillota;Bacilli"\n')
            fh.write(f"genomeB,{os.path.basename(b)},Bacteria;Bacteroidota\n")
        tax = os.path.join(d, "refs.tax")
        build(read_rows(csv_path), out, resolve=lambda f: os.path.join(d, f), tax_path=tax)
        lines = open(tax).read().splitlines()
        assert lines[2] == "#levels: Taxonomy_1 Taxonomy_2 Taxonomy_3 Genome", lines
        assert lines[3:] == [
            "genomeA|0|ctg1\tBacteria;Bacillota;Bacilli;genomeA",
            "genomeA|1|ctg2\tBacteria;Bacillota;Bacilli;genomeA",
            "genomeB|0|ctg1\tBacteria;Bacteroidota;unclassified;genomeB",
        ], lines

        # Taxon rows: to the taxa TSV, not the FASTA; taxa-only writes no FASTA.
        with open(csv_path, "w") as fh:
            fh.write("genome_id,fasta_path,taxon\n")
            fh.write(f"genomeA,{os.path.basename(a)},\n")
            fh.write('bfrag,,"Bacteria;Bacteroidota;Bacteroides;Bacteroides fragilis"\n')
        taxa = os.path.join(d, "panel_taxa.tsv")
        assert build(read_rows(csv_path), out, resolve=lambda f: os.path.join(d, f),
                     tax_path=os.path.join(d, "mixed.tax"), taxa_path=taxa) == 2
        assert [ln[1:].split("|")[0] for ln in open(out) if ln.startswith(">")] == ["genomeA"] * 2
        assert open(taxa).read() == ("id\ttaxon\n"
                                     "bfrag\tBacteria;Bacteroidota;Bacteroides;Bacteroides fragilis\n")
        assert not os.path.exists(os.path.join(d, "mixed.tax")), "genome rows have no taxonomy"
        with open(csv_path, "w") as fh:
            fh.write("genome_id,fasta_path,taxon\nbfrag,,Bacteroides fragilis\n")
        only = os.path.join(d, "taxa_only.fasta")
        assert build(read_rows(csv_path), only, taxa_path=taxa) == 0
        assert not os.path.exists(only), "taxa-only writes no FASTA"
    print("build_sr_refs self-check ok")


if __name__ == "__main__":
    import sys

    if "--selfcheck" in sys.argv:
        _selfcheck()
    else:
        raise SystemExit(main())
