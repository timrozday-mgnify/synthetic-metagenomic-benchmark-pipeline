#!/usr/bin/env python3
"""Build the mapseq DB reference FASTA + taxonomy from a named sequence collection.

Ports the 16S staging + `.tax` synthesis from the standalone `build_profiling_dbs.py`,
but takes an *explicit* per-sequence taxonomy (Kingdom;Genus;Species) from a manifest
instead of deriving it from the id.

Manifest TSV, one row per collection sequence:
    <id> \\t <ssu_fasta> \\t <taxonomy>
Each ssu_fasta may hold >=1 16S records; every record is emitted with a
genome-traceable header `>{id}|{i}|{orig_header}` and tagged with that id's taxonomy.

Outputs:
    --out-fasta    concatenated, header-rewritten 16S records
    --out-tax      mapseq .tax (cutoff/name/levels header + `header\\ttaxonomy` lines)
    --out-headers  JSON sidecar for the OTU step: {fasta_order, header_to_id, id_to_tax}
"""

import argparse
import json
import sys

TAX_CUTOFFS = "0.00:0.08 0.85:0.65 0.95:0.85"  # loose defaults; only 3 ranks used here
TAX_LEVELS = ["Kingdom", "Genus", "Species"]


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


def read_manifest(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) != 3:
                raise ValueError(f"manifest row must be 'id<TAB>ssu<TAB>taxonomy': {line!r}")
            rows.append((parts[0], parts[1], parts[2]))
    return rows


def build(manifest, out_fasta, out_tax, out_headers):
    rows = read_manifest(manifest)
    fasta_order, header_to_id, id_to_tax = [], {}, {}
    with open(out_fasta, "w") as fa, open(out_tax, "w") as tx:
        tx.write(f"#cutoff: {TAX_CUTOFFS}\n#name: custom\n#levels: {' '.join(TAX_LEVELS)}\n")
        for gid, ssu, taxonomy in rows:
            if gid in id_to_tax and id_to_tax[gid] != taxonomy:
                raise ValueError(f"conflicting taxonomy for id {gid}: {id_to_tax[gid]!r} vs {taxonomy!r}")
            id_to_tax[gid] = taxonomy
            for i, (orig, seq) in enumerate(parse_fasta(ssu)):
                header = f"{gid}|{i}|{orig}"
                fa.write(f">{header}\n{seq}\n")
                tx.write(f"{header}\t{taxonomy}\n")
                fasta_order.append(header)
                header_to_id[header] = gid
    if not fasta_order:
        raise SystemExit("build_mapseq_refs: no 16S records found across the collection's sequences")
    with open(out_headers, "w") as fh:
        json.dump({"fasta_order": fasta_order, "header_to_id": header_to_id, "id_to_tax": id_to_tax}, fh)


def _selfcheck():
    import os
    import tempfile

    d = tempfile.mkdtemp()
    ssu = os.path.join(d, "a.fasta")
    with open(ssu, "w") as fh:
        fh.write(">seqX\nACGT\n>seqY\nTTTT\n")
    man = os.path.join(d, "m.tsv")
    with open(man, "w") as fh:
        fh.write(f"bacteroides_fragilis\t{ssu}\tBacteria;Bacteroides;fragilis\n")
    fa, tx, hd = (os.path.join(d, n) for n in ("o.fasta", "o.tax", "o.json"))
    build(man, fa, tx, hd)

    faout = open(fa).read()
    assert ">bacteroides_fragilis|0|seqX" in faout and ">bacteroides_fragilis|1|seqY" in faout, faout
    txout = open(tx).read().splitlines()
    assert txout[0].startswith("#cutoff:") and txout[2].startswith("#levels:"), txout
    assert txout[3] == "bacteroides_fragilis|0|seqX\tBacteria;Bacteroides;fragilis", txout
    meta = json.load(open(hd))
    assert meta["header_to_id"]["bacteroides_fragilis|1|seqY"] == "bacteroides_fragilis"
    assert meta["id_to_tax"]["bacteroides_fragilis"] == "Bacteria;Bacteroides;fragilis"
    assert meta["fasta_order"] == ["bacteroides_fragilis|0|seqX", "bacteroides_fragilis|1|seqY"]
    print("build_mapseq_refs self-check ok")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", required=True)
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-tax", required=True)
    p.add_argument("--out-headers", required=True)
    a = p.parse_args()
    build(a.manifest, a.out_fasta, a.out_tax, a.out_headers)


if __name__ == "__main__":
    if "--selfcheck" in sys.argv:
        _selfcheck()
    else:
        main()
