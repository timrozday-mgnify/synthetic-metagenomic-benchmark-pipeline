#!/usr/bin/env python3
"""Majority-vote each mapseq cluster to one taxonomy string -> mapseq `.otu`.

Ports `parse_mscluster` + `build_mapseq_otu` from the standalone `build_profiling_dbs.py`.
Reads the `.mscluster` mapseq emits (one line per cluster:
`<cluster_id> <member_seq_index> <member_seq_index> ...`) and the `headers.json`
sidecar from `build_mapseq_refs.py`, and writes one row per cluster:

    <cluster_id> \\t <taxonomy> \\t <cluster_id>

The 3rd column is a stable placeholder taxid (real SILVA `.otu` files carry an NCBI
taxid there; a custom DB has none).
"""

import argparse
import json
import sys


def parse_mscluster(path):
    """Yield (cluster_id, [member seq_index, ...]) from a mapseq .mscluster file."""
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            cluster_i, *member_is = (int(x) for x in parts)
            yield cluster_i, member_is


def build_otu(mscluster, headers, out_otu):
    meta = json.load(open(headers))
    fasta_order = meta["fasta_order"]
    header_to_id = meta["header_to_id"]
    id_to_tax = meta["id_to_tax"]
    clusters = {c: members for c, members in parse_mscluster(mscluster)}
    with open(out_otu, "w") as fh:
        for cluster_i in sorted(clusters):
            ids = [header_to_id[fasta_order[i]] for i in clusters[cluster_i]]
            winner = max(set(ids), key=ids.count)
            fh.write(f"{cluster_i}\t{id_to_tax[winner]}\t{cluster_i}\n")


def _selfcheck():
    import os
    import tempfile

    d = tempfile.mkdtemp()
    ms = os.path.join(d, "c.mscluster")
    with open(ms, "w") as fh:
        fh.write("0 0 1 2\n1 3\n")
    hd = os.path.join(d, "h.json")
    with open(hd, "w") as fh:
        json.dump(
            {
                "fasta_order": ["h0", "h1", "h2", "h3"],
                "header_to_id": {"h0": "frag", "h1": "frag", "h2": "bolt", "h3": "ramo"},
                "id_to_tax": {
                    "frag": "Bacteria;Bacteroides;fragilis",
                    "bolt": "Bacteria;Clostridium;bolteae",
                    "ramo": "Bacteria;Clostridium;ramosum",
                },
            },
            fh,
        )
    out = os.path.join(d, "o.otu")
    build_otu(ms, hd, out)
    rows = [line.split("\t") for line in open(out).read().splitlines()]
    # cluster 0 = [frag, frag, bolt] -> majority frag; cluster 1 = [ramo].
    assert rows[0] == ["0", "Bacteria;Bacteroides;fragilis", "0"], rows
    assert rows[1] == ["1", "Bacteria;Clostridium;ramosum", "1"], rows
    print("build_mapseq_otu self-check ok")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mscluster", required=True)
    p.add_argument("--headers", required=True)
    p.add_argument("--out-otu", required=True)
    a = p.parse_args()
    build_otu(a.mscluster, a.headers, a.out_otu)


if __name__ == "__main__":
    if "--selfcheck" in sys.argv:
        _selfcheck()
    else:
        main()
