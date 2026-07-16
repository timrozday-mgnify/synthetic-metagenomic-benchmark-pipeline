#!/usr/bin/env python3
"""Sort/index a genome-blender ground-truth BAM and derive the ground-truth
profile: target (from the input abundances) and realized (from the BAM) relative
abundance per genome_id.

BAM references are named ``{genome_id}:{contig_id}`` (genome-blender convention),
so realized counts are the idxstats mapped reads aggregated by the part before
the first ':'.
"""

import argparse
import csv
from collections import defaultdict

import pysam


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam", required=True, help="Unsorted ground-truth BAM.")
    p.add_argument("--genomes-csv", required=True, help="genome_id,fasta_path,abundance CSV.")
    p.add_argument("--sorted-bam", required=True, help="Output coordinate-sorted BAM.")
    p.add_argument("--truth-tsv", required=True, help="Output ground-truth profile TSV.")
    p.add_argument(
        "--keep-names",
        help="Optional file of read names (one per line); only these reads are "
        "kept in the sorted BAM and counted for the realized profile.",
    )
    return p.parse_args()


def filter_bam_by_names(bam: str, names_file: str, out_bam: str) -> str:
    """Write only reads whose query_name is listed in names_file; return out_bam."""
    keep = {ln.strip() for ln in open(names_file) if ln.strip()}
    with pysam.AlignmentFile(bam, "rb", check_sq=False) as inb:
        with pysam.AlignmentFile(out_bam, "wb", template=inb) as out:
            for read in inb:
                if read.query_name in keep:
                    out.write(read)
    return out_bam


def main() -> int:
    args = parse_args()

    src_bam = args.bam
    if args.keep_names:
        src_bam = filter_bam_by_names(args.bam, args.keep_names, args.sorted_bam + ".filtered.bam")

    pysam.sort("-o", args.sorted_bam, src_bam)
    pysam.index(args.sorted_bam)

    # Realized read counts per genome_id, from idxstats (ref\tlen\tmapped\tunmapped).
    realized: dict[str, int] = defaultdict(int)
    for line in pysam.idxstats(args.sorted_bam).splitlines():
        ref, _length, mapped, _unmapped = line.split("\t")
        if ref == "*":
            continue
        genome_id = ref.split(":", 1)[0]
        realized[genome_id] += int(mapped)

    # Target abundances from the input CSV (normalised to sum to 1).
    target_raw: dict[str, float] = {}
    with open(args.genomes_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            target_raw[row["genome_id"]] = target_raw.get(row["genome_id"], 0.0) + float(
                row["abundance"]
            )
    target_total = sum(target_raw.values()) or 1.0
    realized_total = sum(realized.values()) or 1

    genome_ids = sorted(set(target_raw) | set(realized))
    with open(args.truth_tsv, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            ["genome_id", "target_rel_abundance", "realized_n_reads", "realized_rel_abundance"]
        )
        for gid in genome_ids:
            n = realized.get(gid, 0)
            writer.writerow(
                [
                    gid,
                    f"{target_raw.get(gid, 0.0) / target_total:.6f}",
                    n,
                    f"{n / realized_total:.6f}",
                ]
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
