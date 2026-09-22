#!/usr/bin/env python3
"""Score the Phase 4 arms against the truth, one CSV row per arm.

    python score.py <results_dir> > scores.csv

  tv          species total variation over the panel, both sides renormalised over the
              panel's species (every genome here is its own species, so for the SR arms
              this is also the genome TV)
  background  predicted share outside the panel: SR's out-of-panel share, or AAP reads
              whose label is not one of the panel's SILVA species
  mass_absent predicted share on panel species the sample does not contain

AAP's labels are SILVA 138.1 names. species.tsv maps each genome to its 138.1 species;
a blank means 138.1 has no species label for it, so arm 1 cannot reach that genome.
"""
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SAMPLE = "T1.515-YF-806BR"   # the run id carries the primer pair
SR_ARMS = ["sim_merged", "sim_pairs", "align", "raw"]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def score(pred, background, truth):
    """TV etc. of {species: share} against the truth, both renormalised over the panel."""
    p_total, t_total = sum(pred.values()) or 1.0, sum(truth.values()) or 1.0
    species = set(truth)
    p = {s: pred.get(s, 0.0) / p_total for s in species}
    t = {s: truth[s] / t_total for s in species}
    return {"tv": 0.5 * sum(abs(p[s] - t[s]) for s in species),
            "background": background,
            "mass_absent": sum(p[s] for s in species if t[s] == 0)}


def main():
    results = Path(sys.argv[1])
    sample = results / SAMPLE
    silva = {r["genome_id"]: r["silva_species"] for r in read_tsv(HERE / "species.tsv")}

    truth_rows = read_tsv(next(sample.glob("*.truth.tsv")))
    realized = {r["genome_id"]: float(r["realized_rel_abundance"]) for r in truth_rows}
    truth = {g: realized.get(g, 0.0) for g in silva}   # absent panel genomes count as 0

    writer = csv.DictWriter(sys.stdout, ["arm", "tv", "background", "mass_absent"])
    writer.writeheader()

    # Arm 1: AAP's own counts per SILVA label. Genome ids stand in for species below.
    by_label = {s: g for g, s in silva.items() if s}
    counts, other = {}, 0.0
    for tsv in (sample / "profiling" / "aap" / SAMPLE / "taxonomy-summary").glob("*/*.tsv"):
        for line in tsv.read_text().splitlines():
            if line.startswith("#"):
                continue
            _, n, lineage, *_ = line.split("\t")
            g = by_label.get(lineage.rsplit(";s__", 1)[-1]) if ";s__" in lineage else None
            if g:
                counts[g] = counts.get(g, 0.0) + float(n)
            else:
                other += float(n)
    total = sum(counts.values()) + other
    writer.writerow({"arm": "aap", **score(counts, other / total if total else 0.0, truth)})

    for arm in SR_ARMS:
        profile = sample / f"{SAMPLE}.{arm}.sr_profile.tsv"
        if not profile.exists():
            print(f"missing {profile}", file=sys.stderr)
            continue
        pred = {r["genome_id"]: float(r["predicted_rel_abundance"]) for r in read_tsv(profile)}
        sidecar = profile.with_name(f"{SAMPLE}.{arm}.sr_background.tsv")
        background = float(sidecar.read_text().split()[1]) if sidecar.exists() else 0.0
        writer.writerow({"arm": arm, **score(pred, background, truth)})


if __name__ == "__main__":
    main()
