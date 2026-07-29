#!/usr/bin/env python3
"""Normalise a superresolution `inferred_composition.csv` into the same
genome_id-keyed predicted profile that `normalize_sylph_profile.py` emits, so a
superresolution run is directly comparable to sylph/aap against `truth.tsv`.

Both superresolution pipelines emit one row per genome with `genome_id`,
`observed_rel_abundance` (naive per-reference counts) and `inferred_mean` (the
posterior mean of the Bayesian inversion), plus `inferred_lo`/`inferred_hi`.
We take `inferred_mean` — the pipelines' actual estimate — and renormalise it to
sum to 1 over the reported genomes.
"""

import argparse
import csv


def normalize(rows: list[dict[str, str]]) -> list[tuple[str, float]]:
    """[(genome_id, rel_abundance)] from inferred_mean, renormalised to sum to 1."""
    abund: dict[str, float] = {}
    for row in rows:
        gid = (row.get("genome_id") or "").strip()
        if not gid:
            continue
        value = row.get("inferred_mean")
        if value is None or str(value).strip() == "":
            raise SystemExit(f"normalize_sr_profile: no inferred_mean for genome '{gid}'")
        abund[gid] = abund.get(gid, 0.0) + float(value)
    total = sum(abund.values()) or 1.0
    return [(gid, abund[gid] / total) for gid in sorted(abund)]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--composition", required=True, help="<id>.inferred_composition.csv")
    p.add_argument("--output", required=True, help="Output normalised profile TSV.")
    args = p.parse_args()

    with open(args.composition, newline="") as fh:
        rows = list(csv.DictReader(fh))

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["genome_id", "predicted_rel_abundance", "predicted_tax_rel_abundance"])
        for gid, rel in normalize(rows):
            # ponytail: superresolution has no separate taxonomic-abundance estimate,
            # so both columns carry the same value (the truth-comparable contract needs
            # the column). The naive `observed_rel_abundance` baseline stays available
            # in the raw CSV published under profiling/sr/.
            writer.writerow([gid, f"{rel:.6f}", f"{rel:.6f}"])
    return 0


def _selfcheck() -> None:
    rows = [
        {"sample": "S1", "genome_id": "genomeA", "observed_rel_abundance": "0.5",
         "inferred_mean": "0.6", "inferred_lo": "0.5", "inferred_hi": "0.7"},
        {"sample": "S1", "genome_id": "genomeB", "observed_rel_abundance": "0.5",
         "inferred_mean": "0.2", "inferred_lo": "0.1", "inferred_hi": "0.3"},
    ]
    out = dict(normalize(rows))
    assert set(out) == {"genomeA", "genomeB"}, out
    assert abs(sum(out.values()) - 1.0) < 1e-9, out
    assert abs(out["genomeA"] - 0.75) < 1e-9, out  # 0.6 / (0.6 + 0.2)
    # Already-normalised input is a no-op, and rows are sorted by genome_id.
    same = normalize([{"genome_id": "z", "inferred_mean": "0.25"},
                      {"genome_id": "a", "inferred_mean": "0.75"}])
    assert [g for g, _ in same] == ["a", "z"], same
    assert abs(same[0][1] - 0.75) < 1e-9, same
    print("normalize_sr_profile self-check ok")


if __name__ == "__main__":
    import sys

    if "--selfcheck" in sys.argv:
        _selfcheck()
    else:
        raise SystemExit(main())
