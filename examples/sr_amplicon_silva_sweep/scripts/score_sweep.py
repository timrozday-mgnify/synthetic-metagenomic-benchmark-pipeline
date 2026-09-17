#!/usr/bin/env python3
"""Score every grid point per genus: one CSV row per benchmark dir x method.

    python scripts/score_sweep.py [results_dir] [config.yaml] > scores.csv

Truth is `truth.tsv` rolled up to genus through `panel[].taxonomy`. A prediction is the
run's `inferred_v4_groups.csv`: each V4 group's mass goes to the genus of its `lca`, and a
group whose `lca` stops above genus (members from several genera) goes to `unresolved`.

  tv             total variation over genera plus `unresolved` (truth has none there)
  max_abs_error  the largest per-genus error, and `worst_genus` it is on
  mass_absent    predicted share on genera the community does not contain
  unresolved     predicted share SILVA cannot place at genus

`mapseq_only` is the first grid point's `observed_rel_abundance`: MAPseq's labels with no
inference, which no inference knob changes.
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import silva_sweep as ss  # noqa: E402

HERE = Path(__file__).resolve().parent.parent
COLUMNS = ["sample", "depth", "method", "tv", "max_abs_error", "worst_genus",
           "mass_absent", "unresolved"]
UNRESOLVED = "unresolved"


def truth_genera(path, lineage):
    """{genus: share} from a truth.tsv, genomes rolled up through {genome_id: lineage}."""
    out = defaultdict(float)
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[ss.genus(lineage[row["genome_id"]])] += float(row["realized_rel_abundance"])
    return dict(out)


def predicted_genera(path, column):
    """{genus | 'unresolved': share} from an inferred_v4_groups.csv."""
    out = defaultdict(float)
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            out[ss.genus(row["lca"]) or UNRESOLVED] += float(row[column])
    return dict(out)


def scores(pred, truth):
    """Accuracy of one genus profile against the truth, both {genus: share}."""
    genera = set(pred) | set(truth)
    p_total, t_total = sum(pred.values()) or 1.0, sum(truth.values()) or 1.0
    p = {g: pred.get(g, 0.0) / p_total for g in genera}
    t = {g: truth.get(g, 0.0) / t_total for g in genera}
    err = {g: abs(p[g] - t[g]) for g in genera}
    worst = max(err, key=lambda g: err[g])
    return {"tv": 0.5 * sum(err.values()), "max_abs_error": err[worst],
            "worst_genus": worst.rsplit(";", 1)[-1],
            "mass_absent": sum(p[g] for g in genera if t[g] == 0 and g != UNRESOLVED),
            "unresolved": p.get(UNRESOLVED, 0.0)}


def main():
    results_dir = (Path(sys.argv[1]).resolve() if len(sys.argv) > 1
                   else HERE.parent.parent / "results" / "sr_amplicon_silva_sweep")
    cfg = ss.load_config(sys.argv[2] if len(sys.argv) > 2 else HERE / "config.yaml")
    lineage = {m["id"]: m["taxonomy"] for m in cfg["panel"]}
    names = [s["name"] for s in ss.settings(cfg)]

    writer = csv.DictWriter(sys.stdout, COLUMNS, extrasaction="ignore")
    writer.writeheader()
    skipped = 0
    for sample, directory, depth, _pair in ss.benchmark_dirs(cfg, results_dir):
        truths = sorted(directory.glob("*.truth.tsv"))
        if not truths:
            skipped += 1
            continue
        truth = truth_genera(truths[0], lineage)
        run_id = f"{sample}.sub{depth}" if depth else sample
        groups = {name: directory / "profiling" / "sr" / f"{run_id}.{name}.inferred_v4_groups.csv"
                  for name in names}
        profiles = {name: (path, "inferred_mean") for name, path in groups.items()}
        profiles["mapseq_only"] = (groups[names[0]], "observed_rel_abundance")
        for method, (path, column) in profiles.items():
            if path.exists():
                row = {"sample": sample, "depth": depth or "full", "method": method,
                       **scores(predicted_genera(path, column), truth)}
                writer.writerow({k: f"{v:.5f}" if isinstance(v, float) else v
                                 for k, v in row.items()})
    if skipped:
        print(f"NOTE: {skipped} benchmark dir(s) under {results_dir} have no truth.tsv yet",
              file=sys.stderr)


def _selfcheck():
    import tempfile

    fam = "Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae"
    bac, strep = f"{fam};Bacteroides", "Bacteria;Bacillota;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus"
    with tempfile.TemporaryDirectory() as d:
        truth_tsv, groups_csv = Path(d) / "S.truth.tsv", Path(d) / "S.p.inferred_v4_groups.csv"
        truth_tsv.write_text("genome_id\trealized_rel_abundance\n"
                             "bf\t0.25\nbt\t0.25\nss\t0.5\n")
        groups_csv.write_text(
            "sample,v4_group_id,lca,observed_rel_abundance,inferred_mean\n"
            f"S,g1,{bac};unclassified,0.4,0.5\n"
            f"S,g2,{fam},0.2,0.1\n"                         # stops above genus
            f"S,g3,{strep};unclassified,0.3,0.3\n"
            "S,g4,Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;"
            "Enterobacteriaceae;Escherichia-Shigella,0.1,0.1\n")
        truth = truth_genera(truth_tsv, {"bf": bac, "bt": bac, "ss": strep})
        assert truth == {bac: 0.5, strep: 0.5}, truth
        pred = predicted_genera(groups_csv, "inferred_mean")
        assert abs(pred[UNRESOLVED] - 0.1) < 1e-9 and abs(pred[bac] - 0.5) < 1e-9, pred
        s = scores(pred, truth)
    # |0.5-0.5| + |0.3-0.5| + 0.1 unresolved + 0.1 absent = 0.4, halved.
    assert abs(s["tv"] - 0.2) < 1e-9, s
    assert s["worst_genus"] == "Streptococcus" and abs(s["max_abs_error"] - 0.2) < 1e-9, s
    assert abs(s["mass_absent"] - 0.1) < 1e-9 and abs(s["unresolved"] - 0.1) < 1e-9, s
    assert scores({bac: 2.0}, {bac: 1.0})["tv"] == 0, "both sides renormalised"
    print("score_sweep selfcheck OK")


if __name__ == "__main__":
    if "--selfcheck" in sys.argv:
        _selfcheck()
    else:
        main()
