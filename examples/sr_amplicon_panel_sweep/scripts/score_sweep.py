#!/usr/bin/env python3
"""Score every profile of the sweep against its truth: one CSV row per benchmark dir x method.

    python scripts/score_sweep.py [results_dir] [config.yaml] > scores.csv

  tv                  genome total variation over the panel genomes, both sides
                      renormalised without `background`
  max_abs_error       the largest per-genome error, and `worst_genome` it is on
  mass_absent         predicted share on genomes the community does not contain
  background          reinterpretation's out-of-panel share (0 for the custom arm)
  strain_split_error  |predicted - true| share of `score.strain_pair`'s second genome
                      within the pair; blank when the pair holds <= 1% of the truth

`custom.mapseq_only` is the custom arm's raw per-genome read share: MAPseq against the
panel, no inference. It needs no setting of its own.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import panel_sweep as ps  # noqa: E402

HERE = Path(__file__).resolve().parent.parent
COLUMNS = ["sample", "depth", "method", "tv", "max_abs_error", "worst_genome",
           "mass_absent", "background", "strain_split_error"]


def read_shares(path, column):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t" if path.suffix == ".tsv" else ",")
        shares = {row["genome_id"]: float(row[column]) for row in reader}
    # The profile is over genomes; the unexplained share sits in its sidecar.
    sidecar = path.with_name(path.name.replace(".sr_profile.tsv", ".sr_background.tsv"))
    if sidecar != path and sidecar.exists():
        shares["background"] = float(sidecar.read_text().split()[1])
    return shares


def scores(pred, truth, pair):
    """Accuracy of one profile against the truth, both {genome_id: share}."""
    background = pred.get("background", 0.0)
    pred = {g: v for g, v in pred.items() if g != "background"}
    genomes = set(pred) | set(truth)
    p_total, t_total = sum(pred.values()) or 1.0, sum(truth.values()) or 1.0
    p = {g: pred.get(g, 0.0) / p_total for g in genomes}
    t = {g: truth.get(g, 0.0) / t_total for g in genomes}
    err = {g: abs(p[g] - t[g]) for g in genomes}
    worst = max(err, key=lambda g: err[g])
    both_p = p.get(pair[0], 0.0) + p.get(pair[1], 0.0)
    both_t = t.get(pair[0], 0.0) + t.get(pair[1], 0.0)
    split = (abs(p.get(pair[1], 0.0) / both_p - t.get(pair[1], 0.0) / both_t)
             if both_t > 0.01 and both_p > 0 else "")
    return {"tv": 0.5 * sum(err.values()), "max_abs_error": err[worst], "worst_genome": worst,
            "mass_absent": sum(p[g] for g in genomes if t[g] == 0),
            "background": background, "strain_split_error": split}


def main():
    results_dir = (Path(sys.argv[1]).resolve() if len(sys.argv) > 1
                   else HERE.parent.parent / "results" / "sr_amplicon_panel_sweep")
    cfg = ps.load_config(sys.argv[2] if len(sys.argv) > 2 else HERE / "config.yaml")
    pair = cfg["score"]["strain_pair"]
    names = [s["name"] for arm in ps.ARMS for s in ps.settings(cfg, arm)]
    first_custom = ps.settings(cfg, "custom")[0]["name"]

    writer = csv.DictWriter(sys.stdout, COLUMNS, extrasaction="ignore")
    writer.writeheader()
    skipped = 0
    for sample, run_id, directory, depth, _pair in ps.benchmark_dirs(cfg, results_dir):
        truths = sorted(directory.glob("*.truth.tsv"))
        if not truths:
            skipped += 1
            continue
        truth = read_shares(truths[0], "realized_rel_abundance")
        profiles = {name: (directory / f"{run_id}.{name}.sr_profile.tsv",
                           "predicted_rel_abundance") for name in names}
        # Observed shares do not depend on the inference knobs, so any custom point has them.
        profiles["custom.mapseq_only"] = (
            directory / "profiling" / "sr" / f"{run_id}.{first_custom}.inferred_composition.csv",
            "observed_rel_abundance")
        for method, (path, column) in profiles.items():
            if path.exists():
                row = {"sample": sample, "depth": depth or "full", "method": method,
                       **scores(read_shares(path, column), truth, pair)}
                writer.writerow({k: f"{v:.5f}" if isinstance(v, float) else v
                                 for k, v in row.items()})
    if skipped:
        print(f"NOTE: {skipped} benchmark dir(s) under {results_dir} have no truth.tsv yet",
              file=sys.stderr)


def _selfcheck():
    pair = ["bu", "bu2"]
    truth = {"a": 0.5, "bu": 0.0, "bu2": 0.5, "gone": 0.0}
    perfect = scores({"a": 0.25, "bu2": 0.25, "background": 0.5}, truth, pair)
    assert perfect["tv"] == 0 and perfect["background"] == 0.5, perfect
    assert perfect["strain_split_error"] == 0, perfect
    off = scores({"a": 0.4, "bu": 0.2, "bu2": 0.2, "gone": 0.2}, truth, pair)
    assert abs(off["tv"] - 0.4) < 1e-9 and abs(off["mass_absent"] - 0.4) < 1e-9, off
    assert abs(off["strain_split_error"] - 0.5) < 1e-9 and off["worst_genome"] == "bu2", off
    assert scores({"a": 1.0}, {"a": 1.0}, pair)["strain_split_error"] == "", "pair absent"
    print("score_sweep selfcheck OK")


if __name__ == "__main__":
    if "--selfcheck" in sys.argv:
        _selfcheck()
    else:
        main()
