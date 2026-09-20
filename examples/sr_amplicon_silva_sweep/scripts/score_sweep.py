#!/usr/bin/env python3
"""Score every profiling arm and grid point per genus: one CSV row per benchmark dir x
method.

    python scripts/score_sweep.py [results_dir] [config.yaml] > scores.csv

Genus is the one space all four arms can be compared in - two of them report panel
genomes, one reports SILVA V4 groups and one reports SILVA lineages - so everything is
rolled up to it. Truth is `truth.tsv` through `panel[].taxonomy`. Predictions:

  silva        `inferred_v4_groups.csv` - each V4 group's mass goes to the genus of its
               `lca`; a group whose `lca` stops above genus goes to `unresolved`.
  custom,      `inferred_composition.csv` - each panel genome's mass goes to its own
  silva_panel  lineage's genus. `silva_panel`'s `background` bucket (reads the panel
               cannot explain) goes to `unresolved`.
  aap          the amplicon-analysis-pipeline's krona table - a count per lineage,
               summed per genus; a lineage stopping above genus goes to `unresolved`.

  tv             total variation over genera plus `unresolved` (truth has none there)
  max_abs_error  the largest per-genus error, and `worst_genus` it is on
  mass_absent    predicted share on genera the community does not contain
  unresolved     predicted share that arm cannot place at genus

`mapseq_only` is the first `silva` grid point's `observed_rel_abundance`: MAPseq's labels
against SILVA with no inference, which no inference knob changes.

The `error_model` column names the arm the READS came from (trained / flat / naive), so a
method's row is only comparable with another method's at the same error model.
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import silva_sweep as ss  # noqa: E402

HERE = Path(__file__).resolve().parent.parent
COLUMNS = ["sample", "error_model", "depth", "arm", "method", "tv", "max_abs_error",
           "worst_genus", "mass_absent", "unresolved"]
UNRESOLVED = "unresolved"
# superresolution's bucket for reads a `panel:` reinterpretation cannot explain. It is
# not a genus, and the truth has no mass there, so it scores like `unresolved`.
BACKGROUND = "background"


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


def genome_genera(path, lineage, column="inferred_mean"):
    """{genus | 'unresolved': share} from an inferred_composition.csv, whose rows are
    panel genomes. `background` is superresolution's out-of-panel bucket."""
    out = defaultdict(float)
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            gid = (row.get("genome_id") or "").strip()
            if not gid:
                continue
            key = UNRESOLVED if gid == BACKGROUND else ss.genus(lineage[gid])
            out[key or UNRESOLVED] += float(row[column])
    return dict(out)


def krona_genera(path):
    """{genus | 'unresolved': share} from an amplicon-analysis-pipeline krona table:
    `<count>\t<rank>\t<rank>...`, one line per assigned lineage.

    Lenient about the leading rank on purpose - MAPseq writes the domain first, but a
    krona table built through Root/unclassified prefixes is common enough that dropping
    them is cheaper than pinning one writer's layout.
    """
    out = defaultdict(float)
    for line in Path(path).read_text().splitlines():
        fields = [f.strip() for f in line.split("\t")]
        if len(fields) < 2:
            continue
        try:
            count = float(fields[0])
        except ValueError:                  # a header line, if the writer emits one
            continue
        ranks = [r for r in fields[1:] if r and r.lower() not in ("root", "unclassified")]
        out[ss.genus(";".join(ranks)) or UNRESOLVED] += count
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
    for sample, directory, depth, _pair, em in ss.benchmark_dirs(cfg, results_dir):
        truths = sorted(directory.glob("*.truth.tsv"))
        if not truths:
            skipped += 1
            continue
        truth = truth_genera(truths[0], lineage)
        run_id = f"{sample}.sub{depth}" if depth else sample
        sr = directory / "profiling" / "sr"

        # arm -> {method: callable() -> {genus: share}}, only for files that exist.
        methods = []
        groups = {name: sr / f"{run_id}.{name}.inferred_v4_groups.csv" for name in names}
        for name, path in groups.items():
            methods.append(("silva", name, path, lambda p=path: predicted_genera(p, "inferred_mean")))
        first = groups[names[0]]
        methods.append(("silva", "mapseq_only", first,
                        lambda p=first: predicted_genera(p, "observed_rel_abundance")))
        for arm in ("custom", "silva_panel"):
            for entry in ss.arm_settings(cfg, arm):
                path = sr / f"{run_id}.{entry['name']}.inferred_composition.csv"
                methods.append((arm, entry["name"], path,
                                lambda p=path: genome_genera(p, lineage)))
        # AAP publishes its whole output tree; the krona table is the one file that is a
        # count per lineage, which is all the genus roll-up needs.
        for path in sorted(directory.glob(f"profiling/aap/{run_id}/taxonomy-summary/**/*krona.txt")):
            methods.append(("aap", "aap", path, lambda p=path: krona_genera(p)))

        for arm, method, path, predict in methods:
            if not path.exists():
                continue
            row = {"sample": sample, "error_model": em["name"] or "default",
                   "depth": depth or "full", "arm": arm, "method": method,
                   **scores(predict(), truth)}
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

    # The genome-space arms (custom / silva_panel): panel ids rolled up through the same
    # lineages, with `background` landing on `unresolved`.
    with tempfile.TemporaryDirectory() as d:
        comp = Path(d) / "S.custom.inferred_composition.csv"
        comp.write_text("sample,genome_id,observed_rel_abundance,inferred_mean\n"
                        "S,bf,0.2,0.3\nS,bt,0.2,0.2\nS,ss,0.5,0.4\n"
                        f"S,{BACKGROUND},0.1,0.1\n")
        g = genome_genera(comp, {"bf": bac, "bt": bac, "ss": strep})
    assert abs(g[bac] - 0.5) < 1e-9 and abs(g[UNRESOLVED] - 0.1) < 1e-9, g

    # The aap arm: a krona count-per-lineage table. A leading Root rank is dropped, and a
    # lineage that stops above genus is unresolved.
    with tempfile.TemporaryDirectory() as d:
        krona = Path(d) / "S_SILVA-SSU.krona.txt"
        krona.write_text("\t".join(["50"] + bac.split(";")) + "\n"
                         + "\t".join(["30", "Root"] + strep.split(";")) + "\n"
                         + "\t".join(["20"] + fam.split(";")) + "\n")
        k = krona_genera(krona)
    assert k == {bac: 50.0, strep: 30.0, UNRESOLVED: 20.0}, k
    assert abs(scores(k, {bac: 0.5, strep: 0.3})["unresolved"] - 0.2) < 1e-9
    print("score_sweep selfcheck OK")


if __name__ == "__main__":
    if "--selfcheck" in sys.argv:
        _selfcheck()
    else:
        main()
