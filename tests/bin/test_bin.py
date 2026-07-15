"""Unit tests for the bin/ helper scripts (pure-python logic).

Run: pytest tests/bin/test_bin.py
build_truth.py needs pysam + a BAM and is covered by the pipeline e2e test.
"""

import json
import subprocess
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"


def run(script, *args, **kw):
    return subprocess.run(
        [sys.executable, str(BIN / script), *args],
        capture_output=True,
        text=True,
        **kw,
    )


def test_build_model_config_slugs_and_ids():
    out = run("build_model_config.py", "AdditiveContext(5), AdditiveContext(7)")
    assert out.returncode == 0, out.stderr
    cfg = json.loads(out.stdout)
    assert [m["components"] for m in cfg["models"]] == [
        "AdditiveContext(5)",
        "AdditiveContext(7)",
    ]
    assert cfg["models"][0]["id"] == "m0_additivecontext5"


def test_build_model_config_rejects_empty():
    assert run("build_model_config.py", "  ").returncode != 0


def test_pick_best_model_min_aic(tmp_path):
    csv = tmp_path / "aic.csv"
    csv.write_text(
        "model_id,inference,aic\n"
        "m0,maximum_likelihood,120.5\n"
        "m1,maximum_likelihood,90.2\n"      # winner
        "m1,variational_inference,\n"
        "m2,maximum_likelihood,105.0\n"
    )
    out = run("pick_best_model.py", str(csv), "hq-illumina")
    assert out.returncode == 0, out.stderr
    assert out.stdout == "m1"


def test_rewrite_genomes_csv_basenames(tmp_path):
    src = tmp_path / "in.csv"
    src.write_text(
        "genome_id,fasta_path,abundance\n"
        "genomeA,/abs/path/to/genomeA.fasta,0.7\n"
        "genomeB,tests/data/genomeB.fasta,0.3\n"
    )
    dst = tmp_path / "out.csv"
    out = run("rewrite_genomes_csv.py", str(src), str(dst))
    assert out.returncode == 0, out.stderr
    lines = dst.read_text().splitlines()
    assert lines[0] == "genome_id,fasta_path,abundance"
    assert lines[1] == "genomeA,genomeA.fasta,0.7"
    assert lines[2] == "genomeB,genomeB.fasta,0.3"
