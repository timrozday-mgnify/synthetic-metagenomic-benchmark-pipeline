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


def test_build_mapseq_refs_and_otu(tmp_path):
    # Two collection sequences, each a 1-record 16S FASTA + explicit taxonomy.
    ssu_a = tmp_path / "a.fasta"
    ssu_a.write_text(">seqA\nACGTACGT\n")
    ssu_b = tmp_path / "b.fasta"
    ssu_b.write_text(">seqB\nTTTTGGGG\n")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        f"frag\t{ssu_a}\tBacteria;Bacteroides;fragilis\n"
        f"ramo\t{ssu_b}\tBacteria;Clostridium;ramosum\n"
    )
    fasta, tax, headers = (tmp_path / f"o.{e}" for e in ("fasta", "tax", "json"))
    out = run(
        "build_mapseq_refs.py",
        "--manifest", str(manifest),
        "--out-fasta", str(fasta),
        "--out-tax", str(tax),
        "--out-headers", str(headers),
    )
    assert out.returncode == 0, out.stderr
    # Headers rewritten `>{id}|{i}|{orig}`; tax carries the explicit strings.
    assert ">frag|0|seqA" in fasta.read_text() and ">ramo|0|seqB" in fasta.read_text()
    tax_lines = tax.read_text().splitlines()
    assert tax_lines[0].startswith("#cutoff:") and tax_lines[2].startswith("#levels:")
    assert "frag|0|seqA\tBacteria;Bacteroides;fragilis" in tax_lines
    meta = json.loads(headers.read_text())
    assert meta["header_to_id"] == {"frag|0|seqA": "frag", "ramo|0|seqB": "ramo"}

    # Majority-vote OTU: cluster 0 = both frag seqs + one ramo -> frag; cluster 1 = ramo.
    meta["fasta_order"] = ["frag|0|seqA", "x", "ramo|0|seqB"]
    meta["header_to_id"] = {"frag|0|seqA": "frag", "x": "frag", "ramo|0|seqB": "ramo"}
    headers.write_text(json.dumps(meta))
    mscluster = tmp_path / "c.mscluster"
    mscluster.write_text("0 0 1 2\n1 2\n")
    otu = tmp_path / "o.otu"
    out = run(
        "build_mapseq_otu.py",
        "--mscluster", str(mscluster),
        "--headers", str(headers),
        "--out-otu", str(otu),
    )
    assert out.returncode == 0, out.stderr
    rows = [ln.split("\t") for ln in otu.read_text().splitlines()]
    assert rows[0] == ["0", "Bacteria;Bacteroides;fragilis", "0"], rows
    assert rows[1] == ["1", "Bacteria;Clostridium;ramosum", "1"], rows


def test_build_mapseq_selfchecks():
    assert run("build_mapseq_refs.py", "--selfcheck").returncode == 0
    assert run("build_mapseq_otu.py", "--selfcheck").returncode == 0


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
