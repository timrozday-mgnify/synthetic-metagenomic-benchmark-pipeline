#!/usr/bin/env python3
"""Generate tiny deterministic test fixtures for the pipeline:

  - two small reference genomes (genomeA/genomeB .fasta)
  - a high-coverage "natural" paired-end FASTQ to train the error model on
    (reference-free skiver dump needs recurring k-mers -> high coverage)
  - two genomes CSVs (different abundances) + a samplesheet

Run from the repo root:  python tests/data/generate_fixtures.py
Fixtures are committed so tests need no generation step.
"""

import gzip
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
BASES = "ACGT"
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def random_genome(rng: random.Random, length: int) -> str:
    return "".join(rng.choice(BASES) for _ in range(length))


def write_fasta(path: str, name: str, seq: str, width: int = 70) -> None:
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i : i + width] + "\n")


def mutate(rng: random.Random, base: str, rate: float) -> str:
    if rng.random() < rate:
        return rng.choice([b for b in BASES if b != base])
    return base


def simulate_pairs(rng, genomes, weights, n_pairs, read_len, frag_len, err):
    """Yield (r1_seq, r1_qual, r2_seq, r2_qual) high-coverage paired reads."""
    for _ in range(n_pairs):
        g = rng.choices(genomes, weights=weights, k=1)[0]
        start = rng.randint(0, len(g) - frag_len)
        frag = g[start : start + frag_len]
        r1 = "".join(mutate(rng, b, err) for b in frag[:read_len])
        r2 = "".join(mutate(rng, b, err) for b in revcomp(frag)[:read_len])
        q1 = "".join("I" if x == y else "5" for x, y in zip(r1, frag[:read_len]))
        rc = revcomp(frag)[:read_len]
        q2 = "".join("I" if x == y else "5" for x, y in zip(r2, rc))
        yield r1, q1, r2, q2


def write_fastq_gz(path: str, records) -> None:
    with gzip.open(path, "wt") as fh:
        for i, (seq, qual) in enumerate(records):
            fh.write(f"@read{i}\n{seq}\n+\n{qual}\n")


def main() -> None:
    rng = random.Random(1234)
    genomeA = random_genome(rng, 3000)
    genomeB = random_genome(rng, 2500)
    write_fasta(os.path.join(HERE, "genomeA.fasta"), "genomeA_contig1", genomeA)
    write_fasta(os.path.join(HERE, "genomeB.fasta"), "genomeB_contig1", genomeB)

    pairs = list(
        simulate_pairs(
            rng,
            [genomeA, genomeB],
            weights=[0.6, 0.4],
            n_pairs=6000,
            read_len=150,
            frag_len=300,
            err=0.01,
        )
    )
    write_fastq_gz(
        os.path.join(HERE, "natural_R1.fastq.gz"), [(r1, q1) for r1, q1, _, _ in pairs]
    )
    write_fastq_gz(
        os.path.join(HERE, "natural_R2.fastq.gz"), [(r2, q2) for _, _, r2, q2 in pairs]
    )

    with open(os.path.join(HERE, "genomes_S1.csv"), "w") as fh:
        fh.write("genome_id,fasta_path,abundance\n")
        fh.write("genomeA,tests/data/genomeA.fasta,0.7\n")
        fh.write("genomeB,tests/data/genomeB.fasta,0.3\n")

    with open(os.path.join(HERE, "genomes_S2.csv"), "w") as fh:
        fh.write("genome_id,fasta_path,abundance\n")
        fh.write("genomeA,tests/data/genomeA.fasta,0.4\n")
        fh.write("genomeB,tests/data/genomeB.fasta,0.6\n")

    print("Wrote fixtures to", HERE)


if __name__ == "__main__":
    main()
