#!/usr/bin/env python3
"""Write an amplicon-analysis-pipeline `-c` config exposing one pipeline-built mapseq DB.

Emits the `params.mapseq_databases { <name> { ... } }` block the nested EBI AAP run reads,
pointing at the absolute paths of a mapseq DB quartet built by BUILD_DATABASES. This is the
in-pipeline equivalent of the `aap.config` the standalone `build_profiling_dbs.py` wrote.
"""

import argparse
import os


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--name", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--tax", required=True)
    p.add_argument("--otu", required=True)
    p.add_argument("--mscluster", required=True)
    # Rfam rRNA-detection DBs (params.rrnas_rfam_*). Absolute host paths; AAP fails
    # with "file() ... cannot be empty" if these are unset, so emit them when given.
    p.add_argument("--rfam-covariance-model", dest="rfam_cm")
    p.add_argument("--rfam-claninfo", dest="rfam_claninfo")
    p.add_argument("--output", required=True)
    a = p.parse_args()

    ap = os.path.abspath
    rfam_lines = ""
    if a.rfam_cm:
        rfam_lines += f"    rrnas_rfam_covariance_model = '{a.rfam_cm}'\n"
    if a.rfam_claninfo:
        rfam_lines += f"    rrnas_rfam_claninfo = '{a.rfam_claninfo}'\n"
    with open(a.output, "w") as fh:
        fh.write(
            "params {\n"
            "    mapseq_databases {\n"
            f"        {a.name} {{\n"
            f"            fasta = '{ap(a.fasta)}'\n"
            f"            tax = '{ap(a.tax)}'\n"
            f"            otu = '{ap(a.otu)}'\n"
            f"            mscluster = '{ap(a.mscluster)}'\n"
            f"            label = '{a.name}'\n"
            "            run_otu = true\n"
            "            run_asv = false\n"
            "        }\n"
            "    }\n"
            f"{rfam_lines}"
            "}\n"
        )


if __name__ == "__main__":
    main()
