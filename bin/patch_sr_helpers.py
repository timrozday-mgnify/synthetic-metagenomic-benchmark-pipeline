#!/usr/bin/env python3
"""Make a pulled superresolution pipeline's Python helpers runnable under noexec.

Nextflow puts a pipeline's ``bin/`` on PATH, but on an HPC work filesystem
mounted ``noexec`` a helper cannot be exec'd: the task dies with exit 126,
``Permission denied``. Rewrite every bare ``helper.py`` invocation in the pulled
repo's ``*.nf`` into ``python "${projectDir}/bin/helper.py"``.

Keyed on helper NAME rather than on module path, so a process added upstream is
covered without another fix here. The verify pass uses a wider net than the
rewrite: any surviving mention of a helper outside a ``bin/`` path is a call
shape this patch does not handle, and is reported as an error rather than
silently shipped to the cluster.

Takes the checkout, or an asset root to find it under: Nextflow's asset layout is
version-dependent (``$NXF_ASSETS/<org>/<name>`` before 26.x, then
``$NXF_ASSETS/.repos/<org>/<name>/clones/<rev>``), so the launcher asks for the
resolved path rather than assuming one. The resolved checkout is printed on stdout.

Usage: patch_sr_helpers.py <repo-or-assets-dir> [--repo org/name]
"""

import argparse
import re
import sys
from pathlib import Path


def patch_repo(repo: Path) -> list[str]:
    """Rewrite helper calls in ``repo``; return a list of problem descriptions."""
    helpers = sorted(p.name for p in (repo / "bin").glob("*.py"))
    if not helpers:
        return [f"no bin/*.py helpers in {repo}"]
    nf_files = [p for p in repo.rglob("*.nf") if ".git" not in p.parts]
    if not nf_files:
        return [f"no *.nf files in {repo}"]

    problems = []
    for nf in nf_files:
        text = original = nf.read_text()
        for name in helpers:
            # A helper starts a command line, a command substitution, or a pipe stage.
            bare = re.compile(
                rf"(^\s*|\$\(|\|\s*){re.escape(name)}(?=[\s)]|$)", re.MULTILINE
            )
            text = bare.sub(rf'\1python "${{projectDir}}/bin/{name}"', text)
        if text != original:
            nf.write_text(text)
        for lineno, line in enumerate(text.splitlines(), 1):
            if line.lstrip().startswith("//"):
                continue
            for name in helpers:
                if name in line and f"bin/{name}" not in line:
                    problems.append(f"{nf}:{lineno}: unpatched call to {name}: {line.strip()}")
    return problems


def find_checkout(root: Path, name: str | None = None) -> Path:
    """Resolve the pipeline checkout at, or somewhere under, ``root``."""
    if (root / "bin").is_dir() and (root / "main.nf").is_file():
        return root
    found = [
        nf.parent
        for nf in root.rglob("main.nf")
        if (nf.parent / "bin").is_dir() and (nf.parent / "nextflow.config").is_file()
        and (name is None or name in str(nf.parent))
    ]
    if not found:
        raise SystemExit(f"patch_sr_helpers.py: no pipeline checkout under {root}")
    # The pipeline root is the shallowest match; anything deeper is a vendored one.
    return min(found, key=lambda p: len(p.parts))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("root", type=Path, help="pipeline checkout, or asset dir holding it")
    ap.add_argument("--repo", help="org/name, to disambiguate a multi-repo asset dir")
    args = ap.parse_args()

    checkout = find_checkout(args.root, args.repo)
    problems = patch_repo(checkout)
    for problem in problems:
        print(problem, file=sys.stderr)
    if problems:
        return 1
    print(checkout)
    return 0


if __name__ == "__main__":
    sys.exit(main())
