#!/usr/bin/env python
"""Error-model-driven sub-species composition inference for amplicon runs.

Two stages, run from the reports Taskfile (basicpython env):

``mismapping``
    Estimate the reference-to-reference mis-mapping matrix ``M`` **analytically
    from the trained skiver error model** (no read simulation / re-mapping), plus
    the deterministic genome->reference translation table ``T``. For each pair of
    V4 reference amplicons (A, B) we model mapseq as a max-likelihood assigner and
    form the per-read log-likelihood ratio ``LLR(A,B)=logP(r|A)-logP(r|B)`` for a
    read ``r`` drawn from A's error model. Positions are independent, so
    ``E[LLR]`` (= per-column KL of the emission distributions) and ``Var[LLR]`` are
    closed-form; ``P(A mis-assigned to B) ~= Phi(-E[LLR]/sd[LLR])``.

``infer``
    Bayesian inference (Pyro: NUTS / VI / MLE) on the latent *true genome*
    composition ``theta``. Generative model:
    ``theta ~ Dirichlet(alpha)`` -> ``r_true = theta @ T`` (per reference) ->
    ``r_obs = r_true @ M`` (per reference, observed) -> likelihood on the observed
    per-reference mseq abundance. Started from the observed composition.

Both stages carry a ``--demo`` self-check that needs no run data.

Note: inference uses **Pyro**, not NumPyro, because loading the skiver error model
(`lib.error_application` -> `lib.context_error_models`) already imports pyro, so a
NumPyro+JAX stack would be redundant weight in the same env.

Run:
    subspecies_infer.py mismapping --db-fasta db.fasta --model-pt m.pt -o out/
    subspecies_infer.py infer --run-dir RUN --mismap-dir out/ --mode vi -o out/
    subspecies_infer.py mismapping --demo   # self-check
    subspecies_infer.py infer --demo        # self-check
"""
from __future__ import annotations

import argparse
import gzip
import importlib.util
import json
import logging
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger("subspecies_infer")

_HERE = Path(__file__).resolve().parent
# skiver error-model library (ErrorModel.load, _logits_for_reference, _masked_probs).
_SKIVER_LIB = _HERE.parent.parent / "vendor" / "skiver" / "scripts"

# 515-YF / 806BR (V4). Overridable on the CLI.
DEFAULT_FWD_PRIMER = "GTGYCAGCMGCCGCGGTAA"
DEFAULT_REV_PRIMER = "GGACTACNVGGGTWTCTAAT"

_IUPAC = {
    "A": set("A"), "C": set("C"), "G": set("G"), "T": set("T"),
    "R": set("AG"), "Y": set("CT"), "S": set("GC"), "W": set("AT"),
    "K": set("GT"), "M": set("AC"), "B": set("CGT"), "D": set("AGT"),
    "H": set("ACT"), "V": set("ACG"), "N": set("ACGT"),
}
_COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


# ── FASTA / headers / primers ────────────────────────────────────────────────


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Return ``[(header, sequence)]``; header is the full id line (sans '>')."""
    opener = gzip.open if str(path).endswith(".gz") else open
    records: list[tuple[str, str]] = []
    header, chunks = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header, chunks = line[1:].strip().split()[0], []
            else:
                chunks.append(line.strip().upper())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def genome_of_header(header: str) -> str:
    """DB entry header ``genome|index|orig`` -> genome id (before first '|')."""
    return header.split("|", 1)[0]


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _find_primer(seq: str, primer: str, max_mismatch: int, start: int = 0) -> int | None:
    """Earliest index >= ``start`` where ``primer`` (IUPAC) matches ``seq`` within
    ``max_mismatch`` mismatches, else ``None``."""
    lp = len(primer)
    for i in range(start, len(seq) - lp + 1):
        mm = 0
        for j in range(lp):
            if seq[i + j] not in _IUPAC.get(primer[j], set(seq[i + j])):
                mm += 1
                if mm > max_mismatch:
                    break
        if mm <= max_mismatch:
            return i
    return None


def extract_v4(seq: str, fwd: str, rev: str, max_mismatch: int) -> str | None:
    """In-silico PCR: return the amplicon between the forward primer and the
    reverse-complement of the reverse primer, or ``None`` if not amplifiable.

    ponytail: exact-window IUPAC scan with a mismatch budget; AmpliconHunter may
    amplify a few edge cases this misses, which just drops them from the estimate.
    """
    f = _find_primer(seq, fwd, max_mismatch)
    if f is None:
        return None
    amp_start = f + len(fwd)
    rc = revcomp(rev)
    r = _find_primer(seq, rc, max_mismatch, start=amp_start)
    if r is None:
        return None
    amplicon = seq[amp_start:r]
    return amplicon or None


# ── Error-model emission distributions & mis-mapping matrix ───────────────────


def _load_error_model(model_pt: Path, use_vi: bool):
    if str(_SKIVER_LIB) not in sys.path:
        sys.path.insert(0, str(_SKIVER_LIB))
    from lib.error_application import ErrorModel  # noqa: E402  (needs sys.path)

    return ErrorModel.load(model_pt, use_vi=use_vi)


def emission_distribution(model, seq: str) -> np.ndarray:
    """Per-position emitted-base distribution ``[L, 4]`` for a reference amplicon.

    Uses the error model's masked per-position error-type probabilities: emitting
    the true base is the match probability, emitting another base is that
    substitution probability. Insertion/deletion mass is dropped and the four
    base probabilities renormalised (indels are handled at alignment gaps).
    """
    if str(_SKIVER_LIB) not in sys.path:
        sys.path.insert(0, str(_SKIVER_LIB))
    from lib.error_application import _masked_probs  # noqa: E402
    from lib.error_application import _CHAR_TO_IDX, _ERR_SUB_START

    raw = np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)
    ref_idx = _CHAR_TO_IDX[raw].astype(np.int64)
    probs = _masked_probs(model._logits_for_reference(ref_idx, True), ref_idx)  # [L,10]
    emit = np.empty((ref_idx.shape[0], 4), dtype=np.float64)
    for b in range(4):
        emit[:, b] = np.where(ref_idx == b, probs[:, 0], probs[:, _ERR_SUB_START + b])
    emit /= emit.sum(axis=1, keepdims=True)
    return emit


def _align_columns(a: str, b: str) -> list[tuple[int | None, int | None]]:
    """Align ``a`` (query) to ``b`` (target) with edlib; return aligned columns as
    ``(ai, bj)`` where a gap is ``None``."""
    import edlib

    res = edlib.align(a, b, task="path")
    cols: list[tuple[int | None, int | None]] = []
    ai = bj = 0
    for length, op in _parse_cigar(res["cigar"]):
        for _ in range(length):
            if op in "=X":
                cols.append((ai, bj)); ai += 1; bj += 1
            elif op == "I":  # base in query (a), gap in target (b)
                cols.append((ai, None)); ai += 1
            else:            # 'D': base in target (b), gap in query (a)
                cols.append((None, bj)); bj += 1
    return cols


def _parse_cigar(cigar: str) -> list[tuple[int, str]]:
    out, num = [], ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            out.append((int(num), ch)); num = ""
    return out


def _llr_moments(emit_a: np.ndarray, emit_b: np.ndarray, seq_a: str, seq_b: str,
                 gap_penalty: float) -> tuple[float, float]:
    """Mean and variance of ``LLR(A,B)`` for a read drawn from A, aligned A->B.

    Match/mismatch columns contribute the KL divergence (mean) and second-moment
    of the per-base log-ratio under A's emission distribution. Gap columns add a
    fixed ``gap_penalty`` favouring A (near-identical 16S rarely differ by indels).
    """
    mean = var = 0.0
    for ai, bj in _align_columns(seq_a, seq_b):
        if ai is None or bj is None:
            mean += gap_penalty  # deterministic; no variance contribution
            continue
        pa = emit_a[ai]
        pb = np.clip(emit_b[bj], 1e-9, None)
        logratio = np.log(np.clip(pa, 1e-9, None)) - np.log(pb)  # [4]
        e = float((pa * logratio).sum())               # KL(pa||pb) >= 0
        e2 = float((pa * logratio ** 2).sum())
        mean += e
        var += max(e2 - e * e, 0.0)
    return mean, var


def build_mismapping_matrix(emits: list[np.ndarray], seqs: list[str],
                            gap_penalty: float, progress: bool = True) -> np.ndarray:
    """Row-stochastic mis-mapping matrix ``M[a,b]=P(assign a->b)``.

    Off-diagonal weight ``Phi(-E[LLR]/sd[LLR])`` (probability a read from A scores
    higher under B than under A); the row is then normalised with a unit diagonal.
    """
    from scipy.stats import norm

    n = len(seqs)
    M = np.eye(n, dtype=np.float64)
    rows = _progress(range(n), total=n, desc=f"mismapping {n}x{n} refs",
                     enabled=progress, unit="ref")
    for a in rows:
        for b in range(n):
            if a == b:
                continue
            mean, var = _llr_moments(emits[a], emits[b], seqs[a], seqs[b], gap_penalty)
            sd = np.sqrt(var)
            if sd < 1e-9:
                M[a, b] = 1.0 if mean < 1e-9 else 0.0  # identical -> full confusion
            else:
                M[a, b] = float(norm.cdf(-mean / sd))
    M /= M.sum(axis=1, keepdims=True)
    return M


def _progress(iterable, *, total=None, desc="", enabled=True, unit="it", leave=False):
    """tqdm progress bar when enabled, else the plain iterable (tqdm ships with pyro)."""
    if not enabled:
        return iterable
    from tqdm.auto import tqdm
    return tqdm(iterable, total=total, desc=desc, unit=unit, leave=leave)


def stage_mismapping(args) -> None:
    log.info("loading DB fasta %s", args.db_fasta)
    records = read_fasta(args.db_fasta)
    log.info("loading error model %s", args.model_pt)
    model = _load_error_model(args.model_pt, args.use_vi)
    log.info("extracting V4 amplicons (primers %s / %s, <=%d mismatches) from %d entries",
             args.fwd_primer, args.rev_primer, args.primer_mismatches, len(records))

    refseqs: list[str] = []          # DB entry headers with a valid V4
    v4: list[str] = []
    genomes_of: list[str] = []
    idx_rows: list[dict] = []
    for header, seq in records:
        amp = extract_v4(seq, args.fwd_primer, args.rev_primer, args.primer_mismatches)
        amplifiable = amp is not None
        idx_rows.append({"refseq": header, "genome": genome_of_header(header),
                         "v4_len": len(amp) if amp else 0, "amplifiable": amplifiable})
        if amplifiable:
            refseqs.append(header)
            v4.append(amp)
            genomes_of.append(genome_of_header(header))
    if not refseqs:
        raise SystemExit("no reference produced a V4 amplicon; check primers / DB fasta")
    log.info("%d/%d entries amplifiable; computing mis-mapping via '%s'",
             len(refseqs), len(records), args.method)

    emits = [emission_distribution(model, s) for s in v4]
    t0 = time.time()
    if args.method == "simulate":
        M = _mismapping_by_simulation(model, refseqs, v4, args)
    else:
        M = build_mismapping_matrix(emits, v4, args.gap_penalty, progress=args.progress)
    log.info("mis-mapping matrix computed in %.1fs", time.time() - t0)

    genomes = sorted(set(genomes_of))
    T = np.zeros((len(genomes), len(refseqs)), dtype=np.float64)
    g_idx = {g: i for i, g in enumerate(genomes)}
    for j, g in enumerate(genomes_of):
        T[g_idx[g], j] = 1.0
    # Rows sum to 1: the within-genome distribution of a genome's amplicon reads over
    # its 16S copies (uniform). This makes the latent theta the true *read-space* genome
    # composition (r_true genome-marginal = theta), directly comparable to the read-space
    # realized truth and observed profile — no spurious copy-number reweighting.
    T = T / T.sum(axis=1, keepdims=True)

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(M, index=refseqs, columns=refseqs).to_csv(out / "mismapping_matrix.csv")
    pd.DataFrame(T, index=genomes, columns=refseqs).to_csv(out / "translation_table.csv")
    pd.DataFrame(idx_rows).to_csv(out / "refseq_index.csv", index=False)
    print(f"mismapping: {len(refseqs)} amplifiable refs / {len(records)} DB entries, "
          f"{len(genomes)} genomes -> {out}")


def _mismapping_by_simulation(model, refseqs, v4, args) -> np.ndarray:
    """Cross-check M: simulate N errored V4 reads per ref, identity-assign via edlib."""
    import edlib
    from lib.error_application import apply_batch

    rng = np.random.default_rng(args.seed)
    n = len(refseqs)
    M = np.zeros((n, n), dtype=np.float64)
    for a in range(n):
        reads = apply_batch(model, [(f"r{a}_{i}", v4[a], True) for i in range(args.n_sim)],
                            rng, emit_quality=False)
        for rec in reads:
            best, best_d = a, None
            for b in range(n):
                d = edlib.align(rec.sequence, v4[b], task="distance")["editDistance"]
                if best_d is None or d < best_d:
                    best, best_d = b, d
            M[a, best] += 1.0
    M /= M.sum(axis=1, keepdims=True)
    return M


# ── Observed abundances from mseq ─────────────────────────────────────────────


def observed_refseq_counts(cell_dir: Path, mseq_glob: str, refseqs: list[str]) -> Counter:
    """Per-reference observed read counts from a cell's mseq files (dbhit == full
    DB header). Only counts hits to references in our amplifiable set."""
    keep = set(refseqs)
    counts: Counter = Counter()
    for path in cell_dir.glob(mseq_glob):
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2 or not fields[1]:
                    continue
                dbhit = fields[1]
                if dbhit in keep:
                    counts[dbhit] += 1
    return counts


# ── Pyro inference ────────────────────────────────────────────────────────────


def composition_model(M, T, alpha, N, y_obs=None, likelihood="dirichlet_multinomial",
                      use_mismapping=True, s_sigma=0.3, od_loc=1.1, od_scale=1.0):
    """Generative model of the observed per-reference counts.

    ``theta`` is the true **read-space** genome composition; ``T`` rows sum to 1
    (within-genome copy distribution) so ``r_true = theta@T`` has genome-marginal
    ``theta``. ``theta ~ Dirichlet`` -> ``r_true = theta@T``. When ``use_mismapping`` a sampled
    scalar ``s`` scales the mis-mapping matrix (``M_eff = (1-s)I + s*M``,
    ``s~LogNormal`` centred at 1) which then acts on ``r_true`` to give
    ``r_obs = r_true@M_eff``; otherwise ``r_obs = r_true`` (no correction — the
    baseline control). Counts follow a Dirichlet-Multinomial whose concentration is
    ``conc_frac*N*r_obs`` — overdispersion parameterised as a *fraction of N* so the
    prior is sample-size-invariant and centred near the multinomial limit rather than
    discarding read-count precision — or plain Multinomial.
    """
    import pyro
    import pyro.distributions as dist
    import torch

    dt = M.dtype
    G, Sdim = T.shape[0], M.shape[0]
    theta = pyro.sample("theta", dist.Dirichlet(alpha * torch.ones(G, dtype=dt)))

    r_true = theta @ T
    r_true = r_true / r_true.sum()
    if use_mismapping:
        # Mis-mapping scale applied to M before it acts on r_true.
        s = pyro.sample("s", dist.LogNormal(torch.tensor(0.0, dtype=dt),
                                            torch.tensor(s_sigma, dtype=dt)))
        M_eff = (1.0 - s) * torch.eye(Sdim, dtype=dt) + s * M
        # A large s can push a diagonal negative for very confusable refs; keep M_eff
        # a valid stochastic matrix. ponytail: clamp+renorm only bites in s's upper tail.
        M_eff = torch.clamp(M_eff, min=0.0)
        M_eff = M_eff / M_eff.sum(-1, keepdim=True)
        r_obs = r_true @ M_eff
        r_obs = r_obs / r_obs.sum()
    else:
        r_obs = r_true

    if y_obs is None:
        counts, total = None, int(N)
    else:
        counts = torch.round(y_obs * N).to(torch.long)
        total = int(counts.sum())   # obs must sum to total_count (rounding drifts ±1)
    if likelihood == "multinomial":
        pyro.sample("y", dist.Multinomial(total_count=total, probs=r_obs), obs=counts)
    else:
        # Overdispersion as a fraction of N (sample-size-invariant, centred near the
        # multinomial limit): concentration = conc_frac * N * r_obs.
        conc_frac = pyro.sample("conc_frac", dist.LogNormal(torch.tensor(od_loc, dtype=dt),
                                                            torch.tensor(od_scale, dtype=dt)))
        conc = conc_frac * float(N) * r_obs + 1e-6
        pyro.sample("y", dist.DirichletMultinomial(concentration=conc, total_count=total),
                    obs=counts)


def _fit(mode, M, T, alpha, N, y_obs, likelihood, theta_init, args, desc=""):
    """Return (theta_samples[K,G] or None, theta_point[G], diagnostics dict, loss_trace)."""
    import pyro
    import torch
    from pyro.infer.autoguide.initialization import init_to_value

    use_mm = getattr(args, "use_mismapping", True)
    show = getattr(args, "progress", True)
    mk = {"y_obs": y_obs, "likelihood": likelihood, "use_mismapping": use_mm}

    pyro.clear_param_store()
    init = {"theta": theta_init}
    if use_mm:
        init["s"] = torch.tensor(1.0, dtype=M.dtype)
    if likelihood != "multinomial":
        init["conc_frac"] = torch.tensor(3.0, dtype=M.dtype)  # ~exp(od_loc), prior median

    def _summ(samples):
        mean = samples.mean(0)
        lo = torch.quantile(samples, 0.05, dim=0)
        hi = torch.quantile(samples, 0.95, dim=0)
        return mean, lo, hi

    if mode == "nuts":
        from pyro.infer import MCMC, NUTS
        log.debug("%s NUTS: %d warmup + %d samples", desc, args.warmup, args.num_samples)
        kernel = NUTS(composition_model, init_strategy=init_to_value(values=init))
        mcmc = MCMC(kernel, num_samples=args.num_samples, warmup_steps=args.warmup,
                    disable_progbar=not show)
        mcmc.run(M, T, alpha, N, **mk)
        s = mcmc.get_samples()["theta"]
        diag = mcmc.diagnostics().get("theta", {})
        rhat = np.atleast_1d(np.asarray(diag.get("r_hat", np.nan))).astype(float)
        ess = np.atleast_1d(np.asarray(diag.get("n_eff", np.nan))).astype(float)
        mean, lo, hi = _summ(s)
        return s, mean, {"r_hat": rhat.tolist(), "n_eff": ess.tolist(),
                         "max_r_hat": float(np.nanmax(rhat))}, None

    from pyro.infer import SVI, Trace_ELBO, Predictive
    from pyro.infer.autoguide import AutoNormal, AutoDelta
    from pyro.optim import Adam

    guide_cls = AutoDelta if mode == "mle" else AutoNormal
    guide = guide_cls(composition_model, init_loc_fn=init_to_value(values=init))
    svi = SVI(composition_model, guide, Adam({"lr": args.lr}), Trace_ELBO())
    log.debug("%s %s: %d SVI steps (lr=%g)", desc, mode.upper(), args.steps, args.lr)
    losses = []
    bar = _progress(range(args.steps), total=args.steps, desc=f"{desc} {mode}",
                    enabled=show, unit="step")
    every = max(1, args.steps // 20)
    for step in bar:
        loss = float(svi.step(M, T, alpha, N, **mk))
        losses.append(loss)
        if step % every == 0:
            if hasattr(bar, "set_postfix"):
                bar.set_postfix(loss=f"{loss:.1f}")
            log.debug("%s %s step %d/%d loss=%.3f", desc, mode, step, args.steps, loss)
    if mode == "mle":
        point = guide.median()["theta"].detach()
        return None, point, {"final_loss": losses[-1]}, losses
    pred = Predictive(composition_model, guide=guide, num_samples=args.num_samples,
                      return_sites=["theta"])
    s = pred(M, T, alpha, N, **{**mk, "y_obs": None})["theta"].squeeze()
    mean, lo, hi = _summ(s)
    return s, mean, {"final_loss": losses[-1]}, losses


def stage_infer(args) -> None:
    import torch

    _bp = _load_module("benchmark_preprocess", _HERE / "benchmark_preprocess.py")

    mm_dir = args.mismap_dir
    M_df = pd.read_csv(mm_dir / "mismapping_matrix.csv", index_col=0)
    T_df = pd.read_csv(mm_dir / "translation_table.csv", index_col=0)
    refseqs = list(M_df.columns)
    genomes = list(T_df.index)
    M = torch.tensor(M_df.to_numpy(), dtype=torch.float64)
    T = torch.tensor(T_df.loc[genomes, refseqs].to_numpy(), dtype=torch.float64)
    refseq_genome = {r: genome_of_header(r) for r in refseqs}

    pipeline_root = args.run_dir / args.pipeline_dir
    cells = _bp.find_cells(pipeline_root)
    if args.assay:
        # M is amplicon-specific: only apply it to cells from the matching assay.
        cells = [c for c in cells if c["assay"] == args.assay]
    if not cells:
        raise SystemExit(f"no cells under {pipeline_root}"
                         + (f" for assay {args.assay}" if args.assay else ""))

    log.info("inference: %d cell(s), %d genomes, %d references, mode=%s likelihood=%s%s",
             len(cells), len(genomes), len(refseqs), args.mode, args.likelihood,
             "" if args.use_mismapping else " (no mis-mapping correction)")

    comp_rows: list[dict] = []
    diag_rows: list[dict] = []
    loss_rows: list[dict] = []
    for ci, cell in enumerate(cells, 1):
        cname = f"{cell['sample']}.{cell['assay']}/{cell['depth']}"
        counts = observed_refseq_counts(cell["dir"], args.mseq_glob, refseqs)
        total = sum(counts.values())
        if total == 0:
            log.info("[%d/%d] %s: no reads, skipping", ci, len(cells), cname)
            continue
        log.info("[%d/%d] %s: %d reads", ci, len(cells), cname, total)
        y = np.array([counts.get(r, 0) for r in refseqs], dtype=np.float64)
        y_rel = torch.tensor(y / total, dtype=torch.float64)

        # Observed genome composition (collapse dbhit->genome): baseline + init.
        obs_genome = defaultdict(float)
        for r, c in counts.items():
            obs_genome[refseq_genome[r]] += c / total
        theta_obs = np.array([obs_genome.get(g, 0.0) for g in genomes])
        theta_init = torch.tensor(np.clip(theta_obs, 1e-4, None), dtype=torch.float64)
        theta_init = theta_init / theta_init.sum()

        truth = _bp.read_truth(next(iter(cell["dir"].glob(_bp.TRUTH_GLOB))))
        truth_map = truth.set_index("genome_id")["realized_rel_abundance"].to_dict()

        samples, point, diag, losses = _fit(
            args.mode, M, T, args.alpha, float(total), y_rel, args.likelihood,
            theta_init, args, desc=f"[{ci}/{len(cells)}] {cname}")
        inferred = point.numpy()
        lo = hi = [np.nan] * len(genomes)
        if samples is not None:
            lo = torch.quantile(samples, 0.05, dim=0).numpy()
            hi = torch.quantile(samples, 0.95, dim=0).numpy()

        tag = {"sample": cell["sample"], "assay": cell["assay"], "depth": cell["depth"],
               "sweep_x": cell["sweep_x"]}
        for i, g in enumerate(genomes):
            comp_rows.append({**tag, "genome_id": g,
                              "observed_rel_abundance": theta_obs[i],
                              "inferred_mean": float(inferred[i]),
                              "inferred_lo": float(lo[i]), "inferred_hi": float(hi[i]),
                              "truth_rel_abundance": float(truth_map.get(g, 0.0))})
        diag_rows.append({**tag, "mode": args.mode, "likelihood": args.likelihood,
                          "use_mismapping": args.use_mismapping, "n_reads": int(total),
                          **{k: json.dumps(v) for k, v in diag.items()}})
        if losses is not None:
            for step, loss in enumerate(losses):
                loss_rows.append({**tag, "step": step, "loss": loss})

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(comp_rows).to_csv(out / "inferred_composition.csv", index=False)
    pd.DataFrame(diag_rows).to_csv(out / "inference_diagnostics.csv", index=False)
    if loss_rows:
        pd.DataFrame(loss_rows).to_csv(out / "loss_trace.csv", index=False)
    print(f"infer[{args.mode}]: {len(cells)} cells, {len(genomes)} genomes -> {out}")


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ── Demos ─────────────────────────────────────────────────────────────────────


def demo_mismapping() -> None:
    # Primer extraction: concretise the degenerate primers to ACGT (real refs are
    # ACGT), flank a payload, and recover it.
    def _concrete(p):  # first ACGT option per IUPAC code
        return "".join(sorted(_IUPAC[c])[0] for c in p)
    payload = "ACGTACGTAA" * 3
    seq = ("TTT" + _concrete(DEFAULT_FWD_PRIMER) + payload
           + _concrete(revcomp(DEFAULT_REV_PRIMER)) + "GGG")
    assert extract_v4(seq, DEFAULT_FWD_PRIMER, DEFAULT_REV_PRIMER, 2) == payload

    sys.path.insert(0, str(_SKIVER_LIB))
    from lib.error_application import ErrorModel
    from lib.encoding import NUM_ERROR_TYPES
    # Context-free model: match likely, subs less likely, no indels.
    logits = np.full((4, NUM_ERROR_TYPES), -6.0, dtype=np.float32)
    logits[:, 0] = 3.0                        # match
    for r in range(4):
        for b in range(4):
            if b != r:
                logits[r, 1 + b] = -1.0       # substitutions
    model = ErrorModel.from_spec("BaseContext(1)", {"logits": logits})

    a = "ACGTACGTACGTACGT"
    b = a[:8] + "T" + a[9:]                    # 1 diff
    c = "TGCATGCATGCATGCA"                      # very different
    emits = [emission_distribution(model, s) for s in (a, b, c)]
    M = build_mismapping_matrix(emits, [a, b, c], gap_penalty=4.0)
    assert np.allclose(M.sum(1), 1.0), M.sum(1)
    assert (np.diag(M) > 0.5).all(), M
    assert M[0, 1] > M[0, 2], M            # near-identical confused more than distant
    # identical refs -> ~half mass each
    Mid = build_mismapping_matrix([emits[0], emits[0]], [a, a], gap_penalty=4.0)
    assert abs(Mid[0, 1] - 0.5) < 1e-6, Mid
    print("demo mismapping: OK")


def demo_infer() -> None:
    import torch
    torch.manual_seed(0)

    # 3 genomes, genome 1 has two near-identical copies (refseqs 1 & 2) that map to
    # each other; genome 0 -> refseq 0, genome 2 -> refseq 3.
    genomes = ["g0", "g1", "g2"]
    T = torch.tensor([[1, 0, 0, 0],
                      [0, 1, 1, 0],
                      [0, 0, 0, 1]], dtype=torch.float64)
    M = torch.tensor([[0.95, 0.03, 0.01, 0.01],
                      [0.02, 0.55, 0.40, 0.03],   # copies 1<->2 heavily confused
                      [0.02, 0.42, 0.53, 0.03],
                      [0.01, 0.02, 0.02, 0.95]], dtype=torch.float64)
    theta_true = torch.tensor([0.2, 0.5, 0.3], dtype=torch.float64)
    r_true = theta_true @ T
    r_obs = (r_true / r_true.sum()) @ M
    r_obs = r_obs / r_obs.sum()

    # Naive observed genome composition (collapse refseq->genome by membership).
    memb = T.argmax(0)
    obs_genome = torch.zeros(3, dtype=torch.float64)
    for s in range(4):
        obs_genome[memb[s]] += r_obs[s]

    args = argparse.Namespace(mode="vi", num_samples=200, warmup=0, steps=1200, lr=0.05,
                              progress=False)
    _, point, diag, losses = _fit("vi", M, T, 0.5, 5000.0, r_obs, "dirichlet_multinomial",
                                  obs_genome / obs_genome.sum(), args)
    inferred = point.numpy()
    err_naive = float(abs(obs_genome / obs_genome.sum() - theta_true).sum())
    err_inf = float(np.abs(inferred - theta_true.numpy()).sum())
    assert losses[-1] < losses[0], (losses[0], losses[-1])
    assert err_inf < err_naive, (err_inf, err_naive)
    print(f"demo infer: L1 naive={err_naive:.3f} -> inferred={err_inf:.3f} OK")


# ── CLI ───────────────────────────────────────────────────────────────────────


def _add_common(p) -> None:
    p.add_argument("--verbose", "-v", action="store_true", help="DEBUG-level step logging")
    p.add_argument("--no-progress", dest="progress", action="store_false",
                   help="disable tqdm progress bars (e.g. non-interactive runs)")
    p.set_defaults(progress=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    m = sub.add_parser("mismapping", help="estimate M and T from the error model")
    m.add_argument("--db-fasta", type=Path)
    m.add_argument("--model-pt", type=Path)
    m.add_argument("--use-vi", action="store_true", help="use the model's VI posterior mean")
    m.add_argument("--fwd-primer", default=DEFAULT_FWD_PRIMER)
    m.add_argument("--rev-primer", default=DEFAULT_REV_PRIMER)
    m.add_argument("--primer-mismatches", type=int, default=3)
    m.add_argument("--method", choices=["likelihood", "simulate"], default="likelihood")
    m.add_argument("--gap-penalty", type=float, default=4.0,
                   help="LLR added per alignment gap column (likelihood method)")
    m.add_argument("--n-sim", type=int, default=500, help="reads/ref for --method simulate")
    m.add_argument("--seed", type=int, default=0)
    m.add_argument("-o", "--output-dir", type=Path)
    m.add_argument("--demo", action="store_true")
    _add_common(m)

    inf = sub.add_parser("infer", help="Bayesian inference of true genome composition")
    inf.add_argument("--run-dir", type=Path)
    inf.add_argument("--pipeline-dir", default="results/subspecies_v4_sweep")
    inf.add_argument("--mismap-dir", type=Path, help="dir with mismapping_matrix.csv etc.")
    inf.add_argument("--mseq-glob", default="profiling/aap/*/taxonomy-summary/*/*.mseq.gz")
    inf.add_argument("--assay", default=None,
                     help="only infer cells from this assay (M is amplicon-specific), "
                          "e.g. 515-YF-806BR")
    inf.add_argument("--mode", choices=["nuts", "vi", "mle"], default="vi")
    # multinomial is the default: with deep amplicon data the read count is highly
    # informative, and Dirichlet-Multinomial overdispersion tends to absorb the
    # mis-mapping-correction signal (use it only when counts are genuinely overdispersed).
    inf.add_argument("--likelihood", choices=["multinomial", "dirichlet_multinomial"],
                     default="multinomial")
    inf.add_argument("--no-mismapping", dest="use_mismapping", action="store_false",
                     help="baseline control: infer without the mis-mapping correction "
                          "(r_obs = r_true, no s term)")
    inf.set_defaults(use_mismapping=True)
    inf.add_argument("--alpha", type=float, default=0.5, help="Dirichlet prior concentration")
    inf.add_argument("--num-samples", type=int, default=500)
    inf.add_argument("--warmup", type=int, default=500)
    inf.add_argument("--steps", type=int, default=3000, help="SVI steps (vi/mle)")
    inf.add_argument("--lr", type=float, default=0.02)
    inf.add_argument("-o", "--output-dir", type=Path)
    inf.add_argument("--demo", action="store_true")
    _add_common(inf)

    args = ap.parse_args()
    logging.basicConfig(
        level=logging.DEBUG if getattr(args, "verbose", False) else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
    if args.cmd == "mismapping":
        if args.demo:
            return demo_mismapping()
        for req in ("db_fasta", "model_pt", "output_dir"):
            if getattr(args, req) is None:
                ap.error(f"--{req.replace('_', '-')} required")
        stage_mismapping(args)
    else:
        if args.demo:
            return demo_infer()
        for req in ("run_dir", "mismap_dir", "output_dir"):
            if getattr(args, req) is None:
                ap.error(f"--{req.replace('_', '-')} required")
        stage_infer(args)


if __name__ == "__main__":
    main()
