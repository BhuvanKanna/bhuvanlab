# -*- coding: utf-8 -*-
"""
simulate_truncation.py

Simulation harness: does a truncation metric measure truncation, or expression?

Draws synthetic genes with a **known** truncation on a known side, at a range of
expression levels, censors them at the -1 detection floor exactly as
``log2(TPM+1)-1`` does, and scores two families of metric on each:

- **old** -- the 4-parameter fit on the excluded (> -1) values, as in
  ``outputs/``: ``rti``, ``lti``, ``rti_sigma_dist``, ``lti_sigma_dist``.
- **new** -- ``tailmodel.fit_tail_model`` on all values: ``rti_missing``,
  ``lti_missing``, ``rti_lr``, ``lti_lr``.

A metric worth keeping passes three checks, printed at the end:

1. **Expression confound** -- on untruncated Gaussians, a high-expression gene
   and a low-expression gene must be indistinguishable (AUC ~ 0.5). This is the
   synthetic version of the expression-matched control result.
2. **Power** -- at a fixed expression level, truncated genes must outscore
   untruncated ones (AUC well above 0.5).
3. **Skew specificity** -- an untruncated but left-skewed gene should not read
   as truncated. Reported, not required: no tail metric fully separates a short
   tail from a cut one, but it is worth knowing how much it confuses them.

Latent sigma is 1 throughout; every metric here is location/scale-free away
from the floor, so the expression level ``mu`` is the only thing that matters,
and it matters only through how much of the distribution the floor eats.

    python simulate_truncation.py                    # full grid, 200 reps/cell
    python simulate_truncation.py --reps 40 --n 300  # quick look
"""

import argparse
import os
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, skewnorm, truncnorm

from bhuvanfitter import BhuvanFitter
from tailmodel import FLOOR, fit_tail_model

MUS = (-1.5, -1.0, -0.5, 0.0, 0.5, 1.5, 3.0)
NS = (100, 300, 800)
CUTS = (2.0, 1.0, 0.5)                  # sigma from the latent mean
SKEW_A = -4.0                           # skew-normal shape for the skew null
SEED = 20261009
OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "truncation_sim"

# Every metric oriented so that higher = "more truncated".
OLD = {"rti": 1, "lti": 1, "rti_sigma_dist": -1, "lti_sigma_dist": -1}
NEW = {"rti_missing": 1, "lti_missing": 1, "rti_lr": 1, "lti_lr": 1}
RIGHT = ["rti", "rti_sigma_dist", "rti_missing", "rti_lr"]
LEFT = ["lti", "lti_sigma_dist", "lti_missing", "lti_lr"]


def scenarios():
    """(scenario, side, cut) -- the truncations each expression level gets."""
    yield "untruncated", "none", np.inf
    for c in CUTS:
        yield f"right_{c:g}", "right", c
    for c in CUTS:
        yield f"left_{c:g}", "left", c
    yield "leftskew", "none", np.inf


def draw(rng, n, mu, side, cut, scenario):
    if scenario == "leftskew":
        m, v = skewnorm.stats(SKEW_A, moments="mv")
        z = (skewnorm.rvs(SKEW_A, size=n, random_state=rng) - m) / np.sqrt(v)
    else:
        lo = -cut if side == "left" else -np.inf
        hi = cut if side == "right" else np.inf
        z = truncnorm.rvs(lo, hi, size=n, random_state=rng)
    y = mu + z
    return np.where(y <= FLOOR, FLOOR, y)


def score(values):
    row = {k: np.nan for k in (*OLD, *NEW)}
    row["n_censored"] = int(np.sum(values <= FLOOR))
    detected = values[values > FLOOR]
    if detected.size >= 10:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                r = BhuvanFitter(detected, gene_name="sim").fit("fourparam",
                                                                max_nfev=2000)
            row.update({k: r[k] for k in OLD})
        except RuntimeError:
            pass
    t = fit_tail_model(values)
    if t["success"]:
        row.update({k: t[k] for k in NEW})
    return row


def run_cell(args):
    n, mu, scenario, side, cut, reps, seed = args
    rng = np.random.default_rng(seed)
    rows = []
    for rep in range(reps):
        row = score(draw(rng, n, mu, side, cut, scenario))
        row.update(n=n, mu=mu, scenario=scenario, side=side, cut=cut, rep=rep)
        rows.append(row)
    return rows


def auc(pos, neg):
    """P(pos > neg) + 0.5 P(tie), NaNs dropped. NaN if either side is empty."""
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    pos, neg = pos[np.isfinite(pos)], neg[np.isfinite(neg)]
    if not pos.size or not neg.size:
        return np.nan
    ranks = pd.Series(np.concatenate([pos, neg])).rank().to_numpy()
    u = ranks[:pos.size].sum() - pos.size * (pos.size + 1) / 2
    return float(u / (pos.size * neg.size))


def oriented(df, metric):
    return df[metric] * {**OLD, **NEW}[metric]


def summarize(df):
    out = {}
    null = df[df.scenario == "untruncated"]

    # 1. confound: untruncated high-expression vs untruncated low-expression
    rows = []
    for n in sorted(df.n.unique()):
        hi = null[(null.n == n) & (null.mu == 3.0)]
        for lo_mu in (-0.5, 0.0, 0.5):
            lo = null[(null.n == n) & (null.mu == lo_mu)]
            for m in (*RIGHT, *LEFT):
                rows.append({"n": n, "high_mu": 3.0, "low_mu": lo_mu, "metric": m,
                             "auc_low_reads_more_truncated":
                                 auc(oriented(lo, m), oriented(hi, m)),
                             "defined_low": float(lo[m].notna().mean())})
    out["confound"] = pd.DataFrame(rows)

    # 2. power: truncated vs untruncated at the same n and mu
    rows = []
    for (n, mu), cell in df.groupby(["n", "mu"]):
        base = cell[cell.scenario == "untruncated"]
        for scen, side in (("right_2", "right"), ("right_1", "right"),
                           ("right_0.5", "right"), ("left_2", "left"),
                           ("left_1", "left"), ("left_0.5", "left")):
            trunc = cell[cell.scenario == scen]
            for m in (RIGHT if side == "right" else LEFT):
                rows.append({"n": n, "mu": mu, "scenario": scen, "metric": m,
                             "auc": auc(oriented(trunc, m), oriented(base, m)),
                             "defined": float(trunc[m].notna().mean())})
    out["power"] = pd.DataFrame(rows)

    # 3. skew: untruncated left-skewed vs untruncated Gaussian, same n and mu
    rows = []
    for (n, mu), cell in df.groupby(["n", "mu"]):
        base = cell[cell.scenario == "untruncated"]
        sk = cell[cell.scenario == "leftskew"]
        for m in RIGHT:
            rows.append({"n": n, "mu": mu, "metric": m,
                         "auc": auc(oriented(sk, m), oriented(base, m))})
    out["skew"] = pd.DataFrame(rows)

    # 4. recovery: median estimated missing fraction vs the truth
    rows = []
    for scen, side, cut in scenarios():
        if side == "none":
            continue
        m = "rti_missing" if side == "right" else "lti_missing"
        sub = df[df.scenario == scen]
        for (n, mu), cell in sub.groupby(["n", "mu"]):
            rows.append({"n": n, "mu": mu, "scenario": scen,
                         "true_missing": float(norm.sf(cut)),
                         "median_est": float(cell[m].median()),
                         "defined": float(cell[m].notna().mean())})
    out["recovery"] = pd.DataFrame(rows)
    return out


def print_report(s):
    pd.set_option("display.width", 160)
    pd.set_option("display.max_columns", 20)

    print("\n=== 1. EXPRESSION CONFOUND (untruncated genes only) ===")
    print("AUC that a LOW-expression gene reads 'more truncated' than one at mu=3.")
    print("Ideal 0.5. Far from 0.5 = the metric is measuring expression.\n")
    c = s["confound"].pivot_table(index=["metric"], columns=["n", "low_mu"],
                                  values="auc_low_reads_more_truncated")
    print(c.round(2).to_string())

    print("\n=== 2. POWER (truncated vs untruncated, same n and mu) ===")
    print("AUC, higher = better. Rows: metric x truncation; columns: n, mu.\n")
    p = s["power"]
    for scen in ("right_1", "left_1"):
        print(f"-- {scen} (cut 1 sigma from the mean) --")
        print(p[p.scenario == scen].pivot_table(index="metric", columns=["n", "mu"],
                                                values="auc").round(2).to_string())
        print()

    print("=== 3. SKEW SPECIFICITY (left-skewed, untruncated vs Gaussian) ===")
    print("AUC of the skewed gene reading 'right-truncated'. 0.5 = not fooled.\n")
    print(s["skew"].pivot_table(index="metric", columns=["n", "mu"],
                                values="auc").round(2).to_string())

    r = s["recovery"]
    r = r[r.n == r.n.max()]
    print(f"\n=== 4. RECOVERY of the missing fraction (new metric, n={r.n.max()}) ===\n")
    print(r.pivot_table(index=["scenario", "true_missing"], columns="mu",
                        values="median_est").round(3).to_string())


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--n", type=int, nargs="*", default=list(NS))
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--out", type=Path, default=OUT_DIR)
    a = ap.parse_args()

    tasks, k = [], 0
    for n in a.n:
        for mu in MUS:
            for scen, side, cut in scenarios():
                tasks.append((n, mu, scen, side, cut, a.reps, SEED + k))
                k += 1

    t0 = time.time()
    print(f"{len(tasks)} cells x {a.reps} reps on {a.jobs} workers ...", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        for chunk in ex.map(run_cell, tasks):
            rows.extend(chunk)
    df = pd.DataFrame(rows)
    print(f"done in {time.time() - t0:.0f}s", flush=True)

    a.out.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out / "sim_rows.csv", index=False, lineterminator="\n")
    warnings.simplefilter("ignore", RuntimeWarning)    # all-NaN cells: nothing defined there
    s = summarize(df)
    for name, table in s.items():
        table.to_csv(a.out / f"summary_{name}.csv", index=False, lineterminator="\n")
    print_report(s)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
