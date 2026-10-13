# -*- coding: utf-8 -*-
"""
evaluate_utc.py

Tests P1-P4 and S1 from ../utc/PLAN.md, exactly as locked there.

    cd fourparam
    python evaluate_utc.py --collins ../../../collins2022-TableS7-pHaplo-pTriplo.tsv

Reads ../utc/gene_scores.csv (from utc.py). Writes ../utc/evaluation.json and
prints a verdict per test.

AUC is the Mann-Whitney probability that a matched positive scores lower than
a matched negative (ties count 1/2), computed over the matched positive and
negative samples; the bootstrap resamples matched pairs.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
UTC_DIR = HERE.parent / "utc"
GENELISTS = HERE.parent / "genelists"
DOCS = HERE.parent.parent / "docs"

CALIPER = 0.1
BOOT = 2000
SEED = 0


def resolve(name: str) -> set[str]:
    """Gene-list tokens -> unversioned ids, symbol first, then id (as the browser)."""
    by_sym, ids = {}, set()
    for line in (DOCS / "genes.tsv").open(encoding="utf-8"):
        if not line.strip():
            continue
        gid, sym = line.rstrip("\n").split("\t")[:2]
        gid = gid.split(".", 1)[0]
        ids.add(gid)
        by_sym.setdefault(sym.upper(), gid)
    out = set()
    for tok in (GENELISTS / name).open(encoding="utf-8"):
        tok = tok.split("#", 1)[0].strip()
        if not tok:
            continue
        g = by_sym.get(tok.upper())
        if g is None and tok.split(".", 1)[0].upper() in ids:
            g = tok.split(".", 1)[0].upper()
        if g:
            out.add(g)
    return out


def match(pos: pd.DataFrame, neg: pd.DataFrame, seed: int = SEED) -> list[tuple[str, str]]:
    """Greedy 1:1 nearest-neighbor on `expr`, caliper CALIPER, no replacement,
    positives in random order."""
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(pos))
    neg_expr = neg["expr"].to_numpy()
    used = np.zeros(len(neg), bool)
    pairs = []
    for i in order:
        d = np.abs(neg_expr - pos["expr"].iat[i])
        d[used] = np.inf
        j = int(np.argmin(d))
        if d[j] <= CALIPER:
            used[j] = True
            pairs.append((pos["gene"].iat[i], neg["gene"].iat[j]))
    return pairs


def auc_lower(p: np.ndarray, n: np.ndarray) -> float:
    """P(positive < negative), ties 1/2."""
    r = rankdata(np.concatenate([p, n]))
    u = r[:len(p)].sum() - len(p) * (len(p) + 1) / 2     # P(pos > neg) form
    return 1.0 - u / (len(p) * len(n))


def pair_aucs(scores: pd.DataFrame, pairs, metrics, boot=BOOT, seed=SEED) -> dict:
    s = scores.set_index("gene")
    P = {m: s.loc[[a for a, _ in pairs], m].to_numpy() for m in metrics}
    N = {m: s.loc[[b for _, b in pairs], m].to_numpy() for m in metrics}
    point = {m: auc_lower(P[m], N[m]) for m in metrics}
    rng = np.random.default_rng(seed)
    draws = {m: np.empty(boot) for m in metrics}
    for b in range(boot):
        idx = rng.integers(0, len(pairs), len(pairs))
        for m in metrics:
            draws[m][b] = auc_lower(P[m][idx], N[m][idx])
    return {"point": point, "draws": draws}


def ci(v: np.ndarray) -> list[float]:
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]


def partial_spearman(df: pd.DataFrame, a: str, b: str, covs: list[str],
                     boot=BOOT, seed=SEED) -> dict:
    R = {c: rankdata(df[c]) for c in [a, b] + covs}

    def one(idx):
        X = np.column_stack([np.ones(len(idx))] + [R[c][idx] for c in covs])
        ra = R[a][idx] - X @ np.linalg.lstsq(X, R[a][idx], rcond=None)[0]
        rb = R[b][idx] - X @ np.linalg.lstsq(X, R[b][idx], rcond=None)[0]
        return float(np.corrcoef(ra, rb)[0, 1])

    n = len(df)
    point = one(np.arange(n))
    rng = np.random.default_rng(seed)
    draws = np.array([one(rng.integers(0, n, n)) for _ in range(boot)])
    return {"rho": point, "ci": ci(draws), "n": n}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--collins", type=Path, required=True,
                   help="Collins et al. 2022 Table S7 (ensgid, pHaplo, pTriplo)")
    args = p.parse_args()

    scores = pd.read_csv(UTC_DIR / "gene_scores.csv").dropna(subset=["utc", "ltc", "var"])
    pos_ids, neg_ids = resolve("pos_pTriplo.txt"), resolve("neg_dup_tolerant.txt")
    both = pos_ids & neg_ids
    pos = scores[scores.gene.isin(pos_ids - both)].reset_index(drop=True)
    neg = scores[scores.gene.isin(neg_ids - both)].reset_index(drop=True)
    pairs = match(pos, neg)
    res = {"n_scored_genes": len(scores), "n_pos_scored": len(pos),
           "n_neg_scored": len(neg), "n_overlap_dropped": len(both),
           "n_pairs": len(pairs)}
    s = scores.set_index("gene")
    res["matched_expr_auc"] = auc_lower(s.loc[[a for a, _ in pairs], "expr"].to_numpy(),
                                        s.loc[[b for _, b in pairs], "expr"].to_numpy())

    a = pair_aucs(scores, pairs, ["utc", "ltc", "var"])
    res["P1_utc_auc"] = {"auc": a["point"]["utc"], "ci": ci(a["draws"]["utc"])}
    res["P1_pass"] = res["P1_utc_auc"]["ci"][0] > 0.5
    diff = a["draws"]["utc"] - a["draws"]["ltc"]
    res["P2_ltc_auc"] = {"auc": a["point"]["ltc"], "ci": ci(a["draws"]["ltc"])}
    res["P2_utc_minus_ltc"] = {"diff": a["point"]["utc"] - a["point"]["ltc"], "ci": ci(diff)}
    res["P2_pass"] = res["P2_utc_minus_ltc"]["ci"][0] > 0
    res["S1_var_auc"] = {"auc": a["point"]["var"], "ci": ci(a["draws"]["var"])}

    col = pd.read_csv(args.collins, sep="\t", usecols=["ensgid", "pHaplo", "pTriplo"])
    j = scores.merge(col, left_on="gene", right_on="ensgid").dropna(subset=["pHaplo", "pTriplo"])
    res["P3a_utc_pTriplo"] = partial_spearman(j, "utc", "pTriplo", ["pHaplo", "expr"])
    res["P3a_pass"] = res["P3a_utc_pTriplo"]["ci"][1] < 0
    res["P3b_ltc_pHaplo"] = partial_spearman(j, "ltc", "pHaplo", ["pTriplo", "expr"])
    res["P3b_pass"] = res["P3b_ltc_pHaplo"]["ci"][1] < 0

    halves = {}
    for parity in (0, 1):
        keep = lambda df: df[df.gene.str[4:].astype(int) % 2 == parity]
        hp = match(keep(pos).reset_index(drop=True), keep(neg).reset_index(drop=True))
        h = pair_aucs(scores, hp, ["utc"], boot=0)
        halves["even" if parity == 0 else "odd"] = {"n_pairs": len(hp), "auc": h["point"]["utc"]}
    res["P4_halves"] = halves
    res["P4_pass"] = all(v["auc"] > 0.5 for v in halves.values())
    res["supported"] = bool(res["P1_pass"] and res["P2_pass"])

    (UTC_DIR / "evaluation.json").write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
    f = lambda d: f"{d['auc']:.3f} [{d['ci'][0]:.3f}-{d['ci'][1]:.3f}]"
    print(f"scored genes {len(scores):,}; positives {len(pos)}, negatives {len(neg)}, "
          f"overlap dropped {len(both)}; matched pairs {len(pairs)} "
          f"(expression AUC {res['matched_expr_auc']:.3f})")
    print(f"P1  UTC AUC            {f(res['P1_utc_auc'])}  {'PASS' if res['P1_pass'] else 'FAIL'}")
    print(f"P2  LTC AUC            {f(res['P2_ltc_auc'])}")
    d = res["P2_utc_minus_ltc"]
    print(f"    UTC - LTC          {d['diff']:+.3f} [{d['ci'][0]:+.3f}-{d['ci'][1]:+.3f}]  "
          f"{'PASS' if res['P2_pass'] else 'FAIL'}")
    for k, lab in (("P3a_utc_pTriplo", "P3a UTC~pTriplo|pHaplo,expr"),
                   ("P3b_ltc_pHaplo", "P3b LTC~pHaplo|pTriplo,expr")):
        v = res[k]
        print(f"{lab}  rho {v['rho']:+.3f} [{v['ci'][0]:+.3f}-{v['ci'][1]:+.3f}] n={v['n']:,}  "
              f"{'PASS' if res[k[:3] + '_pass'] else 'FAIL'}")
    print("P4  halves             " + ", ".join(f"{k} {v['auc']:.3f} (n={v['n_pairs']})"
                                                for k, v in halves.items())
          + f"  {'PASS' if res['P4_pass'] else 'FAIL'}")
    print(f"S1  VAR AUC            {f(res['S1_var_auc'])}")
    print(f"\nHypothesis supported (P1 and P2): {'YES' if res['supported'] else 'NO'}")


if __name__ == "__main__":
    main()
