# -*- coding: utf-8 -*-
"""
drti_evaluate.py

Phase 6 of ../drti/PLAN.md, exactly as locked there, run once on the frozen
dRTI (../drti/gene_scores.csv from drti_final.py).

    python drti_evaluate.py --collins <Collins et al. 2022 Table S7 .tsv>

- P1: pTriplo >= 0.94 vs duplication-tolerant, 1:1 matched on expression
  (caliper 0.1) and ICC (caliper 0.03); AUC = P(positive has higher dRTI);
  supported if the 95% bootstrap CI lower bound > 0.5.
- P2: partial Spearman(dRTI, pTriplo | pHaplo, expr, ICC) > 0, CI excludes 0.
- P3: P1 on each ENSG-parity half, point estimate > 0.5 in both.
- S1: donor-level variability (log IQR of the `mean` object); AUC =
  P(positive LESS variable), the same orientation as utc S1.

dRTI is supported as a dosage marker only if P1 and P2 hold. Writes
../drti/evaluation.json.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from drti_bench import OUT
from evaluate_utc import BOOT, SEED, auc_lower, ci, partial_spearman, resolve

CAL_EXPR, CAL_ICC = 0.1, 0.03


def match2(pos: pd.DataFrame, neg: pd.DataFrame, seed: int = SEED) -> list[tuple[str, str]]:
    """Greedy 1:1 matching inside both calipers, nearest by
    (d_expr/0.1)^2 + (d_icc/0.03)^2, no replacement, positives in random order."""
    rng = np.random.default_rng(seed)
    ne, ni = neg["expr"].to_numpy(), neg["icc"].to_numpy()
    used = np.zeros(len(neg), bool)
    pairs = []
    for i in rng.permutation(len(pos)):
        de = np.abs(ne - pos["expr"].iat[i])
        di = np.abs(ni - pos["icc"].iat[i])
        d = (de / CAL_EXPR) ** 2 + (di / CAL_ICC) ** 2
        d[used | (de > CAL_EXPR) | (di > CAL_ICC)] = np.inf
        j = int(np.argmin(d))
        if np.isfinite(d[j]):
            used[j] = True
            pairs.append((pos["gene"].iat[i], neg["gene"].iat[j]))
    return pairs


def paired(scores: pd.DataFrame, pairs, col: str, higher: bool, boot: int = BOOT) -> dict:
    """AUC over matched samples, P(pos > neg) if higher else P(pos < neg);
    bootstrap over pairs."""
    s = scores.set_index("gene")
    P = s.loc[[a for a, _ in pairs], col].to_numpy(float)
    N = s.loc[[b for _, b in pairs], col].to_numpy(float)
    f = (lambda p, n: 1.0 - auc_lower(p, n)) if higher else auc_lower
    point = f(P, N)
    if not boot:
        return {"auc": point}
    rng = np.random.default_rng(SEED)
    d = np.empty(boot)
    for b in range(boot):
        idx = rng.integers(0, len(pairs), len(pairs))
        d[b] = f(P[idx], N[idx])
    return {"auc": point, "ci": ci(d)}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--collins", type=Path, required=True)
    args = p.parse_args()

    scores = pd.read_csv(OUT / "gene_scores.csv").dropna(subset=["drti", "expr", "icc"])
    pos_ids, neg_ids = resolve("pos_pTriplo.txt"), resolve("neg_dup_tolerant.txt")
    both = pos_ids & neg_ids
    pos = scores[scores.gene.isin(pos_ids - both)].reset_index(drop=True)
    neg = scores[scores.gene.isin(neg_ids - both)].reset_index(drop=True)
    pairs = match2(pos, neg)
    res = {"n_scored_genes": len(scores), "n_pos_scored": len(pos), "n_neg_scored": len(neg),
           "n_overlap_dropped": len(both), "n_pairs": len(pairs)}
    res["matched_expr_auc"] = paired(scores, pairs, "expr", True, boot=0)["auc"]
    res["matched_icc_auc"] = paired(scores, pairs, "icc", True, boot=0)["auc"]

    res["P1_drti_auc"] = paired(scores, pairs, "drti", True)
    res["P1_pass"] = res["P1_drti_auc"]["ci"][0] > 0.5

    col = pd.read_csv(args.collins, sep="\t", usecols=["ensgid", "pHaplo", "pTriplo"])
    j = scores.merge(col, left_on="gene", right_on="ensgid").dropna(subset=["pHaplo", "pTriplo"])
    res["P2_partial"] = partial_spearman(j, "drti", "pTriplo", ["pHaplo", "expr", "icc"])
    res["P2_pass"] = res["P2_partial"]["ci"][0] > 0

    halves = {}
    for parity in (0, 1):
        keep = lambda df: df[df.gene.str[4:].astype(int) % 2 == parity].reset_index(drop=True)
        hp = match2(keep(pos), keep(neg))
        halves["even" if parity == 0 else "odd"] = {
            "n_pairs": len(hp), "auc": paired(scores, hp, "drti", True, boot=0)["auc"]}
    res["P3_halves"] = halves
    res["P3_pass"] = all(v["auc"] > 0.5 for v in halves.values())

    sv = scores.dropna(subset=["var_logiqr"])
    res["S1_var_auc_pos_less_variable"] = paired(sv, [pq for pq in pairs
                                                      if pq[0] in set(sv.gene) and pq[1] in set(sv.gene)],
                                                 "var_logiqr", False)
    res["supported"] = bool(res["P1_pass"] and res["P2_pass"])

    (OUT / "evaluation.json").write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
    f = lambda d: f"{d['auc']:.3f} [{d['ci'][0]:.3f}-{d['ci'][1]:.3f}]"
    print(f"scored genes {len(scores):,}; positives {len(pos)}, negatives {len(neg)}, "
          f"overlap dropped {len(both)}; matched pairs {len(pairs)} "
          f"(expression AUC {res['matched_expr_auc']:.3f}, ICC AUC {res['matched_icc_auc']:.3f})")
    print(f"P1  dRTI AUC          {f(res['P1_drti_auc'])}  {'PASS' if res['P1_pass'] else 'FAIL'}")
    v = res["P2_partial"]
    print(f"P2  partial rho       {v['rho']:+.3f} [{v['ci'][0]:+.3f}-{v['ci'][1]:+.3f}] "
          f"n={v['n']:,}  {'PASS' if res['P2_pass'] else 'FAIL'}")
    print("P3  halves            " + ", ".join(f"{k} {h['auc']:.3f} (n={h['n_pairs']})"
                                               for k, h in halves.items())
          + f"  {'PASS' if res['P3_pass'] else 'FAIL'}")
    print(f"S1  less variable     {f(res['S1_var_auc_pos_less_variable'])}")
    print(f"\ndRTI supported as a dosage marker (P1 and P2): {'YES' if res['supported'] else 'NO'}")


if __name__ == "__main__":
    main()
