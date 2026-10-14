# -*- coding: utf-8 -*-
"""
drti_bench.py

Phases 3-4 of ../drti/PLAN.md: the semi-synthetic spike-in benchmark and the
locked selection rule. No dosage labels are used; genes in the published
control sets are excluded from the pool.

    python drti_bench.py --prep <drti_prep.py output dir> [--genes 4000] [--jobs 4]

Writes ../drti/bench_raw.csv.gz (every statistic for every gene x condition),
../drti/bench_scores.csv (power / false fire / confound correlation per
candidate) and ../drti/selection.json.
"""
from __future__ import annotations

import os

# One math-library thread per worker: ten workers each spawning a thread pool
# oversubscribe the cores and left them ~25% busy.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
import zlib
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from sklearn.linear_model import LogisticRegression

import drti_core as core
from evaluate_utc import resolve

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "drti"

MIN_TISSUES = 10
MIN_A, MIN_B = 3, 5
CONTROL_SETS = ("pos_pTriplo.txt", "pos_decipher.txt", "neg_dup_tolerant.txt", "neg_olfactory.txt")
TRUNC = {"t05": 0.05, "t10": 0.10, "t20": 0.20}
CONFOUNDS = ("comp_tissue", "comp_systemic", "skew", "outliers", "fewer")
CONDITIONS = ("control",) + tuple(TRUNC) + tuple(f"ctl_{t}" for t in TRUNC) + \
    ("hard05", "ctl_hard05") + CONFOUNDS
KNN = 200
REALISTIC = {"skew_a": 0.25, "sys_shift": 0.9}
REALISTIC_CONDS = ("skew", "comp_systemic")
BATCH = 500
FF_MAX, CC_MAX, STACK_GAIN = 0.05, 0.10, 0.02


# ------------------------------------------------------------------ data
def load_genes(prep: Path, genes: set[str] | None):
    """{gene: (donors, tissues, Z, median_log per tissue)} for the chosen genes."""
    cols = defaultdict(dict)
    meds = defaultdict(dict)
    donor_sets = {}
    for f in sorted(prep.glob("*.npz")):
        d = np.load(f, allow_pickle=True)
        t = str(d["tissue"])
        g_all, z, ml = d["genes"], d["z"], d["median_log"]
        donor_sets[t] = d["donors"]
        sel = np.flatnonzero(np.isin(g_all, list(genes))) if genes is not None else np.arange(len(g_all))
        for j in sel:
            cols[g_all[j]][t] = z[:, j]
            meds[g_all[j]][t] = float(ml[j])
    out = {}
    for g, tz in cols.items():
        tissues = sorted(tz)
        donors = sorted(set().union(*(set(donor_sets[t]) for t in tissues)))
        idx = {d: i for i, d in enumerate(donors)}
        Z = np.full((len(donors), len(tissues)), np.nan, dtype=float)
        for j, t in enumerate(tissues):
            rows = [idx[d] for d in donor_sets[t]]
            Z[rows, j] = tz[t]
        out[g] = (np.array(donors), tissues, Z.astype(np.float32),
                  np.array([meds[g][t] for t in tissues]))
    return out


def tissue_counts(prep: Path) -> pd.Series:
    n = defaultdict(int)
    for f in sorted(prep.glob("*.npz")):
        for g in np.load(f, allow_pickle=True)["genes"]:
            n[g] += 1
    return pd.Series(n)


# ------------------------------------------------------------ conditions
def conditions(Z: np.ndarray, gene: str, skew_a: float = 0.5, sys_shift: float = 1.0):
    """Yield (condition, Z_B restricted to analysis donors) for one gene.
    skew_a / sys_shift are the locked confound strengths; REALISTIC holds the
    realistic-strength values (PLAN.md, deviation 2)."""
    rng = np.random.default_rng(zlib.crc32(gene.encode()))
    perm = rng.permutation(Z.shape[1])
    A, B = perm[: Z.shape[1] // 2], perm[Z.shape[1] // 2:]
    kA = np.sum(np.isfinite(Z[:, A]), 1)
    kB = np.sum(np.isfinite(Z[:, B]), 1)
    ok = (kA >= MIN_A) & (kB >= MIN_B)
    if ok.sum() < 100:          # too few donors for any statistic after removals
        return
    LA = np.nanmean(Z[ok][:, A], 1)
    ZB = Z[ok][:, B]
    n = len(LA)
    yield "control", ZB
    for name, f in TRUNC.items():
        flag = LA > np.quantile(LA, 1 - f)
        yield name, ZB[~flag]
        yield f"ctl_{name}", ZB[rng.permutation(n)[: n - flag.sum()]]
    mB = np.nanmean(ZB, 1)
    flag = mB > np.quantile(mB, 0.95)
    yield "hard05", ZB[~flag]
    yield "ctl_hard05", ZB[rng.permutation(n)[: n - flag.sum()]]

    Zc = ZB.copy()
    for t in range(Zc.shape[1]):
        has = np.flatnonzero(np.isfinite(Zc[:, t]))
        hit = rng.choice(has, size=int(round(0.15 * len(has))), replace=False)
        Zc[hit, t] -= 2.0
    yield "comp_tissue", Zc
    Zc = ZB.copy()
    Zc[rng.random(n) < 0.15] -= sys_shift
    yield "comp_systemic", Zc
    Zc = (1 - np.exp(-skew_a * ZB)) / skew_a
    Zc = (Zc - np.nanmean(Zc, 0)) / np.nanstd(Zc, 0)
    yield "skew", Zc
    Zc = ZB.copy()
    for i in np.flatnonzero(rng.random(n) < 0.01):
        has = np.flatnonzero(np.isfinite(Zc[i]))
        Zc[i, rng.choice(has)] += 4.0
    yield "outliers", Zc
    yield "fewer", ZB[rng.permutation(n)[: int(round(0.7 * n))]]


def bench_gene(args, realistic: bool = False):
    g, (donors, tissues, Z, med) = args
    Z = Z.astype(float)            # shipped as float32; optimizers want float64
    rows = []
    icc = core.icc_of(Z)
    expr = float(np.mean(med))
    for cond, ZB in conditions(Z, g, **(REALISTIC if realistic else {})):
        if realistic and cond not in REALISTIC_CONDS:
            continue
        s = core.all_stats(ZB)
        row = {"gene": g, "condition": cond, "n": s.pop("n"), "icc": icc, "expr": expr}
        row.update({f"{o}.{k}": v for (o, k), v in s.items()})
        rows.append(row)
    return rows


def bench_gene_realistic(args):
    """Only the two confound rows, at realistic strength (same random draws)."""
    return bench_gene(args, realistic=True)


# ---------------------------------------------------------------- scoring
def calibrate_knn(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    """Residual of each statistic from its median among the KNN control genes
    nearest in (expr, icc, log n), scaled by their MAD. Reference = control."""
    ref = df[df.condition == "control"]
    feats = lambda d: np.column_stack([d.expr, d.icc, np.log(d.n)])
    F_ref = feats(ref)
    mu, sd = F_ref.mean(0), F_ref.std(0)
    Fr = (F_ref - mu) / sd
    Fq = (feats(df) - mu) / sd
    out = df[["gene", "condition", "n", "icc", "expr"]].copy()
    # nearest neighbors in batches
    k = min(KNN, len(Fr) - 1)
    nn = np.empty((len(Fq), k), dtype=int)
    for i in range(0, len(Fq), 2000):
        d2 = ((Fq[i:i + 2000, None, :] - Fr[None, :, :]) ** 2).sum(-1)
        nn[i:i + 2000] = np.argpartition(d2, k, axis=1)[:, :k]
    for c in cols:
        r = ref[c].to_numpy()
        nb = r[nn]
        med = np.nanmedian(nb, 1)
        mad = 1.4826 * np.nanmedian(np.abs(nb - med[:, None]), 1)
        out[c] = (df[c].to_numpy() - med) / np.where(mad > 0, mad, np.nan)
    return out


def auc(pos: np.ndarray, neg: np.ndarray) -> float:
    pos, neg = pos[np.isfinite(pos)], neg[np.isfinite(neg)]
    if len(pos) < 20 or len(neg) < 20:
        return np.nan
    r = rankdata(np.concatenate([pos, neg]))
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def score_candidates(tab: pd.DataFrame, cols: list[str], genes_half=None) -> pd.DataFrame:
    """Cross-gene AUCs per candidate column. Genes split 50/50 (seed 2) into
    condition vs its control."""
    if genes_half is not None:
        tab = tab[tab.gene.isin(genes_half)]
    genes = np.array(sorted(tab.gene.unique()))
    rng = np.random.default_rng(2)
    side = dict(zip(genes, rng.random(len(genes)) < 0.5))
    wide = {c: tab.pivot(index="gene", columns="condition", values=c) for c in cols}
    in_cond = pd.Series(side)
    ctl = tab[tab.condition == "control"].set_index("gene")
    rows = []
    for c in cols:
        w = wide[c]
        r = {"candidate": c}
        pairs = {**{t: f"ctl_{t}" for t in list(TRUNC) + ["hard05"]},
                 **{k: "control" for k in CONFOUNDS}}
        for cond, base in pairs.items():
            a = w.loc[in_cond[w.index].to_numpy(), cond].to_numpy()
            b = w.loc[~in_cond[w.index].to_numpy(), base].to_numpy()
            r[cond] = auc(a, b)
        r["power"] = np.nanmean([r["t05"], r["t10"]])
        r["FF"] = max(abs(r[k] - 0.5) for k in CONFOUNDS)
        r["CC"] = max(abs(spearmanr(ctl[c], ctl[k], nan_policy="omit")[0])
                      for k in ("expr", "icc", "n"))
        rows.append(r)
    out = pd.DataFrame(rows)
    out["admissible"] = (out.FF <= FF_MAX) & (out.CC <= CC_MAX)
    return out.sort_values("power", ascending=False)


def stack(tab: pd.DataFrame, feats: list[str]):
    """Logistic stack trained on even-ENSG genes, scored on odd (plan step 3)."""
    parity = tab.gene.str[4:].astype(int) % 2
    lab = tab.condition.map(lambda c: 1 if c in ("t05", "t10")
                            else 0 if c in ("control",) + CONFOUNDS else -1)
    tr = tab[(parity == 0) & (lab >= 0)].dropna(subset=feats)
    X = tr[feats].to_numpy()
    mu, sd = X.mean(0), X.std(0)
    model = LogisticRegression(max_iter=2000).fit((X - mu) / sd, lab[tr.index])
    out = tab.copy()
    Xa = out[feats].to_numpy()
    good = np.isfinite(Xa).all(1)
    out["stacked"] = np.nan
    out.loc[good, "stacked"] = model.decision_function((Xa[good] - mu) / sd)
    odd = sorted(tab.gene[parity == 1].unique())
    return out, model, mu, sd, odd


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--prep", type=Path, required=True)
    p.add_argument("--genes", type=int, default=4000)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--realistic", action="store_true",
                   help="Rescore with skew / comp_systemic at realistic strength "
                        "(PLAN.md deviation 2); needs bench_raw.csv.gz.")
    args = p.parse_args()

    OUT.mkdir(exist_ok=True)
    raw_path = OUT / "bench_raw.csv.gz"
    if raw_path.exists():
        tab = pd.read_csv(raw_path)
        print(f"reusing {raw_path.name}: {tab.gene.nunique():,} genes")
    else:
        excluded = set().union(*(resolve(f) for f in CONTROL_SETS))
        nt = tissue_counts(args.prep)
        pool = sorted(g for g in nt[nt >= MIN_TISSUES].index if g not in excluded)
        rng = np.random.default_rng(1)
        chosen = set(rng.choice(pool, size=min(args.genes, len(pool)), replace=False))
        print(f"pool {len(pool):,} genes with >= {MIN_TISSUES} tissues "
              f"({len(excluded):,} control-set genes excluded); benchmarking {len(chosen):,}")
        # Checkpointed in batches of BATCH genes, so an interrupted run resumes
        # from the last finished batch instead of starting over.
        parts = OUT / "bench_parts"
        parts.mkdir(exist_ok=True)
        order = sorted(chosen)
        batches = [order[i:i + BATCH] for i in range(0, len(order), BATCH)]
        part = lambda b: parts / f"part_{b:02d}.csv.gz"
        todo = [b for b in range(len(batches)) if not part(b).exists()]
        if todo:
            print(f"{len(batches) - len(todo)} of {len(batches)} batches already done", flush=True)
            data = load_genes(args.prep, set().union(*(batches[b] for b in todo)))
            with ProcessPoolExecutor(args.jobs) as ex:
                for b in todo:
                    items = [(g, data.pop(g)) for g in batches[b] if g in data]
                    rows = [row for r in ex.map(bench_gene, items, chunksize=4) for row in r]
                    pd.DataFrame(rows).to_csv(part(b), index=False, float_format="%.6g")
                    print(f"  batch {b + 1}/{len(batches)} done", flush=True)
        tab = pd.concat([pd.read_csv(part(b)) for b in range(len(batches))], ignore_index=True)
        tab.to_csv(raw_path, index=False, float_format="%.6g")

    sfx = ""
    if args.realistic:
        sfx = "_realistic"
        real_path = OUT / "bench_raw_realistic.csv.gz"
        if real_path.exists():
            real = pd.read_csv(real_path)
        else:
            order = sorted(tab.gene.unique())
            rows = []
            with ProcessPoolExecutor(args.jobs) as ex:
                for i in range(0, len(order), 2 * BATCH):    # bounded memory
                    data = load_genes(args.prep, set(order[i:i + 2 * BATCH]))
                    rows += [row for r in ex.map(bench_gene_realistic, data.items(), chunksize=4)
                             for row in r]
                    del data
                    print(f"  realistic: {min(i + 2 * BATCH, len(order)):,} genes", flush=True)
            real = pd.DataFrame(rows)
            real.to_csv(real_path, index=False, float_format="%.6g")
        tab = pd.concat([tab[~tab.condition.isin(REALISTIC_CONDS)], real], ignore_index=True)
        print(f"realistic strengths {REALISTIC}: replaced {len(real):,} rows")

    stat_cols = [c for c in tab.columns if "." in c and not c.endswith("missing_frac")]
    cal = calibrate_knn(tab, stat_cols)
    both = tab[["gene", "condition"]].copy()
    for c in stat_cols:
        both["raw:" + c] = tab[c].to_numpy()
        both["knn:" + c] = cal[c].to_numpy()
    both = pd.concat([both, tab[["n", "icc", "expr"]]], axis=1)
    cand = [c for c in both.columns if c.startswith(("raw:", "knn:"))]
    scores = score_candidates(both, cand)
    scores.to_csv(OUT / f"bench_scores{sfx}.csv", index=False, float_format="%.4f")

    adm = scores[scores.admissible]
    sel = {"n_genes": int(tab.gene.nunique()), "n_candidates": len(cand),
           "n_admissible": int(len(adm))}
    pd.set_option("display.width", 250)
    show = ["candidate", "t05", "t10", "t20", "hard05", "power"] + list(CONFOUNDS) + ["FF", "CC", "admissible"]
    print("\nTop 15 by power (all candidates):")
    print(scores[show].head(15).round(3).to_string(index=False))
    if adm.empty:
        sel["chosen"] = None
        print("\nNo admissible candidate.")
    else:
        best = adm.iloc[0]
        sel["best_single"] = best[show].to_dict()
        print("\nTop 10 admissible:")
        print(adm[show].head(10).round(3).to_string(index=False))
        feats = list(adm.candidate)
        st, model, mu, sd, odd = stack(both, feats)
        s_odd = score_candidates(st, ["stacked", best.candidate], genes_half=odd)
        so = s_odd.set_index("candidate")
        gain = so.loc["stacked", "power"] - so.loc[best.candidate, "power"]
        stack_ok = bool(so.loc["stacked", "admissible"] and gain >= STACK_GAIN)
        sel["stack_odd_half"] = so.reset_index()[show].to_dict("records")
        sel["stack_gain"] = float(gain)
        sel["stack_features"] = feats
        sel["stack_coef"] = dict(zip(feats, model.coef_[0].tolist()))
        sel["stack_intercept"] = float(model.intercept_[0])
        sel["stack_mu"] = dict(zip(feats, mu.tolist()))
        sel["stack_sd"] = dict(zip(feats, sd.tolist()))
        sel["chosen"] = "stacked" if stack_ok else best.candidate
        print("\nStack vs best single on held-out odd-ENSG genes:")
        print(s_odd[show].round(3).to_string(index=False))
        print(f"gain {gain:+.3f} (needs >= {STACK_GAIN}); stacked admissible: "
              f"{bool(so.loc['stacked', 'admissible'])}")
        print(f"\nCHOSEN: {sel['chosen']}")
    (OUT / f"selection{sfx}.json").write_text(json.dumps(sel, indent=2, default=str) + "\n",
                                        encoding="utf-8")


if __name__ == "__main__":
    main()
