# -*- coding: utf-8 -*-
"""
drti_final.py

Phase 5 of ../drti/PLAN.md: the chosen recipe (../drti/selection.json) on every
gene with >= 10 tissues, using all of its tissues (no A/B split).

    python drti_final.py --prep <drti_prep.py output dir> [--selection selection.json]
                         [--jobs 4] [--boot 200] [--boot-genes N]

Writes ../drti/gene_scores.csv: one row per gene with
  drti           the chosen score, calibrated, higher = more truncated
  drti_se        donor-bootstrap SE of drti (blank where not bootstrapped)
  missing_frac   edge_lr's estimated fraction missing above the ceiling
  var_logiqr     log IQR of the donor `mean` object (S1 variability)
  icc, n_donors, n_tissues, expr
  every raw statistic (`<object>.<stat>`), for transparency

Calibration (`knn`) uses all scored genes as the reference: none is spiked, so
all are "control". A stacked recipe standardizes each feature over the scored
genes and applies the benchmark's coefficients.

The bootstrap resamples donors' object values (objects are not rebuilt per
replicate). Statistics with no optimizer are cheap and are bootstrapped for
every gene; if the recipe needs edge_lr, `--boot-genes N` limits the bootstrap
to a fixed random N genes (seed 3), since 200 likelihood fits per gene for
13,000 genes is weeks of CPU.

Processed in checkpointed batches of BATCH genes (../drti/final_parts/).
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

import drti_core as core
from drti_bench import KNN, MIN_TISSUES, OUT, load_genes, tissue_counts

BATCH = 1500
BOOT_SEED = 3


def parse_recipe(sel: dict) -> dict:
    """{'kind': 'single'|'stacked', 'feats': [(cal, obj, stat)], ...}."""
    chosen = sel.get("chosen")
    if not chosen:
        raise SystemExit("selection has no chosen candidate")
    split = lambda c: (c.split(":", 1)[0], *c.split(":", 1)[1].split(".", 1))
    if chosen == "stacked":
        feats = sel["stack_features"]
        return {"kind": "stacked", "names": feats, "feats": [split(c) for c in feats],
                "coef": np.array([sel["stack_coef"][c] for c in feats])}
    return {"kind": "single", "names": [chosen], "feats": [split(chosen)], "coef": np.array([1.0])}


def stat_value(obj_vals: dict, s2: np.ndarray, obj: str, stat: str) -> float:
    v = obj_vals[obj]
    if stat == "edge_lr":
        return core.edge_lr(v, s2)[0]
    v = v[np.isfinite(v)]
    return core.STATS[stat](v) if len(v) >= 50 else np.nan


def score_gene(args):
    (g, (donors, tissues, Z, med)), feats, boot = args
    Z = Z.astype(float)
    o = core.objects(Z)
    row = {"gene": g, "n_tissues": len(tissues), "n_donors": int(o["keep"].sum()),
           "icc": core.icc_of(Z), "expr": float(np.mean(med)),
           "condition": "control", "n": int(o["keep"].sum())}
    for obj in core.OBJECTS:
        v = o[obj][np.isfinite(o[obj])]
        for name, f in core.STATS.items():
            row[f"{obj}.{name}"] = f(v) if len(v) >= 50 else np.nan
    lr, miss = core.edge_lr(o["mean"], o["s2"])
    row["mean.edge_lr"], row["missing_frac"] = lr, miss
    m = o["mean"][np.isfinite(o["mean"])]
    q25, q75 = np.percentile(m, [25, 75]) if len(m) else (np.nan, np.nan)
    row["var_logiqr"] = float(np.log(q75 - q25)) if q75 > q25 else np.nan

    # Bootstrap of the raw feature statistics (calibration is applied later).
    if boot:
        rng = np.random.default_rng(int(g[4:]) + BOOT_SEED)
        n = len(o["mean"])
        draws = np.full((boot, len(feats)), np.nan)
        for b in range(boot):
            idx = rng.integers(0, n, n)
            ov = {k: o[k][idx] for k in core.OBJECTS}
            for j, (_, obj, stat) in enumerate(feats):
                draws[b, j] = stat_value(ov, o["s2"][idx], obj, stat)
        for j, (_, obj, stat) in enumerate(feats):
            row[f"boot:{obj}.{stat}"] = draws[:, j].tolist()
    return row


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--prep", type=Path, required=True)
    p.add_argument("--selection", type=Path, default=OUT / "selection.json")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--boot", type=int, default=200)
    p.add_argument("--boot-genes", type=int, default=None,
                   help="Bootstrap only this many random genes (seed 3).")
    args = p.parse_args()

    sel = json.loads(args.selection.read_text(encoding="utf-8"))
    rec = parse_recipe(sel)
    print(f"recipe: {rec['kind']} {rec['names']}", flush=True)

    nt = tissue_counts(args.prep)
    genes = sorted(nt[nt >= MIN_TISSUES].index)
    boot_set = set(genes)
    if args.boot_genes:
        boot_set = set(np.random.default_rng(BOOT_SEED).choice(genes, args.boot_genes, replace=False))
    print(f"{len(genes):,} genes with >= {MIN_TISSUES} tissues; "
          f"bootstrapping {len(boot_set):,} x {args.boot}", flush=True)

    parts = OUT / "final_parts"
    parts.mkdir(parents=True, exist_ok=True)
    batches = [genes[i:i + BATCH] for i in range(0, len(genes), BATCH)]
    part = lambda b: parts / f"part_{b:02d}.pkl"
    with ProcessPoolExecutor(args.jobs) as ex:
        for b, batch in enumerate(batches):
            if part(b).exists():
                continue
            data = load_genes(args.prep, set(batch))
            items = [((g, data[g]), rec["feats"], args.boot if g in boot_set else 0)
                     for g in batch if g in data]
            rows = list(ex.map(score_gene, items, chunksize=4))
            pd.DataFrame(rows).to_pickle(part(b))
            del data
            print(f"  batch {b + 1}/{len(batches)} done", flush=True)
    tab = pd.concat([pd.read_pickle(part(b)) for b in range(len(batches))], ignore_index=True)

    # Calibrate each feature (and every bootstrap draw with the same shift/scale).
    stat_cols = [c for c in tab.columns if "." in c and not c.startswith("boot:")]
    vals, draws = [], []
    for (c, obj, stat) in rec["feats"]:
        col = f"{obj}.{stat}"
        raw = tab[col].to_numpy(float)
        if c == "knn":
            med, mad = knn_med_mad(tab, col)
        else:
            med, mad = np.zeros(len(raw)), np.ones(len(raw))
        vals.append((raw - med) / mad)
        bc = f"boot:{col}"
        if bc in tab:
            D = np.array([np.asarray(d, float) if isinstance(d, list) else np.full(args.boot, np.nan)
                          for d in tab[bc]])
            draws.append((D - med[:, None]) / mad[:, None])
        else:
            draws.append(None)
    V = np.column_stack(vals)
    if rec["kind"] == "stacked":
        mu, sd = np.nanmean(V, 0), np.nanstd(V, 0)
        drti = ((V - mu) / sd) @ rec["coef"]
        if all(d is not None for d in draws):
            B = sum(((d - mu[j]) / sd[j]) * rec["coef"][j] for j, d in enumerate(draws))
        else:
            B = None
    else:
        drti, B = V[:, 0], draws[0]
    tab["drti"] = drti
    tab["drti_se"] = np.nan
    if B is not None:
        ok = np.isfinite(B).sum(1) >= 2
        tab.loc[ok, "drti_se"] = np.nanstd(B[ok], axis=1, ddof=1)

    keep = ["gene", "drti", "drti_se", "missing_frac", "var_logiqr", "icc", "n_donors",
            "n_tissues", "expr"] + stat_cols
    out = tab[keep].sort_values("gene")
    out.to_csv(OUT / "gene_scores.csv", index=False, float_format="%.6g")
    print(f"wrote gene_scores.csv: {out.drti.notna().sum():,} scored genes; "
          f"median SE {np.nanmedian(out.drti_se):.3f}; SD of drti {np.nanstd(out.drti):.3f}")


def knn_med_mad(tab: pd.DataFrame, col: str) -> tuple[np.ndarray, np.ndarray]:
    """Per gene, the median and 1.4826*MAD of `col` among the KNN genes nearest
    in standardized (expr, icc, log n) -- drti_bench.calibrate_knn with every
    scored gene as the reference."""
    feats = np.column_stack([tab.expr, tab.icc, np.log(tab.n)])
    mu, sd = feats.mean(0), feats.std(0)
    F = (feats - mu) / sd
    k = min(KNN, len(F) - 1)
    r = tab[col].to_numpy(float)
    med = np.empty(len(F))
    mad = np.empty(len(F))
    for i in range(0, len(F), 250):
        d2 = ((F[i:i + 250, None, :] - F[None, :, :]) ** 2).sum(-1)
        nn = np.argpartition(d2, k, axis=1)[:, :k]
        nb = r[nn]
        med[i:i + 250] = np.nanmedian(nb, 1)
        mad[i:i + 250] = 1.4826 * np.nanmedian(np.abs(nb - med[i:i + 250, None]), 1)
    return med, np.where(mad > 0, mad, np.nan)


if __name__ == "__main__":
    main()
