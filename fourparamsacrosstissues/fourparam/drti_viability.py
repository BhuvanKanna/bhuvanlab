# -*- coding: utf-8 -*-
"""
drti_viability.py

Viability checks for donor-level right truncation (../drti/PLAN.md), on the
output of drti_prep.py.

  (a) coverage: donors with >= K tissues expressing each gene
  (b) ICC: share of a gene's (standardized) variance that is donor-shared
      across tissues, with and without covariate + PC removal
  (c) spike-in: split each gene's tissues into halves A and B; remove donors
      whose A-half mean is in the top SPIKE fraction (a truncation acting on a
      latent level, not on the measured values); measure how far UTC on the
      B data moves, against removing the same number of donors at random.
      Single B tissue vs donor-level mean over B tissues.

    python drti_viability.py --prep <dir of .npz> [--jobs 6]
"""
from __future__ import annotations

import argparse
import zlib
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

SPIKE = 0.10
REPS = 60
MIN_TISSUES = 10
MIN_PER_HALF = 3


def utc(v: np.ndarray) -> np.ndarray:
    """UTC along the last axis."""
    q50, q75, q975 = np.percentile(v, [50, 75, 97.5], axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log((q975 - q50) / (q75 - q50))


def load(prep: Path):
    """{gene: {tissue: Series(z by donor)}} for both z and z_raw."""
    by_gene = defaultdict(dict)
    by_gene_raw = defaultdict(dict)
    for f in sorted(prep.glob("*.npz")):
        d = np.load(f, allow_pickle=True)
        t = str(d["tissue"])
        donors, z, zr = d["donors"], d["z"], d["z_raw"]     # read each array once
        for j, g in enumerate(d["genes"]):
            by_gene[g][t] = pd.Series(z[:, j], index=donors)
            by_gene_raw[g][t] = pd.Series(zr[:, j], index=donors)
    return by_gene, by_gene_raw


def frame(tissues: dict) -> pd.DataFrame:
    """donors x tissues, NaN where the donor has no sample."""
    df = pd.DataFrame(tissues)
    return df[~df.index.duplicated()]


def icc(df: pd.DataFrame) -> float:
    k = df.notna().sum(1)
    df = df[k >= 2]
    k = df.notna().sum(1).to_numpy(float)
    if len(df) < 20:
        return np.nan
    vals = df.to_numpy(float)
    means = np.nanmean(vals, 1)
    grand = np.nansum(vals) / k.sum()
    N, D = k.sum(), len(df)
    msb = (k * (means - grand) ** 2).sum() / (D - 1)
    msw = np.nansum((vals - means[:, None]) ** 2) / (N - D)
    k0 = (N - (k ** 2).sum() / N) / (D - 1)
    return float((msb - msw) / (msb + (k0 - 1) * msw))


def spike_effect(x: np.ndarray, flagged: np.ndarray, rng) -> float:
    """Standardized UTC drop from removing `flagged` vs removing as many at random."""
    n_drop = int(flagged.sum())
    if n_drop < 3 or len(x) - n_drop < 60:
        return np.nan
    spiked = utc(x[~flagged])
    keep = len(x) - n_drop
    idx = np.argsort(rng.random((REPS, len(x))), axis=1)[:, :keep]
    ctrl = utc(x[idx])
    sd = ctrl.std(ddof=1)
    return float((ctrl.mean() - spiked) / sd) if sd > 0 else np.nan


def one_gene(args):
    g, tissues, tissues_raw = args
    df = frame(tissues)
    out = {"gene": g, "n_tissues": df.shape[1],
           "icc": icc(df), "icc_raw": icc(frame(tissues_raw))}
    for K in (3, 5, 10):
        out[f"donors_ge{K}"] = int((df.notna().sum(1) >= K).sum())
    if df.shape[1] < MIN_TISSUES:
        return out
    rng = np.random.default_rng(zlib.crc32(g.encode()))
    cols = rng.permutation(df.columns)
    A, B = cols[: len(cols) // 2], cols[len(cols) // 2:]
    kA = df[A].notna().sum(1)
    LA = df[A].mean(1)[kA >= MIN_PER_HALF]
    kB = df.loc[LA.index, B].notna().sum(1)
    mB = df.loc[LA.index, B].mean(1)[kB >= MIN_PER_HALF]
    if len(mB) < 100:
        return out
    cut = LA.quantile(1 - SPIKE)
    out["corr_A_B"] = float(np.corrcoef(LA[mB.index], mB)[0, 1])
    out["eff_donor"] = spike_effect(mB.to_numpy(), (LA[mB.index] > cut).to_numpy(), rng)
    effs = []
    for t in B:
        x = df.loc[LA.index, t].dropna()
        if len(x) >= 100:
            effs.append(spike_effect(x.to_numpy(), (LA[x.index] > cut).to_numpy(), rng))
    out["eff_tissue_median"] = float(np.nanmedian(effs)) if effs else np.nan
    out["n_B_tissues"] = len(effs)
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--prep", type=Path, required=True)
    p.add_argument("--out", type=Path, default=None)
    p.add_argument("--jobs", type=int, default=6)
    args = p.parse_args()

    by_gene, by_gene_raw = load(args.prep)
    jobs = [(g, by_gene[g], by_gene_raw[g]) for g in sorted(by_gene)]
    with ProcessPoolExecutor(args.jobs) as ex:
        res = pd.DataFrame(list(ex.map(one_gene, jobs, chunksize=20)))
    if args.out:
        res.to_csv(args.out, index=False, float_format="%.4g")

    pd.set_option("display.width", 200)
    print(f"genes: {len(res):,}")
    print("\n(a) coverage, median over genes:")
    print(res[["n_tissues", "donors_ge3", "donors_ge5", "donors_ge10"]]
          .describe(percentiles=[.25, .5, .75]).round(1).to_string())
    print("\n(b) ICC (donor-shared share of variance):")
    print(res[["icc_raw", "icc"]].describe(percentiles=[.1, .25, .5, .75, .9]).round(3).to_string())
    s = res.dropna(subset=["eff_donor", "eff_tissue_median"])
    print(f"\n(c) spike-in, {len(s):,} genes with >= {MIN_TISSUES} tissues "
          f"(top {SPIKE:.0%} of A-half removed; effect = UTC drop in SDs of random removal):")
    print(s[["corr_A_B", "eff_tissue_median", "eff_donor"]]
          .describe(percentiles=[.1, .25, .5, .75, .9]).round(3).to_string())
    print(f"\n  donor-level effect > tissue-level effect in "
          f"{(s.eff_donor > s.eff_tissue_median).mean():.1%} of genes")
    print("  detection (effect > 1.645 SD): tissue "
          f"{(s.eff_tissue_median > 1.645).mean():.1%}, donor {(s.eff_donor > 1.645).mean():.1%}")
    s = s.assign(icc_bin=pd.qcut(s.icc, 4))
    print(s.groupby("icc_bin", observed=True)[["corr_A_B", "eff_tissue_median", "eff_donor"]]
          .median().round(3).to_string())


if __name__ == "__main__":
    main()
