#!/usr/bin/env python
"""
generate_fourparam_auto_bins.py

Refit ONE tissue's fourparam table under Igor Pro's three "auto bin size" rules
(Dr. Pierce, 2026-10) and compare each against the fixed 40 bins:

    sturges   bins  = ceil(1 + log2(n))
    scott     width = 3.49 * sd  * n^(-1/3),  bins = ceil((max - min) / width)
    fd        width = 2    * IQR * n^(-1/3),  bins = ceil((max - min) / width)

Scott and Freedman-Diaconis (fd) set a bin WIDTH from each gene's own spread,
so the bin count is per gene, not per tissue; Sturges depends only on n, so it
varies per gene only through n_obs after the <= threshold exclusion. n is the
gene's n_obs (the values actually histogrammed), not the tissue's sample count.

Two fixed arms ride along as references: 40 (the production tables) and the
earlier linear rule round(0.073 * n_samples + 3) from
generate_fourparam_adaptive_bins.py. All arms are refit in this run, so the
comparison is like-for-like. Everything else is generate_fourparam.py
unchanged: same fit, same filters, same MIN_OBS.

    python generate_fourparam_auto_bins.py \
        --input ../data/v11_log2_brain_spinal_cord_cervical_c-1.csv.gz --threshold -1

Writes to results/auto_binning/:
    <stem>_fourparam_<arm>bins<suffix>.csv   one fourparam table per arm, + n_bins, r_squared
    <stem>_auto_bin_counts<suffix>.csv       per-gene bin count under every rule
    <stem>_auto_bin_comparison<suffix>.csv   one summary row per arm, vs 40 bins

Resumable: an arm whose table is already written is read back, not refit.
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import generate_fourparam as G
import normality as N
from bhuvanfitter import BhuvanFitter
from generate_fourparam_adaptive_bins import bins_for as linear_bins_for

OUT_DIR = Path(__file__).resolve().parent.parent / "results" / "auto_binning"

RULES = ("sturges", "scott", "fd")
MIN_BINS = 5        # curve_fit needs more bins than the 4 parameters
COMPARE_MIN_OBS = 30  # degenerate low-n fits are unstable even at fixed bins


def auto_bins(data: np.ndarray, rule: str) -> int:
    """Bin count under one Igor Pro auto rule; 0 when the rule is undefined."""
    n = data.size
    if n < 2:
        return 0
    if rule == "sturges":
        return int(np.ceil(1 + np.log2(n)))
    span = float(data.max() - data.min())
    if rule == "scott":
        width = 3.49 * float(np.std(data, ddof=1)) * n ** (-1 / 3)
    elif rule == "fd":
        q75, q25 = np.percentile(data, [75, 25])
        width = 2 * float(q75 - q25) * n ** (-1 / 3)
    else:
        raise ValueError(rule)
    # A zero IQR (ties) or zero spread leaves the width undefined.
    if width <= 0 or span <= 0:
        return 0
    return max(1, int(np.ceil(span / width)))


_THRESHOLD = None
_MAX_NFEV = 2000


def _init_worker(threshold, max_nfev) -> None:
    global _THRESHOLD, _MAX_NFEV
    _THRESHOLD = threshold
    _MAX_NFEV = max_nfev
    warnings.filterwarnings("ignore")


def _fit_one(item) -> dict:
    gene, values, n_bins = item
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    if _THRESHOLD is not None:
        data = data[data > _THRESHOLD]
    n_obs = int(data.size)
    lo = float(data.min()) if n_obs else np.nan
    hi = float(data.max()) if n_obs else np.nan
    if n_obs < G.MIN_OBS or n_bins < MIN_BINS:
        row = G._failed_row(gene, n_obs, lo=lo, hi=hi)
    else:
        # BINS is read in __init__, and each worker fits one gene at a time.
        BhuvanFitter.BINS = n_bins
        try:
            row = BhuvanFitter(data, gene_name=gene).fit("fourparam", max_nfev=_MAX_NFEV)
        except (RuntimeError, ValueError):
            row = G._failed_row(gene, n_obs, lo=lo, hi=hi)
    row["n_bins"] = int(n_bins)
    row["r_squared"] = (N.r_squared(data, row["y0"], row["A"], row["x0"], row["w"], bins=n_bins)
                        if row["fit_success"] else np.nan)
    return row


def _side(table: pd.DataFrame) -> pd.Series:
    """Which tail carries the truncation: 'R', 'L', or '' (neither / no fit)."""
    return np.where(table["rti"] > 0, "R", np.where(table["lti"] > 0, "L", ""))


def compare(arms: dict[str, pd.DataFrame], ref: str = "40") -> pd.DataFrame:
    base = arms[ref]
    rows = []
    for arm, t in arms.items():
        ok = t["fit_success"].astype(bool)
        both = ok & base["fit_success"].astype(bool) & (t["n_obs"] >= COMPARE_MIN_OBS)
        fitted_bins = t.loc[ok, "n_bins"]
        row = {
            "arm": arm,
            "fit_success": int(ok.sum()),
            "bins_median": fitted_bins.median(),
            "bins_p5": fitted_bins.quantile(0.05),
            "bins_p95": fitted_bins.quantile(0.95),
            "bins_min": fitted_bins.min(),
            "bins_max": fitted_bins.max(),
            "r_squared_median": t.loc[ok, "r_squared"].median(),
            "n_compared": int(both.sum()),
            "frac_rti_gt0": float((t.loc[both, "rti"] > 0).mean()),
            "lti_median": t.loc[both, "lti"].median(),
            "rti_median_where_gt0": t.loc[both & (t["rti"] > 0), "rti"].median(),
        }
        for col in ("rti", "lti", "rti_sigma_dist", "lti_sigma_dist", "x0", "w"):
            row[f"spearman_{col}_vs_{ref}"] = spearmanr(t.loc[both, col], base.loc[both, col],
                                                        nan_policy="omit").statistic
        row[f"frac_side_flip_vs_{ref}"] = float(
            (_side(t.loc[both]) != _side(base.loc[both])).mean())
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--threshold", type=float, default=None)
    p.add_argument("--jobs", type=int, default=8)
    p.add_argument("--max-nfev", type=int, default=2000, dest="max_nfev")
    args = p.parse_args()

    df = G.load_expression(args.input, id_col="Name", drop_cols=["Description"])
    n_samples = df.shape[0]
    linear = linear_bins_for(n_samples)
    tissue_sturges = int(np.ceil(1 + np.log2(n_samples)))
    print(f"{df.shape[1]} genes x {n_samples} samples; tissue-level Sturges = "
          f"{tissue_sturges}, linear rule = {linear}", flush=True)

    counts = {"gene": [], "n_obs": []} | {r: [] for r in RULES}
    for g in df.columns:
        arr = df[g].to_numpy(dtype=float)
        arr = arr[np.isfinite(arr)]
        if args.threshold is not None:
            arr = arr[arr > args.threshold]
        counts["gene"].append(g)
        counts["n_obs"].append(arr.size)
        for r in RULES:
            counts[r].append(auto_bins(arr, r))
    bin_counts = pd.DataFrame(counts)

    plans = {"40": np.full(len(df.columns), 40), str(linear): np.full(len(df.columns), linear)}
    plans |= {r: bin_counts[r].to_numpy() for r in RULES}

    name_map = G.load_name_map(args.input, "Name", "Description")
    stem = G.input_stem(args.input)
    suffix = "" if args.threshold is None else f"_excluded_at_or_below_{args.threshold:g}"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    cols = [c for c in G.COLUMNS if c not in ("hist", "hist_max", "r_squared")]
    cols.insert(cols.index("n_obs"), "n_bins")
    cols.append("r_squared")
    arms = {}
    with mp.Pool(args.jobs, initializer=_init_worker,
                 initargs=(args.threshold, args.max_nfev)) as pool:
        for arm, bins in plans.items():
            label = arm if arm in RULES else f"{arm}bins"
            out = OUT_DIR / f"{stem}_fourparam_{label}{suffix}.csv"
            # Resumable: each arm is ~20 min, so reuse one already written.
            if out.exists():
                arms[arm] = pd.read_csv(out).drop(columns=["genename"])[cols]
                print(f"Reusing {out}", flush=True)
                continue
            items = ((g, df[g].to_numpy(dtype=float), int(b)) for g, b in zip(df.columns, bins))
            records = list(pool.imap(_fit_one, items, chunksize=64))
            table = pd.DataFrame.from_records(records)[cols]
            arms[arm] = table
            G.insert_genename(table, name_map).to_csv(out, index=False, lineterminator="\n")
            print(f"Wrote {out} ({int(table['fit_success'].sum())} fit_success)", flush=True)

    bin_counts = G.insert_genename(bin_counts, name_map)
    bin_counts.to_csv(OUT_DIR / f"{stem}_auto_bin_counts{suffix}.csv", index=False, lineterminator="\n")
    summary = compare(arms)
    summary.insert(1, "n_samples", n_samples)
    summary.to_csv(OUT_DIR / f"{stem}_auto_bin_comparison{suffix}.csv", index=False, lineterminator="\n")
    with pd.option_context("display.width", 250, "display.max_columns", None):
        print(summary.T.to_string())


if __name__ == "__main__":
    main()
