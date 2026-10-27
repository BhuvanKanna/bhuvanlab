#!/usr/bin/env python
"""
generate_fourparam_adaptive_bins.py

Refit ONE tissue's fourparam table with a sample-size-dependent bin count
instead of the fixed 40 (Dr. Pierce, 2026-10): tissues with >= 500 samples keep
40 bins; below that,

    bins = round(M * n_samples + B),   M = 0.073, B = 3

where n_samples is the tissue's donor-sample count (columns of the matrix), not
a gene's n_obs. Everything else is generate_fourparam.py unchanged: same fit,
same filters, same MIN_OBS. Output drops the browser-only hist/hist_max columns
and adds n_bins and r_squared (computed on the same adaptive binning).

    python generate_fourparam_adaptive_bins.py \
        --input ../data/v11_log2_brain_spinal_cord_cervical_c-1.csv.gz --threshold -1
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
from pathlib import Path

import numpy as np
import pandas as pd

import generate_fourparam as G
import normality as N
from bhuvanfitter import BhuvanFitter

OUT_DIR = Path(__file__).resolve().parent.parent / "results" / "adaptive_bins"

M, B = 0.073, 3
FULL_BINS, FULL_BINS_AT = 40, 500


def bins_for(n_samples: int) -> int:
    if n_samples >= FULL_BINS_AT:
        return FULL_BINS
    return int(round(M * n_samples + B))


def _init_worker(n_bins, threshold, max_nfev) -> None:
    BhuvanFitter.BINS = n_bins
    G._init_worker(threshold, max_nfev)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--threshold", type=float, default=None)
    p.add_argument("--jobs", type=int, default=8)
    p.add_argument("--max-nfev", type=int, default=2000, dest="max_nfev")
    p.add_argument("--bins", type=int, default=None,
                   help="Override the formula with a fixed bin count.")
    args = p.parse_args()

    df = G.load_expression(args.input, id_col="Name", drop_cols=["Description"])
    n_samples = df.shape[0]
    n_bins = args.bins or bins_for(n_samples)
    print(f"{df.shape[1]} genes x {n_samples} samples -> {n_bins} bins", flush=True)

    items = ((g, df[g].to_numpy(dtype=float)) for g in df.columns)
    # Spawned workers re-import bhuvanfitter, so BINS must be set in each one.
    with mp.Pool(args.jobs, initializer=_init_worker,
                 initargs=(n_bins, args.threshold, args.max_nfev)) as pool:
        records = list(pool.imap(G._fit_one, items, chunksize=64))
    table = pd.DataFrame.from_records(records, columns=G.COLUMNS)
    table = table.drop(columns=["hist", "hist_max"])

    r2 = []
    for gene, ok, y0, A, x0, w in table[["gene", "fit_success", "y0", "A", "x0", "w"]].itertuples(index=False):
        if not ok:
            r2.append(np.nan)
            continue
        arr = df[gene].to_numpy(dtype=float)
        arr = arr[np.isfinite(arr)]
        if args.threshold is not None:
            arr = arr[arr > args.threshold]
        r2.append(N.r_squared(arr, y0, A, x0, w, bins=n_bins))
    table["r_squared"] = r2
    table.insert(table.columns.get_loc("n_obs"), "n_bins", n_bins)

    table = G.insert_genename(table, G.load_name_map(args.input, "Name", "Description"))
    stem = G.input_stem(args.input)
    suffix = "" if args.threshold is None else f"_excluded_at_or_below_{args.threshold:g}"
    out = OUT_DIR / f"{stem}_fourparam_{n_bins}bins{suffix}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out, index=False, lineterminator="\n")
    print(f"Wrote {out} ({len(table)} rows, {int(table['fit_success'].sum())} fit_success)")


if __name__ == "__main__":
    main()
