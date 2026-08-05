#!/usr/bin/env python
"""
generate_fourparam.py

Generate a **4-parameter Gaussian fourparam table** from a gene-expression
matrix. One row per gene, with the full metric schema (see COLUMNS below).
Self-contained: the only local dependency is ``bhuvanfitter.py`` (the fit
library) sitting next to this file.

Input matrix
------------
A CSV (optionally gzipped, ``.csv.gz``) laid out **genes as rows**:

    <id-col>,<name-col>,<sample_1>,<sample_2>, ...
    ENSG...,SYMBOL,12.3,0.0, ...

``--id-col`` names the gene-identifier column (default ``Name``), ``--name-col``
names the common-gene-name / symbol column (default ``Description``). Every other
column is treated as a sample. The matrix is transposed internally to
samples-as-rows so each gene's values become one column to fit.

Two table types (the difference is one filter)
----------------------------------------------
* **raw** (no ``--threshold``): every finite expression value is kept.
* **excluded** (``--threshold T``, typically ``-1``): every value **at or below
  T** is dropped from a gene's array *before* fitting, so ``n_obs`` shrinks and a
  gene can fall below ``MIN_OBS``. With this repo's ``log2(TPM+1)-1`` data,
  ``T = -1`` drops exactly the zero-TPM (undetected) samples.

Everything else — the fit, the columns, ``MIN_OBS``, the failure handling — is
**identical** between the two; only this pre-fit value filter differs.

Default filters applied to every table (both types), per gene, in order
-----------------------------------------------------------------------
1. Non-finite values (NaN / inf) are dropped.
2. (excluded tables only) values ``<= threshold`` are dropped.
3. If the surviving ``n_obs < MIN_OBS`` (=10), the gene is written as a
   ``fit_success=False`` row with all metrics NaN (it is never skipped).
4. If ``curve_fit`` fails to converge, same ``fit_success=False`` NaN row.

Output
------
Written to ``../outputs/`` (a sibling of this ``fourparam/`` folder), filename
derived from the input stem:

    raw       -> outputs/<stem>_fourparam.csv
    excluded  -> outputs/<stem>_fourparam_excluded_at_or_below_<T>.csv

A ``genename`` column (from ``--name-col``) is inserted immediately after
``gene``. Parallelism (``--jobs``) is order-preserving, so a parallel run is
bit-identical to ``--jobs 1``.

Usage
-----
    # raw table for one tissue
    python generate_fourparam.py --input ../data/v11_log2_liver.csv.gz \
        --id-col Name --name-col Description --jobs 8 --max-nfev 2000

    # excluded (<= -1) table for the same tissue
    python generate_fourparam.py --input ../data/v11_log2_liver.csv.gz \
        --id-col Name --name-col Description --threshold -1 --jobs 8 --max-nfev 2000
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from bhuvanfitter import BhuvanFitter, encode_histogram

warnings.filterwarnings("ignore")  # silence skew/kurt precision-loss noise on near-constant genes

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE.parent / "outputs"

# Minimum finite observations for a gene to be fit (else fit_success=False row).
MIN_OBS = 10

# Column order == the keys BhuvanFitter.fit("fourparam") returns, exactly.
COLUMNS = [
    "gene", "y0", "A", "x0", "w", "sumsquarevalue",
    "ti_fourparam_sigma_dist", "truncationindex",
    "min", "max", "mean", "std", "skew", "kurt",
    "right", "maxheight", "rightheight", "n_obs", "fit_success",
    # Thumbnail histogram for the browser's Shape column. Quantised and lossy —
    # a rendering aid, not analysis data. See bhuvanfitter.encode_histogram.
    "hist", "hist_max",
]


def input_stem(input_path: Path) -> str:
    """`v11_log2_liver.csv.gz` -> `v11_log2_liver` (strip .csv / .csv.gz)."""
    name = input_path.name
    for suffix in (".csv.gz", ".csv", ".gz"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return input_path.stem


def output_csv_for(input_path: Path, threshold) -> Path:
    stem = input_stem(input_path)
    if threshold is None:
        return OUT_DIR / f"{stem}_fourparam.csv"
    return OUT_DIR / f"{stem}_fourparam_excluded_at_or_below_{threshold:g}.csv"


def load_expression(csv_path: Path, id_col: str, drop_cols=()) -> pd.DataFrame:
    """Read a genes-as-rows matrix and return samples-as-rows (genes as columns).
    `id_col` becomes the column index; `drop_cols` (e.g. the name column) are
    discarded so they aren't treated as samples. Reads .csv and .csv.gz."""
    df = pd.read_csv(csv_path)
    df = df.set_index(id_col)
    drop = [c for c in drop_cols if c in df.columns]
    if drop:
        df = df.drop(columns=drop)
    return df.T  # rows = samples, columns = genes


def load_name_map(csv_path: Path, id_col: str, name_col: str) -> pd.Series:
    m = pd.read_csv(csv_path, usecols=[id_col, name_col])
    return m.drop_duplicates(id_col).set_index(id_col)[name_col]


def insert_genename(table: pd.DataFrame, name_map: pd.Series) -> pd.DataFrame:
    table = table.drop(columns=["genename"], errors="ignore")
    table.insert(1, "genename", table["gene"].map(name_map))
    return table


def _failed_row(gene: str, n_obs: int, hist: str = "", hist_max: int = 0) -> dict:
    row = {c: np.nan for c in COLUMNS}
    row["gene"] = gene
    row["n_obs"] = int(n_obs)
    row["fit_success"] = False
    # 0 rather than NaN: one NaN anywhere in the column makes pandas write
    # "12.0" instead of "12" for every other row.
    row["hist"] = hist
    row["hist_max"] = int(hist_max)
    return row


# -- parallel worker ----------------------------------------------------------

_WORKER_THRESHOLD = None       # None = raw; else drop values <= threshold
_WORKER_MAX_NFEV = 10_000


def _init_worker(threshold, max_nfev) -> None:
    global _WORKER_THRESHOLD, _WORKER_MAX_NFEV
    _WORKER_THRESHOLD = threshold
    _WORKER_MAX_NFEV = max_nfev
    warnings.filterwarnings("ignore")


def _fit_one(item) -> dict:
    gene, values = item
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    if _WORKER_THRESHOLD is not None:
        data = data[data > _WORKER_THRESHOLD]   # drop values <= threshold
    n_obs = int(data.size)

    # Computed before and independently of the fit: a gene whose fit failed
    # still gets to show its real distribution in the browser.
    hist, hist_max = encode_histogram(data)

    if n_obs < MIN_OBS:
        return _failed_row(gene, n_obs, hist, hist_max)
    try:
        bf = BhuvanFitter(data, gene_name=gene)
        row = bf.fit("fourparam", max_nfev=_WORKER_MAX_NFEV)
    except (RuntimeError, ValueError):
        return _failed_row(gene, n_obs, hist, hist_max)
    row["hist"] = hist
    row["hist_max"] = hist_max
    return row


def build_table(df: pd.DataFrame, threshold, jobs: int, max_nfev: int,
                limit=None) -> pd.DataFrame:
    genes = list(df.columns)
    if limit:
        genes = genes[:limit]
    items = ((g, df[g].to_numpy(dtype=float)) for g in genes)

    if jobs and jobs != 1:
        n = jobs if jobs > 0 else (os.cpu_count() or 1)
        with mp.Pool(n, initializer=_init_worker,
                     initargs=(threshold, max_nfev)) as pool:
            records = list(pool.imap(_fit_one, items, chunksize=64))
    else:
        _init_worker(threshold, max_nfev)
        records = [_fit_one(it) for it in items]

    return pd.DataFrame.from_records(records, columns=COLUMNS)


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path, required=True,
                   help="Gene matrix CSV/CSV.GZ (genes as rows).")
    p.add_argument("--id-col", type=str, default="Name",
                   help="Gene-identifier column (default: Name).")
    p.add_argument("--name-col", type=str, default="Description",
                   help="Common-gene-name/symbol column, inserted as `genename` "
                        "(default: Description). Also dropped from the samples.")
    p.add_argument("--threshold", type=float, default=None,
                   help="If set, exclude values <= this before fitting "
                        "(excluded table). Omit for the raw table.")
    p.add_argument("--jobs", type=int, default=8, help="Parallel worker processes.")
    p.add_argument("--max-nfev", type=int, default=2000, dest="max_nfev",
                   help="curve_fit evaluation cap (default 2000).")
    p.add_argument("--limit", type=int, default=None,
                   help="Fit only the first N genes (sanity checks).")
    args = p.parse_args()

    drop_cols = [args.name_col] if args.name_col else []
    print(f"Reading {args.input.name} ...", flush=True)
    df = load_expression(args.input, id_col=args.id_col, drop_cols=drop_cols)
    print(f"  {df.shape[1]} genes x {df.shape[0]} samples", flush=True)

    mode = "raw" if args.threshold is None else f"excluded <= {args.threshold:g}"
    print(f"Fitting ({mode}, jobs={args.jobs}, max_nfev={args.max_nfev}) ...",
          flush=True)
    table = build_table(df, args.threshold, args.jobs, args.max_nfev, args.limit)

    if args.name_col:
        table = insert_genename(
            table, load_name_map(args.input, args.id_col, args.name_col))

    out = output_csv_for(args.input, args.threshold)
    out.parent.mkdir(parents=True, exist_ok=True)
    # Pin LF so every table in outputs/ shares one line terminator regardless of
    # platform; build_gene_major.py copies row text verbatim.
    table.to_csv(out, index=False, lineterminator="\n")
    n_ok = int(table["fit_success"].sum())
    print(f"Wrote {out.relative_to(HERE.parent)} "
          f"({len(table)} rows, {n_ok} fit_success)", flush=True)


if __name__ == "__main__":
    main()
