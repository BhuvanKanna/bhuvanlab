# -*- coding: utf-8 -*-
"""
compute_qc.py

Build a QC / distribution-classification table for one tissue matrix: one row
per gene, saying whether the gene's distribution is Gaussian, how it fails if it
is not, and whether its 4-parameter fit is even usable.

This deliberately does **not** touch `outputs/`. Those tables are mirrored into
4,665 gene_major shards, pinned by four hardcoded column lists, and carry a
byte-identity guarantee with extract_genes.py; adding columns there means
regenerating and re-pushing ~852 MB. A QC table joins on `gene` and costs
nothing.

The science is in normality.py -- this module only does I/O, parallelism and
the null calibration that every gene in the tissue shares.

Usage
-----
    # raw table for one tissue
    python compute_qc.py --input ../data/v11_log2_kidney_cortex.csv.gz \
        --id-col Name --name-col Description --jobs 8

    # excluded (<= -1) counterpart
    python compute_qc.py --input ../data/v11_log2_kidney_cortex.csv.gz \
        --id-col Name --name-col Description --threshold -1 --jobs 8

The input matrices are Git LFS pointers by default. Materialise one with:

    git -c lfs.fetchexclude= lfs fetch origin main --include=<path>
    git lfs checkout <path>

(plain `git lfs pull` is a no-op in this repo -- lfs.fetchexclude=* wins.)
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import normality as N
from generate_fourparam import (input_stem, load_expression, load_name_map,
                                output_csv_for)

warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
QC_DIR = HERE.parent / "qc"

COLUMNS = [
    "gene", "genename", "n_obs",
    # stage 0-1: what kind of object is this
    "frac_at_floor", "is_zero_inflated", "n_modes",
    # stage 2-3: is it normal, and how does it fail
    "shapiro_p", "shapiro_bh_reject", "dagostino_p", "anderson_stat",
    "skew", "kurt",
    # stage 4: truncated-Gaussian MLE, scored against the simulated null
    "mle_mu", "mle_sigma", "mle_sigma_dist", "d_aic", "d_aic_threshold",
    "mle_success",
    # fit validity -- these gate everything downstream
    "r_squared", "sigma_span_ratio", "sigma_span_ok", "x0_in_range",
    "n_obs_ok", "support_ok", "qc_pass",
    # stage 5
    "dist_class",
]


def qc_csv_for(input_path: Path, threshold) -> Path:
    stem = input_stem(input_path)
    if threshold is None:
        return QC_DIR / f"{stem}_qc.csv"
    return QC_DIR / f"{stem}_qc_excluded_at_or_below_{threshold:g}.csv"


def load_fit_params(input_path: Path, threshold) -> pd.DataFrame:
    """
    The `y0, A, x0, w, min, max` columns of the matching fourparam table, needed
    for R^2 and the sigma/span gate. Read as str then coerced, matching
    extract_genes.py -- these tables are the source of truth and nothing here
    should round-trip a float through pandas' writer.

    Returns an empty frame if the table has not been generated yet; the QC run
    still produces every distribution column, and only the fit-quality gates go
    NaN.
    """
    path = output_csv_for(input_path, threshold)
    if not path.exists():
        return pd.DataFrame()
    cols = ["gene", "y0", "A", "x0", "w", "min", "max"]
    df = pd.read_csv(path, usecols=cols, dtype=str, keep_default_na=False)
    for c in cols[1:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.set_index("gene")


# -- parallel worker ----------------------------------------------------------

_W_THRESHOLD = None
_W_PARAMS: dict = {}
_W_NULL: dict = {}


def _init_worker(threshold, params: dict, null: dict) -> None:
    global _W_THRESHOLD, _W_PARAMS, _W_NULL
    _W_THRESHOLD = threshold
    _W_PARAMS = params
    _W_NULL = null
    warnings.filterwarnings("ignore")


def _blank_row(gene: str, n_obs: int, extra: dict | None = None) -> dict:
    row = {c: np.nan for c in COLUMNS}
    row["gene"] = gene
    row["n_obs"] = int(n_obs)
    row["dist_class"] = "undetermined"
    for flag in ("is_zero_inflated", "mle_success", "sigma_span_ok",
                 "x0_in_range", "n_obs_ok", "support_ok", "qc_pass",
                 "shapiro_bh_reject"):
        row[flag] = False
    if extra:
        row.update(extra)
    return row


def _qc_one(item) -> dict:
    gene, values = item
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    if _W_THRESHOLD is not None:
        data = data[data > _W_THRESHOLD]
    n_obs = int(data.size)

    support = N.classify_support(data)
    row = _blank_row(gene, n_obs, {
        "frac_at_floor": support["frac_at_floor"],
        "is_zero_inflated": support["is_zero_inflated"],
    })
    if n_obs == 0:
        return row

    # A zero-inflated gene is not a failed Gaussian, it is a different object.
    # Short-circuiting keeps the expensive stages off the ~63% of genes for
    # which their answer would be meaningless anyway.
    if support["is_zero_inflated"]:
        row["dist_class"] = "zero_inflated"
        p = _W_PARAMS.get(gene)
        if p:
            row.update(_gate_fields(data, p, n_obs, support["frac_at_floor"]))
        return row

    tests = N.run_tests(data)
    row.update({k: tests[k] for k in
                ("shapiro_p", "dagostino_p", "anderson_stat", "skew", "kurt")})
    row["n_modes"] = N.count_modes(data)

    trunc = N.fit_truncated(data)
    row.update({k: trunc[k] for k in
                ("mle_mu", "mle_sigma", "mle_sigma_dist", "d_aic", "mle_success")})
    row["d_aic_threshold"] = _W_NULL.get("d_aic_threshold", np.nan)

    p = _W_PARAMS.get(gene)
    if p:
        row.update(_gate_fields(data, p, n_obs, support["frac_at_floor"]))

    row["dist_class"] = N.classify(support, int(row["n_modes"] or 0), tests,
                                   trunc, _W_NULL)
    return row


def _gate_fields(data, p: dict, n_obs: int, frac_at_floor: float) -> dict:
    out = {"r_squared": N.r_squared(data, p["y0"], p["A"], p["x0"], p["w"])}
    out.update(N.qc_gates(p["x0"], p["w"], p["min"], p["max"], n_obs, frac_at_floor))
    return out


def build_table(df: pd.DataFrame, threshold, params: pd.DataFrame, null: dict,
                jobs: int, limit=None) -> pd.DataFrame:
    genes = list(df.columns)
    if limit:
        genes = genes[:limit]

    # Ship the fit parameters to workers as a plain dict: a DataFrame is pickled
    # per chunk, a dict of small dicts is pickled once at pool start.
    pmap = ({} if params.empty
            else {g: r for g, r in params.to_dict(orient="index").items()})

    items = ((g, df[g].to_numpy(dtype=float)) for g in genes)
    if jobs and jobs != 1:
        n = jobs if jobs > 0 else (os.cpu_count() or 1)
        with mp.Pool(n, initializer=_init_worker,
                     initargs=(threshold, pmap, null)) as pool:
            records = list(pool.imap(_qc_one, items, chunksize=32))
    else:
        _init_worker(threshold, pmap, null)
        records = [_qc_one(it) for it in items]

    table = pd.DataFrame.from_records(records, columns=COLUMNS)

    # BH needs every gene's p-value at once, so it cannot live in the worker.
    table["shapiro_bh_reject"] = N.bh_fdr(table["shapiro_p"].to_numpy())
    return table


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path, required=True,
                   help="Tissue matrix, genes as rows (.csv or .csv.gz).")
    p.add_argument("--id-col", type=str, default="Name")
    p.add_argument("--name-col", type=str, default="Description")
    p.add_argument("--threshold", type=float, default=None,
                   help="Drop values <= this before testing (use -1 for the "
                        "excluded table). Omit for raw.")
    p.add_argument("--jobs", type=int, default=8)
    p.add_argument("--limit", type=int, default=None,
                   help="Only the first N genes (smoke tests).")
    p.add_argument("--null-sims", type=int, default=N.NULL_SIMS,
                   help="Simulated normals for the d_aic threshold.")
    p.add_argument("--force", action="store_true",
                   help="Recompute even if the output exists.")
    args = p.parse_args()

    out = qc_csv_for(args.input, args.threshold)
    if out.exists() and not args.force:
        print(f"exists, skipping: {out}   (--force to recompute)")
        return

    if args.input.stat().st_size < 1000:
        raise SystemExit(
            f"{args.input} is {args.input.stat().st_size} bytes -- that is a Git LFS\n"
            f"pointer, not the matrix. Materialise it with:\n"
            f"  git -c lfs.fetchexclude= lfs fetch origin main --include=\"{args.input}\"\n"
            f"  git lfs checkout \"{args.input}\"")

    df = load_expression(args.input, args.id_col, drop_cols=(args.name_col,))
    names = load_name_map(args.input, args.id_col, args.name_col)
    params = load_fit_params(args.input, args.threshold)
    if params.empty:
        print(f"note: {output_csv_for(args.input, args.threshold).name} not found; "
              f"fit-quality columns will be NaN")

    genes = list(df.columns)[: args.limit] if args.limit else list(df.columns)

    # One calibration for the whole tissue. Every gene shares a sample size and
    # d_aic is scale-free, so this is measured once rather than 74,628 times.
    n_typical = int(np.median([np.isfinite(df[g].to_numpy(dtype=float)).sum()
                               for g in genes[:200]]))
    print(f"calibrating null at n={n_typical} ({args.null_sims} sims)...")
    null = N.calibrate_null(n_typical, n_sims=args.null_sims)
    print(f"  d_aic threshold {null['d_aic_threshold']:.3f} "
          f"(median {null['d_aic_median']:.3f}); "
          f"skew noise band [{null['skew_lo']:.3f}, {null['skew_hi']:.3f}]")

    table = build_table(df, args.threshold, params, null, args.jobs, args.limit)
    # genename is already a COLUMNS slot, so this fills it rather than inserting.
    table["genename"] = table["gene"].map(names)

    QC_DIR.mkdir(parents=True, exist_ok=True)
    table.to_csv(out, index=False, lineterminator="\n")
    print(f"wrote {out}  ({len(table):,} genes)")
    counts = table["dist_class"].value_counts()
    for k, v in counts.items():
        print(f"  {k:<16} {v:>7,}  {v/len(table):6.1%}")
    print(f"  {'qc_pass':<16} {int(table.qc_pass.sum()):>7,}  "
          f"{table.qc_pass.mean():6.1%}")


if __name__ == "__main__":
    main()
