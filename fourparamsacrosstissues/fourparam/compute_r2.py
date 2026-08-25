# -*- coding: utf-8 -*-
"""
compute_r2.py

Add `r_squared` to a fourparam table without refitting anything.

`sumsquarevalue` in `outputs/` is an unnormalised residual sum of squares: it
scales with n and with peak height, so it is useless for comparing one gene to
another or one tissue to another. Measured across the 54 excluded tables, the
Spearman correlation between a tissue's median `sumsquarevalue` and its median
`n_obs` is **0.99** -- the column is very nearly a restatement of donor count.
R^2 divides that by the total sum of squares and is comparable everywhere.

**Nothing is refit.** R^2 = 1 - SSR/TSS, and SSR is already in the table as
`sumsquarevalue`; only TSS = sum((count - mean_count)^2) over the 40 bins was
never stored, and it is not recoverable from the published columns (`hist` is
quantised to 1/63 of peak and exists in only 2 of the 54 excluded tables). So
this re-bins each gene from the source matrix -- the cheap half of the pipeline
-- and reuses the stored `y0, A, x0, w` as-is.

That the re-binning reproduces the original is checkable, and `--verify` checks
it: recomputing SSR from the fresh histogram and comparing to the stored
`sumsquarevalue` agrees to ~1e-14 relative, with no `n_obs` disagreement. If a
tissue ever fails that, its histogram is not the one that was fit and its R^2
would be quietly wrong.

The science is in `normality.r_squared`, which the QC tables already use -- this
module only does I/O. Importing it rather than reimplementing the formula is
what keeps `r2/` and `qc/` from ever disagreeing.

Usage
-----
    cd fourparam
    python compute_r2.py --input ../data/v11_log2_kidney_cortex.csv.gz
    python compute_r2.py --input ../data/v11_log2_kidney_cortex.csv.gz --verify
    python compute_r2.py --all                    # every materialised matrix
    python compute_r2.py --all --raw              # the raw tables instead

Writes `r2/<stem>_r2_excluded_at_or_below_-1.csv` (gene, r_squared) at full
float precision, for analysis. `build_r2.py` turns that into the fixed-width
string the browser reads. Skips existing output unless `--force`.

The input matrices are Git LFS pointers by default. Materialise one with:

    git -c lfs.fetchexclude= lfs fetch origin main --include=<path>
    git lfs checkout <path>

(plain `git lfs pull` is a no-op in this repo -- lfs.fetchexclude=* wins.)
"""
from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import normality as N
from bhuvanfitter import _fourparam_gaussian
from generate_fourparam import input_stem, output_csv_for

warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
DATA_DIR = HERE.parent / "data"
R2_DIR = HERE.parent / "r2"

# Rows are parsed in blocks so a 800-donor matrix never lands in memory whole,
# and so the float conversion happens in pandas' C parser rather than per row.
CHUNK_ROWS = 2000

# Below this a file in data/ is a Git LFS pointer, not a matrix.
POINTER_MAX_BYTES = 1000


def r2_csv_for(input_path: Path, threshold) -> Path:
    stem = input_stem(input_path)
    if threshold is None:
        return R2_DIR / f"{stem}_r2.csv"
    return R2_DIR / f"{stem}_r2_excluded_at_or_below_{threshold:g}.csv"


def load_fit_params(input_path: Path, threshold) -> pd.DataFrame:
    """
    The `y0, A, x0, w, n_obs, fit_success, sumsquarevalue` columns of the
    matching fourparam table. Read as str then coerced, matching
    extract_genes.py -- these tables are the source of truth and nothing here
    should round-trip a float through pandas' writer.
    """
    path = output_csv_for(input_path, threshold)
    if not path.exists():
        raise SystemExit(f"missing fourparam table {path} -- generate it first")
    cols = ["gene", "y0", "A", "x0", "w", "n_obs",
            "fit_success", "sumsquarevalue"]
    df = pd.read_csv(path, usecols=cols, dtype=str, keep_default_na=False)
    for c in ("y0", "A", "x0", "w", "n_obs", "sumsquarevalue"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["fit_success"] = df["fit_success"].str.strip().str.lower() == "true"
    return df.set_index("gene")


def _values_for(row: np.ndarray, threshold) -> np.ndarray:
    """The exact array the fit saw: finite values, then the exclusion cut."""
    arr = row[np.isfinite(row)]
    if threshold is not None:
        arr = arr[arr > threshold]
    return arr


def compute(input_path: Path, threshold, id_col: str, verify: bool) -> pd.DataFrame:
    params = load_fit_params(input_path, threshold)
    genes, values = [], []
    checks: list[float] = []
    n_obs_mismatch = 0

    reader = pd.read_csv(input_path, chunksize=CHUNK_ROWS)
    for chunk in reader:
        ids = chunk[id_col].to_numpy()
        sample_cols = [c for c in chunk.columns if _is_sample(c, id_col)]
        mat = chunk[sample_cols].to_numpy(dtype=float, na_value=np.nan)

        for gid, raw in zip(ids, mat):
            if gid not in params.index:
                continue
            p = params.loc[gid]
            genes.append(gid)

            if not bool(p["fit_success"]):
                values.append(np.nan)
                continue

            arr = _values_for(raw, threshold)
            if arr.size != p["n_obs"]:
                n_obs_mismatch += 1
            values.append(N.r_squared(arr, p["y0"], p["A"], p["x0"], p["w"]))

            if verify and arr.size >= 10:
                counts, edges = np.histogram(arr, bins=N.HIST_BINS)
                centres = 0.5 * (edges[:-1] + edges[1:])
                pred = _fourparam_gaussian(centres, p["y0"], p["A"],
                                           p["x0"], p["w"])
                if np.all(np.isfinite(pred)):
                    ssr = float(np.sum((counts.astype(float) - pred) ** 2))
                    ref = float(p["sumsquarevalue"])
                    checks.append(abs(ssr - ref) / max(abs(ref), 1e-12))

    if verify:
        _report_verification(checks, n_obs_mismatch)
    elif n_obs_mismatch:
        print(f"  WARNING: n_obs disagrees on {n_obs_mismatch:,} genes")

    return pd.DataFrame({"gene": genes, "r_squared": values})


def _is_sample(col: str, id_col: str) -> bool:
    """Sample columns are everything but the id and the gene-name column."""
    return col not in (id_col, "Description", "genename")


def _report_verification(checks: list[float], n_obs_mismatch: int) -> None:
    """
    The re-binning is only sound if it reproduces the histogram that was fit.
    Comparing recomputed SSR against the stored `sumsquarevalue` is a direct
    test of that, and a loud failure is much better than a quietly wrong column.
    """
    if not checks:
        print("  verify: no comparable genes")
        return
    err = np.array(checks)
    worst = float(err.max())
    print(f"  verify: SSR max rel err {worst:.2e}, median {np.median(err):.2e}, "
          f"n_obs mismatches {n_obs_mismatch:,}")
    if worst > 1e-6 or n_obs_mismatch:
        print("  WARNING: recomputed histogram does not match the stored fit; "
              "R^2 for this tissue is not trustworthy")


def run_one(input_path: Path, threshold, id_col: str, force: bool,
            verify: bool) -> bool:
    out = r2_csv_for(input_path, threshold)
    if out.exists() and not force:
        print(f"{input_path.name}: exists, skipping ({out.name})")
        return False
    if input_path.stat().st_size < POINTER_MAX_BYTES:
        print(f"{input_path.name}: Git LFS pointer, not a matrix -- skipping")
        return False

    print(f"{input_path.name}:")
    df = compute(input_path, threshold, id_col, verify)
    R2_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    ok = df["r_squared"].notna()
    print(f"  wrote {out.name}: {len(df):,} genes, {int(ok.sum()):,} with R^2, "
          f"median {df.loc[ok, 'r_squared'].median():.3f}")
    return True


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path, default=None,
                   help="One matrix in data/. Omit with --all.")
    p.add_argument("--all", action="store_true",
                   help="Every materialised matrix in data/ (pointers skipped).")
    p.add_argument("--id-col", type=str, default="Name")
    p.add_argument("--raw", action="store_true",
                   help="Do the raw table. Default is the excluded (<= -1) one.")
    p.add_argument("--force", action="store_true",
                   help="Recompute even if the output exists.")
    p.add_argument("--verify", action="store_true",
                   help="Also check recomputed SSR against sumsquarevalue.")
    p.add_argument("--jobs", type=int, default=1,
                   help="Tissues to process in parallel. Each worker streams its "
                        "own matrix in chunks, so memory is per-worker and small; "
                        "the ceiling is disk throughput, not RAM.")
    args = p.parse_args()

    threshold = None if args.raw else -1.0

    if args.all:
        inputs = sorted(DATA_DIR.glob("v11_log2_*.csv.gz"))
    elif args.input:
        inputs = [args.input]
    else:
        raise SystemExit("pass --input <matrix> or --all")

    if args.jobs > 1 and len(inputs) > 1:
        done = _run_parallel(inputs, threshold, args)
    else:
        done = sum(bool(run_one(path, threshold, args.id_col, args.force,
                                args.verify))
                   for path in inputs)
    print(f"\n{done} table(s) written to {R2_DIR}")


def _run_parallel(inputs, threshold, args) -> int:
    """
    One process per tissue. Output is buffered per tissue and printed whole, so
    two workers finishing at once cannot interleave their lines -- which would
    make a verification warning look like it belonged to the wrong tissue.
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

    done = 0
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(_run_one_captured, path, threshold, args.id_col,
                        args.force, args.verify): path
            for path in inputs
        }
        for i, fut in enumerate(as_completed(futures), 1):
            path = futures[fut]
            try:
                wrote, text = fut.result()
            except Exception as exc:                      # noqa: BLE001
                print(f"[{i}/{len(inputs)}] {path.name}: FAILED -- {exc}")
                continue
            print(f"[{i}/{len(inputs)}] {text.rstrip()}")
            done += bool(wrote)
    return done


def _run_one_captured(path, threshold, id_col, force, verify):
    """run_one with its stdout captured, for orderly parallel logging."""
    import contextlib
    import io

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        wrote = run_one(path, threshold, id_col, force, verify)
    return wrote, buf.getvalue()


if __name__ == "__main__":
    main()
