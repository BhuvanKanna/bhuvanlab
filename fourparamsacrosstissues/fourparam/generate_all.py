#!/usr/bin/env python
"""
generate_all.py

Driver that generates **both** fourparam tables (raw and excluded <= -1) for
**every** converted tissue matrix in ``../data/``, writing all of them to
``../outputs/``. With 54 tissue files this produces 108 tables.

For each ``../data/v11_log2_<tissue>.csv.gz`` it runs ``generate_fourparam.py``
twice:

    raw       -> outputs/v11_log2_<tissue>_fourparam.csv
    excluded  -> outputs/v11_log2_<tissue>_fourparam_excluded_at_or_below_-1.csv

Already-present output tables are **skipped** (so an interrupted run resumes by
re-invoking the same command). Genes are fit with ``--id-col Name
--name-col Description`` (a ``genename`` column is inserted after ``gene``).

Usage
-----
    python generate_all.py                 # all tissues, both tables
    python generate_all.py --jobs 8        # tune parallelism
    python generate_all.py --limit 50      # sanity: first 50 genes per table
"""
from __future__ import annotations

import argparse
from pathlib import Path

from generate_fourparam import (
    load_expression, load_name_map, insert_genename, build_table,
    output_csv_for, OUT_DIR,
)

HERE = Path(__file__).resolve().parent
DATA_DIR = HERE.parent / "data"
ID_COL, NAME_COL = "Name", "Description"
THRESHOLDS = [None, -1.0]   # raw, then excluded <= -1


def run_one(input_path: Path, threshold, jobs: int, max_nfev: int, limit) -> None:
    out = output_csv_for(input_path, threshold)
    if out.exists():
        print(f"  skip (exists): {out.name}", flush=True)
        return
    df = load_expression(input_path, id_col=ID_COL, drop_cols=[NAME_COL])
    table = build_table(df, threshold, jobs, max_nfev, limit)
    table = insert_genename(table, load_name_map(input_path, ID_COL, NAME_COL))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    table.to_csv(out, index=False)
    n_ok = int(table["fit_success"].sum())
    print(f"  wrote {out.name} ({len(table)} rows, {n_ok} fit_success)", flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--jobs", type=int, default=8)
    p.add_argument("--max-nfev", type=int, default=2000, dest="max_nfev")
    p.add_argument("--limit", type=int, default=None)
    args = p.parse_args()

    tissues = sorted(DATA_DIR.glob("v11_log2_*.csv.gz"))
    if not tissues:
        raise SystemExit(f"No tissue matrices found in {DATA_DIR} "
                         f"(expected v11_log2_*.csv.gz).")
    print(f"{len(tissues)} tissue matrices -> {2 * len(tissues)} tables", flush=True)

    for i, tissue in enumerate(tissues, 1):
        print(f"[{i}/{len(tissues)}] {tissue.name}", flush=True)
        for threshold in THRESHOLDS:
            run_one(tissue, threshold, args.jobs, args.max_nfev, args.limit)

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
