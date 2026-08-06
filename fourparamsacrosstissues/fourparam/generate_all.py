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
    output_csv_for, OUT_DIR, COLUMNS,
)

HERE = Path(__file__).resolve().parent
DATA_DIR = HERE.parent / "data"
ID_COL, NAME_COL = "Name", "Description"
THRESHOLDS = [None, -1.0]   # raw, then excluded <= -1


CURRENT_HEADER = ["gene", "genename"] + COLUMNS[1:]


def header_is_current(path: Path) -> bool:
    """True if this table's columns match what the generator writes today."""
    with path.open("r", encoding="utf-8") as fh:
        return fh.readline().rstrip("\r\n").split(",") == CURRENT_HEADER


def run_one(input_path: Path, threshold, jobs: int, max_nfev: int, limit,
            refresh: bool = False) -> None:
    out = output_csv_for(input_path, threshold)
    if out.exists():
        # Plain skip is what makes an interrupted run resumable. But after a
        # column is added, "exists" no longer means "current" -- every stale
        # table would be skipped and the run would silently do nothing.
        # --refresh keeps the resumability and fixes that.
        if not refresh or header_is_current(out):
            print(f"  skip (exists): {out.name}", flush=True)
            return
        print(f"  refresh (stale columns): {out.name}", flush=True)
    df = load_expression(input_path, id_col=ID_COL, drop_cols=[NAME_COL])
    table = build_table(df, threshold, jobs, max_nfev, limit)
    table = insert_genename(table, load_name_map(input_path, ID_COL, NAME_COL))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    # Pin LF so every table in outputs/ shares one line terminator regardless of
    # platform; build_gene_major.py copies row text verbatim.
    table.to_csv(out, index=False, lineterminator="\n")
    n_ok = int(table["fit_success"].sum())
    print(f"  wrote {out.name} ({len(table)} rows, {n_ok} fit_success)", flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--jobs", type=int, default=8)
    p.add_argument("--max-nfev", type=int, default=2000, dest="max_nfev")
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--refresh", action="store_true",
                   help="Also regenerate tables whose columns are out of date "
                        "(header != the generator's current COLUMNS). Without "
                        "this, every existing table is skipped, so a run after "
                        "a column change silently does nothing.")
    args = p.parse_args()

    tissues = sorted(DATA_DIR.glob("v11_log2_*.csv.gz"))
    if not tissues:
        raise SystemExit(f"No tissue matrices found in {DATA_DIR} "
                         f"(expected v11_log2_*.csv.gz).")
    print(f"{len(tissues)} tissue matrices -> {2 * len(tissues)} tables", flush=True)

    for i, tissue in enumerate(tissues, 1):
        print(f"[{i}/{len(tissues)}] {tissue.name}", flush=True)
        for threshold in THRESHOLDS:
            run_one(tissue, threshold, args.jobs, args.max_nfev, args.limit,
                    refresh=args.refresh)

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
