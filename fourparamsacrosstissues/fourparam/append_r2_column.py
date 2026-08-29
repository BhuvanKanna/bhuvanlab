# -*- coding: utf-8 -*-
"""
append_r2_column.py

Append `r_squared` as a real last column to the excluded fourparam tables in
`outputs/`, so a table downloaded straight from the repo carries the same number
the browser shows -- without anyone having to join `r2/` by hand.

Nothing is refit and nothing is recomputed: the values come verbatim from the
full-precision `r2/<stem>_r2_excluded_at_or_below_-1.csv` that `compute_r2.py`
already wrote. This step is pure I/O.

Why last, and why only the excluded tables
------------------------------------------
Last, because every existing consumer indexes these tables positionally --
`extract_genes.py`, the browser's `indexTable`, the shard builders. Appending
leaves every byte left of the new column where it was. Only the excluded
tables, because the browser no longer serves the raw ones and `r2/` was only
ever computed for excluded.

The browser still reads R^2 from `r2/`, not from this column
------------------------------------------------------------
`docs/index.html` joins R^2 by gene index (`R2_COL`) rather than reading it out
of the table, because the gene-major shards in `gene_major/` are a 2.4 GB
mirror built from these files and are not being rebuilt. The join is the only
source that works on *both* load routes, so it stays the browser's source and
`r_squared` is kept out of `manifest.columns`. The two are the same numbers
from the same file; this column exists for whoever downloads the CSV directly.

Safety
------
Gene order is asserted row by row against `r2/`, not merely by row count: a
silent off-by-one would attach every R^2 to the wrong gene. Any mismatch aborts
that tissue and leaves its table untouched. Writes go to a temp file next to
the target and are moved into place, so an interrupted run cannot truncate a
17 MB table. Re-running is a no-op on tables that already have the column.

Usage
-----
    cd fourparam
    python append_r2_column.py --all
    python append_r2_column.py --all --dry-run
    python append_r2_column.py --tissue kidney_cortex
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUTS = ROOT / "outputs"
R2_DIR = ROOT / "r2"

SUFFIX = "_fourparam_excluded_at_or_below_-1.csv"
R2_SUFFIX = "_r2_excluded_at_or_below_-1.csv"
COLUMN = "r_squared"


def r2_path_for(table: Path) -> Path:
    """`<stem>_fourparam_excluded...` -> `<stem>_r2_excluded...`."""
    return R2_DIR / (table.name[: -len(SUFFIX)] + R2_SUFFIX)


def load_r2(path: Path) -> list[tuple[str, str]]:
    """(gene, r_squared-as-written) per row, in file order. Blank stays blank."""
    rows = []
    with path.open("r", encoding="utf-8", newline="") as fh:
        header = fh.readline().rstrip("\r\n").split(",")
        if header != ["gene", COLUMN]:
            raise ValueError(f"{path.name}: unexpected header {header}")
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            gene, _, value = line.partition(",")
            rows.append((gene, value))
    return rows


def append(table: Path, dry_run: bool) -> str:
    r2_file = r2_path_for(table)
    if not r2_file.exists():
        return f"skip (no r2/ file: {r2_file.name})"

    with table.open("r", encoding="utf-8", newline="") as fh:
        header = fh.readline().rstrip("\r\n")
    if header.split(",")[-1] == COLUMN:
        return "skip (already has the column)"

    values = load_r2(r2_file)
    tmp = table.with_suffix(".csv.tmp")
    written = 0
    try:
        with table.open("r", encoding="utf-8", newline="") as src, \
                tmp.open("w", encoding="utf-8", newline="") as dst:
            src.readline()                       # header, re-emitted below
            dst.write(f"{header},{COLUMN}\n")
            for i, line in enumerate(src):
                line = line.rstrip("\r\n")
                if not line:
                    continue
                if i >= len(values):
                    raise ValueError(
                        f"{table.name}: table has more rows than r2/ ({i + 1} > {len(values)})")
                gene, value = values[i]
                # Row-by-row, not just a count: an off-by-one here would put
                # every R^2 on the wrong gene, and nothing downstream would
                # ever notice.
                if not line.startswith(gene + ","):
                    raise ValueError(
                        f"{table.name}: row {i + 1} is {line.split(',')[0]!r}, "
                        f"r2/ has {gene!r}")
                dst.write(f"{line},{value}\n")
                written += 1
        if written != len(values):
            raise ValueError(
                f"{table.name}: wrote {written} rows, r2/ has {len(values)}")
        if dry_run:
            tmp.unlink()
            return f"ok (dry run, {written} rows)"
        os.replace(tmp, table)
        return f"ok ({written} rows)"
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all", action="store_true",
                    help="every excluded table in outputs/")
    ap.add_argument("--tissue", action="append", default=[],
                    help="tissue stem, e.g. kidney_cortex (repeatable)")
    ap.add_argument("--dry-run", action="store_true",
                    help="verify the join and row counts, write nothing")
    args = ap.parse_args()

    if args.tissue:
        tables = [OUTPUTS / f"v11_log2_{t}{SUFFIX}" for t in args.tissue]
        missing = [p for p in tables if not p.exists()]
        if missing:
            print("no such table: " + ", ".join(p.name for p in missing), file=sys.stderr)
            return 2
    elif args.all:
        tables = sorted(OUTPUTS.glob(f"*{SUFFIX}"))
    else:
        ap.error("pass --all or --tissue")

    failures = 0
    for table in tables:
        try:
            note = append(table, args.dry_run)
        except Exception as err:                 # keep going; report at the end
            note = f"FAILED: {err}"
            failures += 1
        print(f"{table.name:<62} {note}", flush=True)

    print(f"\n{len(tables)} table(s), {failures} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
