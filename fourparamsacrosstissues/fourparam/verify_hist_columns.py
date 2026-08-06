#!/usr/bin/env python
"""
verify_hist_columns.py

Check the `hist` / `hist_max` columns in generated fourparam tables.

These two columns feed the browser's Shape cell, and every invariant below is
one the browser silently depends on. Run this after regenerating tables and
before pushing them.

    python verify_hist_columns.py                     # every table in outputs/
    python verify_hist_columns.py --tissues uterus,vagina
    python verify_hist_columns.py --table raw

Invariants checked, per table:

1. `hist` and `hist_max` are the last two columns, in that order.
2. Every line splits into exactly len(header) fields, and no double-quote
   appears anywhere. The browser parses these CSVs with a plain split(",") and
   has no quoted-field handling, so one stray comma shifts every column index
   after it -- silently, and for every row.
3. Every non-empty `hist` is exactly HIST_BINS characters drawn from
   HIST_ALPHABET.
4. `hist_max` never serialises as a float (a NaN anywhere in the column makes
   pandas write "12.0" for every other row).
5. **No histogram ships without its bin edges.** Edges are deliberately not
   stored -- they are linspace(min, max, 41) -- so a row carrying `hist` with a
   NaN `min` or `max` cannot be drawn at all.
6. `hist` is empty exactly when n_obs == 0.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from bhuvanfitter import HIST_ALPHABET, HIST_BINS

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE.parent / "outputs"
ALPHABET = set(HIST_ALPHABET)


def check_table(path: Path) -> list[str]:
    """Return a list of problem strings; empty means the table is clean."""
    problems: list[str] = []
    text = path.read_text(encoding="utf-8")
    lines = text.strip().split("\n")
    header = lines[0].split(",")

    if header[-2:] != ["hist", "hist_max"]:
        return [f"last two columns are {header[-2:]}, expected ['hist', 'hist_max']"]

    widths = {len(line.split(",")) for line in lines}
    if widths != {len(header)}:
        problems.append(f"ragged rows: field counts {sorted(widths)}, "
                        f"expected {len(header)} -- a comma leaked into a field")
    if '"' in text:
        problems.append("a double-quote appears in the file; the browser's "
                        "split(\",\") cannot handle quoted fields")

    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    hist, hmax, n_obs = d["hist"], d["hist_max"], d["n_obs"]

    non_empty = hist[hist != ""]
    bad_len = (non_empty.str.len() != HIST_BINS).sum()
    if bad_len:
        problems.append(f"{bad_len} rows whose hist is not {HIST_BINS} chars")

    stray = {c for s in non_empty.head(20000) for c in s} - ALPHABET
    if stray:
        problems.append(f"characters outside the alphabet in hist: {sorted(stray)}")

    if hmax.str.contains(r"\.", regex=True).any():
        problems.append("hist_max serialised as a float somewhere")

    no_edges = (d["min"] == "") | (d["min"] == "nan") | \
               (d["max"] == "") | (d["max"] == "nan")
    orphan = int(((hist.str.len() == HIST_BINS) & no_edges).sum())
    if orphan:
        problems.append(f"{orphan} histograms carry no min/max, so they cannot "
                        f"be drawn at all")

    zero = n_obs == "0"
    if int((zero & (hist != "")).sum()):
        problems.append("a row with n_obs=0 carries a histogram")
    if int((~zero & (hist == "")).sum()):
        problems.append("a row with n_obs>0 is missing its histogram")

    return problems


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tissues", type=str, default=None,
                   help="Comma-separated tissue names (default: all found).")
    p.add_argument("--table", choices=["raw", "excluded", "both"], default="both",
                   help="Which table type to check (default: both).")
    args = p.parse_args()

    paths = sorted(OUT_DIR.glob("*_fourparam*.csv"))
    if args.tissues:
        wanted = [t.strip() for t in args.tissues.split(",") if t.strip()]
        paths = [p_ for p_ in paths if any(f"_{t}_fourparam" in p_.name for t in wanted)]
    if args.table == "raw":
        paths = [p_ for p_ in paths if p_.name.endswith("_fourparam.csv")]
    elif args.table == "excluded":
        paths = [p_ for p_ in paths if not p_.name.endswith("_fourparam.csv")]

    if not paths:
        sys.exit("No matching tables in outputs/.")

    failures = 0
    for path in paths:
        problems = check_table(path)
        if problems:
            failures += 1
            print(f"FAIL  {path.name}", flush=True)
            for problem in problems:
                print(f"        {problem}", flush=True)
        else:
            print(f"ok    {path.name}", flush=True)

    print(f"\n{len(paths) - failures}/{len(paths)} tables clean", flush=True)
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
