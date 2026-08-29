# -*- coding: utf-8 -*-
"""
stage_worm_table.py

Copy the C. elegans fourparam table into ``worm/`` for publication, repairing
the gene names Excel turned into dates on the way.

Why this exists
---------------
Eleven *C. elegans* gene names are also month abbreviations -- ``mar-1``,
``apr-1``, ``jun-1``, ``sep-1``, ``oct-1`` and friends -- and any spreadsheet
that has touched the table has silently rewritten them as dates. In the table
handed over they arrive as ``2025-03-01 00:00:00``. This is the well-known
Excel gene-name problem, and it is not cosmetic: those genes cannot be found by
name, they sort to the top of an alphabetical listing, and ``apr-1`` in
particular is one nobody wants to lose.

The repair is unambiguous and reversible: month number -> abbreviation, day
number -> the suffix, so ``2025-04-01`` is ``apr-1`` and nothing else. Only the
``genename`` column is touched; ``gene`` and ``wormbasegeneid`` are untouched
and remain the join keys, so every rewrite can be checked against WormBase.

The source file is read, never written. Re-running is a no-op on a table whose
names are already clean.

Usage
-----
    cd fourparam
    python stage_worm_table.py --source ../../../worm_fourparam_table_excluded_at_or_below_-1.csv
    python stage_worm_table.py --source <path> --dry-run
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_DEST = HERE.parent / "worm" / "worm_fourparam_excluded_at_or_below_-1.csv"

# ``2025-03-01 00:00:00`` and the bare ``2025-03-01`` some exporters emit.
DATE_RE = re.compile(r"^(\d{4})-(\d{2})-(\d{2})(?:[ T]00:00:00)?$")
MONTHS = {1: "jan", 2: "feb", 3: "mar", 4: "apr", 5: "may", 6: "jun",
          7: "jul", 8: "aug", 9: "sep", 10: "oct", 11: "nov", 12: "dec"}


def unmangle(name: str) -> str | None:
    """``2025-04-01 00:00:00`` -> ``apr-1``. None when the name is already fine."""
    m = DATE_RE.match(name.strip())
    if not m:
        return None
    month, day = int(m.group(2)), int(m.group(3))
    if month not in MONTHS or not 1 <= day <= 31:
        return None
    return f"{MONTHS[month]}-{day}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", type=Path, required=True,
                    help="The worm fourparam table as handed over.")
    ap.add_argument("--dest", type=Path, default=DEFAULT_DEST,
                    help=f"Where to write it (default: {DEFAULT_DEST.name} in worm/).")
    ap.add_argument("--dry-run", action="store_true",
                    help="Report what would be repaired, write nothing.")
    args = ap.parse_args(argv)

    if not args.source.is_file():
        print(f"no such file: {args.source}", file=sys.stderr)
        return 2

    with args.source.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.reader(fh))
    header, body = rows[0], rows[1:]
    if "genename" not in header:
        print(f"no genename column in {args.source.name}", file=sys.stderr)
        return 2
    gi, ni = header.index("genename"), header.index("gene")

    fixed = []
    for row in body:
        better = unmangle(row[gi])
        if better is not None:
            fixed.append((row[ni], row[gi], better))
            row[gi] = better

    for gene, was, now in fixed:
        print(f"  {gene:<24} {was}  ->  {now}")
    print(f"{len(fixed)} gene name(s) repaired, {len(body):,} rows")

    if args.dry_run:
        print("(dry run, nothing written)")
        return 0

    args.dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.dest.with_suffix(".csv.tmp")
    try:
        # LF, matching outputs/, so the published tables all read the same way.
        with tmp.open("w", encoding="utf-8", newline="\n") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(header)
            w.writerows(body)
        tmp.replace(args.dest)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    print(f"Wrote {args.dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
