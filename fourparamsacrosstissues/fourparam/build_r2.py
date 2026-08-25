# -*- coding: utf-8 -*-
"""
build_r2.py

Publish `r_squared` to the browser as a fixed-width field per gene.

Same trade as `build_qc_class.py` and the `hist` column: quantise hard, keep it
a fixed width, and never put a comma in it. Each `r2/` table becomes a string of
exactly `6 * len(genes.tsv)` characters in **genes.tsv order** -- about 437 KB --
and the browser already holds a gene's row index (`gi`), so the lookup is
`text.substr(gi * 6, 6)` with no join and nothing to parse.

Six characters is `%6.3f`: `" 0.599"`, `"-0.012"`, `" 1.000"`. Three decimals is
far finer than the column is ever read to, and the sign slot is load-bearing --
R^2 below zero means the fit is worse than a flat line at the mean count, which
is a real verdict and must not be clipped to 0. Values outside +/-9.999 are
clamped to the bound so the width holds; `--strict` refuses instead, and the run
reports the count either way. A gene with no R^2 is six spaces, which parses to
NaN rather than to a number.

    cd fourparam
    python build_r2.py                       # every r2/ table that exists
    python build_r2.py --tissues kidney_cortex

Writes `docs/r2/<same stem as the CSV>.txt` and patches the `r2` block of
`docs/manifest.json` in place, leaving every other key untouched.

Re-run this after computing any new R^2 table. `build_gui_data.py` rewrites
manifest.json wholesale but carries over the blocks it finds, so the two do not
have to run in a particular order.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
R2_DIR = HERE.parent / "r2"
DOCS = HERE.parent.parent / "docs"
DOCS_R2 = DOCS / "r2"
GENES_TSV = DOCS / "genes.tsv"
MANIFEST = DOCS / "manifest.json"

WIDTH = 6                 # "%6.3f" -> " 0.599", "-0.012"
DECIMALS = 3
MISSING = " " * WIDTH
CLAMP = 9.999             # the largest magnitude %6.3f can hold


def strip_version(gid: str) -> str:
    return gid.split(".", 1)[0]


def load_gene_order() -> list[str]:
    """Unversioned ids in genes.tsv order -- the browser's `gi` indexes this."""
    ids = []
    with GENES_TSV.open(encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            ids.append(strip_version(line.split("\t", 1)[0]))
    return ids


def encode(r2_csv: Path, order: list[str], strict: bool) -> tuple[str, dict]:
    df = pd.read_csv(r2_csv, dtype={"gene": str})
    df["r_squared"] = pd.to_numeric(df["r_squared"], errors="coerce")
    # Versions drift between GTEx releases; genes.tsv and the r2 tables are built
    # from the same annotation today, but matching on the versioned string is the
    # single easiest way to silently produce an all-blank file.
    lut = dict(zip(df["gene"].map(strip_version), df["r_squared"]))

    parts, n_val, n_clamped, n_missing = [], 0, 0, 0
    for gid in order:
        v = lut.get(gid)
        if v is None or not np.isfinite(v):
            parts.append(MISSING)
            n_missing += 1
            continue
        if abs(v) > CLAMP:
            if strict:
                raise SystemExit(
                    f"{r2_csv.name}: {gid} has R^2 {v!r}, outside +/-{CLAMP} "
                    f"-- would not fit {WIDTH} chars")
            v = float(np.clip(v, -CLAMP, CLAMP))
            n_clamped += 1
        parts.append(f"{v:{WIDTH}.{DECIMALS}f}")
        n_val += 1

    text = "".join(parts)
    assert len(text) == WIDTH * len(order), "field width drifted"
    return text, {"values": n_val, "missing": n_missing, "clamped": n_clamped}


def patch_manifest(published: dict) -> None:
    """
    Merge an `r2` block into manifest.json without disturbing anything else.

    build_gui_data.py rewrites this file from scratch, so this must stay a patch
    rather than a template -- otherwise whichever ran last wins and the other's
    output vanishes.
    """
    man = json.loads(MANIFEST.read_text(encoding="utf-8"))
    man["r2"] = {
        "available": bool(published),
        "base_url": "r2",
        # The browser slices by these rather than hardcoding 6, so the field can
        # widen later without a code change on the page.
        "width": WIDTH,
        "decimals": DECIMALS,
        "clamp": CLAMP,
        # Mirrors qc.kind_ids so the r2 block stands alone even when the others
        # are unavailable.
        "kind_ids": {"raw": "raw", "excluded": "excluded_at_or_below_-1"},
        # {tissue: {kind_id: filename}} -- the browser only fetches what is here,
        # so a tissue whose R^2 has not been computed simply shows no value.
        "files": published,
    }
    MANIFEST.write_text(json.dumps(man, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tissues", type=str, default=None,
                   help="Comma-separated subset; default is every r2/ table found.")
    p.add_argument("--strict", action="store_true",
                   help="Fail on an R^2 too large for the field instead of clamping.")
    args = p.parse_args()

    if not GENES_TSV.exists():
        raise SystemExit(f"missing {GENES_TSV} -- run build_gui_data.py first")
    order = load_gene_order()
    print(f"gene order: {len(order):,} ids from {GENES_TSV.name}")

    wanted = set(args.tissues.split(",")) if args.tissues else None
    DOCS_R2.mkdir(parents=True, exist_ok=True)

    published: dict = {}
    total_clamped = 0
    for csv in sorted(R2_DIR.glob("v11_log2_*_r2*.csv")):
        stem = csv.stem                       # v11_log2_<tissue>_r2[_excluded...]
        rest = stem.split("v11_log2_", 1)[1]
        tissue, _, tail = rest.partition("_r2")
        if wanted and tissue not in wanted:
            continue
        kind = "raw" if tail == "" else tail.lstrip("_")

        text, counts = encode(csv, order, args.strict)
        out = DOCS_R2 / f"{stem}.txt"
        out.write_text(text, encoding="ascii", newline="")
        published.setdefault(tissue, {})[kind] = out.name
        total_clamped += counts["clamped"]

        print(f"{out.name:<62} {len(text):>9,} chars  {counts}")

    patch_manifest(published)
    n_files = sum(len(v) for v in published.values())
    print(f"\npatched {MANIFEST.name}: {len(published)} tissue(s), {n_files} file(s)")
    if total_clamped:
        print(f"clamped {total_clamped:,} value(s) to +/-{CLAMP}")


if __name__ == "__main__":
    main()
