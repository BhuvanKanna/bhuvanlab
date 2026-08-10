# -*- coding: utf-8 -*-
"""
build_qc_class.py

Publish `dist_class` to the browser as one character per gene.

The QC tables are ~14 MB each. Fetching one just to colour a single column would
roughly double what the page already downloads per tissue, so instead each table
becomes a string of exactly `len(genes.tsv)` characters in **genes.tsv order** --
about 73 KB, and the browser already holds a gene's row index (`gi`), so the
lookup is `text[gi]` with no join and no parsing.

This is the same trade the `hist` column makes: quantise hard, keep it a fixed
width, and never put a comma in it.

    cd fourparam
    python build_qc_class.py                 # every qc/ table that exists
    python build_qc_class.py --tissues kidney_cortex

Writes `docs/qc/<same stem as the CSV>.txt` and patches the `qc` block of
`docs/manifest.json` in place, leaving every other key untouched.

Re-run this after generating any new QC table. `build_gui_data.py` rewrites
manifest.json wholesale but carries over whatever `qc` block it finds, so the two
no longer have to run in a particular order.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
QC_DIR = HERE.parent / "qc"
DOCS = HERE.parent.parent / "docs"
DOCS_QC = DOCS / "qc"
GENES_TSV = DOCS / "genes.tsv"
MANIFEST = DOCS / "manifest.json"

# One character per class. Lower case throughout so the string stays visually
# uniform; `.` is "no verdict for this gene", which is not the same as
# `undetermined` (the cascade looked and could not decide).
CLASS_CHAR = {
    "normal": "n",
    "right_truncated": "t",
    "right_skewed": "s",
    "non_normal": "x",
    "multimodal": "m",
    "zero_inflated": "z",
    "undetermined": "u",
}
MISSING = "."


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


def encode(qc_csv: Path, order: list[str]) -> tuple[str, dict]:
    df = pd.read_csv(qc_csv, usecols=["gene", "dist_class"],
                     dtype=str, keep_default_na=False)
    # Versions drift between GTEx releases; genes.tsv and the QC tables are built
    # from the same annotation today, but matching on the versioned string is the
    # single easiest way to silently produce an all-dots file.
    lut = dict(zip(df["gene"].map(strip_version), df["dist_class"]))
    chars = [CLASS_CHAR.get(lut.get(g, ""), MISSING) for g in order]
    text = "".join(chars)
    counts = {}
    for c in text:
        counts[c] = counts.get(c, 0) + 1
    return text, counts


def patch_manifest(published: dict) -> None:
    """
    Merge a `qc` block into manifest.json without disturbing anything else.

    build_gui_data.py rewrites this file from scratch, so this must stay a patch
    rather than a template -- otherwise whichever ran last wins and the other's
    output vanishes.
    """
    man = json.loads(MANIFEST.read_text(encoding="utf-8"))
    man["qc"] = {
        "available": bool(published),
        "base_url": "qc",
        "classes": {v: k for k, v in CLASS_CHAR.items()},
        "missing_char": MISSING,
        # The page's toggle says "raw"/"excluded"; the filenames say what the
        # tables say. Mirrors gene_major.kind_ids so the qc block stands alone
        # even when gene_major is unavailable.
        "kind_ids": {"raw": "raw", "excluded": "excluded_at_or_below_-1"},
        # {tissue: {kind_id: filename}} -- the browser only fetches what is here,
        # so a tissue whose QC has not been computed simply shows no class.
        "files": published,
    }
    MANIFEST.write_text(json.dumps(man, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tissues", type=str, default=None,
                   help="Comma-separated subset; default is every qc/ table found.")
    args = p.parse_args()

    if not GENES_TSV.exists():
        raise SystemExit(f"missing {GENES_TSV} -- run build_gui_data.py first")
    order = load_gene_order()
    print(f"gene order: {len(order):,} ids from {GENES_TSV.name}")

    wanted = set(args.tissues.split(",")) if args.tissues else None
    DOCS_QC.mkdir(parents=True, exist_ok=True)

    published: dict = {}
    for csv in sorted(QC_DIR.glob("v11_log2_*_qc*.csv")):
        stem = csv.stem                       # v11_log2_<tissue>_qc[_excluded...]
        rest = stem.split("v11_log2_", 1)[1]
        tissue, _, tail = rest.partition("_qc")
        if wanted and tissue not in wanted:
            continue
        kind = "raw" if tail == "" else tail.lstrip("_")

        text, counts = encode(csv, order)
        out = DOCS_QC / f"{stem}.txt"
        out.write_text(text, encoding="ascii", newline="")
        published.setdefault(tissue, {})[kind] = out.name

        named = {CLASS_CHAR_INV.get(c, c): n for c, n in sorted(counts.items())}
        print(f"{out.name:<58} {len(text):>7,} chars  {named}")

    patch_manifest(published)
    n_files = sum(len(v) for v in published.values())
    print(f"\npatched {MANIFEST.name}: {len(published)} tissue(s), {n_files} file(s)")


CLASS_CHAR_INV = {v: k for k, v in CLASS_CHAR.items()}
CLASS_CHAR_INV[MISSING] = "missing"


if __name__ == "__main__":
    main()
