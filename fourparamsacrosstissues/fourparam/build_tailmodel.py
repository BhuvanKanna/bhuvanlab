# -*- coding: utf-8 -*-
"""
build_tailmodel.py

Publish the tail model's ``rti_missing`` / ``lti_missing`` and ``idklti`` to the
browser, as the page's **newRTI**, **newLTI** and **idkLTI** columns.

Same trade as ``build_r2.py``: one fixed-width record per gene in **genes.tsv
order**, so the browser reads gene ``gi`` as a fixed-offset slice with no join
and nothing to parse. Each record is three ``%7.5f`` fields -- rti_missing,
lti_missing, idklti -- so 21 characters per gene, ~1.6 MB per tissue. Five decimals because the values that
matter are small: under a true Gaussian of N donors ``rti_missing`` is ~1/(N+1),
0.00122 in muscle. Both are fractions in [0, 1], so nothing needs clamping. A
missing value is seven spaces: a failed or not-``usable`` tail fit blanks
rti/lti (a degenerate fit is not a number worth showing), as does an LTI that is
undefined because donors are censored at the floor. ``idklti`` comes from the
4-parameter fit, not the tail model, so it does not depend on ``usable``.

    cd fourparam
    python build_tailmodel.py                          # every tailmodel/ table
    python build_tailmodel.py --tissues muscle_skeletal

Writes ``docs/tailmodel/<same stem as the CSV>.txt`` and patches the
``tailmodel`` block of ``docs/manifest.json`` in place.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SRC_DIR = HERE.parent / "tailmodel"
DOCS = HERE.parent.parent / "docs"
DOCS_OUT = DOCS / "tailmodel"
GENES_TSV = DOCS / "genes.tsv"
MANIFEST = DOCS / "manifest.json"

FIELDS = ["rti_missing", "lti_missing", "idklti"]
GATED = {"rti_missing", "lti_missing"}       # blanked unless the fit is usable
WIDTH = 7                 # "%7.5f" -> "0.00122"
DECIMALS = 5
MISSING = " " * WIDTH


def strip_version(gid: str) -> str:
    return gid.split(".", 1)[0]


def load_gene_order() -> list[str]:
    ids = []
    with GENES_TSV.open(encoding="utf-8") as f:
        for line in f:
            if line.strip():
                ids.append(strip_version(line.split("\t", 1)[0]))
    return ids


def encode(csv: Path, order: list[str]) -> tuple[str, dict]:
    df = pd.read_csv(csv, dtype={"gene": str})
    keys = df["gene"].map(strip_version)
    usable = (df["usable"].astype(str).str.lower() == "true"
              if "usable" in df else pd.Series(True, index=df.index))
    luts = {}
    for f in FIELDS:
        v = (pd.to_numeric(df[f], errors="coerce") if f in df
             else pd.Series(np.nan, index=df.index))
        if f in GATED:
            v = v.where(usable)
        luts[f] = dict(zip(keys, v))

    parts, counts = [], {f: 0 for f in FIELDS}
    for gid in order:
        for f in FIELDS:
            v = luts[f].get(gid)
            if v is None or not np.isfinite(v):
                parts.append(MISSING)
                continue
            if not 0.0 <= v <= 1.0:
                raise SystemExit(f"{csv.name}: {gid} {f}={v!r} is not a fraction")
            parts.append(f"{v:{WIDTH}.{DECIMALS}f}")
            counts[f] += 1
    text = "".join(parts)
    assert len(text) == WIDTH * len(FIELDS) * len(order), "field width drifted"
    return text, counts


def patch_manifest(published: dict) -> None:
    """Merge a ``tailmodel`` block; build_gui_data.py carries it over."""
    man = json.loads(MANIFEST.read_text(encoding="utf-8"))
    man["tailmodel"] = {
        "available": bool(published),
        "base_url": "tailmodel",
        "width": WIDTH,
        "decimals": DECIMALS,
        "fields": FIELDS,
        # {tissue: filename}. Computed from raw rows, so not keyed by table kind.
        "files": published,
    }
    MANIFEST.write_text(json.dumps(man, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tissues", type=str, default=None)
    a = p.parse_args()

    order = load_gene_order()
    wanted = set(a.tissues.split(",")) if a.tissues else None
    DOCS_OUT.mkdir(parents=True, exist_ok=True)

    published = {}
    for csv in sorted(SRC_DIR.glob("v11_log2_*_tailmodel.csv")):
        tissue = csv.stem[len("v11_log2_"):-len("_tailmodel")]
        if wanted and tissue not in wanted:
            continue
        text, counts = encode(csv, order)
        out = DOCS_OUT / f"{csv.stem}.txt"
        out.write_text(text, encoding="ascii", newline="")
        published[tissue] = out.name
        print(f"{out.name:<48} {len(text):>10,} chars  {counts}")

    patch_manifest(published)
    print(f"patched {MANIFEST.name}: {sorted(published)}")


if __name__ == "__main__":
    main()
