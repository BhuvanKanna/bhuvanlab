#!/usr/bin/env python
"""
backfill_hist.py

Add the ``hist`` / ``hist_max`` columns to excluded tables that predate them,
**without refitting**. The histogram is a property of the data alone (40 bins
over ``[min, max]`` of the post-exclusion values), so it is recomputed from the
tissue's matrix in ``data/`` with the generator's own filter and
``bhuvanfitter.encode_histogram`` — the same code path a full regeneration
would take — and spliced in after ``fit_success``, where ``COLUMNS`` puts it.

Each line is split on commas and the two new fields inserted as text after
``fit_success``, so every existing value stays byte-identical (``r_squared``,
if present, stays last) with one exception: tables written before ``hist``
existed left ``min`` / ``max`` blank on a failed row, and those are the
histogram's bin edges. A blank one is filled with ``repr()`` of the data's
min/max, which is the text the current generator writes there.

Before writing, every gene's recomputed ``n_obs`` -- and every ``min`` / ``max``
that is not blank -- must match the table's own. Any mismatch means the matrix is not the one the table
was fitted from, and that tissue is refused rather than given histograms that
disagree with its curve.

    python backfill_hist.py                      # every excluded table missing hist
    python backfill_hist.py --tissues liver,lung
    python backfill_hist.py --check uterus       # recompute a table that HAS hist, compare only
"""
from __future__ import annotations

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from bhuvanfitter import encode_histogram

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT_DIR = ROOT / "outputs"
DATA_DIR = ROOT / "data"
SUFFIX = "_fourparam_excluded_at_or_below_-1.csv"
THRESHOLD = -1.0
SPLICE_AFTER = "fit_success"


def table_path(tissue: str) -> Path:
    return OUT_DIR / f"v11_log2_{tissue}{SUFFIX}"


def histograms(tissue: str) -> dict:
    """gene -> (hist, hist_max, n_obs, min, max), exactly as the generator computes them."""
    df = pd.read_csv(DATA_DIR / f"v11_log2_{tissue}.csv.gz")
    ids = df["Name"].to_numpy()
    vals = df.drop(columns=["Name", "Description"]).to_numpy(dtype=float)
    out = {}
    for gene, row in zip(ids, vals):
        data = row[np.isfinite(row)]
        data = data[data > THRESHOLD]
        h, hm = encode_histogram(data)
        n = int(data.size)
        out[gene] = (h, hm, n,
                     float(data.min()) if n else np.nan,
                     float(data.max()) if n else np.nan)
    return out


def _close(a: str, b: float) -> bool:
    if a == "" or a.lower() == "nan":
        return np.isnan(b)
    return np.isclose(float(a), b, rtol=1e-12, atol=1e-12)


def process(tissue: str, check_only: bool) -> str:
    path = table_path(tissue)
    lines = path.read_text(encoding="utf-8").split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    header = lines[0].split(",")
    has = "hist" in header
    if has and not check_only:
        return f"{tissue}: already has hist, skipped"
    if not has and check_only:
        return f"{tissue}: no hist to check against"

    ci = {c: i for i, c in enumerate(header)}
    k = ci[SPLICE_AFTER] + 1
    H = histograms(tissue)

    new, bad, filled = [], [], 0
    for ln in lines[1:]:
        f = ln.split(",")
        if len(f) != len(header):
            return f"{tissue}: REFUSED — ragged row ({len(f)} fields): {ln[:60]}"
        g = f[ci["gene"]]
        if g not in H:
            return f"{tissue}: REFUSED — gene {g} not in data matrix"
        h, hm, n, lo, hi = H[g]
        if int(float(f[ci["n_obs"]])) != n:
            bad.append(g)
            continue
        # Tables written before `hist` existed leave min/max blank on a failed
        # row; the current generator fills them, because they are the
        # histogram's bin edges. Fill exactly those, with repr() -- the same
        # text pandas writes -- and require every filled one to agree.
        mismatch = False
        for col, val in (("min", lo), ("max", hi)):
            if f[ci[col]] == "" and n:
                if not check_only:
                    f[ci[col]] = repr(val)
                    filled += 1
            elif not _close(f[ci[col]], val):
                mismatch = True
        if mismatch:
            bad.append(g)
            continue
        if check_only:
            if f[ci["hist"]] != h or int(float(f[ci["hist_max"]] or 0)) != hm:
                bad.append(g)
            continue
        new.append(",".join(f[:k] + [h, str(hm)] + f[k:]))

    if len(lines) - 1 != len(H):
        return f"{tissue}: REFUSED — {len(lines) - 1} rows vs {len(H)} genes in matrix"
    if bad:
        return (f"{tissue}: REFUSED — {len(bad)} genes disagree with the table "
                f"(first: {bad[:3]})")
    if check_only:
        return f"{tissue}: OK — all {len(H):,} recomputed histograms match"

    out_header = header[:k] + ["hist", "hist_max"] + header[k:]
    tmp = path.with_suffix(".csv.tmp")
    tmp.write_text("\n".join([",".join(out_header)] + new) + "\n", encoding="utf-8",
                   newline="\n")
    tmp.replace(path)
    return (f"{tissue}: added hist to {len(new):,} rows, "
            f"filled {filled:,} blank min/max")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tissues", help="Comma-separated; default every excluded table.")
    ap.add_argument("--check", metavar="TISSUES",
                    help="Recompute tables that already carry hist and compare only.")
    ap.add_argument("--jobs", type=int, default=4,
                    help="Tissues in parallel (each holds one matrix in memory).")
    a = ap.parse_args()

    check = a.check is not None
    if check:
        tissues = a.check.split(",")
    elif a.tissues:
        tissues = a.tissues.split(",")
    else:
        tissues = sorted(p.name[len("v11_log2_"):-len(SUFFIX)]
                         for p in OUT_DIR.glob(f"v11_log2_*{SUFFIX}"))

    failed = False
    with ProcessPoolExecutor(a.jobs) as ex:
        for msg in ex.map(process, tissues, [check] * len(tissues)):
            print(msg, flush=True)
            failed |= "REFUSED" in msg
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
