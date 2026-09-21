# -*- coding: utf-8 -*-
"""
migrate_rti_lti.py

Migrate existing fourparam tables to the RTI / LTI schema, **without refitting**.

Three things happen per table:

1. ``ti_fourparam_sigma_dist`` -> ``rti_sigma_dist`` and ``truncationindex`` ->
   ``rti``. Pure renames; the numbers are untouched.
2. Four new left-side columns are computed from the stored fit parameters:
   ``lti_sigma_dist``, ``lti``, ``left`` and ``leftheight``.
3. Columns are reordered into the new schema (left/right pairs adjacent), with
   ``r_squared`` -- which is not part of ``COLUMNS`` -- kept last.

Nothing is refit. Every new value is derived from ``y0, A, x0, w, min, max``,
which the table already carries, so the curve is reconstructed rather than
re-estimated. ``--verify`` proves the reconstruction is the right one by
recomputing ``maxheight`` and ``rightheight`` from it and comparing against the
stored columns.

**Carried-over values are copied as text and never parsed.** The byte-identity
guarantee between ``outputs/``, ``gene_major/`` and ``extract_genes.py`` holds
only because nothing re-serialises a float: pandas' writer emits ~16 significant
digits rather than the shortest round-tripping repr, which silently turns
``0.012596832467784065`` into ``0.012596832467784``. Floats are parsed here only
to *compute the new columns*; the parsed values are then thrown away.

Usage
-----
    python migrate_rti_lti.py --all --verify        # every table in outputs/
    python migrate_rti_lti.py --tissues liver       # just one
    python migrate_rti_lti.py --input path.csv --output out.csv
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from bhuvanfitter import _fourparam_gaussian
from generate_fourparam import COLUMNS

# The schema as it appears in a table: COLUMNS with genename inserted at 1.
NEW_COLUMNS = ["gene", "genename"] + COLUMNS[1:]

RENAMES = {"ti_fourparam_sigma_dist": "rti_sigma_dist",
           "truncationindex": "rti"}

# The columns this script derives. Everything else is copied verbatim.
DERIVED = ["lti_sigma_dist", "lti", "left", "leftheight"]

# Must match BhuvanFitter._CURVE_GRID -- the grid the stored maxheight and
# rightheight were measured on. A different resolution would give subtly
# different baselines and the --verify check would (correctly) fail.
CURVE_GRID = 600

# Genes per vectorised block. 600 grid points x 5000 genes is ~24 MB of float64.
BLOCK = 5000

PREFIX = "v11_log2_"


def _curve_geometry(y0, A, x0, w, lo, hi):
    """
    Reconstruct the fitted curve's geometry over ``[lo, hi]`` for many genes.

    Returns ``(baseline, maxheight, rightheight, leftheight)``, each an array.
    The baseline is the curve's minimum over the interval -- the same quantity
    ``BhuvanFitter._curve_baseline`` computes, on the same 600-point grid -- and
    both edge heights are measured above it, exactly as the class does.
    """
    # np.linspace with array endpoints, NOT lo + (hi-lo)*t. The two are
    # mathematically identical and numerically are not, and BhuvanFitter uses
    # np.linspace. For a collapsed fit (w ~ 1e-6) a one-ulp shift in x moves the
    # curve height by percent, so the grid has to be built the same way to
    # reproduce the stored heights.
    x = np.linspace(lo, hi, CURVE_GRID, axis=-1)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        curve = y0[:, None] + A[:, None] * np.exp(
            -((x - x0[:, None]) / w[:, None]) ** 2)
        baseline = np.nanmin(curve, axis=1)
        top = np.nanmax(curve, axis=1)
        f_hi = _fourparam_gaussian(hi, y0, A, x0, w)
        f_lo = _fourparam_gaussian(lo, y0, A, x0, w)
    return baseline, top - baseline, f_hi - baseline, f_lo - baseline


def _ratio(height, maxheight):
    """Edge height over full height, clamped to [0, 1]; NaN when the curve is flat."""
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(maxheight == 0, np.nan, height / maxheight)
    return np.clip(r, 0.0, 1.0)


def _fmt(values, usable):
    """Shortest round-tripping repr per value; empty string where unusable."""
    out = np.full(values.shape, "", dtype=object)
    ok = usable & np.isfinite(values)
    idx = np.flatnonzero(ok)
    out[idx] = [repr(float(v)) for v in values[idx]]
    return out


def migrate_table(src, dest, verify=False):
    """
    Rewrite one table into the RTI/LTI schema. Returns a stats dict.

    ``src`` may already be migrated, in which case nothing is written and the
    returned dict carries ``skipped=True`` -- the script is re-runnable.
    """
    src, dest = Path(src), Path(dest)
    df = pd.read_csv(src, dtype=str, keep_default_na=False)

    if "rti" in df.columns:
        return {"skipped": True, "rows": len(df), "verified": 0,
                "max_rel_err": 0.0}
    missing = [c for c in RENAMES if c not in df.columns]
    if missing:
        raise SystemExit(f"{src.name}: not an old-schema table (missing {missing})")

    df = df.rename(columns=RENAMES)

    num = {c: pd.to_numeric(df[c], errors="coerce").to_numpy(dtype=float)
           for c in ("y0", "A", "x0", "w", "min", "max", "maxheight",
                     "rightheight")}
    # A row is usable when the fit parameters and the interval are all real and
    # the interval is non-degenerate. Failed-fit rows have blank parameters and
    # fall out here, so they receive empty strings rather than a number derived
    # from NaN.
    usable = (np.isfinite(num["y0"]) & np.isfinite(num["A"])
              & np.isfinite(num["x0"]) & np.isfinite(num["w"])
              & np.isfinite(num["min"]) & np.isfinite(num["max"])
              & (num["w"] != 0) & (num["max"] > num["min"]))

    n = len(df)
    baseline = np.full(n, np.nan)
    maxheight = np.full(n, np.nan)
    rightheight = np.full(n, np.nan)
    leftheight = np.full(n, np.nan)
    idx = np.flatnonzero(usable)
    for s in range(0, idx.size, BLOCK):
        b = idx[s:s + BLOCK]
        baseline[b], maxheight[b], rightheight[b], leftheight[b] = _curve_geometry(
            num["y0"][b], num["A"][b], num["x0"][b], num["w"][b],
            num["min"][b], num["max"][b])

    sigma = num["w"] / np.sqrt(2.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        lti_sigma_dist = (num["x0"] - num["min"]) / sigma
    lti = _ratio(leftheight, maxheight)

    stats = {"skipped": False, "rows": n, "verified": 0, "max_rel_err": 0.0}
    if verify:
        # The reconstruction must reproduce the two height columns already in
        # the table. If it does not, the curve being measured is not the curve
        # that was fit and every new value here would be quietly wrong.
        # The error is measured as a fraction of the CURVE's height, not of the
        # value being checked. A ceiling far out in the tail leaves
        # `rightheight` at ~5e-17 on a curve ~1e3 tall -- a floating-point zero,
        # which both sides agree on. Dividing by that value would report a 100%
        # error and condemn a migration that is exactly right.
        for name, got in (("maxheight", maxheight), ("rightheight", rightheight)):
            stored = num[name]
            cmp = usable & np.isfinite(stored) & np.isfinite(got)
            scale = np.maximum(np.abs(maxheight[cmp]), 1e-12)
            err = np.abs(got[cmp] - stored[cmp]) / scale
            if err.size:
                stats["max_rel_err"] = max(stats["max_rel_err"], float(err.max()))
            if name == "maxheight":
                stats["verified"] = int(cmp.sum())

    df["lti_sigma_dist"] = _fmt(lti_sigma_dist, usable)
    df["lti"] = _fmt(lti, usable)
    df["leftheight"] = _fmt(leftheight, usable)
    # `left` is the x_min the metrics were measured against, the mirror of
    # `right`. It is the `min` column, copied as TEXT so it stays byte-identical,
    # and blank on rows with no fit -- exactly how `right` behaves.
    df["left"] = np.where(usable, df["min"].to_numpy(dtype=object), "")

    # Only 4 of the 108 tables carry hist/hist_max, so the schema is applied as
    # an ordering over the columns this table actually has, not as a requirement.
    ordered = [c for c in NEW_COLUMNS if c in df.columns]
    extras = [c for c in df.columns if c not in NEW_COLUMNS]
    out = df[ordered + extras]
    dest.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(dest, index=False, lineterminator="\n")
    return stats


def discover(outputs, tissues=None):
    paths = sorted(outputs.glob(f"{PREFIX}*_fourparam*.csv"))
    if tissues:
        want = {t.strip() for t in tissues}
        paths = [p for p in paths
                 if any(p.name[len(PREFIX):].startswith(t) for t in want)]
    return paths


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all", action="store_true", help="every table in outputs/")
    ap.add_argument("--tissues", help="comma-separated tissue names")
    ap.add_argument("--input", type=Path)
    ap.add_argument("--output", type=Path)
    ap.add_argument("--outputs", type=Path,
                    default=Path(__file__).resolve().parents[1] / "outputs")
    ap.add_argument("--verify", action="store_true",
                    help="recompute maxheight/rightheight and compare to stored")
    args = ap.parse_args(argv)

    if args.input:
        jobs = [(args.input, args.output or args.input)]
    elif args.all or args.tissues:
        tissues = args.tissues.split(",") if args.tissues else None
        jobs = [(p, p) for p in discover(args.outputs, tissues)]
    else:
        ap.error("pass --all, --tissues or --input")

    if not jobs:
        print("nothing to do")
        return 0

    worst, done, skipped = 0.0, 0, 0
    for i, (src, dest) in enumerate(jobs, 1):
        st = migrate_table(src, dest, verify=args.verify)
        if st["skipped"]:
            skipped += 1
            print(f"[{i}/{len(jobs)}] {src.name}: already migrated, skipped")
            continue
        done += 1
        worst = max(worst, st["max_rel_err"])
        note = (f", verified {st['verified']} rows, max rel err {st['max_rel_err']:.2e}"
                if args.verify else "")
        print(f"[{i}/{len(jobs)}] {src.name}: {st['rows']} rows{note}", flush=True)

    print(f"\nmigrated {done}, skipped {skipped}")
    if args.verify:
        print(f"worst relative error against stored heights: {worst:.2e}")
        if worst > 1e-6:
            print("FAIL: reconstruction does not match the stored fit", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
