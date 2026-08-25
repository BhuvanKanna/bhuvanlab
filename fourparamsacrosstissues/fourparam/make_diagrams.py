#!/usr/bin/env python3
"""Render a 5-histogram summary sheet for every fourparam table.

One PNG per table in ``outputs/``: raw sheets to ``diagrams/`` and excluded
(<= -1) sheets to ``excluded_diagrams/``, so 54 tissues x 2 filters =
108 images, each holding five histograms over the ~74,628 genes in that table:

    truncationindex   sumsquarevalue   mean   std   ti_fourparam_sigma_dist

The image file name mirrors the table name exactly (``.csv`` -> ``.png``), so a
sheet and the table it came from sort next to each other.

Only rows with ``fit_success == True`` and a finite value are histogrammed; the
counts that were dropped are printed on the sheet rather than silently omitted.

Axis choices, which are per-metric on purpose
---------------------------------------------
These five columns do not share a shape, so one axis rule cannot serve them all.

``truncationindex``
    Bounded [0, 1] by construction, and violently zero-inflated - in liver
    roughly 60% of genes sit at exactly 0 and p99 is under 0.09. Plotted on its
    full fixed [0, 1] domain with a **log count axis**, because on a linear one
    the entire panel is a single bar at zero and the informative tail is
    invisible.

``ti_fourparam_sigma_dist``
    Spans about six orders of magnitude either side of zero: degenerate fits
    with w ~ 1e-4 push it past 1e5 while the biologically meaningful range is
    single digits. No linear window shows the bulk without hiding ~10% of the
    genes, so this one gets a **symlog x axis** (linear within +/-1, log beyond)
    with matching symlog-spaced bins. Nothing is clipped.

``sumsquarevalue``, ``mean``, ``std``
    Heavy right tails but no sign change and no extreme dynamic range. Each gets
    an adaptive linear window: the narrower of the 6x-IQR fences and the
    p0.5-p99.5 span that hides at most ``MAX_HIDDEN`` of the rows, falling back
    to whichever hides less. The hidden count is stamped on the panel.

Usage::

    python make_diagrams.py                    # all 108, 8 processes
    python make_diagrams.py --tissues liver,lung
    python make_diagrams.py --force            # redo sheets that already exist
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter

HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUTS = HERE.parent / "outputs"
DEFAULT_DIAGRAMS = HERE.parent / "diagrams"
# The two filters are answers to different questions and are almost never read
# side by side, so they get separate directories rather than one folder of 108
# sheets in which the pairs interleave alphabetically.
DEFAULT_EXCLUDED_DIAGRAMS = HERE.parent / "excluded_diagrams"

PREFIX = "v11_log2_"
RAW_SUFFIX = "_fourparam.csv"
EXCLUDED_SUFFIX = "_fourparam_excluded_at_or_below_-1.csv"

# (column, panel title, one-line gloss)
METRICS = [
    ("truncationindex", "Truncation index",
     "f(x_max) / f(peak), baseline removed - higher = more truncated"),
    ("sumsquarevalue", "Sum of squares",
     "residual sum of squares of the 4-parameter fit - lower = better fit"),
    ("mean", "Mean",
     "mean expression across donors, log2(TPM + 1) - 1"),
    ("std", "Standard deviation",
     "spread of expression across donors"),
    ("ti_fourparam_sigma_dist", "Sigma distance",
     "(x_max - x0) / sigma - how many sigma the ceiling sits above the peak"),
]

BINS = 60
MAX_HIDDEN = 0.02        # a linear window may hide at most 2% of the rows
IQR_K = 6.0              # fence width, in IQRs, for the adaptive window
SYMLOG_LINTHRESH = 1.0   # sigma-distance is linear within +/-1, log beyond

INK = "#1F2328"
INK_SOFT = "#57606A"
INK_FAINT = "#8B949E"
BAR = "#4C7DB0"
BAR_EDGE = "#33587F"
RULE = "#C2603A"         # median marker; a line against a fill, plus labelled
SURFACE = "#FFFFFF"
GRID = "#E6E9EC"


# --------------------------------------------------------------------------- #
# windowing
# --------------------------------------------------------------------------- #
def _hidden(values: np.ndarray, lo: float, hi: float) -> float:
    if values.size == 0:
        return 0.0
    return float(((values < lo) | (values > hi)).mean())


def adaptive_window(values: np.ndarray) -> tuple[float, float, int]:
    """A linear x window that shows the bulk without discarding much.

    Considers the 6x-IQR fences and the p0.5-p99.5 span, prefers the narrower of
    the two that hides <= MAX_HIDDEN, and otherwise takes whichever hides less.
    Returns (lo, hi, n_hidden).
    """
    lo_d, hi_d = float(values.min()), float(values.max())
    if not np.isfinite(lo_d) or not np.isfinite(hi_d) or lo_d == hi_d:
        pad = max(abs(lo_d) * 0.05, 1e-6)
        return lo_d - pad, hi_d + pad, 0

    q1, q3 = np.percentile(values, [25, 75])
    iqr = q3 - q1
    candidates = []
    if iqr > 0:
        candidates.append((max(q1 - IQR_K * iqr, lo_d), min(q3 + IQR_K * iqr, hi_d)))
    p_lo, p_hi = np.percentile(values, [0.5, 99.5])
    candidates.append((max(float(p_lo), lo_d), min(float(p_hi), hi_d)))
    candidates.append((lo_d, hi_d))

    scored = [(lo, hi, _hidden(values, lo, hi)) for lo, hi in candidates if hi > lo]
    if not scored:
        return lo_d, hi_d, 0

    ok = [c for c in scored if c[2] <= MAX_HIDDEN]
    lo, hi, frac = min(ok, key=lambda c: c[1] - c[0]) if ok \
        else min(scored, key=lambda c: c[2])
    return lo, hi, int(round(frac * values.size))


def symlog_edges(values: np.ndarray, n_bins: int, linthresh: float) -> np.ndarray:
    """Bin edges that are evenly spaced on a symlog axis, so bars look uniform."""
    lo, hi = float(values.min()), float(values.max())

    def fwd(x):      # value -> symlog position
        s = np.sign(x)
        a = np.abs(x)
        return np.where(a <= linthresh, x / linthresh,
                        s * (1.0 + np.log10(np.maximum(a, linthresh) / linthresh)))

    def inv(p):      # symlog position -> value
        s = np.sign(p)
        a = np.abs(p)
        return np.where(a <= 1.0, p * linthresh,
                        s * linthresh * 10.0 ** (a - 1.0))

    p_lo, p_hi = float(fwd(np.array([lo]))[0]), float(fwd(np.array([hi]))[0])
    if p_hi <= p_lo:
        p_hi = p_lo + 1.0
    edges = inv(np.linspace(p_lo, p_hi, n_bins + 1))
    edges[0], edges[-1] = lo, hi
    return np.unique(edges)


def _count_label(v: float, _pos=None) -> str:
    """Y tick label. Keeps a decimal on the k suffix so 1500/2000/2500 do not
    all collapse to the same '2k'."""
    if v <= 0:
        return ""
    if v >= 1000:
        return f"{v / 1000:.1f}".rstrip("0").rstrip(".") + "k"
    return f"{v:g}"


def _fmt(x: float) -> str:
    """Compact number for annotations: 1234.5 -> 1.23e3, 0.0821 -> 0.0821."""
    if not np.isfinite(x):
        return "n/a"
    if x == 0:
        return "0"
    a = abs(x)
    if a >= 1e4 or a < 1e-3:
        return f"{x:.2e}"
    if a >= 100:
        return f"{x:,.1f}"
    return f"{x:.4g}"


# --------------------------------------------------------------------------- #
# one panel
# --------------------------------------------------------------------------- #
def draw_panel(ax, values: np.ndarray, column: str, title: str, gloss: str,
               n_total: int) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
        ax.spines[side].set_linewidth(1.0)
    ax.tick_params(colors=INK_SOFT, labelsize=8.5, length=3, width=0.8)
    ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)

    ax.text(0.0, 1.115, title, transform=ax.transAxes, color=INK,
            fontsize=12.5, fontweight="600", va="bottom", ha="left")
    ax.text(0.0, 1.038, gloss, transform=ax.transAxes, color=INK_FAINT,
            fontsize=7.8, va="bottom", ha="left")
    ax.set_xlabel(column, color=INK_SOFT, fontsize=8.5, labelpad=5)
    ax.set_ylabel("genes", color=INK_SOFT, fontsize=8.5, labelpad=4)

    if values.size == 0:
        ax.text(0.5, 0.5, "no rows with a finite value", transform=ax.transAxes,
                ha="center", va="center", color=INK_FAINT, fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        return

    note = ""
    log_y = False

    if column == "truncationindex":
        # Bounded [0, 1] and zero-inflated: fixed domain, log counts.
        edges = np.linspace(0.0, 1.0, BINS + 1)
        lo, hi = 0.0, 1.0
        hidden = int(((values < 0.0) | (values > 1.0)).sum())
        log_y = True
        at_zero = int((values == 0).sum())
        note = f"{at_zero:,} exactly 0 ({at_zero / values.size:.0%}) - log counts"
        if hidden:
            note += f" - {hidden:,} outside [0, 1]"
    elif column == "ti_fourparam_sigma_dist":
        # Six orders of magnitude either side of zero: symlog, nothing clipped.
        edges = symlog_edges(values, BINS, SYMLOG_LINTHRESH)
        lo, hi = float(values.min()), float(values.max())
        hidden = 0
        log_y = True
        note = (f"symlog x (linear within +/-{SYMLOG_LINTHRESH:g}) - "
                f"full range, nothing clipped - log counts")
    else:
        lo, hi, hidden = adaptive_window(values)
        edges = np.linspace(lo, hi, BINS + 1)
        note = (f"x clipped to [{_fmt(lo)}, {_fmt(hi)}] - {hidden:,} outside "
                f"({hidden / values.size:.1%})") if hidden else \
               f"x spans the full range [{_fmt(lo)}, {_fmt(hi)}]"

    inside = values[(values >= lo) & (values <= hi)]
    ax.hist(inside, bins=edges, color=BAR, edgecolor=BAR_EDGE,
            linewidth=0.25, zorder=2)

    if column == "ti_fourparam_sigma_dist":
        ax.set_xscale("symlog", linthresh=SYMLOG_LINTHRESH)
    if log_y:
        ax.set_yscale("log")
    ax.set_xlim(lo, hi)

    med = float(np.median(values))
    if lo <= med <= hi:
        ax.axvline(med, color=RULE, linewidth=1.6, zorder=3)
        # Flip the label to the left of the rule when the rule sits in the
        # right-hand third, so it never runs off the panel.
        frac = ax.transLimits.transform((med, 0))[0] if lo < hi else 0.0
        right = frac > 0.66
        ax.annotate(f"median {_fmt(med)}", xy=(med, 1.0),
                    xycoords=("data", "axes fraction"),
                    xytext=(-5 if right else 5, -6), textcoords="offset points",
                    color=RULE, fontsize=8, fontweight="600",
                    ha="right" if right else "left", va="top", zorder=4,
                    bbox=dict(boxstyle="square,pad=0.18", facecolor=SURFACE,
                              edgecolor="none", alpha=0.85))

    ax.yaxis.set_major_formatter(FuncFormatter(_count_label))

    stats = (f"n = {values.size:,} of {n_total:,}    "
             f"mean {_fmt(float(values.mean()))}    "
             f"sd {_fmt(float(values.std(ddof=1))) if values.size > 1 else 'n/a'}")
    ax.text(0.0, -0.235, stats, transform=ax.transAxes, color=INK_SOFT,
            fontsize=8, va="top", ha="left")
    ax.text(0.0, -0.315, note, transform=ax.transAxes, color=INK_FAINT,
            fontsize=7.5, va="top", ha="left")


def draw_info_panel(ax, tissue: str, kind_label: str, n_rows: int,
                    n_fit: int, per_metric: dict[str, int]) -> None:
    """The sixth cell of the 2x3 grid: what this sheet is and what it excluded."""
    ax.set_facecolor(SURFACE)
    ax.axis("off")

    lines = [
        ("Tissue", tissue.replace("_", " ")),
        ("Table", kind_label),
        ("Genes in table", f"{n_rows:,}"),
        ("fit_success = True", f"{n_fit:,}  ({n_fit / n_rows:.1%})" if n_rows else "0"),
    ]
    y = 0.94
    for label, value in lines:
        ax.text(0.0, y, label, color=INK_FAINT, fontsize=8.5, va="top")
        ax.text(0.0, y - 0.058, value, color=INK, fontsize=11.5,
                fontweight="600", va="top")
        y -= 0.155

    ax.text(0.0, y - 0.01, "Finite values histogrammed", color=INK_FAINT,
            fontsize=8.5, va="top")
    y -= 0.065
    for column, _, _ in METRICS:
        n = per_metric.get(column, 0)
        ax.text(0.0, y, column, color=INK_SOFT, fontsize=8, va="top")
        ax.text(1.0, y, f"{n:,}", color=INK_SOFT, fontsize=8, va="top", ha="right")
        y -= 0.052

    ax.text(0.0, y - 0.02,
            "Rows with fit_success = False, or a NaN/inf value for that\n"
            "metric, are excluded from the histogram above.",
            color=INK_FAINT, fontsize=7.5, va="top", linespacing=1.5)


# --------------------------------------------------------------------------- #
# one sheet
# --------------------------------------------------------------------------- #
def render_sheet(table: Path, dest: Path, tissue: str, kind_label: str,
                 dpi: int) -> str:
    df = pd.read_csv(table, usecols=["fit_success"] + [m[0] for m in METRICS],
                     low_memory=False)
    n_rows = len(df)
    ok = df["fit_success"].astype(str).str.lower().eq("true")
    n_fit = int(ok.sum())
    fitted = df[ok]

    series, per_metric = {}, {}
    for column, _, _ in METRICS:
        s = pd.to_numeric(fitted[column], errors="coerce").to_numpy(dtype=float)
        s = s[np.isfinite(s)]
        series[column] = s
        per_metric[column] = int(s.size)

    fig, axes = plt.subplots(2, 3, figsize=(17.5, 10.0), dpi=dpi,
                             facecolor=SURFACE)
    fig.subplots_adjust(left=0.052, right=0.985, top=0.838, bottom=0.088,
                        wspace=0.235, hspace=0.62)

    fig.text(0.052, 0.972, f"{tissue.replace('_', ' ')}",
             color=INK, fontsize=22, fontweight="700", va="top")
    fig.text(0.052, 0.928, f"{kind_label}   -   {table.name}",
             color=INK_SOFT, fontsize=10.5, va="top")
    fig.text(0.985, 0.972,
             "GTEx v11  -  4-parameter Gaussian fit  -  bhuvanlab",
             color=INK_FAINT, fontsize=9, va="top", ha="right")

    flat = axes.ravel()
    for ax, (column, title, gloss) in zip(flat, METRICS):
        draw_panel(ax, series[column], column, title, gloss, n_fit)
    draw_info_panel(flat[5], tissue, kind_label, n_rows, n_fit, per_metric)

    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, facecolor=SURFACE)
    plt.close(fig)
    return f"{dest.name}  ({n_fit:,}/{n_rows:,} fitted)"


def _job(args):
    table, dest, tissue, kind_label, dpi = args
    try:
        return True, render_sheet(table, dest, tissue, kind_label, dpi)
    except Exception as exc:                                   # noqa: BLE001
        return False, f"{table.name}: {type(exc).__name__}: {exc}"


# --------------------------------------------------------------------------- #
def discover(outputs: Path, tissues=None):
    found = []
    for path in sorted(outputs.glob(f"{PREFIX}*{EXCLUDED_SUFFIX}")):
        found.append((path.name[len(PREFIX):-len(EXCLUDED_SUFFIX)],
                      "Excluded <= -1 (drops zero-expression samples)", path))
    excluded_names = {p.name for _, _, p in found}
    for path in sorted(outputs.glob(f"{PREFIX}*{RAW_SUFFIX}")):
        if path.name in excluded_names:
            continue
        found.append((path.name[len(PREFIX):-len(RAW_SUFFIX)],
                      "Raw (every finite value)", path))
    if tissues:
        wanted = {t.strip().lower() for t in tissues}
        found = [f for f in found if f[0].lower() in wanted]
    return sorted(found, key=lambda f: f[2].name)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Render a 5-histogram summary sheet for every fourparam table.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outputs", type=Path, default=DEFAULT_OUTPUTS)
    ap.add_argument("--diagrams", type=Path, default=DEFAULT_DIAGRAMS,
                    help="Destination for raw-table sheets.")
    ap.add_argument("--excluded-diagrams", type=Path,
                    default=DEFAULT_EXCLUDED_DIAGRAMS,
                    help="Destination for excluded (<= -1) sheets.")
    ap.add_argument("--tissues", default=None,
                    help="Comma-separated tissue names (default: all).")
    ap.add_argument("--dpi", type=int, default=130)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--force", action="store_true",
                    help="Re-render sheets that already exist.")
    args = ap.parse_args(argv)

    if not args.outputs.is_dir():
        ap.error(f"outputs directory not found: {args.outputs}")

    tissue_list = args.tissues.split(",") if args.tissues else None
    tables = discover(args.outputs, tissue_list)
    if not tables:
        ap.error(f"no tables matching {PREFIX}*_fourparam*.csv in {args.outputs}")

    args.diagrams.mkdir(parents=True, exist_ok=True)
    args.excluded_diagrams.mkdir(parents=True, exist_ok=True)
    jobs, skipped = [], 0
    for tissue, kind_label, path in tables:
        # Route by the table's own suffix, not by kind_label -- the label is
        # display text and changing its wording must not silently relocate files.
        into = (args.excluded_diagrams if path.name.endswith(EXCLUDED_SUFFIX)
                else args.diagrams)
        dest = into / (path.stem + ".png")
        if dest.exists() and not args.force:
            skipped += 1
            continue
        jobs.append((path, dest, tissue, kind_label, args.dpi))

    n_tissues = len({t for t, _, _ in tables})
    print(f"Tables found : {len(tables)}  ({n_tissues} tissues x 2 filters)",
          file=sys.stderr)
    print(f"Already done : {skipped}" + ("  (use --force to redo)" if skipped else ""),
          file=sys.stderr)
    print(f"To render    : {len(jobs)}  -> {args.diagrams} / "
          f"{args.excluded_diagrams}", file=sys.stderr)
    if not jobs:
        return 0

    failures = []
    with mp.Pool(processes=max(1, min(args.jobs, len(jobs)))) as pool:
        for i, (ok, message) in enumerate(pool.imap_unordered(_job, jobs), 1):
            tag = "ok  " if ok else "FAIL"
            if not ok:
                failures.append(message)
            print(f"  [{i:>3}/{len(jobs)}] {tag} {message}", file=sys.stderr)

    print(f"\nRendered {len(jobs) - len(failures)}/{len(jobs)} sheets "
          f"({len(jobs) * 5 - len(failures) * 5} histograms) -> {args.diagrams}",
          file=sys.stderr)
    for message in failures:
        print(f"  FAILED {message}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
