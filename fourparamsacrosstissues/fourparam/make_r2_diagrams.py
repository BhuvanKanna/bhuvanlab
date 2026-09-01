# -*- coding: utf-8 -*-
"""
make_r2_diagrams.py

One R^2 histogram per tissue -> ``r2_histograms/``, plus a single overview sheet
putting all 54 on shared axes.

**Excluded set only.** ``r2/`` holds only the ``<= -1`` tables, which is also the
default for every analysis in this project.

Why this exists: ``sumsquarevalue`` cannot be compared across genes -- it is
unnormalised, and within a tissue its distribution is very nearly the donor-count
distribution wearing an SSE costume (log10 SSE ~ 1.2 * log10 n_obs, R^2 0.83-0.87).
``r_squared`` divides that out. These sheets are how you check what it left behind.

Read them with the donor count in mind. R^2 is *not* comparable across tissues:
Spearman(donor count, median R^2) = 0.975, so kidney_medulla (n=11) sits near 0.21
and muscle_skeletal (n=818) near 0.91 without either fit being better or worse in
any way you should act on. The subtitle carries the donor count for that reason.

    cd fourparam
    python make_r2_diagrams.py
    python make_r2_diagrams.py --tissues kidney_cortex,liver --force
    python make_r2_diagrams.py --no-overview

Styling matches ``make_diagrams.py`` so the sheets read as siblings; the bar and
median-rule colours were checked for contrast (4.3 / 4.2 against white) and for
CVD separation (worst-case OKLab dE 17.0, target >= 8).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DEFAULT_R2 = HERE.parent / "r2"
DEFAULT_OUTPUTS = HERE.parent / "outputs"
DEFAULT_DEST = HERE.parent / "r2_histograms"

PREFIX = "v11_log2_"
R2_SUFFIX = "_r2_excluded_at_or_below_-1.csv"
TABLE_SUFFIX = "_fourparam_excluded_at_or_below_-1.csv"

# Shared with make_diagrams.py -- these sheets sit beside those.
INK = "#1F2328"
INK_SOFT = "#57606A"
INK_FAINT = "#8B949E"
BAR = "#4C7DB0"
BAR_EDGE = "#33587F"
RULE = "#C2603A"
SURFACE = "#FFFFFF"
GRID = "#E6E9EC"

BINS = 40                 # over [0, 1]; matches the 40-bin fit histogram
X_LIM = (0.0, 1.0)


def tissue_of(path: Path) -> str:
    return path.name[len(PREFIX):-len(R2_SUFFIX)]


def load(tissue: str, r2_dir: Path, outputs: Path):
    """R^2 for converged fits, plus the tissue's donor count.

    Joined against the fourparam table rather than read alone: a gene with
    ``fit_success = False`` has no curve, so its blank R^2 is not a data point
    and must not land in the histogram as a zero.
    """
    r2 = pd.read_csv(r2_dir / f"{PREFIX}{tissue}{R2_SUFFIX}")
    tab = pd.read_csv(outputs / f"{PREFIX}{tissue}{TABLE_SUFFIX}",
                      usecols=["gene", "fit_success", "n_obs"])
    d = tab.merge(r2, on="gene")
    d = d[d["fit_success"] == True]                            # noqa: E712
    d["r_squared"] = pd.to_numeric(d["r_squared"], errors="coerce")
    d = d.replace([np.inf, -np.inf], np.nan).dropna(subset=["r_squared"])
    donors = int(d["n_obs"].max()) if len(d) else 0
    return d["r_squared"].to_numpy(), donors


def _frame(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
        ax.spines[side].set_linewidth(0.8)
    ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)


def _bars(ax, values, width_scale=0.86):
    """Histogram drawn as bars so adjacent fills keep a surface gap."""
    inside = values[(values >= X_LIM[0]) & (values <= X_LIM[1])]
    counts, edges = np.histogram(inside, bins=BINS, range=X_LIM)
    centres = 0.5 * (edges[:-1] + edges[1:])
    w = (edges[1] - edges[0]) * width_scale
    ax.bar(centres, counts, width=w, color=BAR, edgecolor=BAR_EDGE,
           linewidth=0.5, zorder=2)
    return counts, len(values) - len(inside)


def render_one(tissue: str, values, donors: int, dest: Path, dpi: int) -> str:
    # Text is placed in figure coordinates against an explicit axes rectangle
    # rather than as axes-relative offsets: the offsets collide with the tick
    # labels once tight_layout moves the axes, and the amount they move by
    # depends on how wide the y tick labels happen to be for that tissue.
    fig = plt.figure(figsize=(8.0, 4.9))
    fig.patch.set_facecolor(SURFACE)
    ax = fig.add_axes([0.093, 0.255, 0.885, 0.545])
    _frame(ax)

    counts, n_clipped = _bars(ax, values)
    med = float(np.median(values))
    ax.axvline(med, color=RULE, linewidth=1.6, zorder=3)
    # Label on whichever side of the rule has room, so it never runs off frame.
    right = med < 0.72
    ax.annotate(f"median {med:.3f}", xy=(med, counts.max()),
                xytext=(5 if right else -5, 0), textcoords="offset points",
                color=RULE, fontsize=8, fontweight="600", va="top",
                ha="left" if right else "right",
                bbox=dict(boxstyle="square,pad=0.2", facecolor=SURFACE,
                          edgecolor="none", alpha=0.85))

    ax.set_xlim(*X_LIM)
    ax.tick_params(colors=INK_SOFT, labelsize=8.5, length=3, width=0.8)
    ax.set_xlabel("r_squared", color=INK_SOFT, fontsize=8.5, labelpad=4)
    ax.set_ylabel("genes", color=INK_SOFT, fontsize=8.5, labelpad=4)

    fig.text(0.093, 0.935, tissue.replace("_", " "),
             color=INK, fontsize=13, fontweight="600", va="top")
    fig.text(0.093, 0.872,
             f"fit quality of the 4-parameter Gaussian · excluded set "
             f"(values ≤ −1 dropped) · {donors} donors",
             color=INK_FAINT, fontsize=8.5, va="top")

    q1, q3 = np.percentile(values, [25, 75])
    fig.text(0.093, 0.115,
             f"n = {len(values):,} converged fits     median {med:.3f}     "
             f"IQR {q1:.3f}–{q3:.3f}     below 0: {(values < 0).mean():.2%}",
             color=INK_SOFT, fontsize=8.5, va="top")
    note = ("R² = 1 − SSR/TSS over the 40 fit bins. Not comparable across "
            "tissues; median R² tracks donor count at ρ = 0.975.")
    if n_clipped:
        note += f"\n{n_clipped:,} value(s) outside [0, 1] not drawn."
    fig.text(0.093, 0.062, note, color=INK_FAINT, fontsize=7.5,
             va="top", linespacing=1.5)

    fig.savefig(dest, dpi=dpi, facecolor=SURFACE)
    plt.close(fig)
    return f"{dest.name}  n={len(values):,}  median={med:.3f}"


def render_overview(rows, dest: Path, dpi: int) -> str:
    """All tissues on shared axes, ordered by donor count.

    Sorted by donors rather than alphabetically because the gradient *is* the
    finding; alphabetical order hides it.
    """
    rows = sorted(rows, key=lambda r: r[2])
    ncol, nrow = 6, 9
    fig, axes = plt.subplots(nrow, ncol, figsize=(15.5, 16.5),
                             sharex=True, sharey=False)
    fig.patch.set_facecolor(SURFACE)

    for ax, (tissue, values, donors) in zip(axes.ravel(), rows):
        _frame(ax)
        counts, _ = _bars(ax, values, width_scale=0.95)
        med = float(np.median(values))
        ax.axvline(med, color=RULE, linewidth=1.2, zorder=3)
        ax.set_xlim(*X_LIM)
        # Headroom plus a corner opposite the median: the bulk sits under the
        # median, so the far corner is the one place guaranteed to be empty
        # across a set that ranges from mass-at-0 to mass-at-1.
        ax.set_ylim(0, counts.max() * 1.22)
        ax.set_yticks([])
        ax.set_xticks([0, 0.5, 1.0])
        ax.tick_params(colors=INK_SOFT, labelsize=7, length=3, width=0.8)
        ax.set_title(tissue.replace("_", " "), fontsize=8, color=INK,
                     pad=4, loc="left")
        left = med > 0.5
        ax.text(0.03 if left else 0.97, 0.96,
                f"n={donors}\nmed {med:.2f}",
                transform=ax.transAxes, ha="left" if left else "right",
                va="top", fontsize=6.5, color=INK_FAINT, linespacing=1.35)

    for ax in axes.ravel()[len(rows):]:
        ax.set_visible(False)

    fig.suptitle("Fit quality (R²) by tissue, ordered by donor count",
                 fontsize=15, color=INK, x=0.055, ha="left", y=0.985)
    fig.text(0.055, 0.967,
             "GTEx v11, excluded set (values ≤ −1 dropped before fitting) · "
             "converged fits only · shared x-axis, independent y · "
             "rust rule = median",
             fontsize=9, color=INK_FAINT, ha="left")
    fig.text(0.5, 0.011,
             "r_squared:  1 − SSR/TSS over the 40 fit bins. Rises with donor "
             "count (ρ = 0.975), so read down the panels, not across.",
             fontsize=10, color=INK_SOFT, ha="center")

    fig.tight_layout(rect=[0.012, 0.022, 0.995, 0.958])
    fig.savefig(dest, dpi=dpi, facecolor=SURFACE)
    plt.close(fig)
    return f"{dest.name}  {len(rows)} tissues"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Render an R^2 histogram per tissue (excluded set).",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--r2", type=Path, default=DEFAULT_R2)
    ap.add_argument("--outputs", type=Path, default=DEFAULT_OUTPUTS)
    ap.add_argument("--dest", type=Path, default=DEFAULT_DEST)
    ap.add_argument("--tissues", default=None,
                    help="Comma-separated tissue names (default: all).")
    ap.add_argument("--dpi", type=int, default=130)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--no-overview", action="store_true",
                    help="Skip the all-tissue sheet.")
    args = ap.parse_args(argv)

    if not args.r2.is_dir():
        ap.error(f"r2 directory not found: {args.r2} -- run compute_r2.py first")

    tissues = sorted(tissue_of(p) for p in args.r2.glob(f"{PREFIX}*{R2_SUFFIX}"))
    if args.tissues:
        wanted = {t.strip().lower() for t in args.tissues.split(",")}
        tissues = [t for t in tissues if t.lower() in wanted]
    if not tissues:
        ap.error(f"no r2 tables matching {PREFIX}*{R2_SUFFIX} in {args.r2}")

    args.dest.mkdir(parents=True, exist_ok=True)
    print(f"Tissues found : {len(tissues)}  -> {args.dest}", file=sys.stderr)

    rows, made, skipped = [], 0, 0
    for i, tissue in enumerate(tissues, 1):
        values, donors = load(tissue, args.r2, args.outputs)
        if len(values) == 0:
            print(f"  [{i:>3}/{len(tissues)}] SKIP {tissue}: no converged fits",
                  file=sys.stderr)
            continue
        rows.append((tissue, values, donors))
        dest = args.dest / f"{PREFIX}{tissue}_r2_hist.png"
        if dest.exists() and not args.force:
            skipped += 1
            continue
        print(f"  [{i:>3}/{len(tissues)}] ok   "
              f"{render_one(tissue, values, donors, dest, args.dpi)}",
              file=sys.stderr)
        made += 1

    print(f"\nRendered {made} sheet(s)"
          + (f", skipped {skipped} existing (use --force)" if skipped else ""),
          file=sys.stderr)

    if not args.no_overview and rows:
        dest = args.dest / "_overview_r2_by_tissue.png"
        if dest.exists() and not args.force:
            print(f"Overview     : exists, skipping ({dest.name})", file=sys.stderr)
        else:
            print(f"Overview     : {render_overview(rows, dest, args.dpi)}",
                  file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
