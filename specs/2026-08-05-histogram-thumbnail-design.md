# Real histogram thumbnails in the Truncation Browser

**Date:** 2026-08-05
**Status:** approved, not yet implemented
**Touches:** `fourparamsacrosstissues/fourparam/`, `docs/index.html`

## Problem

The `Shape` column in `docs/index.html` draws a *synthetic* miniature: the fitted
4-parameter Gaussian evaluated over a window of the peak's own ±3.4σ, with the
truncation ceiling marked. It shows the model and never the data, and when the
fit is degenerate the picture actively misleads.

The motivating case is **APP in kidney cortex** (raw table):

```
y0 = -358.83   A = 362.35   x0 = 6.6887   w = 22.666   (sigma = 16.03)
min = 4.7252   max = right = 8.6772       n_obs = 104
ti_fourparam_sigma_dist = 0.124           truncationindex = 0.0
```

`truncationindex` is exactly 0 because the height-ratio baseline is the curve's
minimum over `[min, max]`, and with the peak a hair left of that interval's
midpoint the curve's minimum over the interval *is* `f(x_max)` — so
`rightheight = f(x_max) - min(f) = 0` by construction. (96.4% of converged
kidney-cortex rows are exactly 0 for this reason.)

Meanwhile the sparkline picks its window as `[x0 - 3.4*sigma, x0 + 3.4*sigma]`
= **[-52.2, 65.5]**, over which the curve does fall to its baseline, so the
ceiling renders at 99.2% of full height. The table says 0, the picture says ~1,
and both are internally consistent — they normalise over windows that differ by
a factor of ~30.

The TI value is not the problem and is not being changed. The picture is.

## Goal

Replace the synthetic miniature with the gene's **real 40-bin histogram**, with
the fitted 4-parameter Gaussian overlaid on the same axes, so the fit and the
data can be compared at a glance.

## Constraint that shapes everything

The browser has **no access to raw expression values**. It fetches only the
fourparam CSVs. Locally, 50 of the 54 `data/*.csv.gz` matrices are 133-byte Git
LFS pointers (~5.92 GB of payload deliberately not downloaded); only `uterus`,
`vagina`, `whole_blood` and `thyroid` are real files on disk.

So the bin counts must be **precomputed and published**.

A second constraint fixes the encoding: the browser parses CSV with a plain
`line.split(",")` (`docs/index.html:1137`) with no quoted-field handling. Any new
column must contain no comma and no quote, or every column index downstream
shifts.

## Part 1 — Data

### Two new columns

Appended to `COLUMNS` in `generate_fourparam.py` (currently
`generate_fourparam.py:85`), after `fit_success`:

| column | content |
|---|---|
| `hist` | 40 characters, one per bin, each drawn from `A-Z a-z 0-9 + /` representing 0-63. Bin value is `round(count * 63 / hist_max)`. Fixed width, always exactly 40 chars. |
| `hist_max` | the true count in the tallest bin (integer) |

Together these reconstruct the histogram to within 1/63 of peak height — under
one pixel at 44px tall.

**`hist` is a rendering aid, not analysis data.** It is quantised and lossy. This
must be stated in `fourparamsacrosstissues/CLAUDE.md` so nobody downstream reads
it as a source of truth for counts. `n_obs` and `hist_max` remain exact.

The alphabet is the standard base64 set, chosen because it is comma-free,
quote-free, and one character per bin — not because anything is base64-encoded.
Decoding in JS is `ALPHABET.indexOf(ch)` per character.

### Bin edges are derived, not stored

`bhuvanfitter.py:155` calls `np.histogram(arr, bins=40)` with no explicit range,
so edges are `linspace(arr.min(), arr.max(), 41)` — and `min` / `max` are already
columns.

One numpy rule must be mirrored client-side: **when `max == min`, numpy widens
the range to `[min - 0.5, max + 0.5]`.** This is not a rare corner — 8,264 rows
in kidney cortex alone are genes sitting entirely at −1. Missing it renders those
rows as a divide-by-zero.

### Written whenever there is data

`hist` and `hist_max` are written for any row with `n_obs >= 1`, **independent of
`fit_success`**. A gene whose fit failed still shows its real distribution, which
is arguably the case where seeing the data matters most. Rows with `n_obs = 0`
get empty strings for both.

### Cost

- ~46 bytes/row × 8,059,824 rows ≈ **370 MB per mirror**, so ~740 MB across
  `outputs/` and `gene_major/`. The encoding is high-entropy and will not
  compress much in git, so the push grows by roughly that.
- A one-time `git lfs pull` of ~5.92 GB for the 50 missing matrices.
  `run_cluster.py` already fetches one tissue, writes its tables, and deletes its
  copy, so this never needs 5.9 GB resident at once.

### What propagates for free

- `build_gui_data.py:155` derives `manifest.columns` from a reference table, so
  the browser learns the new columns automatically.
- `build_gene_major.py` concatenates row text verbatim, so shards carry them
  without change.
- The working-set CSV export follows `manifest.columns` order, so the
  byte-identity guarantee with `extract_genes.py` holds untouched.

## Part 2 — The cell

### Geometry

Cell grows from 116×34 to **160×44**.

```
lo   = edges[0]  = min          (or min - 0.5 when max == min)
hi   = edges[40] = max          (or max + 0.5 when max == min)
span = hi - lo
x-domain = [lo, hi + 0.12 * span]        // pad right only
y-domain = [0, 1.15 * hist_max]          // count units
```

The drawn domain is `1.12 * span`, so the bars occupy the first `1/1.12` = 89.3%
of the width: 142.9px / 40 bins = **3.57px per bin**. The tallest bar reaches
`1/1.15` = 87% of cell height, leaving headroom above it.

**Padding is right-only, deliberately.** Symmetric padding would read as a tidy
plot; the asymmetry is what says the cap is on the right, which is the entire
subject of this project.

### The curve shares the histogram's y-axis

`f` was fit to bin counts, so it is drawn directly in count units. That is what
makes the overlay meaningful rather than decorative.

It is **clipped, not rescaled**. The histogram is the observation and must always
stay legible; a degenerate fit that shoots off the top or sinks below zero leaves
the frame, and that exit is itself the signal. Implemented as
`.spark { overflow: hidden }` (currently `visible` at `docs/index.html:298`) —
no per-cell `clipPath` IDs to collide across up to 800 rendered rows.

For APP in kidney cortex this renders as bars showing the real right-skewed
pile-up, with the fitted curve a nearly flat teal line across the bottom third.

### Marks, back to front

| | mark | style |
|---|---|---|
| 1 | baseline rule | `--rule`, 1px |
| 2 | histogram | one `<path>` step outline; `--muted` fill @ .30 (dark theme .42); same path stroked 0.8px @ .55 so each bin's step edge stays crisp |
| 3 | fitted curve, observed range | `--curve`, solid 1.6px |
| 4 | fitted curve, past the ceiling | `--curve`, dashed, opacity .32 |
| 5 | ceiling | `--ceiling`, 1.4px, full height |

Mark 4 is the one part of the old sparkline worth keeping — the missing tail. It
now sits over genuinely empty space rather than an invented ±3.4σ window.

The existing translucent `.fill-live` area is removed; the bars are the fill now.

**Colour rationale:** the data is quiet grey ground, the model is the only
saturated thing in the cell, so a fit that does not match its data is immediately
visible. `--ceiling` red stays exclusive to the truncation wall, which is the
page's strongest single signal.

**One path, not 40 rects.** The comment at `docs/index.html:615` already treats
node count as load-bearing; 800 rows × 40 rects would add ~32k nodes. A single
stepped path keeps it at ~5 nodes per cell and gives sharper bin edges than
gapped rectangles.

### Fallbacks

| condition | behaviour |
|---|---|
| `hist` empty, fit params present | draw the curve alone, no bars (covers `n_obs = 0`, and any stale table cached before regeneration). Same x-domain from `min`/`max`; with no `hist_max` the y-domain falls back to the curve's own min/max over that domain |
| `hist_max` absent or 0 | treat as no histogram |
| neither hist nor fit params | empty cell, as today |

### Other browser changes

- `title` attribute on the `<td>` carrying `n_obs`, peak count and the x-range,
  so hovering gives the numbers without building a second UI surface.
- Legend at `docs/index.html:547` gains a bar swatch.
- `.c-shape` width and row height adjust for 160×44.

## Testing and rollout

Ordered so that nothing expensive happens until the cheap things have proven out.

1. **Python round-trip test.** Encode → decode → assert error ≤ 1/63 of peak, over
   real distributions plus the degenerate cases: all values identical, a single
   observation, `n_obs = 0`.
2. **Build the 8 tables for `uterus`, `vagina`, `whole_blood`, `thyroid` first.**
   Those four matrices are already real files on disk — no LFS pull needed. Point
   the browser at them and look at the actual cells.
3. **Static fixture page** rendering four known-hard rows side by side: APP kidney
   cortex (degenerate wide fit), a clean unimodal gene, an all-at-−1 gene, and a
   `fit_success = False` gene.
4. **`build_gene_major.py --verify`** — already byte-compares shards against
   `outputs/`.
5. **Only then** pull the remaining 50 matrices and run the full regeneration.

If 1–4 look right, 5 is mechanical. If they don't, no bandwidth has been spent.

## Explicitly out of scope

- **The `truncationindex` definition does not change.** TI = 0 is the signal being
  looked for; the complaint was about the picture, not the metric.
- No new TI columns, no `ti_valid` flag, no filtering of structurally-zero rows.
- No click-to-enlarge panel. The inline cell at 160×44 is the whole feature.
- `make_diagrams.py` and the `diagrams/` sheets are untouched.
