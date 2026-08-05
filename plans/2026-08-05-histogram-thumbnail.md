# Real Histogram Thumbnails Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the synthetic fitted-curve miniature in the browser's `Shape` column with the gene's real 40-bin histogram, with the fitted 4-parameter Gaussian overlaid on the same count axis.

**Architecture:** Two new CSV columns (`hist`, `hist_max`) carry a quantised 40-bin histogram from the Python pipeline to the browser; `sparkline()` in `docs/index.html` decodes them and draws bars + curve + ceiling in one SVG. Bin edges are not stored — they are derived from the existing `min`/`max` columns.

**Tech Stack:** Python 3 (numpy, pandas, scipy, pytest), vanilla JS/SVG in a single self-contained HTML file.

**Spec:** `specs/2026-08-05-histogram-thumbnail-design.md`

## Global Constraints

- The browser parses CSV with a plain `line.split(",")` (`docs/index.html:1137`). **New columns must contain no comma and no double-quote.**
- `truncationindex` and every other existing metric are **unchanged**. This work touches presentation and adds columns; it changes no science.
- Nothing may re-serialise a float. `extract_genes.py` reads `dtype=str`, `build_gene_major.py` concatenates row text verbatim — the byte-identity guarantee between the browser's CSV export and `extract_genes.py` depends on it.
- All CSV writers pin `lineterminator="\n"`.
- Histogram is always **40 bins**, from `np.histogram(arr, bins=40)` with no explicit range.
- `hist` is a **rendering aid, quantised and lossy**. `n_obs` and `hist_max` stay exact.
- Do **not** run the full 108-table regeneration or `git lfs pull` as part of this plan. Task 3 covers only the four tissues already present as real files on disk.

---

### Task 1: `encode_histogram` in the fit library

**Files:**
- Modify: `fourparamsacrosstissues/fourparam/bhuvanfitter.py`
- Test: `fourparamsacrosstissues/fourparam/tests/test_histogram_encoding.py`

**Interfaces:**
- Consumes: nothing
- Produces: `encode_histogram(data, bins=40) -> tuple[str, int]`, plus module constants `HIST_BINS = 40`, `HIST_LEVELS = 63`, `HIST_ALPHABET` (64 chars)

- [ ] **Step 1: Write the failing test**

Create `fourparamsacrosstissues/fourparam/tests/test_histogram_encoding.py`:

```python
"""The thumbnail histogram encoding is lossy on purpose; these tests pin how
lossy, and pin the degenerate cases that would otherwise divide by zero."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bhuvanfitter import (HIST_ALPHABET, HIST_BINS, HIST_LEVELS,
                          encode_histogram)


def decode(s, hist_max):
    """What the browser does: char -> level -> count."""
    return np.array([HIST_ALPHABET.index(c) for c in s]) * hist_max / HIST_LEVELS


def test_alphabet_is_64_unique_csv_safe_chars():
    assert len(HIST_ALPHABET) == 64
    assert len(set(HIST_ALPHABET)) == 64
    assert "," not in HIST_ALPHABET
    assert '"' not in HIST_ALPHABET
    assert "\n" not in HIST_ALPHABET


def test_encoding_is_always_exactly_40_chars_from_the_alphabet():
    rng = np.random.default_rng(0)
    for _ in range(50):
        data = rng.normal(size=rng.integers(10, 500))
        s, hmax = encode_histogram(data)
        assert len(s) == HIST_BINS
        assert set(s) <= set(HIST_ALPHABET)
        assert hmax > 0


def test_tallest_bin_encodes_to_the_top_level():
    rng = np.random.default_rng(1)
    s, hmax = encode_histogram(rng.normal(size=400))
    assert HIST_ALPHABET.index(max(s, key=HIST_ALPHABET.index)) == HIST_LEVELS


def test_hist_max_is_the_true_peak_count():
    data = np.concatenate([np.zeros(37), np.linspace(1, 2, 9)])
    counts, _ = np.histogram(data, bins=HIST_BINS)
    _, hmax = encode_histogram(data)
    assert hmax == counts.max()


def test_round_trip_error_is_within_one_level_of_peak():
    rng = np.random.default_rng(2)
    for _ in range(50):
        data = rng.normal(size=rng.integers(20, 1000))
        counts, _ = np.histogram(data, bins=HIST_BINS)
        s, hmax = encode_histogram(data)
        assert np.max(np.abs(decode(s, hmax) - counts)) <= hmax / HIST_LEVELS


def test_all_values_identical_gives_one_spike():
    s, hmax = encode_histogram(np.full(64, -1.0))
    assert hmax == 64
    assert len(s) == HIST_BINS
    assert sum(1 for c in s if HIST_ALPHABET.index(c) > 0) == 1


def test_single_observation():
    s, hmax = encode_histogram([3.5])
    assert hmax == 1
    assert sum(1 for c in s if HIST_ALPHABET.index(c) > 0) == 1


def test_no_finite_data_returns_empty():
    assert encode_histogram([]) == ("", 0)
    assert encode_histogram([np.nan, np.inf, -np.inf]) == ("", 0)


def test_non_finite_values_are_dropped_not_counted():
    s_clean, max_clean = encode_histogram([1.0, 2.0, 3.0, 2.0])
    s_dirty, max_dirty = encode_histogram([1.0, 2.0, np.nan, 3.0, 2.0, np.inf])
    assert (s_clean, max_clean) == (s_dirty, max_dirty)


def test_matches_the_fitter_own_histogram():
    """encode_histogram must bin identically to BhuvanFitter, or the browser
    draws bars that do not line up with the curve fitted to them."""
    from bhuvanfitter import BhuvanFitter
    rng = np.random.default_rng(3)
    data = rng.normal(size=300)
    bf = BhuvanFitter(data, gene_name="X")
    s, hmax = encode_histogram(data)
    assert hmax == int(bf.hist_counts.max())
    assert np.allclose(decode(s, hmax), bf.hist_counts, atol=hmax / HIST_LEVELS)
```

- [ ] **Step 2: Run it and watch it fail**

```bash
cd fourparamsacrosstissues/fourparam && python -m pytest tests/test_histogram_encoding.py -q
```
Expected: collection error, `ImportError: cannot import name 'HIST_ALPHABET'`.

- [ ] **Step 3: Implement**

In `bhuvanfitter.py`, after the imports and before `_fourparam_gaussian`:

```python
# -- thumbnail histogram encoding ---------------------------------------------
# One character per bin, so the column is fixed width and contains no comma or
# quote: the browser parses these CSVs with a plain split(","), and any quoted
# field would shift every column index after it.
HIST_BINS = 40
HIST_LEVELS = 63
HIST_ALPHABET = ("ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                 "abcdefghijklmnopqrstuvwxyz"
                 "0123456789+/")


def encode_histogram(data, bins: int = HIST_BINS):
    """Quantise a gene's 40-bin histogram to ``(hist, hist_max)``.

    ``hist`` is exactly ``bins`` characters from ``HIST_ALPHABET``, each holding
    ``round(count * HIST_LEVELS / hist_max)``; ``hist_max`` is the true count in
    the tallest bin. ``("", 0)`` when there is no finite data.

    Binning matches ``BhuvanFitter`` exactly (``np.histogram(arr, bins=40)``,
    no explicit range) so the bars line up with the curve fitted to them.

    This is **a rendering aid, not analysis data** — it is lossy to within
    1/63 of the peak. ``n_obs`` and ``hist_max`` remain exact.
    """
    arr = np.asarray(data, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return "", 0
    counts, _ = np.histogram(arr, bins=bins)
    hist_max = int(counts.max())
    if hist_max <= 0:
        return "", 0
    levels = np.rint(counts * (HIST_LEVELS / hist_max)).astype(int)
    return "".join(HIST_ALPHABET[v] for v in levels), hist_max
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
cd fourparamsacrosstissues/fourparam && python -m pytest tests/test_histogram_encoding.py -q
```
Expected: 10 passed.

- [ ] **Step 5: Commit**

```bash
git add fourparamsacrosstissues/fourparam/bhuvanfitter.py fourparamsacrosstissues/fourparam/tests/
git commit -m "Add encode_histogram: 40-bin thumbnail encoding for the browser"
```

---

### Task 2: Emit `hist` / `hist_max` from the table generator

**Files:**
- Modify: `fourparamsacrosstissues/fourparam/generate_fourparam.py` (`COLUMNS` at :85, `_failed_row` at :132, `_fit_one` at :153)
- Test: `fourparamsacrosstissues/fourparam/tests/test_generate_columns.py`

**Interfaces:**
- Consumes: `encode_histogram(data) -> (str, int)` from Task 1
- Produces: `COLUMNS` ending `..., "n_obs", "fit_success", "hist", "hist_max"`; `_failed_row(gene, n_obs, hist="", hist_max=0)`

**Why `hist_max = 0` rather than NaN on empty rows:** a NaN anywhere in the column forces pandas to write `12.0` instead of `12` for every other row. `0` keeps it an integer column, and the browser already treats `hist_max = 0` as "no histogram".

- [ ] **Step 1: Write the failing test**

Create `fourparamsacrosstissues/fourparam/tests/test_generate_columns.py`:

```python
"""The histogram columns must survive the paths that produce a row, including
the two that never construct a BhuvanFitter."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import generate_fourparam as gf
from bhuvanfitter import HIST_BINS


def test_hist_columns_are_last_and_named():
    assert gf.COLUMNS[-2:] == ["hist", "hist_max"]


def test_failed_row_carries_the_histogram_it_was_given():
    row = gf._failed_row("G", 4, "A" * HIST_BINS, 7)
    assert row["fit_success"] is False
    assert row["hist"] == "A" * HIST_BINS
    assert row["hist_max"] == 7


def test_failed_row_defaults_to_empty_histogram():
    row = gf._failed_row("G", 0)
    assert row["hist"] == ""
    assert row["hist_max"] == 0


def test_below_min_obs_still_gets_a_histogram():
    """A gene that never reaches the fit still shows its real distribution."""
    gf._init_worker(None, 2000)
    row = gf._fit_one(("G", np.array([1.0, 2.0, 3.0])))
    assert row["fit_success"] is False
    assert row["n_obs"] == 3
    assert len(row["hist"]) == HIST_BINS
    assert row["hist_max"] == 1


def test_successful_fit_gets_a_histogram():
    gf._init_worker(None, 2000)
    rng = np.random.default_rng(4)
    row = gf._fit_one(("G", rng.normal(size=300)))
    assert row["fit_success"] is True
    assert len(row["hist"]) == HIST_BINS
    assert row["hist_max"] > 0


def test_zero_observations_gives_empty_histogram():
    gf._init_worker(None, 2000)
    row = gf._fit_one(("G", np.array([np.nan, np.nan])))
    assert row["n_obs"] == 0
    assert row["hist"] == ""
    assert row["hist_max"] == 0


def test_threshold_is_applied_before_binning():
    """The excluded table's histogram must describe the excluded data."""
    values = np.concatenate([np.full(50, -1.0), np.linspace(0, 3, 50)])
    gf._init_worker(None, 2000)
    raw = gf._fit_one(("G", values.copy()))
    gf._init_worker(-1.0, 2000)
    excl = gf._fit_one(("G", values.copy()))
    assert raw["hist"] != excl["hist"]
    assert raw["hist_max"] == 50
    assert excl["n_obs"] == 50


def test_hist_max_column_stays_integer_in_the_written_csv(tmp_path):
    """A NaN in the column would make pandas write 12.0 for every row."""
    gf._init_worker(None, 2000)
    rng = np.random.default_rng(5)
    rows = [gf._fit_one(("A", rng.normal(size=200))),
            gf._fit_one(("B", np.array([np.nan])))]
    table = pd.DataFrame.from_records(rows, columns=gf.COLUMNS)
    out = tmp_path / "t.csv"
    table.to_csv(out, index=False, lineterminator="\n")
    text = out.read_text()
    assert ".0," not in text.split("\n")[1].split(",")[-1] + ","
    for line in text.strip().split("\n")[1:]:
        assert "." not in line.split(",")[-1]


def test_no_commas_or_quotes_in_the_hist_column():
    gf._init_worker(None, 2000)
    rng = np.random.default_rng(6)
    row = gf._fit_one(("G", rng.normal(size=200)))
    assert "," not in row["hist"]
    assert '"' not in row["hist"]
```

- [ ] **Step 2: Run it and watch it fail**

```bash
cd fourparamsacrosstissues/fourparam && python -m pytest tests/test_generate_columns.py -q
```
Expected: `test_hist_columns_are_last_and_named` fails — `COLUMNS[-2:] == ["n_obs", "fit_success"]`.

- [ ] **Step 3: Implement**

In `generate_fourparam.py`, change the import to pull in the encoder:

```python
from bhuvanfitter import BhuvanFitter, encode_histogram
```

Extend `COLUMNS` (:85):

```python
COLUMNS = [
    "gene", "y0", "A", "x0", "w", "sumsquarevalue",
    "ti_fourparam_sigma_dist", "truncationindex",
    "min", "max", "mean", "std", "skew", "kurt",
    "right", "maxheight", "rightheight", "n_obs", "fit_success",
    # Thumbnail histogram for the browser's Shape column. Quantised and lossy —
    # a rendering aid, not analysis data. See bhuvanfitter.encode_histogram.
    "hist", "hist_max",
]
```

Replace `_failed_row` (:132):

```python
def _failed_row(gene: str, n_obs: int, hist: str = "", hist_max: int = 0) -> dict:
    row = {c: np.nan for c in COLUMNS}
    row["gene"] = gene
    row["n_obs"] = int(n_obs)
    row["fit_success"] = False
    # 0 rather than NaN: one NaN in the column makes pandas write "12.0"
    # instead of "12" for every other row.
    row["hist"] = hist
    row["hist_max"] = int(hist_max)
    return row
```

Replace `_fit_one` (:153):

```python
def _fit_one(item) -> dict:
    gene, values = item
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    if _WORKER_THRESHOLD is not None:
        data = data[data > _WORKER_THRESHOLD]   # drop values <= threshold
    n_obs = int(data.size)

    # Computed before and independently of the fit: a gene whose fit failed
    # still gets to show its real distribution in the browser.
    hist, hist_max = encode_histogram(data)

    if n_obs < MIN_OBS:
        return _failed_row(gene, n_obs, hist, hist_max)
    try:
        bf = BhuvanFitter(data, gene_name=gene)
        row = bf.fit("fourparam", max_nfev=_WORKER_MAX_NFEV)
    except (RuntimeError, ValueError):
        return _failed_row(gene, n_obs, hist, hist_max)
    row["hist"] = hist
    row["hist_max"] = hist_max
    return row
```

- [ ] **Step 4: Run both test files**

```bash
cd fourparamsacrosstissues/fourparam && python -m pytest tests/ -q
```
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add fourparamsacrosstissues/fourparam/generate_fourparam.py fourparamsacrosstissues/fourparam/tests/
git commit -m "Emit hist/hist_max columns, written even when the fit fails"
```

---

### Task 3: Regenerate the four tissues that are real files on disk

**Files:**
- Modify: `fourparamsacrosstissues/outputs/v11_log2_{uterus,vagina,whole_blood,thyroid}_fourparam[_excluded_at_or_below_-1].csv` (8 tables)

**Interfaces:**
- Consumes: the generator from Task 2
- Produces: 8 regenerated tables carrying `hist`/`hist_max`, for Task 4 to render against

`uterus` is the smallest of the four (52 MB) — do it first and check it before spending time on `thyroid` (236 MB). The other 50 matrices are LFS pointers and are **out of scope**.

- [ ] **Step 1: Confirm the four matrices are real files, not pointers**

```bash
cd fourparamsacrosstissues/data && python -c "
import glob
for f in sorted(glob.glob('*.csv.gz')):
    if not open(f,'rb').read(40).startswith(b'version https://git-lfs'):
        print('REAL', f)
"
```
Expected: exactly `uterus`, `vagina`, `whole_blood`, `thyroid`.

- [ ] **Step 2: Regenerate uterus, both tables**

```bash
cd fourparamsacrosstissues/fourparam
python generate_fourparam.py --input ../data/v11_log2_uterus.csv.gz \
    --id-col Name --name-col Description --jobs 8 --force
python generate_fourparam.py --input ../data/v11_log2_uterus.csv.gz \
    --id-col Name --name-col Description --threshold -1 --jobs 8 --force
```

(If `generate_fourparam.py` has no `--force`, delete the two target files first.)

- [ ] **Step 3: Verify the new columns landed correctly**

```bash
cd fourparamsacrosstissues/outputs && python -c "
import pandas as pd
d = pd.read_csv('v11_log2_uterus_fourparam.csv', dtype=str, keep_default_na=False)
print(list(d.columns)[-2:])
assert list(d.columns)[-2:] == ['hist','hist_max']
h = d.hist[d.hist != '']
assert (h.str.len() == 40).all(), 'ragged hist column'
assert not h.str.contains(',').any()
assert not d.hist_max.str.contains(r'\.').any(), 'hist_max went float'
print('rows', len(d), 'with hist', len(h))
print(d[d.genename=='APP'][['genename','n_obs','hist_max','hist']].to_string())
"
```
Expected: `['hist', 'hist_max']`, no assertion failure, an APP row with a 40-char `hist`.

- [ ] **Step 4: Regenerate the other three tissues**

```bash
cd fourparamsacrosstissues/fourparam
for t in vagina whole_blood thyroid; do
  python generate_fourparam.py --input ../data/v11_log2_$t.csv.gz \
      --id-col Name --name-col Description --jobs 8 --force
  python generate_fourparam.py --input ../data/v11_log2_$t.csv.gz \
      --id-col Name --name-col Description --threshold -1 --jobs 8 --force
done
```

- [ ] **Step 5: Re-run the Step 3 check against all eight tables, then commit**

```bash
git add fourparamsacrosstissues/outputs/
git commit -m "Regenerate uterus/vagina/whole_blood/thyroid tables with hist columns"
```

---

### Task 4: Draw the real histogram in the browser

**Files:**
- Modify: `docs/index.html` — CSS at :297-303, `sparkline()` at :959-1014, the `Shape` cell at :1297, the legend at :547

**Interfaces:**
- Consumes: `hist` (40 chars) and `hist_max` (int) columns from Task 2
- Produces: `decodeHistogram(s) -> number[]|null` (fractions of `hist_max`), `sparkline(row, colIdx) -> string` (SVG markup or `""`), `shapeCell(fields, colIdx) -> string`

- [ ] **Step 1: Replace the sparkline CSS block (:297-303)**

```css
/* the signature: per-gene histogram with its fitted curve overlaid.
   overflow:hidden is load-bearing — a degenerate fit is drawn clipped rather
   than rescaled, so the bars stay legible and the curve leaving the frame is
   itself the signal. Clipping this way needs no per-cell clipPath id, which
   would have to be unique across up to RENDER_LIMIT rows. */
.spark { display: block; width: 160px; height: 44px; overflow: hidden; }
.spark .bars {
  fill: var(--muted); fill-opacity: .30;
  stroke: var(--muted); stroke-opacity: .55; stroke-width: .8;
  stroke-linejoin: miter;
}
.spark .curve-live { fill: none; stroke: var(--curve); stroke-width: 1.6; }
.spark .curve-censored { fill: none; stroke: var(--curve); stroke-width: 1.6; stroke-dasharray: 2 2.5; opacity: .32; }
.spark .ceiling { stroke: var(--ceiling); stroke-width: 1.4; }
.spark .base { stroke: var(--rule); stroke-width: 1; }
@media (prefers-color-scheme: dark) { .spark .bars { fill-opacity: .42; } }
:root[data-theme="dark"] .spark .bars { fill-opacity: .42; }
:root[data-theme="light"] .spark .bars { fill-opacity: .30; }
```

The `.fill-live` rule is deleted — the bars are the fill now.

- [ ] **Step 2: Replace `sparkline()` (:956-1014) with the decoder plus the new renderer**

```js
/* =====================================================================
   The signature: this row's real 40-bin histogram, with its fitted
   4-parameter Gaussian overlaid on the same count axis.

   Bin edges are not transported: numpy bins with no explicit range, so they
   are linspace(min, max, 41) and both ends are already columns. The one rule
   that has to be mirrored is numpy's zero-width fallback, or every gene
   sitting entirely at -1 divides by zero.
   ===================================================================== */
const HIST_ALPHABET =
  "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
const HIST_BINS = 40;
const HIST_LEVELS = 63;
const SPARK_W = 160, SPARK_H = 44;
const SPARK_PAD_FRAC = 0.12;   // dead space past the ceiling, right side only
const SPARK_HEADROOM = 1.15;   // y-axis top, as a multiple of hist_max

function decodeHistogram(s) {
  if (typeof s !== "string" || s.length !== HIST_BINS) return null;
  const out = new Array(HIST_BINS);
  for (let i = 0; i < HIST_BINS; i++) {
    const v = HIST_ALPHABET.indexOf(s[i]);
    if (v < 0) return null;
    out[i] = v / HIST_LEVELS;          // fraction of hist_max
  }
  return out;
}

function sparkline(row, colIdx) {
  const num = (c) => {
    const i = colIdx[c];
    if (i === undefined) return NaN;
    const v = Number(row[i]);
    return Number.isFinite(v) ? v : NaN;
  };
  const txt = (c) => {
    const i = colIdx[c];
    return i === undefined || row[i] == null ? "" : String(row[i]).trim();
  };

  const y0 = num("y0"), A = num("A"), x0 = num("x0"), w = num("w");
  const lo0 = num("min"), hi0 = num("max"), histMax = num("hist_max");
  const bars = decodeHistogram(txt("hist"));

  const haveCurve = [y0, A, x0, w].every(Number.isFinite) && w > 0;
  const haveBars = !!bars && Number.isFinite(histMax) && histMax > 0 &&
                   Number.isFinite(lo0) && Number.isFinite(hi0);
  if (!haveBars && !haveCurve) return "";

  // -- x domain ------------------------------------------------------------
  let lo, hi;
  if (Number.isFinite(lo0) && Number.isFinite(hi0)) {
    if (hi0 > lo0) { lo = lo0; hi = hi0; }
    else { lo = lo0 - 0.5; hi = hi0 + 0.5; }   // numpy's zero-width fallback
  } else {
    const sigma = w / Math.SQRT2;
    lo = x0 - 3 * sigma; hi = x0 + 3 * sigma;
  }
  const a = lo, b = hi + SPARK_PAD_FRAC * (hi - lo);
  if (!(b > a)) return "";

  const W = SPARK_W, H = SPARK_H, base = H - 1.5, usable = H - 3;
  const sx = (x) => ((x - a) / (b - a)) * W;
  const f = (x) => y0 + A * Math.exp(-Math.pow((x - x0) / w, 2));

  // -- y domain, in count units -------------------------------------------
  let yLo = 0, yHi;
  if (haveBars) {
    yHi = SPARK_HEADROOM * histMax;
  } else {
    let mn = Infinity, mx = -Infinity;
    for (let i = 0; i <= 96; i++) {
      const v = f(a + (b - a) * i / 96);
      if (v < mn) mn = v;
      if (v > mx) mx = v;
    }
    yLo = mn; yHi = mx > mn ? mx : mn + 1;
  }
  const sy = (v) => base - ((v - yLo) / (yHi - yLo)) * usable;

  // -- bars: one stepped path, not 40 rects (node count, see RENDER_LIMIT) --
  let barsPath = "";
  if (haveBars) {
    const binW = (hi - lo) / HIST_BINS;
    let d = `M${sx(lo).toFixed(1)} ${base.toFixed(1)}`;
    for (let i = 0; i < HIST_BINS; i++) {
      const y = sy(bars[i] * histMax).toFixed(1);
      d += `L${sx(lo + i * binW).toFixed(1)} ${y}` +
           `L${sx(lo + (i + 1) * binW).toFixed(1)} ${y}`;
    }
    d += `L${sx(hi).toFixed(1)} ${base.toFixed(1)}Z`;
    barsPath = `<path class="bars" d="${d}"/>`;
  }

  // -- curve, split at the ceiling ----------------------------------------
  const cut = Number.isFinite(num("right")) ? num("right") : hi;
  let livePath = "", censPath = "";
  if (haveCurve) {
    const N = 96, live = [], cens = [];
    for (let i = 0; i <= N; i++) {
      const x = a + (b - a) * i / N;
      (x <= cut ? live : cens).push([sx(x), sy(f(x))]);
    }
    if (live.length && cens.length) {
      const j = [sx(cut), sy(f(cut))];
      live.push(j); cens.unshift(j);
    }
    const d = (pts) => pts.map((p, i) =>
      `${i ? "L" : "M"}${p[0].toFixed(1)} ${p[1].toFixed(1)}`).join("");
    if (cens.length) censPath = `<path class="curve-censored" d="${d(cens)}"/>`;
    if (live.length) livePath = `<path class="curve-live" d="${d(live)}"/>`;
  }

  const cx = sx(cut).toFixed(1);
  const nObs = num("n_obs");
  const label = haveBars
    ? `40-bin histogram, n = ${Number.isFinite(nObs) ? nObs : "?"}, ` +
      `peak bin ${histMax}, range ${lo.toFixed(2)} to ${hi.toFixed(2)}`
    : `fitted curve only, no histogram available`;

  return `<svg class="spark" role="img" viewBox="0 0 ${W} ${H}" preserveAspectRatio="none">
    <title>${escapeHtml(label)}</title>
    <line class="base" x1="0" y1="${base}" x2="${W}" y2="${base}"/>
    ${barsPath}
    ${censPath}
    ${livePath}
    <line class="ceiling" x1="${cx}" y1="1" x2="${cx}" y2="${base}"/>
  </svg>`;
}

/* A failed fit still has a histogram worth seeing, so the cell no longer keys
   off fit_success — it falls back only when there is nothing at all to draw. */
function shapeCell(fields, colIdx) {
  return sparkline(fields, colIdx) || '<span class="nil">—</span>';
}
```

- [ ] **Step 3: Use `shapeCell` in the row renderer (:1297)**

Replace:

```js
      <td class="c-shape">${ok ? sparkline(fields, loadedColIdx)
        : '<span class="nil">no fit</span>'}</td>
```

with:

```js
      <td class="c-shape">${shapeCell(fields, loadedColIdx)}</td>
```

- [ ] **Step 4: Add a bar swatch to the legend (:547)**

Add before the existing `fitted curve` key, and add the CSS rule next to the other `.swatch` rules:

```html
        <span class="key"><span class="swatch b"></span> observed histogram</span>
```

```css
.legend .swatch.b { background: var(--muted); opacity: .45; height: 9px; width: 13px; }
```

- [ ] **Step 5: Build the fixture page and look at it**

```bash
python - <<'PY'
import pandas as pd, pathlib, json
out = pathlib.Path("fourparamsacrosstissues/outputs")
d = pd.read_csv(out / "v11_log2_uterus_fourparam.csv", dtype=str, keep_default_na=False)
pick = {}
pick["clean unimodal"] = d[(d.fit_success=="True") & (d.n_obs.astype(int)>100)].iloc[0]
pick["all at -1"]      = d[d["min"] == d["max"]].iloc[0]
pick["failed fit"]     = d[d.fit_success=="False"].iloc[0]
pick["degenerate wide"] = d[(d.fit_success=="True") &
    (d.w.replace("","0").astype(float) > 5*(d["max"].replace("","0").astype(float)
     - d["min"].replace("","0").astype(float)))].iloc[0]
print(json.dumps({k: v.to_dict() for k, v in pick.items()}, indent=1)[:1500])
PY
```

Then open `docs/index.html` in a browser, load a gene in `uterus` from the four
categories, and confirm: bars visible, curve overlaid in count units, ceiling at
the right edge with dead space past it, and the all-at-−1 gene renders a single
spike rather than a blank cell.

- [ ] **Step 6: Commit**

```bash
git add docs/index.html
git commit -m "Draw the real histogram with the fitted curve overlaid in the Shape column"
```

---

### Task 5: Rebuild the browser's static inputs and update the docs

**Files:**
- Modify: `docs/manifest.json` (generated), `docs/genes.tsv` (generated)
- Modify: `fourparamsacrosstissues/CLAUDE.md`

**Interfaces:**
- Consumes: regenerated tables from Task 3
- Produces: a `manifest.columns` array ending `"hist", "hist_max"`

**Note:** `build_gui_data.py` reads a *reference* table for the column list. It must be pointed at one of the four regenerated tissues, not a stale one, or the manifest will not list the new columns. `build_gene_major.py` is **not** run here — the shards cannot be rebuilt consistently until all 54 tissues carry the columns, so the browser stays on the tissue-major route for these four until the full regeneration happens.

- [ ] **Step 1: Rebuild the manifest against a regenerated table**

```bash
cd fourparamsacrosstissues/fourparam && python build_gui_data.py
```

- [ ] **Step 2: Verify the manifest picked up the columns**

```bash
python -c "
import json; m = json.load(open('docs/manifest.json'))
assert m['columns'][-2:] == ['hist','hist_max'], m['columns'][-4:]
print(len(m['columns']), 'columns')"
```
Expected: `22 columns`.

- [ ] **Step 3: Update `fourparamsacrosstissues/CLAUDE.md`**

In the "Generated statistics" table, append two rows:

```markdown
| `hist` | 40-char thumbnail histogram, one character per bin from `A-Za-z0-9+/` (0-63), each `round(count * 63 / hist_max)`. **A rendering aid for the browser, not analysis data** — quantised and lossy to 1/63 of peak. Empty when `n_obs = 0`. Bin edges are not stored: they are `linspace(min, max, 41)`, and numpy widens a zero-width range to `[min-0.5, max+0.5]`. |
| `hist_max` | exact count in the tallest bin; `0` when there is no histogram (never NaN, which would make pandas write every other row's value as a float) |
```

In "The browser GUI" section, replace the `Shape` description: the cell shows the
row's real 40-bin histogram with the fitted curve overlaid on the same count
axis, 160×44, padded 12% past the ceiling on the right only. The curve is clipped
rather than rescaled, so a degenerate fit visibly leaves the frame.

Add to "Things that will bite you if you don't know them":

```markdown
- **The `hist` column must never contain a comma or a quote.** The browser parses
  these CSVs with a plain `split(",")` and has no quoted-field handling, so a
  single comma there shifts every column index after it. That is why the
  histogram is one character per bin rather than comma-separated counts.
```

- [ ] **Step 4: Commit**

```bash
git add docs/manifest.json docs/genes.tsv fourparamsacrosstissues/CLAUDE.md
git commit -m "Rebuild GUI manifest for hist columns; document them"
```

---

## Deferred: the full regeneration

Out of scope for this plan, to be run once Tasks 1-5 are confirmed good:

1. `git lfs pull` the remaining 50 matrices (~5.92 GB), or let `run_cluster.py` fetch-and-delete per tissue.
2. `python generate_all.py` — 108 tables.
3. `python build_gene_major.py` then `--verify`.
4. `python build_gui_data.py`.
5. Push (~740 MB of new data across both mirrors) — **confirm before pushing**, per `fourparamsacrosstissues/CLAUDE.md`.
