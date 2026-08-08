# CLAUDE.md

This file tells Claude Code what this repository is and what to do in it. Assume
no prior context — everything you need is here or in `fourparam/`.

## What this project is

For each of the **54 GTEx v11 human tissues**, we characterise the shape of every
gene's expression-level distribution across the population of donors. Per gene we
fit a **4-parameter Gaussian** to its 40-bin expression histogram and compute a
set of shape / truncation metrics (see "Generated statistics" below). The headline
metric is the **truncation index**: how strongly the right tail of a gene's
distribution is "capped", which is a candidate signal for genes whose
over-expression is deleterious (individuals expressing past a ceiling are censored
from the healthy population, truncating the distribution).

The expression matrices have already been downloaded and transformed (see
`data/`). **Your job is to generate the fourparam tables** (see "What to do").

## Layout

```
CLAUDE.md                     <- this file
fourparam/                    <- all the code (kept separate from the data)
  bhuvanfitter.py             <- the fit library (4-param Gaussian, KDE, truncated MLE + metrics); single source of truth
  normality.py                <- distribution-classification cascade + calibrated thresholds
  compute_qc.py               <- generate ONE qc table from ONE matrix (raw or excluded)
  generate_fourparam.py       <- generate ONE fourparam table from ONE matrix (raw or excluded)
  generate_all.py             <- driver: both tables for every tissue -> 108 tables
  run_cluster.py              <- same 108 tables, spread over several machines
  convert_gct_to_log2.py      <- (reference) raw GTEx GCT -> log2(TPM+1)-1 CSV
  download_data.py            <- (reference) how data/ was produced
  extract_genes.py            <- pull a gene set out of the tables -> one tidy CSV
  build_gene_major.py         <- re-orient outputs/ into gene_major/ shards
  build_gui_data.py           <- regenerate the browser GUI's static inputs
  make_diagrams.py            <- 5-histogram summary sheet per table -> diagrams/
  verify_hist_columns.py      <- check hist/hist_max invariants in outputs/
  tests/                      <- pytest suite: `python -m pytest tests/ -q`
data/                         <- the 54 input matrices (nothing else)
  v11_log2_<tissue>.csv.gz    <- one per tissue, already log2(TPM+1)-1 transformed
outputs/                      <- the generated tables go here (starts empty)
  v11_log2_<tissue>_fourparam.csv                              <- raw
  v11_log2_<tissue>_fourparam_excluded_at_or_below_-1.csv      <- excluded <= -1
qc/                           <- distribution class + fit-validity, joined on `gene`
  v11_log2_<tissue>_qc[_excluded_at_or_below_-1].csv
gene_major/                   <- the same rows re-oriented for the browser
  shard_NNNN.csv              <- 16 genes x all 54 tissues x both filters
diagrams/                     <- one 5-histogram summary sheet per table
  v11_log2_<tissue>_fourparam[_excluded_at_or_below_-1].png
genelists/                    <- reusable gene sets, one token per line
results/                      <- extracted gene subsets (output of extract_genes.py)
```

At the **repo root** (one level up) `docs/` holds the GitHub Pages GUI — see
"The browser GUI" below.

Requires Python with `numpy`, `pandas`, `scipy` (and `matplotlib`, imported by
`bhuvanfitter.py`).

## The data (`data/`)

Each `data/v11_log2_<tissue>.csv.gz` is a **genes-as-rows** matrix:

    Name,Description,<sample_1>,<sample_2>, ...
    ENSG00000223972.6,DDX11L1,-1.0,0.585, ...

- `Name` = versioned Ensembl gene id (the gene identifier).
- `Description` = common gene name / symbol.
- Every other column is one donor sample.

Values are already **`log2(TPM + 1) - 1`** transformed (linear TPM → log2). Under
this transform **TPM = 0 maps to exactly `-1`** (the floor, no clamping) and every
value is `>= -1`. Integer TPMs 0,1,2,3,… → −1, 0, 0.585, 1, …

## What to do (produces 108 tables)

Generate **two** fourparam tables **per tissue** — 54 tissues × 2 = **108 tables**
in `outputs/`:

```bash
cd fourparam
python generate_all.py            # both tables for all 54 tissues (skips any already written)
```

`generate_all.py` is resumable — it skips tables that already exist, so if it is
interrupted just run it again. Tune `--jobs N` for parallelism (default 8);
`--max-nfev` caps the curve-fit effort (default 2000).

### Running it across several machines (much faster)

The work is ~74.6k curve fits per table and ~7.5M in total, which is hours on
one box. `run_cluster.py` spreads the 54 tissues over several machines:

```bash
cd fourparam
python3 run_cluster.py                    # all configured nodes
python3 run_cluster.py --dry-run          # show the plan, run nothing
python3 run_cluster.py --nodes local,kup  # a subset
```

It reimplements none of the science — each node runs `generate_fourparam.py`
unchanged, so the tables match a local run. Nodes pull from a shared work queue,
so machines with different core counts and speeds each do as much as they can
and a slow or dead node just does less rather than stalling the run. Per tissue
a node receives the one `.csv.gz`, writes both tables, sends them back, then
deletes its copy. Like `generate_all.py` it is resumable, so re-running it picks
up whatever is missing.

Nodes are listed in the `NODES` table at the top of the script (an
`~/.ssh/config` host alias plus the worker count to use); `local` means this
machine, with no copying. Each remote node needs a `~/fourparam-venv` holding
`numpy`, `pandas`, `scipy` and `matplotlib` — a preflight check imports
`bhuvanfitter` and runs a real fit on every node, and drops any node that fails
instead of letting it die mid-run.

To fetch a tissue matrix that is missing from `data/` (the converted
`log2(TPM+1)−1` CSV is derived from GTEx's raw GCT):

```bash
cd fourparam
python download_data.py --tissues thyroid,whole_blood   # omit --tissues for all 54
```

To fit one tissue by hand:

```bash
cd fourparam
# raw table
python generate_fourparam.py --input ../data/v11_log2_liver.csv.gz \
    --id-col Name --name-col Description --jobs 8
# excluded (<= -1) table
python generate_fourparam.py --input ../data/v11_log2_liver.csv.gz \
    --id-col Name --name-col Description --threshold -1 --jobs 8
```

### Check the histogram columns after regenerating

```bash
cd fourparam
python verify_hist_columns.py                      # every table in outputs/
python verify_hist_columns.py --tissues uterus     # just one
```

Exits non-zero on any problem. It checks the things the browser depends on and
cannot defend itself against: uniform field counts (one stray comma in `hist`
shifts every column index after it), fixed 40-char `hist`, integer `hist_max`,
`hist` empty exactly when `n_obs = 0`, and — the one that has already bitten
once — **no histogram shipped without its `min`/`max` bin edges**.

### The two table types differ by exactly one filter

- **raw** (`v11_log2_<tissue>_fourparam.csv`): keeps **every finite** expression
  value for each gene.
- **excluded ≤ −1** (`..._fourparam_excluded_at_or_below_-1.csv`): additionally
  drops every value **at or below −1** from a gene's array *before* fitting.
  Because `log2(TPM+1)−1` hits `−1` only at TPM = 0, this removes exactly the
  **zero-expression (undetected) samples**, so a gene's `n_obs` shrinks and it may
  fall below the minimum.

Nothing else changes between the two — same fit, same columns, same handling.

### Default filters (applied to BOTH table types)

These are the standard filters and are **not** specific to the excluded table.
Per gene, in order:

1. **Drop non-finite values** (NaN / inf).
2. *(excluded table only — this is the one exclusion-specific step)* drop values
   **≤ −1**.
3. **Minimum observations:** if the surviving count `n_obs < 10` (`MIN_OBS`), the
   gene is still written, as a row with `fit_success = False` and all metrics
   `NaN`. Genes are never silently dropped.
4. **Fit convergence:** if `curve_fit` fails to converge, same
   `fit_success = False` / `NaN` row.

So the **only** thing that makes an "excluded ≤ −1" table different from a "raw"
table is step 2. Everything else — the finite filter, the `MIN_OBS = 10` floor,
the failure rows, the columns — is default behaviour shared by both.

> Note on downstream analysis filters: some analyses additionally keep only
> `fit_success == True`, `0 < truncationindex < 1`, and raise the floor to
> `n_obs >= 30`. Those are **analysis-time** filters you apply when *using* a
> table — they are deliberately **not** baked into the generated tables, which
> keep every gene (with `fit_success` flags) so nothing is thrown away up front.

## Generated statistics (columns of every table)

One row per gene. Columns, in order:

| column | meaning |
|---|---|
| `gene` | Ensembl gene id (`Name`) |
| `genename` | common gene name / symbol (`Description`) |
| `y0` | fitted Gaussian baseline offset |
| `A` | fitted amplitude (peak height above baseline) |
| `x0` | fitted peak centre |
| `w` | fitted width (`w = sigma * sqrt(2)`) |
| `sumsquarevalue` | residual sum of squares of the fit (lower = better) |
| `ti_fourparam_sigma_dist` | `(x_max − x0)/(w/√2)` — how many σ the ceiling sits above the peak; **lower = more truncated** |
| `truncationindex` | **height-ratio truncation index**, `f(x_max)/f(peak)` with the curve's interval-minimum subtracted from both; **bounded [0, 1]**; higher = more truncated (0 = ceiling at curve min, 1 = ceiling at peak) |
| `min`, `max` | min / max of the values used for the fit |
| `mean`, `std` | mean and **sample** std (ddof=1) of those values |
| `skew`, `kurt` | skewness and **Fisher excess** kurtosis (normal → 0) of those values |
| `right` | the truncation ceiling `x_max` used by the metrics (the data max) |
| `maxheight`, `rightheight` | curve height at the peak / at `x_max`, above the curve's interval-minimum baseline (`truncationindex = rightheight/maxheight`) |
| `n_obs` | number of finite values fit (after any exclusion) |
| `fit_success` | `True` if the fit converged with `n_obs >= 10`, else `False` (metrics `NaN`) |
| `hist` | 40-character thumbnail histogram, one character per bin from `A-Za-z0-9+/` (= 0–63), each `round(count * 63 / hist_max)`. **Not analysis data** — quantised and lossy to 1/63 of the peak. Written for a browser thumbnail that no longer exists: the Shape cell now draws the fitted curve, and `outputs/` never carried these two columns in the first place, so nothing currently reads them. Empty when `n_obs = 0`. Written whenever `n_obs >= 1`, **independent of `fit_success`**, so a gene whose fit failed still shows its real distribution. Bin edges are not stored: they are `linspace(min, max, 41)`, and numpy widens a zero-width range to `[min-0.5, max+0.5]` |
| `hist_max` | exact count in the tallest bin; `0` when there is no histogram (never `NaN` — one `NaN` in the column makes pandas write every other row's value as a float) |

Fitting details (in `bhuvanfitter.py`, the single source of truth): the histogram
is always **40 bins**; the curve is fit by ordinary least squares (Trust Region
Reflective) minimising the residual sum of squares; the truncation ceiling
`x_max` defaults to the gene's observed max. The `mean/std/skew/kurt` are plain
summary statistics of the (post-exclusion) values and are reported regardless of
whether the fit converged (NaN on a `fit_success=False` row).

## Extracting a gene set (`extract_genes.py`)

Pulls the rows for a gene set out of all 108 tables into one tidy CSV, one row
per `gene × tissue × table`, carrying every statistic column:

```bash
cd fourparam
python extract_genes.py --genes APP SNCA "ALDH*" -o ../results/genes.csv
python extract_genes.py --genes-file ../genelists/adh_aldh_plus.txt -o ../results/out.csv
python extract_genes.py --genes "ADH*" --tissues liver,lung --table raw -o out.csv
```

Tokens may be gene symbols, Ensembl ids, or globs (matched against symbols first,
then ids). **Ensembl versions are always stripped before matching** — GTEx bumps
the `.NN` suffix between releases, so an id list written against an older
annotation (`ENSG00000142192.20`) still resolves against these tables (`.22`).
Matching the full versioned string would silently return nothing, which is the
single easiest way to get a wrong answer here.

Each query token is reported as resolved (with what it hit) or `NOT FOUND`, so a
typo or an absent symbol is loud rather than quietly empty.

## Distribution sheets (`diagrams/`, `make_diagrams.py`)

One PNG per table — 54 tissues × 2 filters = **108 sheets, 540 histograms** —
covering `truncationindex`, `sumsquarevalue`, `mean`, `std`, and
`ti_fourparam_sigma_dist` over the genes in that table, plus a panel stating what
was excluded. File names mirror the table names, so a sheet sorts next to its
source.

```bash
cd fourparam
python make_diagrams.py                 # all 108, skips existing
python make_diagrams.py --tissues liver --force
```

Only `fit_success == True` rows with a finite value are histogrammed, and the
dropped counts are printed on the sheet rather than silently omitted.

**The axis rules are per-metric on purpose** — these five columns do not share a
shape, so one rule cannot serve them all. Do not "simplify" them to a common
linear axis; each choice is load-bearing:

- `truncationindex` is bounded [0, 1] and violently zero-inflated (~91% of liver
  genes sit at exactly 0). Fixed [0, 1] domain with a **log count axis**; on a
  linear one the panel is a single bar and the tail is invisible.
- `ti_fourparam_sigma_dist` spans ~6 orders of magnitude either side of zero
  (degenerate fits with `w ≈ 1e-4` reach 1e5 while the meaningful range is single
  digits). **Symlog x**, linear within ±1, with symlog-spaced bins. Nothing is
  clipped — no linear window shows the bulk without hiding ~10% of the genes.
- `sumsquarevalue`, `mean`, `std` get an adaptive linear window: the narrower of
  the 6×IQR fences and the p0.5–p99.5 span that hides ≤ 2%, falling back to
  whichever hides less. The hidden count is stamped on the panel.

## Distribution QC (`qc/`, `compute_qc.py`, `normality.py`)

Answers two questions the fourparam tables cannot: **is this gene's distribution
even Gaussian**, and **is this particular fit usable**. One row per gene, joined
to `outputs/` on `gene`.

```bash
cd fourparam
python compute_qc.py --input ../data/v11_log2_kidney_cortex.csv.gz \
    --id-col Name --name-col Description --jobs 8              # raw
python compute_qc.py --input ../data/v11_log2_kidney_cortex.csv.gz \
    --id-col Name --name-col Description --threshold -1 --jobs 8  # excluded
```

Skips existing output unless `--force`. ~16 MB and a few minutes per table.

**This deliberately does not add columns to `outputs/`.** Those tables are
mirrored into 4,665 gene_major shards, pinned by four hardcoded column lists,
and carry a byte-identity guarantee with `extract_genes.py`. A separate table
joined on `gene` costs none of that.

### The cascade (`normality.py`)

One p-value is not an answer, because "not normal" covers zero-inflation,
bimodality, right skew and right-truncation, and only the last is this project's
hypothesis. Each stage only sees what survives the one before:

| stage | function | rules out |
|---|---|---|
| 0 | `classify_support` | ≥20% of donors at the −1 floor — a spike plus a smear, not a failed Gaussian |

| 1 | `count_modes` | >1 KDE mode (reuses `gene_peaks`) — wants a mixture, not a truncated Gaussian |
| 2 | `run_tests` | Shapiro–Wilk, D'Agostino–Pearson, Anderson–Darling |
| 3 | skew direction | right-truncation removes the **upper** tail, so it makes skew **negative** |
| 4 | `fit_truncated` + `calibrate_null` | truncated-Gaussian MLE vs plain Gaussian, against a simulated null |
| 5 | `classify` | one `dist_class` label |

**Kolmogorov–Smirnov is deliberately absent** — against parameters estimated from
the same sample its null distribution is wrong, and Lilliefors is still weaker
than Shapiro–Wilk here.

**`FLOOR` is −0.75, but the excluded table cuts at −1.0, and that gap is
deliberate.** Values in `(−1, −0.75]` survive the exclusion filter and still
count toward `frac_at_floor`, so 59.6% of genes stay `zero_inflated` in the
excluded kidney-cortex table. That is correct, not a leak: the band maps to
`0 < TPM ≤ 0.19`, which is trace expression, and a gene made up of it is no more
Gaussian than one made of exact zeros. It matches `has_minus_one_peak` in the
notebook. Genes where *every* donor was at TPM=0 come back `undetermined`
(n_obs = 0 after filtering), which is 14.1% of the excluded table.

### `calibrate_null` is the load-bearing piece — do not replace it with a constant

The truncation point is fixed at the **observed maximum**, which is
data-dependent and therefore biased toward the truncated model. Measured against
2,000 samples drawn from a true Gaussian at n=104:

```
dAIC > 2      fires on 41.4% of genuinely normal data   <-- unusable
dAIC > 5.87   is the calibrated 95th percentile          <-- correct at n=104
skew noise    +/-0.37 at n=104                           <-- |skew| below this is nothing
```

So `d_aic` and `skew` are **always** judged against `calibrate_null(n)`, never a
literal. `d_aic` is a likelihood ratio between two location-scale families and
`skew` is standardised, so both are scale-free — one calibration per sample size
serves every gene in a tissue, and it is `lru_cache`d for exactly that reason.
`tests/test_normality.py` pins both the ~5% flag rate on held-out null data and
the fact that the naive threshold is wrong.

**Power is the real limit.** At n=104 a ceiling at 0.5σ or 1σ is detectable; at
2σ or 3σ it is not, because too little mass is missing to see. Read a null result
in a small tissue as "underpowered", not "no truncation". Kidney cortex has 104
donors and ranks 49th of 54; muscle skeletal has 818.

### QC gates are absolute, not percentile-based

A percentile rule flags a fixed fraction however good the fits are. These are
derived quantities with meaning, so they get real thresholds:

| gate | keep | catches (kidney cortex, raw) |
|---|---|---|
| `sigma_span_ratio = (abs(w)/√2)/(max−min)` | 0.05 – 1.0 | **~25%** collapsed spikes below, **3.3%** APP-mode above |
| `x0_in_range` | `x0` ∈ `[min, max]` | **23.5%** — peak never observed |
| `n_obs` | ≥ 30 | thin excluded-table rows |
| `frac_at_floor` | < 0.20 | **63%** zero-inflated |

Two opposite degeneracies both report `fit_success = True`: `w → 0` collapses the
Gaussian onto one bin (DDX11L1: `w ≈ 1e-6`, `sigma_dist = −17,676`), and `w → ∞`
fits a shallow slice off the top of an enormous one (APP: `w = 22.7`, σ = 16.0
against a 3.95-unit span, `sigma_dist = 0.124` reading as maximal truncation).
Neither is biology.

**`r_squared` replaces `sumsquarevalue` for any cross-gene comparison.** The
stored SSR is unnormalised — it scales with n and peak height, so p50 = 58 says
nothing about whether one fit is good. Expect low values regardless: 40 bins over
~104 donors is ~2.6 donors per bin, and real genes land around R² 0.2–0.4.

### What the first two tissues actually showed — read this before trusting a candidate list

Judge the `right_truncated` rate against the **joint** null of the compound rule
(skew below the band **and** `d_aic` past the threshold). The two criteria are
positively correlated — left skew is exactly what makes the truncated model fit
better — so the joint rate is nowhere near `0.05 x 0.05`. Measured by simulation:
**2.00% at n=104**, **0.67% at n=818**.

| tissue | n | observed (of testable) | joint null | vs chance |
|---|---|---|---|---|
| kidney_cortex | 104 | 0.59% (111/18,819) | 2.00% | **0.29x** |
| muscle_skeletal | 818 | 2.08% (5,000-gene sample) | 0.67% | **3.10x** |

**Power was the limiting factor, and the plan's suspicion was right.** At n=104
there is no excess at all. At n=818 there is a real one, three times chance.
Do not read a null result in a small tissue as absence of truncation — 20 of the
54 tissues have n < 200.

**But the excess is very likely compositional, not lethality.** The candidates
are each tissue's own highly-expressed identity genes, and the two lists do not
overlap at all:

```
kidney_cortex    NPHS2 NPHS1 WT1 SLC22A6 SLC22A8 SLC5A2 SLC6A19 SLC13A3 ANPEP + 11 MT-*
muscle_skeletal  ATP2A1 TNNC1 CASQ2 MYOM3 FHL1 SPEG PPP1R3C RCAN2 TMEM143
```

Both sit at the **92nd-93rd percentile of expression** in their tissue, and MT
genes are **60x enriched** among kidney-cortex candidates (10% vs 0.17%).

Donor-to-donor variation in tissue composition and RNA quality produces exactly
this: a biopsy with less pure cortex has less podocyte and tubule signal, a
degraded sample has less mitochondrial content, and either drags a subset of
donors downward into a long **left** tail. Sample heterogeneity and a biological
expression ceiling are indistinguishable to this statistic.

A dosage-lethality signal should be at least partly **shared** across tissues,
since dosage-sensitive genes tend to be broadly dosage-sensitive. Getting each
tissue's own marker genes back instead is the signature of an artefact.

**Before treating any of these as biology**, regress the candidates against GTEx
sample covariates (RIN, ischemic time, `SMTSISCH`) or a cell-type deconvolution,
and check whether the left tail survives. Those covariates are not in this repo —
they come from the GTEx sample attributes file.

### `dist_class` values

`zero_inflated` · `multimodal` · `right_truncated` · `right_skewed` ·
`non_normal` · `normal` · `undetermined`

`right_truncated` requires **both** skew below the null band and `d_aic` past the
calibrated threshold. Either alone is common noise.

## The browser GUI (`docs/`, GitHub Pages)

`docs/index.html` is a self-contained page published at
<https://bhuvankanna.github.io/bhuvanlab/>. It has two tabs:

- **Gene selection** — type-ahead multi-select over genes *and* tissues (both
  boxes work the same way; tissue matching is underscore-insensitive so "whole
  blood" finds `whole_blood`), a raw / excluded ≤ −1 toggle, and an explicit
  **Load tables** button. Results are one row per gene × tissue, each carrying a
  checkbox and **two views of that gene's fitted curve** — one in real log2
  units, one on a fixed sigma axis (see "The two shape columns" below).
- **Working set** — rows ticked on the first tab, accumulated across any number
  of separate queries, filterable, and exported as one CSV.

### The two shape columns

Each row draws its fitted 4-parameter Gaussian **twice**, 160 x 44px each. The
two answer different questions and neither replaces the other, which is why both
ship. The headers name the axis rather than leaving it to the legend.

| column | axis | answers |
|---|---|---|
| **Fitted curve** | real log2 units, `[min, max]` + 12% right pad | what shape was fit to this gene, on the scale the `x0` / `w` / `min` / `max` columns are printed in |
| **Curve coverage** | fixed −4σ … +4σ, same for every gene | how much of that curve the data actually spans, and where the ceiling sits |

Shared ink: solid + filled is observed, dashed is censored past the ceiling,
rust vertical is the ceiling `x_max`. Coverage adds a dotted curve below `min`
(modelled, never observed) and a muted dashed line at the data floor `min`.

**Fitted curve** — `fittedCurveCell()`, reads `y0`, `A`, `x0`, `w`, `min`,
`max`, `right`.

- **The pad is 12% past `max`, on the right only.** `x_max` *is* the data max,
  so without it the ceiling lands exactly on the frame border and carries no
  information. The asymmetry is the message: the cap is on the right. Do not
  tidy it into symmetric padding.
- **y is anchored at zero; only the top is data-derived.** `outputs/` carries no
  count axis to borrow, so some self-scaling is unavoidable — anchoring at zero
  is what stops a near-flat fit from being stretched to fill the cell.
- **Both scales are per-gene**, so two of these cells side by side are *not*
  comparable. That is the honest cost of real units, and precisely what the
  column next to it exists to fix.

**Curve coverage** — `shapeCell()`, reads `ti_fourparam_sigma_dist`, `w`, `x0`,
`min`.

- **The bell is universal; only the window is per-gene.** `((x − x0)/w)²` is
  exactly `z²/2` for `z = (x − x0)/σ`, `σ = w/√2`. So in z units every gene's
  fitted curve is the same `exp(−z²/2)` once `y0` and `A` divide out — there is
  no per-gene shape to draw. The content is `z_min = (min − x0)/σ` to
  `z_max = (max − x0)/σ`, and `z_max` **is** `ti_fourparam_sigma_dist`.
- **Marking `z_min` is what makes a degenerate fit visible, and is the whole
  reason this is not just a canonical bell with a cut.** APP in kidney cortex
  fits `w = 22.7` (σ = 16.0) to a 3.95-unit data span, so its entire observed
  range is a **0.25σ sliver** at the apex. `ti_fourparam_sigma_dist` reads
  0.124 — apparently maximal truncation — and `truncationindex` reads 0, and
  neither means anything. Here that row is a thin spike, not a plausible bell.
- **The axis is fixed and never data-derived**, so cells are comparable by
  construction and nothing is rescaled to its own range. Both edges clamp to the
  frame: a gene at 40σ and one at 4σ both read as "no visible truncation", which
  is exactly true at this scale.
- **A window is only drawn when it is one.** An inverted or non-finite `z_min`
  (`w <= 0`, or a table lacking `min`/`x0`) degrades to the bare cut rather than
  inventing an extent.

Neither cell keys off `fit_success` — a converged-but-degenerate fit is exactly
what these pictures are for. Each falls back to `—` only when its own inputs are
missing or non-finite.

> Two earlier designs and why they went:
>
> - A **synthetic sparkline** over the peak's own ±3.4σ, which never showed the
>   data. For APP that window is [−52, 66], rendering the ceiling at 99% of full
>   height while `truncationindex` is 0 — table and picture normalising over
>   windows ~30× apart. See `specs/2026-08-05-histogram-thumbnail-design.md`.
> - The **real 40-bin histogram** from `hist` / `hist_max` with the curve
>   overlaid on the count axis. Correct in principle, but `outputs/` does not
>   carry those two columns, so in production the cell rendered a dash for every
>   row. `hist` and `hist_max` are now in the skip set in `visibleStatCols()` —
>   nothing draws them and an encoded 40-char blob is not a readable table
>   column. They remain in the CSV export, which mirrors the source table rather
>   than this view.

Three things about it are load-bearing and easy to undo by accident:

1. **Loading is explicit, and picks an orientation.** The same numbers exist
   twice — see "Two orientations" below. The page prices the query before
   running it, fetches sequentially, and can be cancelled mid-run (keeping the
   rows already gathered). Do not make it auto-fetch on selection change: the
   tissue-major route is still chosen for many-genes-in-one-tissue queries, and
   there it is the half-gigabyte accident this design exists to prevent.
2. **The working set stores values, not references.** Each entry keeps all 18
   statistic columns, so it survives a reload without re-fetching. It is keyed
   `tissue|kind|unversioned-id`, which is why the same gene can appear for many
   tissues and for both raw and excluded at once. It is persisted to
   `localStorage` under `tb-working-set-v1`; if the quota is exceeded the page
   falls back to memory-only and says so rather than failing silently.
3. **Its CSV is byte-identical to `extract_genes.py`'s** — header `tissue, table,
   gene, genename` then the table's own columns in `manifest.json` order, and
   values passed through as text. Exports always carry the full statistic set
   regardless of the compact/all column toggle, and always cover every loaded row
   rather than the rendered subset (rendering is capped at 800 rows). Keep these
   in sync if columns change.

   > This is only true because **nothing re-serialises a float**. `extract_genes.py`
   > reads with `dtype=str, keep_default_na=False`; `build_gene_major.py` builds
   > rows by string concatenation. Parsing to float64 and writing back is *not*
   > lossless: pandas' CSV writer emits ~16 significant digits rather than the
   > shortest round-tripping repr, which silently turned `0.012596832467784065`
   > into `0.012596832467784`. If you make either tool go through floats, this
   > guarantee quietly dies.

**It stores no copy of the tables.** It fetches the real CSV from
`raw.githubusercontent.com` at runtime, so it cannot drift out of sync with
`outputs/`. `raw.githubusercontent.com` serves `Access-Control-Allow-Origin: *`,
so the cross-origin fetch works.

### Two orientations: `outputs/` and `gene_major/`

The same 8,059,824 rows are published twice, because the two shapes answer
opposite questions and a browser query is nearly always the second kind:

| | file | one fetch gets you | cost |
|---|---|---|---|
| tissue-major | `outputs/v11_log2_<tissue>_fourparam*.csv` | every gene, **one** tissue | ~8 MB |
| gene-major | `gene_major/shard_NNNN.csv` | 16 genes, **every** tissue and both filters | ~220 KB |

`shard_NNNN.csv` carries the `tissue,table,` prefix columns that `outputs/` rows
lack — that is exactly the `extract_genes.py` header, so a shard and a CLI
extract concatenate without reconciling anything. The GUI strips those two
fields on read so both routes produce identical row objects downstream.

**Genes are sharded in symbol order, not id order.** That is deliberate: a family
lands in adjacent shards, so all 27 `ALDH*` cost 2 fetches rather than 27. The
whole `adh_aldh_plus` list (38 genes × 54 tissues × 2 = 2,052 rows) is 8 shards,
~1.8 MB, against 432 MB tissue-major.

`docs/genes.tsv` is `id \t symbol \t shard`; the third column is that map. The
GUI picks whichever route fetches fewer bytes, counting only what is not already
cached, and labels the choice in the status line. Rebuild after regenerating any
table:

```bash
cd fourparam
python build_gene_major.py            # ~80 s, 4,665 shards, 2.46 GB
python build_gene_major.py --verify   # byte-compares shards against outputs/
```

Line endings matter here. `build_gene_major.py` copies row text verbatim, so it
normalises CRLF on read; the writers now pin `lineterminator="\n"`. The eight
tables for thyroid/uterus/vagina/whole_blood were written with CRLF by an earlier
run and have been normalised in place.

Pages publishes **only `docs/`** (branch `main`, folder `/docs`). That matters:
the Pages site size limit is 1 GB and `outputs/` alone is 2.05 GB, so publishing
the whole repo would fail.

Regenerate its two static inputs after adding tissues or changing columns:

```bash
cd fourparam
python build_gui_data.py     # writes ../../docs/manifest.json and ../../docs/genes.tsv
```

`genes.tsv` is built from a single table because the gene set is byte-identical
across all 108; `build_gui_data.py` verifies that and aborts if it ever stops
being true.

**`--reference <tissue>` picks the table that supplies `manifest.columns`.**
It defaults to the alphabetically first raw table, which is wrong mid-migration:
if you have added a column and regenerated only some tissues, the default table
may not carry it yet and you publish a manifest that silently omits it. The
manifest drives the browser's CSV export, so an omitted column is an export that
quietly loses data.

### Four places hardcode the column list — they must agree

| where | what it is |
|---|---|
| `generate_fourparam.COLUMNS` | what gets written to `outputs/` |
| `build_gene_major.SHARD_HEADER` | what the gene-major mirror expects **and** emits |
| `extract_genes.STAT_COLUMNS` | what the CLI extract emits |
| `docs/manifest.json` `"columns"` | what the browser's CSV export follows |

Drift is silent and expensive: `build_gene_major.py` rejects every table with
"unexpected header" (it fails safe, but the whole run dies), or the browser
export loses columns, or the byte-identity guarantee with `extract_genes.py`
quietly stops holding. `tests/test_column_lists_agree.py` pins all four —
run `python -m pytest tests/ -q` after changing any column.

## Keep this file current

If you change the code, columns, filters, or layout, update this file in the same
change so it always describes the repo accurately.

## Push changes to GitHub as we work

This folder is the `fourparamsacrosstissues/` subdirectory of the GitHub repo
**`BhuvanKanna/bhuvanlab`**. The git root is the **parent** directory
(`bhuvanlab-main/bhuvanlab-main/`), not this folder — run git from there, or with
`git -C ..`.

**After completing any change to files here, commit and push it** without waiting
to be asked:

```bash
cd ..                      # the git root
git add -A                 # stages only altered/new files
git commit -m "<what changed and why>"
git push origin main
```

Write a real commit message describing the change. Batch one logical change per
commit rather than committing every intermediate edit. If a task is still
in progress and the tree is in a broken state, finish it first, then push.

### Things that will bite you if you don't know them

- **The `hist` column must never contain a comma or a quote.** The browser parses
  these CSVs with a plain `split(",")` and has no quoted-field handling, so a
  single comma there silently shifts every column index after it. That is exactly
  why the histogram is one character per bin rather than comma-separated counts.
- **`hist` is quantised and lossy.** It is a rendering aid for the browser
  thumbnail, never a source of truth for counts. `n_obs` and `hist_max` are exact;
  use those.
- **Always write `df["hist"]`, never `df.hist`.** `DataFrame.hist` is pandas' own
  histogram-plotting method, so attribute access silently returns a bound method
  instead of the column and fails somewhere later with a confusing error. The
  generator is unaffected (it builds rows as dicts keyed by `COLUMNS`), but any
  analysis code reading these tables will hit this.
- **`data/*.csv.gz` are Git LFS pointers, not real files.** The repo is
  configured `--skip-smudge` with `lfs.fetchexclude=*`, so those 134-byte
  pointers are correct and must stay that way. ~5.7 GB of LFS payload is
  deliberately not downloaded. Never "fix" them. To get one on purpose:
  `git lfs pull --include="fourparamsacrosstissues/data/v11_log2_liver.csv.gz"`.
  `download_data.py` recognises pointers and leaves them alone; without that
  guard its gzip-integrity resume check would call every pointer corrupt and
  re-download the lot.
- **This is a blobless partial clone** (`--filter=blob:none`), so `.git` is a few
  hundred KB rather than 813 MB. Old file contents are fetched on demand; that is
  intended. `git fetch --refetch origin main` backfills if ever needed.
- **`outputs/` is plain git, not LFS, on purpose.** These CSVs compress 2.46× in
  git (2.05 GB → ~852 MB); LFS stores blobs uncompressed and would cost more.
- **A full 108-table regeneration is a ~852 MB push.** That is slow on a home
  connection and permanently adds that much to repo history. Before pushing a
  regeneration of many tissues, say so and confirm — don't push it silently. A
  single-tissue regeneration (~8 MB compressed) is fine to push normally.
- **Never commit while a generation run is still writing** — a partially written
  CSV looks like a valid file. Wait for the run to finish.
