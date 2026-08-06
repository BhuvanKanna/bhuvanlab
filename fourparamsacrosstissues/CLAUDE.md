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
  bhuvanfitter.py             <- the fit library (4-param Gaussian + metrics); single source of truth
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
| `hist` | 40-character thumbnail histogram, one character per bin from `A-Za-z0-9+/` (= 0–63), each `round(count * 63 / hist_max)`. **A rendering aid for the browser, not analysis data** — quantised and lossy to 1/63 of the peak. Empty when `n_obs = 0`. Written whenever `n_obs >= 1`, **independent of `fit_success`**, so a gene whose fit failed still shows its real distribution. Bin edges are not stored: they are `linspace(min, max, 41)`, and numpy widens a zero-width range to `[min-0.5, max+0.5]` |
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

## The browser GUI (`docs/`, GitHub Pages)

`docs/index.html` is a self-contained page published at
<https://bhuvankanna.github.io/bhuvanlab/>. It has two tabs:

- **Gene selection** — type-ahead multi-select over genes *and* tissues (both
  boxes work the same way; tissue matching is underscore-insensitive so "whole
  blood" finds `whole_blood`), a raw / excluded ≤ −1 toggle, and an explicit
  **Load tables** button. Results are one row per gene × tissue, each carrying a
  checkbox and that gene's **real 40-bin histogram with the fitted curve
  overlaid** and the truncation ceiling marked (see "The Shape cell" below).
- **Working set** — rows ticked on the first tab, accumulated across any number
  of separate queries, filterable, and exported as one CSV.

### The Shape cell

Each row draws that gene's **real 40-bin histogram** (from `hist` / `hist_max`)
with its fitted 4-parameter Gaussian overlaid **on the same count axis** — the
curve was fit to bin counts, which is what makes the overlay a comparison rather
than a decoration. 160 × 44px; grey bars are the observation, the teal curve is
the model, and the rust line is the ceiling.

- **The x-domain is padded 12% past `max`, on the right only.** `x_max` *is* the
  data max, so without that pad the ceiling would coincide with the panel's right
  border and carry no information. The asymmetry is the point: it says the cap is
  on the right. Do not "tidy" it into symmetric padding.
- **The curve is clipped, not rescaled.** The y-axis is fixed to
  `[0, 1.15 * hist_max]`, so a degenerate fit visibly leaves the frame instead of
  squashing the bars to nothing. That exit is the signal. Clipping is
  `.spark { overflow: hidden }` — deliberately not a `clipPath`, whose id would
  have to stay unique across up to `RENDER_LIMIT` rows.
- **The bars are one stepped `<path>`, not 40 `<rect>`s.** At 800 rendered rows,
  40 rects apiece would add ~32k DOM nodes.
- **The cell no longer keys off `fit_success`.** A failed fit still has a
  histogram worth seeing; the cell falls back to `—` only when there is neither a
  histogram nor fit parameters.

> Why the earlier synthetic sparkline was replaced: it drew the fitted curve over
> the peak's own ±3.4σ and never the data. For APP in kidney cortex (σ = 16 for a
> 4-unit data span) that window is [−52, 66], which renders the ceiling at 99% of
> full height while `truncationindex` is exactly 0 — the table and the picture
> normalising over windows ~30× apart. See
> `specs/2026-08-05-histogram-thumbnail-design.md`.

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
