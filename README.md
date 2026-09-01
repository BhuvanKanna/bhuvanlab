# bhuvanlab

Distribution-shape and truncation-index analysis of human gene expression across
the **54 GTEx v11 tissues**.

## → [Open the Truncation Browser](https://bhuvankanna.github.io/bhuvanlab/)

**[bhuvankanna.github.io/bhuvanlab](https://bhuvankanna.github.io/bhuvanlab/)** —
search genes *and* tissues by type-ahead, load the tables, then tick rows into a
**working set** that persists across queries and exports as one spreadsheet.
Every row draws its fitted curve twice. **Fitted curve** puts it in real log2
units over the gene's own observed range, so you see what was actually fit next
to the `x0` / `w` / `min` / `max` columns. **Curve coverage** puts every gene on
one fixed −4σ … +4σ ruler and fills in the slice the data spans, so the cells are
comparable to each other: a wide band is a well-determined fit, a sliver is a fit
whose σ dwarfs its own data. Both mark the truncation ceiling `x_max`.

Two tabs:

- **Gene selection** — pick any number of genes and any number of tissues (up to
  all 54), choose raw or excluded ≤ −1, and press **Load tables**. Results are
  one row per gene × tissue, each with a checkbox.
- **Working set** — everything ticked so far, gathered across as many separate
  queries as you like, filterable and exportable to CSV.

**One button per control-gene set.** Under the gene box, a **Control sets** row
loads any of the four rosters the control-gene analysis scores against —
`pTriplo > 0.94` (1,385 genes) and `DECIPHER dominant` (146) as positives,
`Duplication-tolerant` (630) and `Olfactory receptors` (384) as negatives. A
click is exactly a paste: the same resolver, the same "N not recognised"
report, and the sets **add** to whatever is already chosen, so you can put
positives and negatives on screen together and sort by `truncationindex`. The
lists live in [`genelists/`](fourparamsacrosstissues/genelists/) and are the
same files `extract_genes.py --genes-file` reads.

**Click any gene name** in either table to open that gene on its own page: its
donor histogram in the tissue you were looking at — 40 real bins, drawn as soon
as the page opens, because clicking already named the gene, the tissue and the
filter, so the fetch is a **~5 KB** slice of `hist_major/` rather than the ~8 MB
table those columns live in — with the fitted Gaussian, a
moment-matched normal, `x₀`, `x_max` and the σ marks each on a checkbox; its mean
expression and its truncation index as ranked bar charts across all 54 tissues,
each with the rank spelled out in words; what the gene actually is; and any
associated phenotype — starting with the curated set where OMIM or G2P name
**over-expression itself** as the disease mechanism, which is the hypothesis this
whole repository exists to test. Every panel is one gene-major fetch, and Back
(the button, `Escape`, or the browser's own) returns to your table untouched.

- Type-ahead over all **74,628** genes (`ALDH2`, `ENSG00000111275`) and all 54
  tissues (`whole blood` matches `whole_blood`)
- **Ensembl versions are ignored**, so ids from an older annotation still resolve
- The working set survives a reload, and mixes tissues *and* raw/excluded rows
- Its CSV is **byte-identical** to what `extract_genes.py` writes — same 22
  columns, same digits — so GUI and CLI output are genuinely interchangeable
- Tables are read live from this repository, so the page always matches `outputs/`

> Loading is deliberately an explicit button, not automatic. The page states the
> cost before you commit to it, streams the files one at a time, and lets you
> cancel mid-run while keeping whatever already arrived.

### Why it is fast

The same 8,059,824 rows are published in **two orientations**, and each query
takes whichever is cheaper:

| | one fetch gets you | size |
|---|---|---|
| `outputs/` — tissue-major | every gene, **one** tissue | ~8 MB |
| `gene_major/` — gene-major | 16 genes, **every** tissue and both filters | ~220 KB |

Asking for 38 genes across all 54 tissues used to mean downloading 54 whole
tissue tables — **432 MB** — to surface 2,052 rows. Gene-major answers it in
**8 fetches, ~1.8 MB**, because genes are sharded in *symbol* order, so a family
shares shards: all 27 `ALDH*` cost 2 fetches, not 27. Tissue-major is still there
and still wins when you want many genes from a single tissue.

---

## What the analysis does

For each gene in each tissue, expression across the donor population is binned
into a 40-bin histogram and fit with a 4-parameter Gaussian:

```
y = y0 + A · exp(−((x − x0) / w)²)        w = σ√2
```

A gene whose over-expression is not tolerated should be missing its right tail:
individuals expressing past some ceiling are absent from a healthy cohort, so the
observed distribution is **right-truncated**. The **truncation index** quantifies
that:

| metric | meaning |
|---|---|
| `truncationindex` | `f(x_max)/f(peak)`, baseline removed → bounded **[0, 1]**. Higher = ceiling sits closer to the peak = more truncated. |
| `ti_fourparam_sigma_dist` | `(x_max − x0)/σ` — how many σ the ceiling sits above the peak. **Lower = more truncated.** |

Both are reported per gene alongside the fit parameters and summary statistics.
The full column list is in
[`fourparamsacrosstissues/CLAUDE.md`](fourparamsacrosstissues/CLAUDE.md).

## Distribution QC

Neither truncation metric means anything if the gene was never Gaussian or the
fit was degenerate, and both are common: in kidney cortex **63%** of genes are
zero-inflated, **23.5%** of converged fits place the peak `x0` outside the
observed data, and **96.4%** of `truncationindex` values are exactly 0.

[`qc/`](fourparamsacrosstissues/qc/) carries one row per gene classifying the
distribution (`zero_inflated`, `multimodal`, `right_truncated`, `right_skewed`,
`non_normal`, `normal`) and flagging whether the fit is usable at all. It joins
to `outputs/` on `gene`.

The truncation call is made against a **simulated null**, not a textbook cutoff:
because the ceiling is fixed at the observed maximum, a `ΔAIC > 2` rule fires on
41% of data drawn from a true Gaussian. See
[`fourparamsacrosstissues/CLAUDE.md`](fourparamsacrosstissues/CLAUDE.md#distribution-qc-qc-compute_qcpy-normalitypy).

## Distribution sheets

[`fourparamsacrosstissues/diagrams/`](fourparamsacrosstissues/diagrams/) holds one
PNG per table — **108 sheets, 540 histograms** — covering `truncationindex`,
`sumsquarevalue`, `mean`, `std`, and `ti_fourparam_sigma_dist` across the genes in
that tissue, plus a panel stating what was excluded. File names mirror the table
names.

Only converged fits with a finite value are histogrammed, and the dropped counts
are printed on the sheet. Two of the five axes are not linear, for reasons that
show up immediately in the data: `truncationindex` is bounded [0, 1] and ~91% of
liver genes sit at exactly 0, so it gets a log count axis; `ti_fourparam_sigma_dist`
spans six orders of magnitude either side of zero once degenerate fits are
included, so it gets a symlog axis and nothing is clipped.

## Extracting genes at the command line

`extract_genes.py` pulls a gene set out of all 108 tables into one tidy CSV
(`tissue × table × gene`, every statistic column):

```bash
cd fourparamsacrosstissues/fourparam

# symbols, ids, and globs all work; ids may carry any version suffix
python extract_genes.py --genes APP SNCA "ALDH*" -o ../results/genes.csv

# or from a file, one token per line
python extract_genes.py --genes-file ../genelists/adh_aldh_plus.txt \
    -o ../results/adh_aldh_plus_all_tissues.csv

# narrow it down
python extract_genes.py --genes "ADH*" --tissues liver,lung --table raw -o out.csv
```

Every query token is reported as resolved or unresolved, so a typo is loud rather
than silently returning nothing.

### Worked example, already run

[`genelists/adh_aldh_plus.txt`](fourparamsacrosstissues/genelists/adh_aldh_plus.txt)
holds APP, SNCA, PCSK9, SOX9, SERPINA1, the six requested alcohol dehydrogenases,
and `ALDH*`. That resolves to **38 genes**, giving 38 × 54 tissues × 2 tables =
**4,104 rows** in
[`results/adh_aldh_plus_all_tissues.csv`](fourparamsacrosstissues/results/adh_aldh_plus_all_tissues.csv)
(4,007 with a converged fit).

> **Note on Ensembl versions.** The requested ids came from an older GTEx
> annotation and do **not** match these tables verbatim — APP was given as
> `ENSG00000142192.20` but is `.22` here; SNCA `.15` → `.17`; PCSK9 `.10` → `.13`.
> Matching on the full versioned string would have silently returned zero rows, so
> versions are always stripped before matching.

## Layout

```
README.md                     <- this file
docs/                         <- the GitHub Pages site (the Truncation Browser)
  index.html                  <- the whole GUI, self-contained
  manifest.json               <- tissue list + column list + gene-set index
  genes.tsv                   <- gene index for type-ahead
  genelists/                  <- the four control sets, one button each
  overexpression_phenotypes.tsv <- 51 curated over-expression driver genes
  qc/                         <- one distribution class per gene, per tissue
fourparamsacrosstissues/
  CLAUDE.md                   <- full spec: columns, filters, how tables are made
  data/                       <- 54 input matrices (Git LFS pointers, see below)
  fourparam/                  <- the code
    bhuvanfitter.py           <- the fit library; single source of truth
    generate_fourparam.py     <- one table from one matrix
    generate_all.py           <- all 108 tables
    run_cluster.py            <- the same, spread over several machines
    extract_genes.py          <- pull a gene set out of the tables
    build_gene_major.py       <- re-orient outputs/ into gene_major/ shards
    build_gui_data.py         <- regenerate docs/manifest.json, genes.tsv, genelists/
    make_diagrams.py          <- render the distribution sheets
  genelists/                  <- reusable gene sets; the control sets are published
                                 to docs/genelists/ as browser buttons
  outputs/                    <- the 108 generated tables (tissue-major)
  gene_major/                 <- the same rows, gene-major: 4,665 shards
  hist_major/                 <- just the histograms, same sharding: ~5 KB each
  diagrams/                   <- 108 distribution sheets, 5 histograms each
  results/                    <- extracted gene subsets
```

Requires Python with `numpy`, `pandas`, `scipy`, `matplotlib`.

## Regenerating

```bash
cd fourparamsacrosstissues/fourparam
python generate_all.py          # all 108 tables; resumable, skips existing
python build_gene_major.py      # re-orient them into gene_major/ (~80 s)
python build_hist_major.py      # mirror the histograms into hist_major/
python build_gui_data.py --reference uterus   # manifest + gene index + gene sets
python make_diagrams.py         # 108 distribution sheets -> diagrams/
```

`build_gui_data.py` verifies that the gene set is identical across tables before
writing the index — the GUI depends on that being true.

## Notes on the data

- **`data/*.csv.gz` are Git LFS pointers.** The 54 matrices are ~5.7 GB of LFS
  payload and are not downloaded by default. Fetch one deliberately:
  ```bash
  git lfs pull --include="fourparamsacrosstissues/data/v11_log2_liver.csv.gz"
  ```
- Values are `log2(TPM + 1) − 1`, so **TPM = 0 maps to exactly `−1`** and every
  value is `≥ −1`. That is what the "excluded ≤ −1" table drops: the
  zero-expression samples, and nothing else.
- `outputs/` and `gene_major/` are tracked as ordinary git rather than LFS on
  purpose — these CSVs compress ~2.5× in git (2.05 GB → ~852 MB, and 2.46 GB →
  ~1.0 GB), whereas LFS stores blobs uncompressed.
- `gene_major/` is a **re-orientation, not new data**. It holds exactly the same
  8,059,824 rows as `outputs/`, byte for byte;
  `python build_gene_major.py --verify` re-derives a sample and checks that.
  Regenerate it whenever you regenerate a table, or the browser will serve stale
  numbers for that tissue.
- Nothing in this pipeline re-serialises a float. `extract_genes.py` reads with
  `dtype=str` and `build_gene_major.py` concatenates row text, because parsing to
  float64 and writing back is *not* lossless: pandas' CSV writer emits ~16
  significant digits rather than the shortest round-tripping repr, which turned
  `0.012596832467784065` into `0.012596832467784`. That is why the browser's
  export and the CLI's agree exactly.
- Tables keep **every** gene, including failures (`fit_success = False`, metrics
  `NaN`). Analysis-time filters such as `fit_success == True`,
  `0 < truncationindex < 1`, and `n_obs >= 30` are deliberately *not* baked in.
  The GUI shows unfiltered rows, so degenerate fits are visible rather than
  hidden — e.g. the pseudogene ALDH7A1P2 in liver reports σ-dist ≈ 292 with
  `w ≈ 0.0007`.
