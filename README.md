# bhuvanlab

Distribution-shape and truncation-index analysis of human gene expression across
the **54 GTEx v11 tissues**.

## → [Open the Truncation Browser](https://bhuvankanna.github.io/bhuvanlab/)

**[bhuvankanna.github.io/bhuvanlab](https://bhuvankanna.github.io/bhuvanlab/)** —
search genes *and* tissues by type-ahead, load the tables, then tick rows into a
**working set** that persists across queries and exports as one spreadsheet.
Every row draws that gene's own fitted curve with its truncation ceiling marked.

Two tabs:

- **Gene selection** — pick any number of genes and any number of tissues (up to
  all 54), choose raw or excluded ≤ −1, and press **Load tables**. Results are
  one row per gene × tissue, each with a checkbox.
- **Working set** — everything ticked so far, gathered across as many separate
  queries as you like, filterable and exportable to CSV.

- Type-ahead over all **74,628** genes (`ALDH2`, `ENSG00000111275`) and all 54
  tissues (`whole blood` matches `whole_blood`)
- **Ensembl versions are ignored**, so ids from an older annotation still resolve
- The working set survives a reload, and mixes tissues *and* raw/excluded rows
- Its CSV has the **same 22 columns** `extract_genes.py` writes, so GUI and CLI
  output are interchangeable
- Tables are read live from this repository, so the page always matches `outputs/`

> Loading is deliberately an explicit button, not automatic: a table is ~8 MB
> over the wire, so "every tissue" is a ~432 MB errand. The page states the cost
> before you commit to it, streams the tables one at a time, and lets you cancel
> mid-run while keeping whatever already arrived.

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
  manifest.json               <- tissue list + column list
  genes.tsv                   <- gene index for type-ahead
fourparamsacrosstissues/
  CLAUDE.md                   <- full spec: columns, filters, how tables are made
  data/                       <- 54 input matrices (Git LFS pointers, see below)
  fourparam/                  <- the code
    bhuvanfitter.py           <- the fit library; single source of truth
    generate_fourparam.py     <- one table from one matrix
    generate_all.py           <- all 108 tables
    run_cluster.py            <- the same, spread over several machines
    extract_genes.py          <- pull a gene set out of the tables
    build_gui_data.py         <- regenerate docs/manifest.json + docs/genes.tsv
  genelists/                  <- reusable gene sets
  outputs/                    <- the 108 generated tables
  results/                    <- extracted gene subsets
```

Requires Python with `numpy`, `pandas`, `scipy`, `matplotlib`.

## Regenerating

```bash
cd fourparamsacrosstissues/fourparam
python generate_all.py          # all 108 tables; resumable, skips existing
python build_gui_data.py        # refresh the GUI's manifest + gene index
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
- `outputs/` is tracked as ordinary git rather than LFS on purpose — these CSVs
  compress 2.46× in git (2.05 GB → ~852 MB), whereas LFS stores blobs
  uncompressed.
- Tables keep **every** gene, including failures (`fit_success = False`, metrics
  `NaN`). Analysis-time filters such as `fit_success == True`,
  `0 < truncationindex < 1`, and `n_obs >= 30` are deliberately *not* baked in.
  The GUI shows unfiltered rows, so degenerate fits are visible rather than
  hidden — e.g. the pseudogene ALDH7A1P2 in liver reports σ-dist ≈ 292 with
  `w ≈ 0.0007`.
