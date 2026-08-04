#!/usr/bin/env python3
"""Extract fourparam statistics for a set of genes across tissues.

Pulls the rows matching a list of genes out of the generated fourparam tables in
``outputs/`` and writes them as one tidy (long-format) CSV: one row per
gene x tissue x table, carrying every statistic column.

Genes may be given as
  * a versioned Ensembl id     ENSG00000142192.22
  * an unversioned Ensembl id  ENSG00000142192      <- version is ignored
  * a gene symbol              APP
  * a glob pattern             ALDH*, ADH1?, ENSG000001421*

**Versions are always ignored when matching.** GTEx bumps the ``.NN`` suffix
between releases, so a gene list written against an older annotation (e.g.
``ENSG00000142192.20``) still resolves against these tables (which carry
``.22``). Matching on the full versioned string would silently return nothing.

Every query token is reported as resolved or unresolved, so a typo or a symbol
that is absent from the annotation is loud rather than silently empty.

Examples
--------
Explicit genes plus every ALDH, all 54 tissues, both tables::

    python extract_genes.py --genes APP SNCA PCSK9 "ALDH*" -o ../results/genes.csv

From a file, one token per line, raw tables only::

    python extract_genes.py --genes-file mygenes.txt --table raw -o out.csv

A couple of tissues only::

    python extract_genes.py --genes "ADH*" --tissues liver,lung -o out.csv
"""

from __future__ import annotations

import argparse
import fnmatch
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

# outputs/ lives next to fourparam/
DEFAULT_OUTPUTS = Path(__file__).resolve().parent.parent / "outputs"

RAW_SUFFIX = "_fourparam.csv"
EXCLUDED_SUFFIX = "_fourparam_excluded_at_or_below_-1.csv"
PREFIX = "v11_log2_"

# Statistic columns, in table order. gene/genename are handled separately.
STAT_COLUMNS = [
    "y0", "A", "x0", "w", "sumsquarevalue",
    "ti_fourparam_sigma_dist", "truncationindex",
    "min", "max", "mean", "std", "skew", "kurt",
    "right", "maxheight", "rightheight", "n_obs", "fit_success",
]


def strip_version(gene_id: str) -> str:
    """ENSG00000142192.22 -> ENSG00000142192 (leaves non-Ensembl strings alone)."""
    return re.sub(r"\.\d+$", "", gene_id.strip())


def discover_tables(outputs: Path, tissues=None, table="both"):
    """Find the (tissue, kind, path) triples to read."""
    found = []
    for path in sorted(outputs.glob(f"{PREFIX}*{RAW_SUFFIX}")):
        # The excluded tables also end in _fourparam.csv? No -- distinct suffix.
        tissue = path.name[len(PREFIX):-len(RAW_SUFFIX)]
        found.append((tissue, "raw", path))
    for path in sorted(outputs.glob(f"{PREFIX}*{EXCLUDED_SUFFIX}")):
        tissue = path.name[len(PREFIX):-len(EXCLUDED_SUFFIX)]
        found.append((tissue, "excluded_at_or_below_-1", path))

    if table == "raw":
        found = [t for t in found if t[1] == "raw"]
    elif table == "excluded":
        found = [t for t in found if t[1] != "raw"]

    if tissues:
        wanted = {t.strip().lower() for t in tissues}
        found = [t for t in found if t[0].lower() in wanted]

    return sorted(found)


def build_gene_lookup(reference_table: Path):
    """Read the gene/genename columns once.

    Safe to do from a single table because the gene set is identical across all
    of them (every gene is written to every table, even when the fit fails).
    """
    df = pd.read_csv(reference_table, usecols=["gene", "genename"], dtype=str)
    df["gene"] = df["gene"].fillna("")
    df["genename"] = df["genename"].fillna("")
    df["gene_nover"] = df["gene"].map(strip_version)
    return df


def resolve_genes(lookup: pd.DataFrame, tokens):
    """Map query tokens onto rows of the annotation.

    Returns (matched_versioned_ids, report) where report maps each token to the
    list of gene ids it resolved to (empty list = unresolved).
    """
    by_nover = {}
    for gid_nover, gid in zip(lookup["gene_nover"], lookup["gene"]):
        by_nover.setdefault(gid_nover.upper(), []).append(gid)

    by_symbol = {}
    for sym, gid in zip(lookup["genename"], lookup["gene"]):
        if sym:
            by_symbol.setdefault(sym.upper(), []).append(gid)

    all_symbols = list(by_symbol.keys())
    all_novers = list(by_nover.keys())

    matched, report = [], {}
    for token in tokens:
        raw = token.strip()
        if not raw:
            continue
        key = strip_version(raw).upper()
        hits = []

        if any(ch in raw for ch in "*?["):
            # Glob: try symbols first, then unversioned Ensembl ids.
            pattern = raw.upper()
            for sym in fnmatch.filter(all_symbols, pattern):
                hits.extend(by_symbol[sym])
            for nov in fnmatch.filter(all_novers, strip_version(pattern)):
                hits.extend(by_nover[nov])
        else:
            hits.extend(by_nover.get(key, []))
            if not hits:
                hits.extend(by_symbol.get(key, []))

        # de-duplicate, keep order
        seen, unique_hits = set(), []
        for gid in hits:
            if gid not in seen:
                seen.add(gid)
                unique_hits.append(gid)

        report[raw] = unique_hits
        matched.extend(unique_hits)

    seen, unique = set(), []
    for gid in matched:
        if gid not in seen:
            seen.add(gid)
            unique.append(gid)
    return unique, report


def _read_one(job):
    """Worker: pull the wanted genes out of one table.

    Read as **text**, never as floats. Parsing to float64 and writing back is not
    lossless here: pandas' CSV writer emits ~16 significant digits rather than
    the shortest round-tripping repr, so ``0.012596832467784065`` in the table
    came back out as ``0.012596832467784``. Keeping every field a string makes
    the extract byte-identical to ``outputs/`` - and therefore to what the
    browser exports, which reads the same text.
    """
    tissue, kind, path, wanted = job
    try:
        df = pd.read_csv(path, dtype=str, keep_default_na=False, na_values=[])
    except Exception as exc:  # noqa: BLE001 - report and carry on
        return tissue, kind, None, f"{path.name}: {exc}"

    sub = df[df["gene"].isin(wanted)].copy()
    sub.insert(0, "tissue", tissue)
    sub.insert(1, "table", kind)
    return tissue, kind, sub, None


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Extract fourparam statistics for a set of genes across tissues.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--genes", nargs="+", default=[],
                    help="Gene symbols, Ensembl ids (version ignored), or globs.")
    ap.add_argument("--genes-file", type=Path,
                    help="File with one gene token per line ('#' comments allowed).")
    ap.add_argument("--outputs", type=Path, default=DEFAULT_OUTPUTS,
                    help="Directory holding the fourparam tables (default: ../outputs).")
    ap.add_argument("--tissues", default=None,
                    help="Comma-separated tissue names. Default: all 50.")
    ap.add_argument("--table", choices=["both", "raw", "excluded"], default="both",
                    help="Which table type to read (default: both).")
    ap.add_argument("-o", "--out", type=Path, required=True,
                    help="Destination CSV (long format).")
    ap.add_argument("--jobs", type=int, default=8, help="Parallel readers (default 8).")
    args = ap.parse_args(argv)

    tokens = list(args.genes)
    if args.genes_file:
        for line in args.genes_file.read_text(encoding="utf-8").splitlines():
            line = line.split("#", 1)[0].strip()
            if line:
                tokens.append(line)
    if not tokens:
        ap.error("no genes given; use --genes and/or --genes-file")

    if not args.outputs.is_dir():
        ap.error(f"outputs directory not found: {args.outputs}")

    tissue_list = args.tissues.split(",") if args.tissues else None
    tables = discover_tables(args.outputs, tissue_list, args.table)
    if not tables:
        ap.error("no matching tables found in " + str(args.outputs))

    print(f"Tables to read : {len(tables)}", file=sys.stderr)

    # ---- resolve the gene list against the annotation -----------------------
    lookup = build_gene_lookup(tables[0][2])
    wanted, report = resolve_genes(lookup, tokens)

    symbol_of = dict(zip(lookup["gene"], lookup["genename"]))

    print(f"Query tokens   : {len(report)}", file=sys.stderr)
    unresolved = [tok for tok, hits in report.items() if not hits]
    for token, hits in report.items():
        if hits:
            shown = ", ".join(f"{symbol_of.get(g, '?')} ({g})" for g in hits[:6])
            if len(hits) > 6:
                shown += f", +{len(hits) - 6} more"
            print(f"  {token:<22} -> {len(hits):>3}  {shown}", file=sys.stderr)
    for token in unresolved:
        print(f"  {token:<22} -> NOT FOUND", file=sys.stderr)

    if not wanted:
        print("\nNothing resolved; no output written.", file=sys.stderr)
        return 1
    print(f"Genes resolved : {len(wanted)}\n", file=sys.stderr)

    # ---- read the tables ----------------------------------------------------
    wanted_set = set(wanted)
    jobs = [(t, k, p, wanted_set) for t, k, p in tables]
    frames, errors = [], []

    with ProcessPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        for i, (tissue, kind, sub, err) in enumerate(pool.map(_read_one, jobs), 1):
            if err:
                errors.append(err)
            elif sub is not None and not sub.empty:
                frames.append(sub)
            print(f"\r  read {i}/{len(jobs)} tables", end="", file=sys.stderr)
    print(file=sys.stderr)

    for err in errors:
        print(f"  WARNING {err}", file=sys.stderr)

    if not frames:
        print("No rows matched in any table; no output written.", file=sys.stderr)
        return 1

    out = pd.concat(frames, ignore_index=True)
    out = out.sort_values(["genename", "gene", "tissue", "table"], kind="stable")

    cols = ["tissue", "table", "gene", "genename"] + \
           [c for c in STAT_COLUMNS if c in out.columns]
    out = out[cols]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    # Pin LF. On Windows an unpinned writer emits CRLF, which would make this
    # file differ from outputs/ (and from the browser's export) by a stray \r.
    out.to_csv(args.out, index=False, lineterminator="\n")

    n_fit = int(out["fit_success"].astype(str).str.lower().eq("true").sum()) \
        if "fit_success" in out.columns else -1

    print(f"Wrote {len(out):,} rows x {len(cols)} cols -> {args.out}", file=sys.stderr)
    print(f"  genes {out['gene'].nunique()} | tissues {out['tissue'].nunique()} "
          f"| tables {out['table'].nunique()}", file=sys.stderr)
    if n_fit >= 0:
        print(f"  fit_success=True rows: {n_fit:,} / {len(out):,}", file=sys.stderr)
    if unresolved:
        print(f"  unresolved tokens: {', '.join(unresolved)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
