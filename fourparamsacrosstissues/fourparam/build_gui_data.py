#!/usr/bin/env python3
"""Generate the small static files the browser GUI needs.

The GUI does NOT get its own copy of the fourparam tables. It fetches the real
tables straight from raw.githubusercontent.com at runtime, so it can never drift
out of sync with ``outputs/`` and the repo does not grow by ~540 MB of duplicated
data. All this script produces is:

  docs/manifest.json  tissue list, column list, and the two data base URLs
  docs/genes.tsv      "<versioned ensembl id>\\t<symbol>\\t<shard>" per line

``genes.tsv`` is built from a single table because the gene set is byte-identical
across all 108 of them (every gene is written to every table, even when its fit
fails). The script asserts that rather than trusting it.

The third column is the **gene-major shard** that holds that gene's rows for
every tissue, as laid out by ``build_gene_major.py``. The GUI uses it to answer
a few-genes-many-tissues query with a handful of ~220 KB fetches instead of one
~8 MB tissue table per tissue. It is recomputed here from the same sort rule
rather than read back from ``gene_major/``, and then checked against the shards
actually on disk, so the index and the shards cannot silently disagree.

Run after adding tissues or changing the table columns::

    python build_gui_data.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from build_gene_major import (EXCLUDED_KIND, RAW_KIND, SHARD_SIZE, gene_order)

REPO_SUBDIR = "fourparamsacrosstissues"
PREFIX = "v11_log2_"
RAW_SUFFIX = "_fourparam.csv"
EXCLUDED_SUFFIX = "_fourparam_excluded_at_or_below_-1.csv"

HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUTS = HERE.parent / "outputs"
DEFAULT_DOCS = HERE.parent.parent / "docs"      # repo-root/docs

RAW_BASE = "https://raw.githubusercontent.com/BhuvanKanna/bhuvanlab/main"
DATA_BASE = f"{RAW_BASE}/{REPO_SUBDIR}/outputs"
GENE_MAJOR_BASE = f"{RAW_BASE}/{REPO_SUBDIR}/gene_major"

DEFAULT_GENE_MAJOR = HERE.parent / "gene_major"


def gene_column_digest(path: Path) -> str:
    """md5 of the gene column, used to prove the gene set matches across tables."""
    col = pd.read_csv(path, usecols=["gene"], dtype=str)["gene"].fillna("")
    return hashlib.md5("\n".join(col).encode("utf-8")).hexdigest()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outputs", type=Path, default=DEFAULT_OUTPUTS)
    ap.add_argument("--docs", type=Path, default=DEFAULT_DOCS)
    ap.add_argument("--gene-major", type=Path, default=DEFAULT_GENE_MAJOR,
                    help="Directory of gene-major shards (default: ../gene_major).")
    ap.add_argument("--check-tables", type=int, default=4,
                    help="How many tables to verify share the gene set (0 = all).")
    args = ap.parse_args(argv)

    if not args.outputs.is_dir():
        ap.error(f"outputs directory not found: {args.outputs}")

    raw_tables = sorted(args.outputs.glob(f"{PREFIX}*{RAW_SUFFIX}"))
    exc_tables = sorted(args.outputs.glob(f"{PREFIX}*{EXCLUDED_SUFFIX}"))
    if not raw_tables:
        ap.error(f"no tables matching {PREFIX}*{RAW_SUFFIX} in {args.outputs}")

    tissues = [p.name[len(PREFIX):-len(RAW_SUFFIX)] for p in raw_tables]
    exc_tissues = {p.name[len(PREFIX):-len(EXCLUDED_SUFFIX)] for p in exc_tables}

    missing_exc = [t for t in tissues if t not in exc_tissues]
    if missing_exc:
        print(f"  WARNING: {len(missing_exc)} tissue(s) have no excluded table: "
              f"{', '.join(missing_exc[:5])}", file=sys.stderr)

    print(f"Tissues        : {len(tissues)}", file=sys.stderr)
    print(f"Raw tables     : {len(raw_tables)}", file=sys.stderr)
    print(f"Excluded tables: {len(exc_tables)}", file=sys.stderr)

    # ---- verify the gene set really is shared -------------------------------
    to_check = raw_tables + exc_tables
    if args.check_tables > 0:
        step = max(1, len(to_check) // args.check_tables)
        to_check = to_check[::step][:args.check_tables]
    print(f"\nVerifying gene set across {len(to_check)} table(s) ...", file=sys.stderr)

    digests = {}
    for path in to_check:
        digests.setdefault(gene_column_digest(path), []).append(path.name)
    if len(digests) != 1:
        print("  ERROR: tables do not share an identical gene column:", file=sys.stderr)
        for digest, names in digests.items():
            print(f"    {digest[:12]}  {len(names)} table(s), e.g. {names[0]}",
                  file=sys.stderr)
        print("  The GUI relies on one shared gene index; aborting.", file=sys.stderr)
        return 1
    print(f"  OK - single gene set ({next(iter(digests))[:12]})", file=sys.stderr)

    # ---- gene index + gene-major shard map ----------------------------------
    reference = raw_tables[0]
    df = pd.read_csv(reference, usecols=["gene", "genename"], dtype=str)
    df["gene"] = df["gene"].fillna("")
    df["genename"] = df["genename"].fillna("")

    # Recompute the shard assignment from build_gene_major's own sort rule, so
    # there is one source of truth for it.
    ordered = gene_order(reference)
    shard_of_gene = {gid: i // SHARD_SIZE for i, gid in enumerate(ordered)}
    n_shards = (len(ordered) + SHARD_SIZE - 1) // SHARD_SIZE

    on_disk = sorted(args.gene_major.glob("shard_*.csv")) \
        if args.gene_major.is_dir() else []
    if not on_disk:
        print(f"\n  WARNING: no shards in {args.gene_major} - run build_gene_major.py.\n"
              f"           The GUI will fall back to tissue-major tables for every "
              f"query.", file=sys.stderr)
        have_gene_major = False
    elif len(on_disk) != n_shards:
        print(f"\n  ERROR: {args.gene_major} holds {len(on_disk):,} shards but the "
              f"index expects {n_shards:,}.", file=sys.stderr)
        print(f"         Re-run build_gene_major.py; aborting rather than "
              f"publishing an index that points at missing files.", file=sys.stderr)
        return 1
    else:
        have_gene_major = True
        print(f"\nGene-major     : {len(on_disk):,} shards, "
              f"{sum(p.stat().st_size for p in on_disk) / 1e9:.2f} GB",
              file=sys.stderr)

    args.docs.mkdir(parents=True, exist_ok=True)
    genes_path = args.docs / "genes.tsv"
    with genes_path.open("w", encoding="utf-8", newline="\n") as fh:
        for gid, sym in zip(df["gene"], df["genename"]):
            fh.write(f"{gid}\t{sym}\t{shard_of_gene[gid]}\n")

    size_kb = genes_path.stat().st_size / 1024
    print(f"Wrote {genes_path.name}: {len(df):,} genes, {size_kb:,.0f} KB "
          f"(id, symbol, shard)", file=sys.stderr)

    # ---- manifest -----------------------------------------------------------
    columns = list(pd.read_csv(reference, nrows=0).columns)

    manifest = {
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "data_base_url": DATA_BASE,
        "file_prefix": PREFIX,
        "table_kinds": [
            {"id": "raw",
             "label": "Raw (every finite value)",
             "suffix": RAW_SUFFIX},
            {"id": "excluded",
             "label": "Excluded ≤ −1 (drops zero-expression samples)",
             "suffix": EXCLUDED_SUFFIX},
        ],
        "tissues": tissues,
        "columns": columns,
        "n_genes": int(len(df)),
        # Gene-major mirror of the same numbers: one shard covers SHARD_SIZE
        # genes across every tissue and both filters.
        "gene_major": {
            "available": have_gene_major,
            "base_url": GENE_MAJOR_BASE,
            "shard_size": SHARD_SIZE,
            "n_shards": n_shards,
            "file_pattern": "shard_{shard:04d}.csv",
            "kind_ids": {"raw": RAW_KIND, "excluded": EXCLUDED_KIND},
        },
    }
    manifest_path = args.docs / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {manifest_path.name}: {len(tissues)} tissues, "
          f"{len(columns)} columns", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
