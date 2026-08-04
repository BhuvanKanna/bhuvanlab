#!/usr/bin/env python3
"""Re-orient the fourparam tables from tissue-major to gene-major.

``outputs/`` is tissue-major: one file per tissue x filter, every gene inside.
That is the right shape for generating and for whole-tissue work, but it is the
wrong shape for the browser, where a query is "these few genes, across many
tissues". Answering that from tissue-major files means downloading a whole
~8 MB table per tissue - 54 tissues is a 432 MB errand to surface a few hundred
rows.

This script writes the transpose into ``gene_major/``: each shard holds every
row (all 54 tissues x both filters = 108 rows per gene) for a contiguous block
of ``SHARD_SIZE`` genes. Looking up one gene across every tissue becomes a
single ~220 KB fetch instead of 432 MB.

Exactness
---------
Shard rows are built by **text concatenation, never by parsing**: a table row
after ``gene,genename`` is already exactly ``extract_genes.py``'s STAT_COLUMNS
in order, so a shard row is literally ``f"{tissue},{kind},{original_line}"``.
No float is ever re-formatted, so the browser's CSV export stays byte-identical
to what the CLI writes. ``--verify`` re-checks that against the source tables.

Gene order
----------
Genes are sorted by symbol (then by id, and symbol-less genes last), so a family
lands in adjacent shards: every ``ALDH*`` is 2-3 fetches, not 19 scattered ones.
The resulting shard id per gene is published in ``docs/genes.tsv`` by
``build_gui_data.py``.

Memory
------
Two passes with bounded memory, rather than holding 2.2 GB of transpose at once:

  pass 1  each table is streamed once and split into SUPER_BUCKETS temp files
  pass 2  each temp bucket (~35 MB) is loaded, grouped, and written as shards

Usage::

    python build_gene_major.py                 # full build, 10 processes
    python build_gene_major.py --verify        # check shards against outputs/
    python build_gene_major.py --limit-tables 4   # quick smoke test
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import shutil
import sys
import time
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUTS = HERE.parent / "outputs"
DEFAULT_GENE_MAJOR = HERE.parent / "gene_major"

PREFIX = "v11_log2_"
RAW_SUFFIX = "_fourparam.csv"
EXCLUDED_SUFFIX = "_fourparam_excluded_at_or_below_-1.csv"

RAW_KIND = "raw"
EXCLUDED_KIND = "excluded_at_or_below_-1"      # matches extract_genes.py

SHARD_SIZE = 16          # genes per shard -> ~4,665 shards, ~220 KB gzipped each
SUPER_BUCKETS = 64       # pass-1 fan-out; bounds pass-2 memory

# The shard header. Identical to extract_genes.py's output header, so a shard
# and a CLI extract can be concatenated without reconciling columns.
SHARD_HEADER = ("tissue,table,gene,genename,y0,A,x0,w,sumsquarevalue,"
                "ti_fourparam_sigma_dist,truncationindex,min,max,mean,std,"
                "skew,kurt,right,maxheight,rightheight,n_obs,fit_success\n")

TABLE_HEADER = SHARD_HEADER[len("tissue,table,"):]


def discover(outputs: Path):
    """(tissue, kind, path) for every table, excluded first so names are unique."""
    found, seen = [], set()
    for path in sorted(outputs.glob(f"{PREFIX}*{EXCLUDED_SUFFIX}")):
        found.append((path.name[len(PREFIX):-len(EXCLUDED_SUFFIX)], EXCLUDED_KIND, path))
        seen.add(path.name)
    for path in sorted(outputs.glob(f"{PREFIX}*{RAW_SUFFIX}")):
        if path.name not in seen:
            found.append((path.name[len(PREFIX):-len(RAW_SUFFIX)], RAW_KIND, path))
    return sorted(found, key=lambda f: (f[0], f[1]))


def gene_order(reference: Path):
    """Sort genes by symbol so families cluster; return the ordered id list."""
    df = pd.read_csv(reference, usecols=["gene", "genename"], dtype=str)
    df["gene"] = df["gene"].fillna("")
    df["genename"] = df["genename"].fillna("")
    # "~" sorts after every letter, so symbol-less genes land at the end.
    df["_key"] = df["genename"].str.upper().where(df["genename"] != "", "~")
    df = df.sort_values(["_key", "gene"], kind="stable")
    return list(df["gene"])


# --------------------------------------------------------------------------- #
# pass 1 - split each table into SUPER_BUCKETS temp files
# --------------------------------------------------------------------------- #
def _split_table(job):
    tissue, kind, path, tmp_dir, super_of_gene, n_super = job
    prefix = f"{tissue},{kind},"
    handles = [(tmp_dir / f"{path.stem}_b{b:02d}.part").open(
        "w", encoding="utf-8", newline="") for b in range(n_super)]
    n_rows, unknown = 0, 0
    try:
        with path.open("r", encoding="utf-8", newline="") as fh:
            header = fh.readline()
            if header.replace("\r\n", "\n") != TABLE_HEADER:
                return False, f"{path.name}: unexpected header", 0
            for line in fh:
                if not line.strip():
                    continue
                # Normalise the terminator. A table written on Windows by an
                # older pandas carries CRLF, and copying that verbatim would
                # leave shards with mixed line endings and a stray \r glued to
                # the last field.
                line = line.rstrip("\r\n") + "\n"
                gene = line[:line.index(",")]
                bucket = super_of_gene.get(gene)
                if bucket is None:
                    unknown += 1
                    continue
                handles[bucket].write(prefix + line)
                n_rows += 1
    finally:
        for fh in handles:
            fh.close()
    if unknown:
        return False, f"{path.name}: {unknown} rows with a gene absent from the index", n_rows
    return True, f"{path.name}: {n_rows:,} rows", n_rows


# --------------------------------------------------------------------------- #
# pass 2 - turn one temp bucket into its shards
# --------------------------------------------------------------------------- #
def _emit_bucket(job):
    bucket, tmp_dir, dest_dir, shard_of_gene, rank_of_gene, first_shard, last_shard = job
    rows: dict[int, list[tuple]] = {s: [] for s in range(first_shard, last_shard + 1)}
    for part in sorted(tmp_dir.glob(f"*_b{bucket:02d}.part")):
        with part.open("r", encoding="utf-8", newline="") as fh:
            for line in fh:
                if not line.strip():
                    continue
                tissue, kind, rest = line.split(",", 2)
                gene = rest[:rest.index(",")]
                shard = shard_of_gene[gene]
                # sort key: gene position, then tissue, then filter
                rows[shard].append((rank_of_gene[gene], tissue, kind, line))

    written, n_rows = 0, 0
    for shard, items in rows.items():
        if not items:
            continue
        items.sort(key=lambda r: (r[0], r[1], r[2]))
        dest = dest_dir / f"shard_{shard:04d}.csv"
        with dest.open("w", encoding="utf-8", newline="") as fh:
            fh.write(SHARD_HEADER)
            for item in items:
                fh.write(item[3])
        written += 1
        n_rows += len(items)
    return bucket, written, n_rows


# --------------------------------------------------------------------------- #
def verify(outputs: Path, dest_dir: Path, shard_ids: list[int], n_tables: int) -> int:
    """Re-read a few shards and confirm every row is byte-identical to source."""
    tables = {(t, k): p for t, k, p in discover(outputs)}
    problems = 0
    for shard in shard_ids:
        path = dest_dir / f"shard_{shard:04d}.csv"
        if not path.exists():
            print(f"  shard {shard:04d}: MISSING", file=sys.stderr)
            problems += 1
            continue
        lines = path.read_text(encoding="utf-8").splitlines()
        if lines[0] + "\n" != SHARD_HEADER:
            print(f"  shard {shard:04d}: header mismatch", file=sys.stderr)
            problems += 1
        by_table: dict[tuple, dict] = {}
        for line in lines[1:]:
            tissue, kind, rest = line.split(",", 2)
            by_table.setdefault((tissue, kind), {})[rest[:rest.index(",")]] = rest

        genes = {g for rows in by_table.values() for g in rows}
        if len(by_table) != n_tables:
            print(f"  shard {shard:04d}: {len(by_table)} tables, expected {n_tables}",
                  file=sys.stderr)
            problems += 1

        checked = 0
        for (tissue, kind), rows in by_table.items():
            src = tables.get((tissue, kind))
            if src is None:
                print(f"  shard {shard:04d}: unknown table {tissue}/{kind}",
                      file=sys.stderr)
                problems += 1
                continue
            with src.open("r", encoding="utf-8", newline="") as fh:
                fh.readline()
                for line in fh:
                    gene = line[:line.index(",")]
                    if gene in rows:
                        if rows[gene] != line.rstrip("\r\n"):
                            print(f"  shard {shard:04d}: {tissue}/{kind}/{gene} DIFFERS",
                                  file=sys.stderr)
                            problems += 1
                        checked += 1
        print(f"  shard {shard:04d}: {len(genes)} genes x {len(by_table)} tables, "
              f"{checked:,} rows byte-compared -> "
              f"{'OK' if problems == 0 else 'PROBLEMS'}", file=sys.stderr)
    return problems


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Re-orient the fourparam tables from tissue-major to gene-major.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outputs", type=Path, default=DEFAULT_OUTPUTS)
    ap.add_argument("--dest", type=Path, default=DEFAULT_GENE_MAJOR)
    ap.add_argument("--tmp", type=Path, default=None,
                    help="Scratch dir for pass 1 (default: <dest>/.tmp).")
    ap.add_argument("--jobs", type=int, default=10)
    ap.add_argument("--shard-size", type=int, default=SHARD_SIZE)
    ap.add_argument("--limit-tables", type=int, default=0,
                    help="Smoke test: use only the first N tables.")
    ap.add_argument("--keep-tmp", action="store_true")
    ap.add_argument("--verify", action="store_true",
                    help="Only verify existing shards against outputs/.")
    args = ap.parse_args(argv)

    if not args.outputs.is_dir():
        ap.error(f"outputs directory not found: {args.outputs}")

    tables = discover(args.outputs)
    if args.limit_tables:
        tables = tables[:args.limit_tables]
    if not tables:
        ap.error(f"no fourparam tables in {args.outputs}")

    genes = gene_order(tables[0][2])
    n_shards = (len(genes) + args.shard_size - 1) // args.shard_size
    shard_of_gene = {g: i // args.shard_size for i, g in enumerate(genes)}
    rank_of_gene = {g: i for i, g in enumerate(genes)}
    per_super = (n_shards + SUPER_BUCKETS - 1) // SUPER_BUCKETS
    super_of_gene = {g: shard_of_gene[g] // per_super for g in genes}

    print(f"Tables      : {len(tables)}", file=sys.stderr)
    print(f"Genes       : {len(genes):,}", file=sys.stderr)
    print(f"Shards      : {n_shards:,}  ({args.shard_size} genes each)", file=sys.stderr)

    if args.verify:
        probe = [0, n_shards // 3, n_shards // 2, n_shards - 1]
        print(f"\nVerifying shards {probe} against {args.outputs} ...", file=sys.stderr)
        bad = verify(args.outputs, args.dest, probe, len(tables))
        print(f"\n{'VERIFY FAILED' if bad else 'VERIFY OK - every compared row is byte-identical'}",
              file=sys.stderr)
        return 1 if bad else 0

    tmp_dir = args.tmp or (args.dest / ".tmp")
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True)
    args.dest.mkdir(parents=True, exist_ok=True)
    for stale in args.dest.glob("shard_*.csv"):
        stale.unlink()

    started = time.time()

    # ---- pass 1 -------------------------------------------------------------
    print(f"\nPass 1/2 - splitting {len(tables)} tables into {SUPER_BUCKETS} buckets",
          file=sys.stderr)
    jobs = [(t, k, p, tmp_dir, super_of_gene, SUPER_BUCKETS) for t, k, p in tables]
    total_rows, failures = 0, []
    with mp.Pool(processes=max(1, min(args.jobs, len(jobs)))) as pool:
        for i, (ok, message, n) in enumerate(pool.imap_unordered(_split_table, jobs), 1):
            total_rows += n
            if not ok:
                failures.append(message)
            print(f"\r  {i}/{len(jobs)} tables, {total_rows:,} rows", end="",
                  file=sys.stderr)
    print(file=sys.stderr)
    for message in failures:
        print(f"  ERROR {message}", file=sys.stderr)
    if failures:
        return 1

    expected = len(genes) * len(tables)
    if total_rows != expected:
        print(f"  ERROR row count {total_rows:,} != expected {expected:,}",
              file=sys.stderr)
        return 1
    print(f"  {total_rows:,} rows staged ({time.time() - started:.0f}s)", file=sys.stderr)

    # ---- pass 2 -------------------------------------------------------------
    print(f"\nPass 2/2 - writing {n_shards:,} shards -> {args.dest}", file=sys.stderr)
    jobs = []
    for bucket in range(SUPER_BUCKETS):
        first = bucket * per_super
        last = min(first + per_super - 1, n_shards - 1)
        if first <= last:
            jobs.append((bucket, tmp_dir, args.dest, shard_of_gene, rank_of_gene,
                         first, last))
    n_written, n_rows_out = 0, 0
    with mp.Pool(processes=max(1, min(args.jobs, len(jobs)))) as pool:
        for i, (bucket, written, n) in enumerate(pool.imap_unordered(_emit_bucket, jobs), 1):
            n_written += written
            n_rows_out += n
            print(f"\r  {i}/{len(jobs)} buckets, {n_written:,} shards, "
                  f"{n_rows_out:,} rows", end="", file=sys.stderr)
    print(file=sys.stderr)

    if not args.keep_tmp:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    size = sum(p.stat().st_size for p in args.dest.glob("shard_*.csv"))
    print(f"\nWrote {n_written:,} shards, {n_rows_out:,} rows, {size / 1e9:.2f} GB "
          f"in {time.time() - started:.0f}s", file=sys.stderr)
    if n_written != n_shards or n_rows_out != expected:
        print(f"  ERROR expected {n_shards:,} shards / {expected:,} rows",
              file=sys.stderr)
        return 1
    print(f"\nNow run:  python build_gui_data.py    "
          f"(publishes the gene -> shard index)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
