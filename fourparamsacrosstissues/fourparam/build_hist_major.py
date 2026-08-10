#!/usr/bin/env python3
"""Publish the thumbnail histograms on their own, sharded by gene.

The browser's gene page needs exactly one thing the gene-major shards cannot
give it: the 40-bin ``hist`` / ``hist_max`` pair for one gene in one tissue under
one filter. Clicking a gene names all three, so the fetch should be that one
row -- but ``hist`` lives only in the tissue-major tables, and reaching into one
of those costs ~8 MB to read 43 bytes.

So the two columns are mirrored here on their own, sharded with **the same rule
and the same shard numbers** as ``gene_major/``. A gene page therefore fetches
``hist_major/shard_NNNN.csv`` -- about 4.5 KB -- alongside the gene-major shard
it was already fetching, and draws real bars with no large download anywhere in
the path.

Why a sidecar rather than adding the columns to ``gene_major/``:

* ``gene_major/`` is 2.46 GB and every shard would be rewritten to carry two
  columns that are empty for 104 of the 108 tables -- a ~1 GB push to publish
  ~21 MB of actual data.
* ``SHARD_HEADER`` is one of four hardcoded column lists that must agree, and it
  underwrites a byte-identity guarantee with ``extract_genes.py``. This adds
  nothing to that contract.
* Tables without ``hist`` contribute no rows at all, so the sidecar is sized by
  what has actually been generated rather than by the gene set.

``min`` / ``max`` are deliberately not repeated here: they are the histogram's
bin edges, but the gene page already holds them from the gene-major shard.

Run after regenerating any table that carries ``hist``::

    python build_hist_major.py
    python build_hist_major.py --verify
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from build_gene_major import (EXCLUDED_KIND, EXCLUDED_SUFFIX, PREFIX, RAW_KIND,
                              RAW_SUFFIX, SHARD_SIZE, gene_order)

HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUTS = HERE.parent / "outputs"
DEFAULT_OUT = HERE.parent / "hist_major"

HEADER = "tissue,table,gene,hist,hist_max"
FILE_PATTERN = "shard_{shard:04d}.csv"


def tables_with_hist(outputs: Path):
    """(path, tissue, kind) for every table that actually carries the columns."""
    found = []
    for suffix, kind in ((RAW_SUFFIX, RAW_KIND), (EXCLUDED_SUFFIX, EXCLUDED_KIND)):
        for path in sorted(outputs.glob(f"{PREFIX}*{suffix}")):
            with path.open("r", encoding="utf-8") as fh:
                cols = fh.readline().rstrip("\r\n").split(",")
            if "hist" in cols and "hist_max" in cols:
                found.append((path, path.name[len(PREFIX):-len(suffix)], kind))
    return found


def build(outputs: Path, out_dir: Path, reference: Path) -> int:
    tables = tables_with_hist(outputs)
    if not tables:
        print(f"No table in {outputs} carries hist/hist_max -- nothing to "
              f"publish. Regenerate at least one table first.", file=sys.stderr)
        return 1

    print(f"Tables with hist : {len(tables)}", file=sys.stderr)
    for _, tissue, kind in tables:
        print(f"  {tissue} / {kind}", file=sys.stderr)

    shard_of = {gid: i // SHARD_SIZE for i, gid in enumerate(gene_order(reference))}
    n_shards = (len(shard_of) + SHARD_SIZE - 1) // SHARD_SIZE

    # shard -> list of row strings. ~21 MB total at 4 tables, so held in memory
    # rather than streamed; a full 108-table set would want the bucketed
    # two-pass approach build_gene_major.py uses.
    buckets: dict[int, list[str]] = {}
    for path, tissue, kind in tables:
        # dtype=str and no NA handling: `hist` is text and `hist_max` is an
        # integer that must not become "11.0". Nothing here re-serialises.
        df = pd.read_csv(path, usecols=["gene", "hist", "hist_max"],
                         dtype=str, keep_default_na=False)
        kept = 0
        for gene, hist, hist_max in zip(df["gene"], df["hist"], df["hist_max"]):
            # No histogram is a real state (n_obs = 0), not an error. Skipping
            # keeps the sidecar to rows that can actually be drawn.
            if not hist:
                continue
            shard = shard_of.get(gene)
            if shard is None:
                continue
            buckets.setdefault(shard, []).append(
                f"{tissue},{kind},{gene},{hist},{hist_max}")
            kept += 1
        print(f"  {path.name}: {kept:,} histograms", file=sys.stderr)

    out_dir.mkdir(parents=True, exist_ok=True)
    written = total = 0
    for shard in range(n_shards):
        rows = buckets.get(shard)
        if not rows:
            continue
        target = out_dir / FILE_PATTERN.format(shard=shard)
        target.write_text(HEADER + "\n" + "\n".join(rows) + "\n",
                          encoding="utf-8", newline="\n")
        written += 1
        total += len(rows)

    size = sum(p.stat().st_size for p in out_dir.glob("shard_*.csv"))
    print(f"\nWrote {written:,} shards, {total:,} rows, {size / 1e6:.1f} MB "
          f"({size / max(written, 1) / 1024:.1f} KB per shard)", file=sys.stderr)
    return 0


def verify(outputs: Path, out_dir: Path) -> int:
    """Re-read a sample of shards and check each row against its source table."""
    tables = {(t, k): p for p, t, k in tables_with_hist(outputs)}
    shards = sorted(out_dir.glob("shard_*.csv"))
    if not shards:
        print(f"No shards in {out_dir}", file=sys.stderr)
        return 1

    sample = shards[::max(1, len(shards) // 20)][:20]
    cache: dict[tuple, dict] = {}
    checked = bad = 0
    for path in sample:
        for line in path.read_text(encoding="utf-8").split("\n")[1:]:
            if not line:
                continue
            tissue, kind, gene, hist, hist_max = line.split(",")
            key = (tissue, kind)
            if key not in cache:
                df = pd.read_csv(tables[key], usecols=["gene", "hist", "hist_max"],
                                 dtype=str, keep_default_na=False)
                cache[key] = {g: (h, m) for g, h, m
                              in zip(df["gene"], df["hist"], df["hist_max"])}
            want = cache[key].get(gene)
            checked += 1
            if want != (hist, hist_max):
                bad += 1
                print(f"  MISMATCH {gene} {tissue}/{kind}: {want} != "
                      f"{(hist, hist_max)}", file=sys.stderr)

    print(f"Checked {checked:,} rows across {len(sample)} shards: "
          f"{'OK - all match' if not bad else f'{bad} MISMATCHES'}",
          file=sys.stderr)
    return 1 if bad else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outputs", type=Path, default=DEFAULT_OUTPUTS)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--reference", type=str, default=None,
                    help="Tissue whose raw table fixes the gene order "
                         "(default: first alphabetically, same as gene_major).")
    ap.add_argument("--verify", action="store_true",
                    help="Byte-compare a sample of shards against outputs/.")
    args = ap.parse_args(argv)

    if not args.outputs.is_dir():
        ap.error(f"outputs directory not found: {args.outputs}")
    if args.verify:
        return verify(args.outputs, args.out)

    raw = sorted(args.outputs.glob(f"{PREFIX}*{RAW_SUFFIX}"))
    if not raw:
        ap.error(f"no raw tables in {args.outputs}")
    reference = raw[0]
    if args.reference:
        match = [p for p in raw if p.name == f"{PREFIX}{args.reference}{RAW_SUFFIX}"]
        if not match:
            ap.error(f"no raw table for --reference {args.reference!r}")
        reference = match[0]
    return build(args.outputs, args.out, reference)


if __name__ == "__main__":
    raise SystemExit(main())
