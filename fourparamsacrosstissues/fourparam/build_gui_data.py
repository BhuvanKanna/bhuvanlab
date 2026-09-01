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
DEFAULT_HIST_MAJOR = HERE.parent / "hist_major"
HIST_MAJOR_BASE = f"{RAW_BASE}/{REPO_SUBDIR}/hist_major"

# C. elegans. Its own directory rather than outputs/, so the tissue pipeline's
# `*_fourparam*.csv` globs cannot sweep a second organism into a GTEx run.
DEFAULT_WORM = HERE.parent / "worm"
DEFAULT_GENELISTS = HERE.parent / "genelists"
WORM_BASE = f"{RAW_BASE}/{REPO_SUBDIR}/worm"
WORM_FILE = "worm_fourparam_excluded_at_or_below_-1.csv"


# One curated row per (gene, disorder) where OMIM/G2P name **over-expression**
# as the driving mechanism. Small enough to publish verbatim; the browser's gene
# page shows it as that gene's phenotype panel.
PHENOTYPES_FILE = "overexpression_phenotypes.tsv"

# Manifest blocks owned by another script, which patches manifest.json in place.
# This script rewrites the file wholesale and must carry these forward verbatim.
PATCHED_BLOCKS = ("qc", "r2")


# Reusable gene sets published beside the page, one token per line, so the
# browser can load a whole roster with one button instead of a paste.
#
# `genelists/` is the source of truth and the same files feed
# `extract_genes.py --genes-file`, so a set is identical at the CLI and in the
# GUI. They are copied into docs/ rather than linked because GitHub Pages only
# serves what is under docs/.
#
# `polarity` is what the control sets are FOR: positive sets should come out
# truncated if the hypothesis holds, negative sets should not. The browser
# colours the buttons by it.
GENESETS_DIR = "genelists"
GENESETS = [
    {"id": "pos_pTriplo",
     "file": "pos_pTriplo.txt",
     "label": "pTriplo > 0.94",
     "polarity": "positive",
     "note": "triplication-sensitive (Collins et al. 2022)"},
    {"id": "pos_decipher",
     "file": "pos_decipher.txt",
     "label": "DECIPHER dominant",
     "polarity": "positive",
     "note": "dominant-mechanism disease genes"},
    {"id": "neg_dup_tolerant",
     "file": "neg_dup_tolerant.txt",
     "label": "Duplication-tolerant",
     "polarity": "negative",
     "note": "tolerated as an extra copy in healthy people"},
    {"id": "neg_olfactory",
     "file": "neg_olfactory.txt",
     "label": "Olfactory receptors",
     "polarity": "negative",
     "note": "the cleanest negative control there is"},
]

# `adh_aldh_plus.txt` is deliberately NOT in that list. It carries `#` comments
# and an unexpanded `ALDH*` glob, and the page already has its own button that
# expands that glob against the live index at runtime. Publishing it here as a
# literal roster would resolve to 12 genes instead of 38.


def geneset_block(src_dir: Path, docs: Path, df) -> dict:
    """Copy each set into docs/ and count how many tokens the index resolves.

    Resolution mirrors the browser and extract_genes.py: symbol first, then
    Ensembl id with the version suffix stripped. The count is published so a
    button can state its size before anyone clicks it, and so a set that
    silently stops resolving is loud on the next rebuild.

    Globs are deliberately NOT expanded. A published set is an exact roster,
    and the browser reports an unexpanded `ALDH*` as unrecognised rather than
    quietly turning it into 27 chips.
    """
    by_symbol = {s.upper() for s in df["genename"] if s}
    by_id = {g.split(".")[0].upper() for g in df["gene"] if g}

    out_dir = docs / GENESETS_DIR
    sets = []
    for spec in GENESETS:
        src = src_dir / spec["file"]
        if not src.is_file():
            print(f"  WARNING: gene set {spec['id']} not found at {src} - "
                  f"skipped", file=sys.stderr)
            continue
        # Strip `#` comments, whole-line and trailing, so a source file written
        # for extract_genes.py publishes as clean one-token-per-line.
        tokens = []
        for line in src.read_text(encoding="utf-8").splitlines():
            tok = line.split("#", 1)[0].strip()
            if tok:
                tokens.append(tok)
        resolved = sum(1 for t in tokens
                       if t.upper() in by_symbol
                       or t.split(".")[0].upper() in by_id)

        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / spec["file"]).write_text("\n".join(tokens) + "\n",
                                            encoding="utf-8", newline="\n")

        sets.append({**spec, "n_tokens": len(tokens), "n_genes": resolved})
        miss = len(tokens) - resolved
        print(f"  {spec['id']:18s} {len(tokens):5d} tokens -> {resolved:5d} "
              f"in GTEx v11" + (f", {miss} not found" if miss else ""),
              file=sys.stderr)

    return {"available": bool(sets), "base_url": GENESETS_DIR, "sets": sets}


def gene_column_digest(path: Path) -> str:
    """md5 of the gene column, used to prove the gene set matches across tables."""
    col = pd.read_csv(path, usecols=["gene"], dtype=str)["gene"].fillna("")
    return hashlib.md5("\n".join(col).encode("utf-8")).hexdigest()


def header_of(path: Path) -> list[str]:
    """First line of a CSV, split on commas. Reads one line, not the file."""
    with path.open("r", encoding="utf-8") as fh:
        return fh.readline().rstrip("\r\n").split(",")


def hist_availability(raw_tables, exc_tables) -> dict:
    """
    Which tables actually carry ``hist`` / ``hist_max``.

    They were added after the first full generation, so most of ``outputs/`` is
    still without them and a regeneration is ~852 MB of push. The browser's gene
    page draws the real 40-bin histogram when they are there and says so plainly
    when they are not — but it can only make that call *before* committing to an
    ~8 MB fetch if the manifest tells it. Reads one line per table.
    """
    def carries(path: Path) -> bool:
        cols = header_of(path)
        return "hist" in cols and "hist_max" in cols

    exc = [p.name[len(PREFIX):-len(EXCLUDED_SUFFIX)] for p in exc_tables if carries(p)]
    return {"excluded": exc}


def worm_block(worm_dir: Path) -> dict:
    """Describe the C. elegans table for the browser's Worm Data tab.

    A second organism rather than a 55th tissue: its identifiers are WormBase,
    `genes.tsv` does not contain them, and there is only ever one table -- so it
    is published as its own block and read by its own panel rather than being
    squeezed into ``tissues``. Deliberately kept out of ``outputs/``, where a
    glob for ``*_fourparam*.csv`` would sweep it into the tissue pipeline.

    ``available: False`` (rather than a missing key) when the file is absent, so
    the page can say "no worm table published" instead of failing to boot.
    """
    path = worm_dir / WORM_FILE
    if not path.is_file():
        print(f"  (no {WORM_FILE} in {worm_dir} - the Worm Data tab will be empty)",
              file=sys.stderr)
        return {"available": False}

    with path.open(encoding="utf-8", newline="") as fh:
        header = fh.readline().rstrip("\r\n").split(",")
        n_genes = sum(1 for line in fh if line.strip())
    print(f"Worm table     : {n_genes:,} genes, {len(header)} columns",
          file=sys.stderr)
    if "r_squared" not in header:
        print("  (no r_squared column yet - run compute_r2.py against the worm "
              "matrix, then append_r2_column.py)", file=sys.stderr)
    return {
        "available": True,
        "base_url": WORM_BASE,
        "file": WORM_FILE,
        "n_genes": n_genes,
        "columns": header,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outputs", type=Path, default=DEFAULT_OUTPUTS)
    ap.add_argument("--docs", type=Path, default=DEFAULT_DOCS)
    ap.add_argument("--gene-major", type=Path, default=DEFAULT_GENE_MAJOR,
                    help="Directory of gene-major shards (default: ../gene_major).")
    ap.add_argument("--hist-major", type=Path, default=DEFAULT_HIST_MAJOR,
                    help="Directory of histogram shards (default: ../hist_major).")
    ap.add_argument("--worm", type=Path, default=DEFAULT_WORM,
                    help="Directory holding the C. elegans table (default: ../worm).")
    ap.add_argument("--genelists", type=Path, default=DEFAULT_GENELISTS,
                    help="Directory of reusable gene sets (default: ../genelists). "
                         "Each set named in GENESETS is copied into docs/ and "
                         "indexed in the manifest.")
    ap.add_argument("--reference", type=str, default=None,
                   help="Tissue whose raw table supplies manifest.columns "
                        "(default: first alphabetically). During a partial "
                        "column migration the alphabetically-first table may "
                        "not carry the new columns yet, which would publish a "
                        "manifest that silently omits them.")
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
    if args.reference:
        match = [p for p in raw_tables
                 if p.name == f"{PREFIX}{args.reference}{RAW_SUFFIX}"]
        if not match:
            print(f"  ERROR: no raw table for --reference {args.reference!r}.",
                  file=sys.stderr)
            return 1
        reference = match[0]
    print(f"  Reference table: {reference.name}", file=sys.stderr)
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

    print("\nGene sets      : published to "
          f"{(args.docs / GENESETS_DIR).name}/", file=sys.stderr)
    genesets = geneset_block(args.genelists, args.docs, df)

    # ---- manifest -----------------------------------------------------------
    columns = list(pd.read_csv(reference, nrows=0).columns)
    # `r_squared` is a real last column of every excluded table (appended by
    # append_r2_column.py) but must NOT appear here: the browser joins R^2 from
    # r2/ by gene index, which is the only source that also covers the
    # gene-major route, and listing it here would render and export it twice.
    columns = [c for c in columns if c != "r_squared"]

    hist_tissues = hist_availability(raw_tables, exc_tables)
    n_hist = len(hist_tissues["excluded"])
    print(f"Histogram cols : {n_hist}/{len(exc_tables)} excluded tables "
          f"carry hist/hist_max", file=sys.stderr)
    if n_hist == 0:
        print("  (the gene page will draw fitted curves only, with no bars)",
              file=sys.stderr)

    # The sidecar is what makes a single gene's histogram a ~5 KB fetch instead
    # of an ~8 MB one, so the gene page prefers it and only falls back to a whole
    # tissue table when it is absent.
    hist_shards = sorted(args.hist_major.glob("shard_*.csv")) \
        if args.hist_major.is_dir() else []
    hist_tissues["gene_major"] = {
        "available": bool(hist_shards),
        "base_url": HIST_MAJOR_BASE,
        "file_pattern": "shard_{shard:04d}.csv",
        "n_shards": len(hist_shards),
        "kind_ids": {"excluded": EXCLUDED_KIND},
    }
    if hist_shards:
        hist_size = sum(p.stat().st_size for p in hist_shards)
        print(f"Histogram shards: {len(hist_shards):,}, {hist_size / 1e6:.1f} MB "
              f"({hist_size / len(hist_shards) / 1024:.1f} KB each)", file=sys.stderr)
    elif n_hist:
        print(f"  WARNING: {n_hist} table(s) carry hist but {args.hist_major} is "
              f"empty - run build_hist_major.py, or the gene page falls back to "
              f"fetching whole ~8 MB tissue tables.", file=sys.stderr)

    if not (args.docs / PHENOTYPES_FILE).is_file():
        print(f"  WARNING: {PHENOTYPES_FILE} not in {args.docs} - the gene page's "
              f"phenotype panel will be empty.", file=sys.stderr)

    manifest = {
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "data_base_url": DATA_BASE,
        "file_prefix": PREFIX,
        # Excluded only. The raw tables are still generated and still in
        # outputs/, but they are not published to the page: the zero-expression
        # spike drags the fit, so every analysis works off the excluded set and
        # offering the choice only offered a worse answer. Re-listing "raw" here
        # is all it would take to bring the toggle back -- index.html builds the
        # control from this list.
        "table_kinds": [
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
            "kind_ids": {"excluded": EXCLUDED_KIND},
        },
        # Per-table, because hist/hist_max landed after the first full run.
        "hist": hist_tissues,
        # Curated over-expression phenotypes, published beside the page itself.
        "phenotypes": {
            "available": (args.docs / PHENOTYPES_FILE).is_file(),
            "file": PHENOTYPES_FILE,
        },
        # Reusable gene sets, one button each in the browser. See geneset_block().
        "genesets": genesets,
        # C. elegans, read by its own tab. See worm_block().
        "worm": worm_block(args.worm),
    }
    manifest_path = args.docs / "manifest.json"

    # `qc` and `r2` are owned by build_qc_class.py and build_r2.py, which patch
    # this file in place. This script rewrites it wholesale, so without carrying
    # those blocks forward the published indexes vanish whenever the manifest is
    # rebuilt — and the only way back is to re-run the owning script, which is
    # not safe mid-run (a live compute_qc.py / compute_r2.py leaves half-written
    # tables behind). Preserve them instead.
    #
    # Any future side-car that patches this file belongs in this tuple. A block
    # that is not listed here is silently dropped on the next rebuild.
    if manifest_path.is_file():
        try:
            previous = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            previous = {}
        for key in PATCHED_BLOCKS:
            block = previous.get(key)
            if isinstance(block, dict):
                manifest[key] = block
                n = len(block.get("files", {}))
                print(f"Carried over   : {key} block, {n} tissue(s)",
                      file=sys.stderr)

    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {manifest_path.name}: {len(tissues)} tissues, "
          f"{len(columns)} columns", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
