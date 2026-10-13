# -*- coding: utf-8 -*-
"""
build_gene_type.py

Publish each gene's GENCODE `gene_type` to the browser as one character per gene.

GTEx v11 is annotated with GENCODE v47 -- every versioned id in genes.tsv
matches a v47 `gene` record exactly -- so that release is the default. The
gene set is identical across tissues, so this is one file, not one per tissue:
a string of exactly `len(genes.tsv)` characters in **genes.tsv order** (~75 KB),
read by the browser as `text[gi]`, the same trade build_qc_class.py makes.

The browser shows two columns from it: the raw `gene_type` and a coarse
`gene_group` (protein_coding / lncRNA / small_ncRNA / pseudogene / IG_TR /
other), so "protein-coding only" or "lncRNA only" is a single filter.

    cd fourparam
    python build_gene_type.py                          # downloads GENCODE v47
    python build_gene_type.py --gtf gencode.v47.annotation.gtf.gz

Writes `docs/genetype/gencode_v47.txt` and patches the `gene_type` block of
`docs/manifest.json` in place. `build_gui_data.py` carries that block over, but
re-run this whenever genes.tsv is rebuilt, since the file is in its order.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import string
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
DOCS = HERE.parent.parent / "docs"
DOCS_GT = DOCS / "genetype"
GENES_TSV = DOCS / "genes.tsv"
MANIFEST = DOCS / "manifest.json"

GENCODE_RELEASE = 47
GENCODE_URL = ("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
               "release_{v}/gencode.v{v}.annotation.gtf.gz")

# `.` is "not in the annotation", never a real type.
MISSING = "."
# Characters handed out to types in sorted order. 62 is far more than the ~40
# types GENCODE uses, and none of them is MISSING.
ALPHABET = string.ascii_letters + string.digits

SMALL_NCRNA = {"miRNA", "misc_RNA", "snRNA", "snoRNA", "scaRNA", "rRNA",
               "Mt_tRNA", "Mt_rRNA", "ribozyme", "vault_RNA", "sRNA", "scRNA"}

_GENE_ID = re.compile(r'gene_id "([^"]+)"')
_GENE_TYPE = re.compile(r'gene_type "([^"]+)"')


def strip_version(gid: str) -> str:
    return gid.split(".", 1)[0]


def group_of(gene_type: str) -> str:
    """Coarse class. Pseudogenes are tested first because IG_V_pseudogene and
    friends are pseudogenes, not immune genes."""
    if gene_type == "protein_coding":
        return "protein_coding"
    if gene_type == "lncRNA":
        return "lncRNA"
    if gene_type in SMALL_NCRNA:
        return "small_ncRNA"
    if gene_type.endswith("pseudogene"):
        return "pseudogene"
    if gene_type.startswith(("IG_", "TR_")):
        return "IG_TR"
    return "other"                            # TEC, artifact


def load_gene_order() -> list[str]:
    """Unversioned ids in genes.tsv order -- the browser's `gi` indexes this."""
    ids = []
    with GENES_TSV.open(encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            ids.append(strip_version(line.split("\t", 1)[0]))
    return ids


def read_gtf_gene_types(lines) -> dict[str, str]:
    """{unversioned gene id: gene_type} from the `gene` records of a GTF."""
    out = {}
    for line in lines:
        if line.startswith("#"):
            continue
        parts = line.split("\t", 8)
        if len(parts) < 9 or parts[2] != "gene":
            continue
        gid = _GENE_ID.search(parts[8])
        gt = _GENE_TYPE.search(parts[8])
        if not gid or not gt:
            continue
        # PAR_Y copies carry a suffix on the same ENSG; the first record wins
        # and the types agree anyway.
        out.setdefault(strip_version(gid.group(1)), gt.group(1))
    return out


def encode(order: list[str], lut: dict[str, str]) -> tuple[str, dict[str, str]]:
    """The per-gene string and its {char: gene_type} legend."""
    types = sorted({lut[g] for g in order if g in lut})
    if len(types) > len(ALPHABET):
        raise SystemExit(f"{len(types)} gene types; the alphabet holds {len(ALPHABET)}")
    char_of = {t: ALPHABET[i] for i, t in enumerate(types)}
    text = "".join(char_of[lut[g]] if g in lut else MISSING for g in order)
    return text, {c: t for t, c in char_of.items()}


def patch_manifest(filename: str, classes: dict[str, str], release: int) -> None:
    """Merge a `gene_type` block into manifest.json without disturbing anything
    else -- build_gui_data.py rewrites the file and carries this block over."""
    man = json.loads(MANIFEST.read_text(encoding="utf-8"))
    man["gene_type"] = {
        "available": True,
        "base_url": "genetype",
        "file": filename,
        "source": f"GENCODE v{release} gene_type",
        "missing_char": MISSING,
        "classes": classes,
        "groups": {t: group_of(t) for t in sorted(classes.values())},
    }
    MANIFEST.write_text(json.dumps(man, indent=2) + "\n", encoding="utf-8")


def open_gtf(path: str | None, release: int):
    if path:
        return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") \
            else open(path, encoding="utf-8")
    url = GENCODE_URL.format(v=release)
    print(f"downloading {url}", file=sys.stderr)
    return gzip.open(urllib.request.urlopen(url), "rt", encoding="utf-8")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gtf", type=str, default=None,
                   help="Local GENCODE annotation GTF (.gtf or .gtf.gz); "
                        "default downloads the release below.")
    p.add_argument("--release", type=int, default=GENCODE_RELEASE)
    args = p.parse_args()

    if not GENES_TSV.exists():
        raise SystemExit(f"missing {GENES_TSV} -- run build_gui_data.py first")
    order = load_gene_order()
    print(f"gene order: {len(order):,} ids from {GENES_TSV.name}")

    with open_gtf(args.gtf, args.release) as f:
        lut = read_gtf_gene_types(f)
    print(f"annotation: {len(lut):,} gene records")

    text, classes = encode(order, lut)
    n_missing = text.count(MISSING)
    if n_missing > len(order) // 100:
        raise SystemExit(f"{n_missing:,} genes absent from GENCODE v{args.release} "
                         "-- wrong release for this genes.tsv?")

    DOCS_GT.mkdir(parents=True, exist_ok=True)
    out = DOCS_GT / f"gencode_v{args.release}.txt"
    out.write_text(text, encoding="ascii", newline="")
    patch_manifest(out.name, classes, args.release)

    counts: dict[str, int] = {}
    for c in text:
        g = "missing" if c == MISSING else group_of(classes[c])
        counts[g] = counts.get(g, 0) + 1
    print(f"{out.name}: {len(text):,} chars, {len(classes)} types")
    for g, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"  {g:<15} {n:>7,}")
    print(f"patched {MANIFEST.name}: gene_type block")


if __name__ == "__main__":
    main()
