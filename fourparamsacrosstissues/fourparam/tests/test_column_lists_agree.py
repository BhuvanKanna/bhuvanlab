"""Four places hardcode the table's column list. They must agree.

generate_fourparam.COLUMNS       -- what gets written to outputs/
build_gene_major.SHARD_HEADER    -- what the gene-major mirror expects and emits
extract_genes.STAT_COLUMNS       -- what the CLI extract emits
docs/manifest.json "columns"     -- what the browser's CSV export follows

Drift between them is silent and expensive: build_gene_major rejects every
table with "unexpected header", or the browser export loses columns, or the
byte-identity guarantee between the export and extract_genes.py dies.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import build_gene_major as bgm
import extract_genes as eg
import generate_fourparam as gf

REPO = Path(__file__).resolve().parents[3]


def test_shard_header_matches_generated_columns():
    """SHARD_HEADER is the generated columns with gene/genename in place and a
    tissue,table prefix."""
    shard_cols = bgm.SHARD_HEADER.strip().split(",")
    assert shard_cols[:2] == ["tissue", "table"]
    assert shard_cols[2:4] == ["gene", "genename"]
    # generate_fourparam.COLUMNS has no genename (insert_genename adds it at 1)
    expected = ["gene", "genename"] + gf.COLUMNS[1:]
    assert shard_cols[2:] == expected


def test_table_header_matches_generated_columns():
    assert bgm.TABLE_HEADER.strip().split(",") == ["gene", "genename"] + gf.COLUMNS[1:]


def test_extract_genes_stat_columns_match():
    assert eg.STAT_COLUMNS == gf.COLUMNS[1:]


def test_manifest_columns_match():
    manifest = REPO / "docs" / "manifest.json"
    if not manifest.exists():
        return  # manifest is generated; skip when absent
    cols = json.loads(manifest.read_text(encoding="utf-8"))["columns"]
    assert cols == ["gene", "genename"] + gf.COLUMNS[1:], (
        "docs/manifest.json is stale -- re-run "
        "build_gui_data.py --reference <a regenerated tissue>")
