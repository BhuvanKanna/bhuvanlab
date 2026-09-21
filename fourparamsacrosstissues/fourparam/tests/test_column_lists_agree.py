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


def table_columns():
    """The table schema: COLUMNS with genename inserted at 1 (insert_genename)."""
    return ["gene", "genename"] + gf.COLUMNS[1:]


def test_shard_header_is_the_schema_without_the_histogram_columns():
    """SHARD_HEADER is the table schema, minus `hist`/`hist_max`, behind a
    tissue,table prefix.

    The mirror deliberately does not carry the histogram: it is 2.4 GB, the
    columns are empty for 104 of the 108 tables, and `hist_major/` publishes
    them separately for ~24 MB. Putting them in SHARD_HEADER (which happened
    when they were added to COLUMNS) makes every table fail the header check,
    because only uterus and vagina actually have those columns -- so the mirror
    silently could not be rebuilt at all.
    """
    shard_cols = bgm.SHARD_HEADER.strip().split(",")
    assert shard_cols[:2] == ["tissue", "table"]
    expected = [c for c in table_columns() if c not in bgm.SHARD_EXCLUDED]
    assert shard_cols[2:] == expected
    assert "hist" not in shard_cols and "hist_max" not in shard_cols


def test_table_header_matches_the_shard_body():
    assert bgm.TABLE_HEADER.strip().split(",") == \
        [c for c in table_columns() if c not in bgm.SHARD_EXCLUDED]


def test_every_real_table_starts_with_the_shard_schema():
    """Whatever a table carries past `fit_success` -- hist, r_squared, both or
    neither -- it must start with exactly the columns the mirror emits, or the
    truncation that drops the extras would cut in the wrong place."""
    outputs = REPO / "fourparamsacrosstissues" / "outputs"
    tables = sorted(outputs.glob("v11_log2_*_fourparam*.csv"))
    if not tables:
        return
    prefix = bgm.TABLE_HEADER.strip().split(",")
    for path in tables:
        with path.open(encoding="utf-8") as fh:
            cols = fh.readline().strip().split(",")
        assert cols[:len(prefix)] == prefix, f"{path.name} diverges from the schema"


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
