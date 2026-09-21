"""The gene-major mirror must tolerate a table carrying extra trailing columns.

The 54 excluded tables end with `r_squared`, which is deliberately not part of
SHARD_HEADER (it is joined from r2/ by gene index so the gene-major route gets
it too). Before this, `_split_table` demanded an exact header match and rejected
every excluded table -- which is why the mirror was last built from the
pre-r_squared tables and silently drifted a schema behind.

Extra columns are dropped rather than carried, so the emitted shard is exactly
SHARD_HEADER and byte-identity with extract_genes.py still holds. No field in
these tables may contain a comma (the browser parses them with a plain
split(",")), so truncating at the Nth comma is safe and stays string-only -- no
float is ever re-serialised.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import build_gene_major as bgm


def _row(gene, n_fields, tail=None):
    fields = [gene] + [f"v{i}" for i in range(1, n_fields)]
    return ",".join(fields + ([tail] if tail else []))


def split(tmp_path, header, rows, gene_index=("ENSG1",)):
    src = tmp_path / "v11_log2_test_fourparam.csv"
    src.write_text(header + "\n" + "\n".join(rows) + "\n",
                   encoding="utf-8", newline="")
    out = tmp_path / "parts"
    out.mkdir()
    ok, msg, n = bgm._split_table(
        ("test", "excluded", src, out, {g: 0 for g in gene_index}, 1))
    return ok, msg, n, out


def test_exact_header_still_works(tmp_path):
    header = bgm.TABLE_HEADER.strip()
    n = len(header.split(","))
    ok, msg, rows, out = split(tmp_path, header, [_row("ENSG1", n)])
    assert ok, msg
    assert rows == 1


def test_trailing_r_squared_is_accepted_and_dropped(tmp_path):
    header = bgm.TABLE_HEADER.strip() + ",r_squared"
    n = len(bgm.TABLE_HEADER.strip().split(","))
    ok, msg, rows, out = split(tmp_path, header, [_row("ENSG1", n, "0.87")])
    assert ok, msg
    assert rows == 1
    written = next(out.glob("*.part")).read_text(encoding="utf-8")
    assert written.endswith("\n")
    fields = written.rstrip("\n").split(",")
    # tissue,table prefix + exactly the shard's own columns, and no r_squared.
    assert len(fields) == 2 + n
    assert "0.87" not in fields


def test_a_genuinely_wrong_header_is_still_rejected(tmp_path):
    header = "gene,genename,something_else"
    ok, msg, rows, out = split(tmp_path, header, ["ENSG1,x,y"])
    assert not ok
    assert "unexpected header" in msg
