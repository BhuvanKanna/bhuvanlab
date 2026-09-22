"""`verify()` must compare against the same thing `_split_table()` emits.

The mirror carries the schema only: whatever a table holds past `fit_success`
(`hist`/`hist_max`, `r_squared`, both or neither) is dropped on the way in. So
the byte-comparison has to drop it too. Comparing a 24-field shard row against
a 27-field source row reports every row as DIFFERS -- a real failure signal,
but on the check rather than on the data.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import build_gene_major as bgm

SCHEMA = bgm.TABLE_HEADER.strip().split(",")


def _fields(gene):
    return [gene] + [f"v{i}" for i in range(1, len(SCHEMA))]


def _write_pair(tmp_path, extra_cols, extra_vals):
    """One source table carrying extra trailing columns, and the shard for it."""
    outputs = tmp_path / "outputs"
    outputs.mkdir()
    dest = tmp_path / "gene_major"
    dest.mkdir()

    gene = "ENSG00000000001.1"
    fields = _fields(gene)
    (outputs / "v11_log2_test_fourparam.csv").write_text(
        ",".join(SCHEMA + extra_cols) + "\n"
        + ",".join(fields + extra_vals) + "\n",
        encoding="utf-8", newline="")

    (dest / "shard_0000.csv").write_text(
        bgm.SHARD_HEADER + "test,raw," + ",".join(fields) + "\n",
        encoding="utf-8", newline="")
    return outputs, dest


def test_verify_passes_when_the_table_carries_r_squared(tmp_path):
    outputs, dest = _write_pair(tmp_path, ["r_squared"], ["0.87"])
    assert bgm.verify(outputs, dest, [0], 1) == 0


def test_verify_passes_when_the_table_carries_hist_and_r_squared(tmp_path):
    outputs, dest = _write_pair(tmp_path, ["hist", "hist_max", "r_squared"],
                                ["A" * 40, "210", "0.87"])
    assert bgm.verify(outputs, dest, [0], 1) == 0


def test_verify_passes_on_a_table_with_no_extra_columns(tmp_path):
    outputs, dest = _write_pair(tmp_path, [], [])
    assert bgm.verify(outputs, dest, [0], 1) == 0


def test_verify_still_catches_a_genuinely_wrong_row(tmp_path):
    """The check has to be able to fail, or dropping columns just hid the bug."""
    outputs, dest = _write_pair(tmp_path, ["r_squared"], ["0.87"])
    shard = dest / "shard_0000.csv"
    text = shard.read_text(encoding="utf-8")
    shard.write_text(text.replace(",v5,", ",CORRUPTED,"), encoding="utf-8",
                     newline="")
    assert bgm.verify(outputs, dest, [0], 1) > 0
