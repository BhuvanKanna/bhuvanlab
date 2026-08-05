"""The histogram columns must survive the paths that produce a row, including
the two that never construct a BhuvanFitter."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import generate_fourparam as gf
from bhuvanfitter import HIST_BINS


def test_hist_columns_are_last_and_named():
    assert gf.COLUMNS[-2:] == ["hist", "hist_max"]


def test_failed_row_carries_the_histogram_it_was_given():
    row = gf._failed_row("G", 4, "A" * HIST_BINS, 7)
    assert row["fit_success"] is False
    assert row["hist"] == "A" * HIST_BINS
    assert row["hist_max"] == 7


def test_failed_row_defaults_to_empty_histogram():
    row = gf._failed_row("G", 0)
    assert row["hist"] == ""
    assert row["hist_max"] == 0


def test_a_histogram_never_ships_without_its_bin_edges():
    """THE contract the browser depends on: bin edges are not stored, they are
    linspace(min, max, 41). A row carrying `hist` but NaN min/max cannot be
    drawn at all, which is exactly what failed-fit rows used to do."""
    gf._init_worker(None, 2000)
    rng = np.random.default_rng(9)
    rows = [gf._fit_one(("ok", rng.normal(size=250))),        # converges
            gf._fit_one(("small", np.array([1.0, 2.0, 3.0]))),  # below MIN_OBS
            gf._fit_one(("flat", np.full(40, -1.0))),           # all identical
            gf._fit_one(("none", np.array([np.nan])))]          # no data
    for row in rows:
        if row["hist"]:
            assert np.isfinite(row["min"]), row["gene"]
            assert np.isfinite(row["max"]), row["gene"]


def test_failed_row_carries_min_and_max_when_there_is_data():
    gf._init_worker(None, 2000)
    row = gf._fit_one(("G", np.array([1.5, 2.0, 9.25])))
    assert row["fit_success"] is False
    assert row["min"] == 1.5
    assert row["max"] == 9.25


def test_failed_row_with_no_data_has_nan_min_and_max():
    gf._init_worker(None, 2000)
    row = gf._fit_one(("G", np.array([np.nan])))
    assert row["hist"] == ""
    assert np.isnan(row["min"]) and np.isnan(row["max"])


def test_below_min_obs_still_gets_a_histogram():
    """A gene that never reaches the fit still shows its real distribution."""
    gf._init_worker(None, 2000)
    row = gf._fit_one(("G", np.array([1.0, 2.0, 3.0])))
    assert row["fit_success"] is False
    assert row["n_obs"] == 3
    assert len(row["hist"]) == HIST_BINS
    assert row["hist_max"] == 1


def test_successful_fit_gets_a_histogram():
    gf._init_worker(None, 2000)
    rng = np.random.default_rng(4)
    row = gf._fit_one(("G", rng.normal(size=300)))
    assert row["fit_success"] is True
    assert len(row["hist"]) == HIST_BINS
    assert row["hist_max"] > 0


def test_zero_observations_gives_empty_histogram():
    gf._init_worker(None, 2000)
    row = gf._fit_one(("G", np.array([np.nan, np.nan])))
    assert row["n_obs"] == 0
    assert row["hist"] == ""
    assert row["hist_max"] == 0


def test_threshold_is_applied_before_binning():
    """The excluded table's histogram must describe the excluded data."""
    values = np.concatenate([np.full(50, -1.0), np.linspace(0, 3, 50)])
    gf._init_worker(None, 2000)
    raw = gf._fit_one(("G", values.copy()))
    gf._init_worker(-1.0, 2000)
    excl = gf._fit_one(("G", values.copy()))
    assert raw["hist"] != excl["hist"]
    assert raw["hist_max"] == 50
    assert excl["n_obs"] == 50


def test_hist_max_column_stays_integer_in_the_written_csv(tmp_path):
    """A NaN in the column would make pandas write 12.0 for every row."""
    gf._init_worker(None, 2000)
    rng = np.random.default_rng(5)
    rows = [gf._fit_one(("A", rng.normal(size=200))),
            gf._fit_one(("B", np.array([np.nan])))]
    table = pd.DataFrame.from_records(rows, columns=gf.COLUMNS)
    out = tmp_path / "t.csv"
    table.to_csv(out, index=False, lineterminator="\n")
    text = out.read_text()
    for line in text.strip().split("\n")[1:]:
        assert "." not in line.split(",")[-1]


def test_no_commas_or_quotes_in_the_hist_column():
    gf._init_worker(None, 2000)
    rng = np.random.default_rng(6)
    row = gf._fit_one(("G", rng.normal(size=200)))
    assert "," not in row["hist"]
    assert '"' not in row["hist"]


def test_written_csv_has_a_uniform_field_count(tmp_path):
    """The browser splits on "," with no quoted-field handling, so a single
    stray comma or quote would shift every column index after it. Ragged rows
    are the failure mode that would cause, so assert against them directly."""
    gf._init_worker(None, 2000)
    rng = np.random.default_rng(7)
    rows = [gf._fit_one(("A", rng.normal(size=250))),          # clean fit
            gf._fit_one(("B", np.array([1.0, 2.0, 3.0]))),     # below MIN_OBS
            gf._fit_one(("C", np.full(40, -1.0))),             # all identical
            gf._fit_one(("D", np.array([np.nan])))]            # no data
    table = pd.DataFrame.from_records(rows, columns=gf.COLUMNS)
    out = tmp_path / "t.csv"
    table.to_csv(out, index=False, lineterminator="\n")

    lines = out.read_text().strip().split("\n")
    assert {len(line.split(",")) for line in lines} == {len(gf.COLUMNS)}
    assert '"' not in out.read_text()


def test_hist_column_is_reached_by_key_not_attribute():
    """DataFrame.hist is pandas' plotting method, so df.hist silently returns a
    bound method instead of the column. Pin that df["hist"] is the real one."""
    gf._init_worker(None, 2000)
    rng = np.random.default_rng(8)
    table = pd.DataFrame.from_records(
        [gf._fit_one(("A", rng.normal(size=200)))], columns=gf.COLUMNS)
    assert callable(table.hist)
    assert len(table["hist"].iloc[0]) == HIST_BINS
