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
