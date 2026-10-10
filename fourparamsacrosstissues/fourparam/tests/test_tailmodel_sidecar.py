"""The browser encoding of newRTI / newLTI (``build_tailmodel.py``).

The page reads gene ``gi``'s field ``f`` as a fixed-offset slice. One record of
the wrong width shifts every later gene onto its neighbor's values, which
nothing downstream can detect, so the geometry is pinned here.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import build_tailmodel as B


def _decode(text, gi, f):
    w, k = B.WIDTH, len(B.FIELDS)
    field = text[(gi * k + f) * w:(gi * k + f + 1) * w]
    return float(field) if field.strip() else np.nan


def test_round_trip_in_gene_order_with_blanks(tmp_path):
    csv = tmp_path / "v11_log2_x_tailmodel.csv"
    pd.DataFrame({
        "gene": ["ENSG2.4", "ENSG1.9"],                   # versions differ, order differs
        "rti_missing": [0.00122, 0.5],
        "lti_missing": [np.nan, 0.123456],
        "usable": [True, True],
        "idklti": [0.0, 0.4],
    }).to_csv(csv, index=False)
    order = ["ENSG1", "ENSG3", "ENSG2"]                   # ENSG3 absent from the table

    text, counts = B.encode(csv, order)

    assert len(text) == B.WIDTH * len(B.FIELDS) * len(order)
    assert _decode(text, 0, 0) == pytest.approx(0.5)
    assert _decode(text, 0, 1) == pytest.approx(0.12346)
    assert np.isnan(_decode(text, 1, 0)) and np.isnan(_decode(text, 1, 1))
    assert _decode(text, 2, 0) == pytest.approx(0.00122)
    assert np.isnan(_decode(text, 2, 1))                  # censored: unknown, not 0
    assert _decode(text, 0, 2) == pytest.approx(0.4)
    assert _decode(text, 2, 2) == 0.0                     # zeroed, not blank
    assert counts == {"rti_missing": 2, "lti_missing": 1, "idklti": 2}


def test_unusable_fit_blanks_tail_fields_but_not_idklti(tmp_path):
    csv = tmp_path / "v11_log2_x_tailmodel.csv"
    pd.DataFrame({"gene": ["ENSG1"], "rti_missing": [1.0], "lti_missing": [np.nan],
                  "usable": [False], "idklti": [0.0]}).to_csv(csv, index=False)
    text, _ = B.encode(csv, ["ENSG1"])
    assert np.isnan(_decode(text, 0, 0))
    assert _decode(text, 0, 2) == 0.0


def test_refuses_a_value_that_is_not_a_fraction(tmp_path):
    csv = tmp_path / "v11_log2_x_tailmodel.csv"
    pd.DataFrame({"gene": ["ENSG1"], "rti_missing": [1.5],
                  "lti_missing": [0.1]}).to_csv(csv, index=False)
    with pytest.raises(SystemExit):
        B.encode(csv, ["ENSG1"])
