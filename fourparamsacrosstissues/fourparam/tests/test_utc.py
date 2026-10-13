"""UTC metric, expression adjustment and evaluation helpers (see ../utc/PLAN.md)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import evaluate_utc as ev
import utc


def test_gaussian_utc_matches_quantile_ratio():
    """For a Gaussian, UTC = ln(z.975 / z.75) = ln(1.960 / 0.674)."""
    rng = np.random.default_rng(0)
    m = utc.tail_metrics(rng.standard_normal((1, 200_000)))
    assert m["utc"][0] == pytest.approx(np.log(1.95996 / 0.67449), abs=0.01)
    assert m["ltc"][0] == pytest.approx(np.log(1.95996 / 0.67449), abs=0.01)


def test_utc_is_shift_and_scale_invariant():
    rng = np.random.default_rng(1)
    x = rng.standard_normal((3, 500))
    a, b = utc.tail_metrics(x), utc.tail_metrics(3.0 * x + 7.0)
    np.testing.assert_allclose(a["utc"], b["utc"])
    np.testing.assert_allclose(a["ltc"], b["ltc"])


def test_flat_upper_half_is_missing_not_infinite():
    x = np.array([[0, 1, 2, 3, 3, 3, 3, 3, 3, 3.0]])
    assert np.isnan(utc.tail_metrics(x)["utc"][0])


def test_lower_contamination_leaves_utc_alone():
    rng = np.random.default_rng(2)
    clean = rng.standard_normal((1, 200_000))
    dirty = clean - 2.0 * (rng.random((1, 200_000)) < 0.15)
    assert abs(utc.tail_metrics(dirty)["utc"][0] - utc.tail_metrics(clean)["utc"][0]) < 0.05


def test_residual_z_removes_an_expression_trend():
    rng = np.random.default_rng(3)
    q50 = rng.uniform(0, 10, 5000)
    v = 0.3 * q50 + rng.standard_normal(5000)
    z = utc.residual_z(q50, v)
    assert abs(np.corrcoef(q50, z)[0, 1]) < 0.05
    assert np.median(np.abs(z)) == pytest.approx(0.6745, abs=0.05)


def test_auc_direction_and_ties():
    assert ev.auc_lower(np.array([0.0, 1.0]), np.array([2.0, 3.0])) == 1.0
    assert ev.auc_lower(np.array([2.0, 3.0]), np.array([0.0, 1.0])) == 0.0
    assert ev.auc_lower(np.array([1.0]), np.array([1.0])) == 0.5


def test_match_respects_caliper_and_no_replacement():
    pos = pd.DataFrame({"gene": ["P1", "P2", "P3"], "expr": [1.0, 1.02, 5.0]})
    neg = pd.DataFrame({"gene": ["N1", "N2"], "expr": [1.05, 9.0]})
    pairs = ev.match(pos, neg)
    assert len(pairs) == 1 and pairs[0][1] == "N1"
