"""The censored-left / truncated-right tail model (``tailmodel.py``).

The detection floor (TPM = 0 -> -1) *censors*: those donors exist and we know
how many there are. A lethality ceiling *truncates*: those donors are absent and
we do not know how many. The model treats each side the way it actually is, and
reports truncation as the estimated fraction of the population that is missing.
"""
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import norm, truncnorm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tailmodel import FLOOR, fit_tail_model


def latent(rng, n, mu=0.0, sigma=1.0, lo_z=-np.inf, hi_z=np.inf):
    """n draws from N(mu, sigma) truncated to [mu + lo_z*sigma, mu + hi_z*sigma]."""
    return truncnorm.rvs(lo_z, hi_z, loc=mu, scale=sigma, size=n, random_state=rng)


def censor(y, floor=FLOOR):
    return np.where(y <= floor, floor, y)


def test_censored_values_are_counted_not_dropped():
    rng = np.random.default_rng(1)
    x = censor(latent(rng, 2000, mu=-0.5))
    r = fit_tail_model(x)
    assert r["n_total"] == 2000
    assert r["n_censored"] == int(np.sum(x <= FLOOR))
    assert r["left_mode"] == "censored"


def test_recovers_mean_under_heavy_censoring():
    """Dropping the floor (the excluded tables) inflates the mean; censoring does not."""
    rng = np.random.default_rng(2)
    x = censor(latent(rng, 5000, mu=-0.8, sigma=1.0))
    r = fit_tail_model(x)
    assert r["success"]
    assert r["mu"] == pytest.approx(-0.8, abs=0.06)
    assert r["sigma"] == pytest.approx(1.0, abs=0.06)
    assert x[x > FLOOR].mean() > -0.8 + 0.5      # the bias the excluded tables carry


def test_right_missing_fraction_recovered():
    rng = np.random.default_rng(3)
    x = latent(rng, 5000, mu=2.0, hi_z=0.5)       # ceiling 0.5 sigma above the mean
    r = fit_tail_model(x)
    assert r["rti_missing"] == pytest.approx(1 - norm.cdf(0.5), abs=0.05)


def test_untruncated_right_tail_reads_near_zero():
    rng = np.random.default_rng(4)
    x = latent(rng, 5000, mu=2.0)
    r = fit_tail_model(x)
    assert r["rti_missing"] < 0.01


def test_lti_undefined_when_floor_is_censoring():
    rng = np.random.default_rng(5)
    x = censor(latent(rng, 1000, mu=0.0))
    r = fit_tail_model(x)
    assert np.isnan(r["lti_missing"])


def test_left_missing_fraction_recovered_without_censoring():
    rng = np.random.default_rng(6)
    x = latent(rng, 5000, mu=3.0, lo_z=-0.5)      # floor 0.5 sigma below, far from -1
    r = fit_tail_model(x)
    assert r["left_mode"] == "truncated"
    assert r["lti_missing"] == pytest.approx(norm.cdf(-0.5), abs=0.05)


def test_location_and_scale_invariance():
    """Away from the detection floor the missing fractions are shape, not level."""
    rng = np.random.default_rng(7)
    z = latent(rng, 800, hi_z=1.0)
    a = fit_tail_model(5.0 + z)
    b = fit_tail_model(9.0 + 2.5 * z)
    assert a["rti_missing"] == pytest.approx(b["rti_missing"], abs=1e-4)


def test_too_few_observations_fail_softly():
    r = fit_tail_model([0.1, 0.2, 0.3])
    assert not r["success"]
    assert np.isnan(r["rti_missing"])
