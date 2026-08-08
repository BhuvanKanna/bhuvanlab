"""The normality cascade is only trustworthy if its threshold is calibrated, so
the load-bearing test here is the null self-check: data that IS Gaussian must be
flagged at the nominal rate and no more. Everything else pins the QC gates and
the degenerate cases."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import normality as N
from bhuvanfitter import BhuvanFitter, _fit_mle_truncated


# -- truncated MLE ------------------------------------------------------------

def test_mle_recovers_untruncated_parameters():
    """The point of the MLE fit: estimate the distribution BEFORE the ceiling.
    The naive sample moments of truncated data are biased and must not match."""
    rng = np.random.default_rng(3)
    x = rng.normal(5.0, 2.0, 20000)
    x = x[x <= 6.0]

    mu, sigma, _, ok = _fit_mle_truncated(x, 6.0)
    assert ok
    assert mu == pytest.approx(5.0, abs=0.1)
    assert sigma == pytest.approx(2.0, abs=0.1)

    # The bias being corrected is large, not marginal.
    assert x.mean() < 4.5
    assert x.std(ddof=1) < 1.7


def test_mle_sigma_dist_matches_the_definition():
    rng = np.random.default_rng(4)
    x = rng.normal(0.0, 1.0, 8000)
    x = x[x <= 1.0]
    r = BhuvanFitter(x, gene_name="g", x_max=1.0).fit("mle")
    assert r["mle_sigma_dist"] == pytest.approx(
        (1.0 - r["mle_mu"]) / r["mle_sigma"], rel=1e-9)


def test_mle_survives_degenerate_input():
    """Constant and near-empty input must return a failed row, never raise."""
    for bad in (np.zeros(50), np.array([1.0]), np.array([])):
        out = N.fit_truncated(bad)
        assert out["mle_success"] is False
        assert not np.isfinite(out["d_aic"])


# -- the calibration, which is the whole argument -----------------------------

@pytest.mark.parametrize("n", [104, 300])
def test_null_calibration_flags_the_nominal_rate(n):
    """Feed the calibrated threshold data drawn from a true Gaussian, using a
    different seed than the calibration itself. It must flag ~5%, not 41%."""
    null = N.calibrate_null(n, n_sims=400)
    rng = np.random.default_rng(9000 + n)
    flags = sum(N.fit_truncated(rng.normal(0, 1, n))["d_aic"] > null["d_aic_threshold"]
                for _ in range(300))
    assert 0.01 <= flags / 300 <= 0.12, f"{flags/300:.1%} at n={n}, expected ~5%"


def test_naive_threshold_of_two_is_badly_miscalibrated():
    """Pins the reason calibrate_null exists at all. If this ever starts passing
    at ~5%, the data-dependent truncation point has stopped biasing the
    comparison and the calibration machinery could be retired."""
    null = N.calibrate_null(104, n_sims=400)
    assert null["d_aic_threshold"] > 3.5
    rng = np.random.default_rng(77)
    naive = sum(N.fit_truncated(rng.normal(0, 1, 104))["d_aic"] > 2.0
                for _ in range(300))
    assert naive / 300 > 0.20


def test_null_thresholds_are_scale_free():
    """One calibration serves every gene in a tissue only if d_aic does not care
    about the gene's own mean and variance."""
    null = N.calibrate_null(104, n_sims=300)
    rng = np.random.default_rng(5)
    base = rng.normal(0, 1, 104)
    a = N.fit_truncated(base)["d_aic"]
    b = N.fit_truncated(base * 37.0 + 1000.0)["d_aic"]
    assert a == pytest.approx(b, rel=1e-6)
    assert null["skew_lo"] < 0 < null["skew_hi"]


def test_skew_band_tightens_as_n_grows():
    small = N.calibrate_null(104, n_sims=300)
    large = N.calibrate_null(600, n_sims=300)
    assert (large["skew_hi"] - large["skew_lo"]) < (small["skew_hi"] - small["skew_lo"])


# -- positive control ---------------------------------------------------------

def test_strong_truncation_is_detected_and_weak_truncation_is_not():
    """Truncation close to the mean removes enough mass to see at n=104; a
    ceiling at 3 sigma removes ~0.1% and is genuinely undetectable. Reporting
    the second as truncated would be a false positive, so this pins both ends."""
    null = N.calibrate_null(104, n_sims=400)
    rng = np.random.default_rng(21)

    def sample(cut):
        x = rng.normal(0, 1, 200000)
        return x[x <= cut][:104]

    assert N.fit_truncated(sample(0.5), x_max=0.5)["d_aic"] > null["d_aic_threshold"]
    assert N.fit_truncated(sample(3.0), x_max=3.0)["d_aic"] < null["d_aic_threshold"]


def test_classify_labels_the_obvious_cases():
    null = N.calibrate_null(104, n_sims=400)
    rng = np.random.default_rng(31)

    def label(x, x_max=None):
        s = N.classify_support(x)
        return N.classify(s, N.count_modes(x), N.run_tests(x),
                          N.fit_truncated(x, x_max=x_max), null)

    assert label(rng.normal(0, 1, 104)) == "normal"

    zero_inflated = np.concatenate([np.full(60, -1.0), rng.normal(3, 1, 44)])
    assert label(zero_inflated) == "zero_inflated"

    bimodal = np.concatenate([rng.normal(0, 0.3, 52), rng.normal(6, 0.3, 52)])
    assert label(bimodal) == "multimodal"

    assert label(rng.exponential(1.0, 104)) == "right_skewed"


# -- QC gates -----------------------------------------------------------------

def test_sigma_span_gate_catches_both_degeneracies():
    """The two failure modes are opposite and both pass fit_success=True."""
    # APP in kidney cortex, excluded table: sigma 16.0 against a 3.95 span.
    wide = N.qc_gates(x0=6.689, w=22.666, lo=4.725, hi=8.677,
                      n_obs=104, frac_at_floor=0.0)
    assert wide["sigma_span_ratio"] > 1.0
    assert not wide["sigma_span_ok"] and not wide["qc_pass"]

    # DDX11L1: w driven to ~1e-6, the collapsed-spike mode.
    spike = N.qc_gates(x0=-0.988, w=1.0001e-06, lo=-1.0, hi=-0.563,
                       n_obs=104, frac_at_floor=0.0)
    assert spike["sigma_span_ratio"] < 0.05
    assert not spike["sigma_span_ok"]

    # SLC34A1: sigma 1.493 against a 8.55 span, a fit worth believing.
    good = N.qc_gates(x0=4.925, w=2.111, lo=-0.963, hi=7.590,
                      n_obs=104, frac_at_floor=0.0)
    assert good["sigma_span_ok"] and good["x0_in_range"] and good["qc_pass"]


def test_x0_outside_the_data_fails_the_gate():
    g = N.qc_gates(x0=99.0, w=1.0, lo=0.0, hi=10.0, n_obs=104, frac_at_floor=0.0)
    assert not g["x0_in_range"] and not g["qc_pass"]


def test_zero_inflation_and_thin_samples_fail_the_gate():
    assert not N.qc_gates(1.0, 1.0, 0.0, 5.0, 104, frac_at_floor=0.55)["qc_pass"]
    assert not N.qc_gates(1.0, 1.0, 0.0, 5.0, 12, frac_at_floor=0.0)["qc_pass"]


def test_r_squared_beats_raw_ssr_for_comparability():
    """A good fit and a bad fit on the same shape must be ordered correctly, and
    R^2 must be bounded above by 1 where raw SSR is unbounded."""
    rng = np.random.default_rng(41)
    data = rng.normal(0, 1, 4000)
    counts, edges = np.histogram(data, bins=40)
    centres = 0.5 * (edges[:-1] + edges[1:])
    peak = float(counts.max())

    good = N.r_squared(data, 0.0, peak, 0.0, np.sqrt(2.0))
    flat = N.r_squared(data, float(counts.mean()), 0.0, 0.0, 1.0)
    assert good > 0.8
    assert good <= 1.0
    assert flat < good
    assert np.isnan(N.r_squared(data, np.nan, 1.0, 0.0, 1.0))


# -- support and multiple testing ---------------------------------------------

def test_classify_support_counts_the_floor():
    x = np.concatenate([np.full(30, -1.0), np.linspace(0, 5, 70)])
    s = N.classify_support(x)
    assert s["frac_at_floor"] == pytest.approx(0.30)
    assert s["is_zero_inflated"]
    assert not N.classify_support(np.linspace(0, 5, 100))["is_zero_inflated"]


def test_bh_fdr_is_stricter_than_raw_alpha_and_ignores_nan():
    rng = np.random.default_rng(51)
    p = rng.uniform(size=10000)          # pure null
    rejected = N.bh_fdr(p, q=0.05)
    assert rejected.sum() < (p < 0.05).sum()
    assert rejected.sum() <= 5

    mixed = np.array([1e-12, 1e-10, np.nan, 0.9])
    out = N.bh_fdr(mixed, q=0.05)
    assert out[0] and out[1] and not out[2] and not out[3]


def test_bh_fdr_recovers_strong_signal():
    p = np.concatenate([np.full(50, 1e-8), np.random.default_rng(6).uniform(size=950)])
    assert N.bh_fdr(p, q=0.05).sum() >= 50


def test_robust_z_is_not_moved_by_an_extreme_outlier():
    """The exact failure that breaks mean/SD rules on ti_fourparam_sigma_dist,
    where degenerate fits reach 1e5."""
    x = np.concatenate([np.random.default_rng(7).normal(0, 1, 500), [1e5]])
    z = N.robust_z(x)
    assert abs(z[-1]) > 3.5
    assert np.nanmax(np.abs(z[:-1])) < 20

    # Over half the values identical (as for truncationindex) -> MAD is 0.
    assert np.all(np.isnan(N.robust_z(np.concatenate([np.zeros(96), np.arange(4)]))))
