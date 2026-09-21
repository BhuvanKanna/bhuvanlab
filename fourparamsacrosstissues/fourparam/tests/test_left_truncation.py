"""The left truncation index (LTI), the mirror of the right one (RTI).

RTI asks how high the fitted curve still is at the observed ceiling ``x_max``.
LTI asks the same question at the observed floor ``x_min``. Both are read
against the same denominator (``maxheight``) and the same baseline (the curve's
minimum over the histogram interval), so the two are directly comparable.

The baseline is what couples them: it is the lower of the two edge heights, so
whichever edge sits *farther* from the peak defines the baseline and zeroes its
own side. Exactly one of LTI / RTI is non-zero for any bell-shaped fit. That is
pinned below because it is surprising and easy to "fix" by accident.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bhuvanfitter import BhuvanFitter


def make_fitter(values, name="TEST"):
    bf = BhuvanFitter(values, gene_name=name)
    bf.fit("fourparam")
    return bf


@pytest.fixture
def skewed():
    """A right-skewed sample: peak near the floor, long tail to the right.

    This is the shape real expression data takes, so it is the case the metrics
    are actually read on. The floor sits close to the peak and the ceiling far
    from it, which makes RTI the zeroed side and LTI the informative one.
    """
    rng = np.random.default_rng(0)
    return rng.lognormal(mean=0.0, sigma=0.6, size=4000)


def test_lti_sigma_dist_counts_sigmas_from_peak_down_to_the_floor(skewed):
    """(x0 - x_min) / sigma, the exact mirror of rti_sigma_dist's (x_max - x0)."""
    bf = make_fitter(skewed)
    sigma = bf.fourparam_w / np.sqrt(2.0)
    expected = (bf.fourparam_x0 - bf.min()) / sigma
    assert bf.lti_sigma_dist == pytest.approx(expected)


def test_leftheight_is_the_curve_at_the_floor_above_the_baseline(skewed):
    """Mirror of rightheight: f(x_min) - min(f), not f(x_min) alone."""
    bf = make_fitter(skewed)
    baseline = bf.fourparam_function(
        np.linspace(bf.hist_edges[0], bf.hist_edges[-1], bf._CURVE_GRID)).min()
    assert bf.leftheight == pytest.approx(
        float(bf.fourparam_function(bf.min()) - baseline))


def test_lti_is_leftheight_over_maxheight_and_is_bounded(skewed):
    bf = make_fitter(skewed)
    assert bf.lti == pytest.approx(bf.leftheight / bf.maxheight)
    assert 0.0 <= bf.lti <= 1.0


def test_exactly_one_of_lti_and_rti_is_zero(skewed):
    """The shared baseline zeroes whichever edge is farther from the peak."""
    bf = make_fitter(skewed)
    assert min(bf.lti, bf.rti) == pytest.approx(0.0, abs=1e-12)
    assert max(bf.lti, bf.rti) > 0.0


def test_mirroring_the_data_swaps_lti_and_rti(skewed):
    """Negating every value reflects the distribution, so the floor becomes the
    ceiling. LTI of the reflection must equal RTI of the original."""
    bf = make_fitter(skewed)
    mirrored = make_fitter(-np.asarray(skewed), name="TEST_MIRROR")
    assert mirrored.lti == pytest.approx(bf.rti, rel=1e-6, abs=1e-9)
    assert mirrored.lti_sigma_dist == pytest.approx(
        bf.rti_sigma_dist, rel=1e-6, abs=1e-9)


def test_rti_is_the_new_name_for_truncationindex(skewed):
    """The rename keeps the old attribute working so nothing downstream breaks."""
    bf = make_fitter(skewed)
    assert bf.rti == bf.truncationindex
    assert bf.rti_sigma_dist == bf.ti_fourparam_sigma_dist


def test_left_metrics_require_the_fit(skewed):
    bf = BhuvanFitter(skewed, gene_name="TEST")
    for attr in ("lti", "lti_sigma_dist", "leftheight"):
        with pytest.raises(RuntimeError):
            getattr(bf, attr)


def test_fourparam_fit_reports_the_left_metrics(skewed):
    """The dict _fit_fourparam returns is what the table rows are built from."""
    bf = BhuvanFitter(skewed, gene_name="TEST")
    result = bf.fit("fourparam")
    for key in ("lti", "lti_sigma_dist", "leftheight", "left",
                "rti", "rti_sigma_dist"):
        assert key in result, f"missing {key}"
    assert result["left"] == pytest.approx(bf.min())
    assert result["lti"] == pytest.approx(bf.lti)
