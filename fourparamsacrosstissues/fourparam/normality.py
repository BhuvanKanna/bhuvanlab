# -*- coding: utf-8 -*-
"""
normality.py

Decide whether one gene's expression distribution is Gaussian, and if it is not,
say **how** it fails. "Not normal" covers zero-inflation, bimodality, right skew
and right-truncation, and only the last is the hypothesis this project is about,
so a single p-value is not a useful answer.

The cascade, in order (each stage only sees what survives the one before):

    0  classify_support   fraction at the -1 floor. >= 20% and the gene is a
                          spike plus a smear; no unimodal continuous model
                          applies and a normality test just rediscovers the
                          spike.
    1  count_modes        KDE modes via bhuvanfitter.gene_peaks. > 1 and the
                          right model is a mixture, not a truncated Gaussian.
    2  run_tests          Shapiro-Wilk (primary), D'Agostino-Pearson (says how
                          it fails), Anderson-Darling (tail-sensitive, and
                          truncation is a tail phenomenon).
    3  skew direction     right-truncation removes the UPPER tail, so it makes
                          skew NEGATIVE. Judged against sampling noise at this
                          n, not against zero.
    4  d_aic vs null      truncated-Gaussian MLE against a plain Gaussian,
                          compared to a simulated null (see calibrate_null).
    5  classify           one verdict string.

Why stage 4 needs a simulated null
----------------------------------
The truncation point is fixed at the observed maximum, which is data-dependent
and therefore biased toward the truncated model. At n=104 the textbook
``d_aic > 2`` rule fires on **41.4%** of samples drawn from a true Gaussian.
``calibrate_null(n)`` measures that bias directly by running the same comparison
on simulated normals, and returns the 95th percentile as the threshold. It is
scale-free -- d_aic and skew are both invariant to location and scale -- so one
calibration per sample size serves every gene in a tissue.

Public API
----------
classify_support(values)                  -- floor fraction / zero-inflation
count_modes(values)                       -- KDE mode count
run_tests(values)                         -- SW / DP / AD / skew / kurt
fit_truncated(values, x_max)              -- truncated-Gaussian MLE + d_aic
calibrate_null(n)                         -- simulated thresholds at this n
r_squared(values, y0, A, x0, w)           -- normalised fit error
qc_gates(...)                             -- fit-validity flags
classify(...)                             -- the verdict string
bh_fdr(pvals, q)                          -- Benjamini-Hochberg rejections
robust_z(values)                          -- median/MAD z-scores
"""

from functools import lru_cache

import numpy as np
from scipy.stats import shapiro, normaltest, anderson, skew as _skew, kurtosis as _kurt

from bhuvanfitter import (BhuvanFitter, gene_peaks, _fourparam_gaussian,
                          _fit_mle_truncated, HIST_BINS)

# The log2(TPM+1)-1 transform sends TPM=0 to exactly -1, so this floor is the
# "not detected" sentinel rather than a measurement. -0.75 keeps the same margin
# the notebook's has_minus_one_peak uses.
FLOOR = -0.75
FLOOR_FRAC = 0.20

# Shapiro-Wilk is unreliable below ~20 points and scipy warns above 5000.
MIN_TEST_OBS = 20
MAX_SHAPIRO_OBS = 5000

# Null calibration. 2000 draws puts the Monte-Carlo error on a 95th percentile
# at roughly +/-0.2 d_aic units, which is far finer than the decision needs.
NULL_SIMS = 2000
NULL_SEED = 20260808
NULL_Q = 0.95


def _finite(values):
    arr = np.asarray(values, dtype=float)
    return arr[np.isfinite(arr)]


# -- stage 0: support ---------------------------------------------------------

def classify_support(values, floor=FLOOR, frac_threshold=FLOOR_FRAC) -> dict:
    """
    Fraction of a gene's values sitting on the not-detected floor.

    Returns ``{n_obs, frac_at_floor, is_zero_inflated}``. A zero-inflated gene is
    not a failed Gaussian, it is a different kind of object: the -1 spike is a
    detection artefact, not a measurement, so fitting any continuous unimodal
    density to it describes the assay rather than the biology.
    """
    arr = _finite(values)
    n = int(arr.size)
    frac = float((arr <= floor).mean()) if n else float("nan")
    return {
        "n_obs": n,
        "frac_at_floor": frac,
        "is_zero_inflated": bool(n and frac >= frac_threshold),
    }


# -- stage 1: modes -----------------------------------------------------------

def count_modes(values, **kwargs) -> int:
    """
    Number of KDE modes, via ``bhuvanfitter.gene_peaks`` -- bin-independent, so
    it cannot be talked into a second peak by an unlucky bin edge. Returns 0 for
    degenerate input, which ``classify`` treats as un-assessable rather than
    unimodal.
    """
    return len(gene_peaks(values, **kwargs))


# -- stage 2: normality tests -------------------------------------------------

def run_tests(values) -> dict:
    """
    Shapiro-Wilk, D'Agostino-Pearson, Anderson-Darling, skew, excess kurtosis.

    Kolmogorov-Smirnov is deliberately absent: against parameters estimated from
    the same sample its null distribution is wrong, and the Lilliefors correction
    that fixes it is still less powerful than Shapiro-Wilk here.

    Anderson-Darling is reported as a statistic, not a p-value, because scipy
    returns fixed critical values for it. It earns its place by weighting the
    tails, which is where truncation lives.

    All test fields are NaN when n < MIN_TEST_OBS or the values have no spread.
    """
    arr = _finite(values)
    out = {"shapiro_p": np.nan, "dagostino_p": np.nan, "anderson_stat": np.nan,
           "skew": np.nan, "kurt": np.nan}
    if arr.size < MIN_TEST_OBS or float(np.std(arr)) < 1e-12:
        return out

    out["skew"] = float(_skew(arr))
    out["kurt"] = float(_kurt(arr))
    try:
        out["shapiro_p"] = float(shapiro(arr[:MAX_SHAPIRO_OBS]).pvalue)
    except ValueError:
        pass
    try:
        out["dagostino_p"] = float(normaltest(arr).pvalue)
    except ValueError:
        pass
    try:
        out["anderson_stat"] = float(anderson(arr, "norm").statistic)
    except (ValueError, ZeroDivisionError):
        pass
    return out


# -- stage 4: truncated fit ---------------------------------------------------

def fit_truncated(values, x_max=None) -> dict:
    """
    Truncated-Gaussian MLE plus its score against a plain Gaussian.

    Thin wrapper over ``BhuvanFitter.fit("mle")`` so this module owns one entry
    point per stage. Returns NaNs rather than raising on degenerate input.
    """
    arr = _finite(values)
    if arr.size < 2 or float(np.std(arr)) < 1e-12:
        return {"mle_mu": np.nan, "mle_sigma": np.nan, "mle_sigma_dist": np.nan,
                "d_aic": np.nan, "mle_success": False}
    try:
        r = BhuvanFitter(arr, gene_name="_", x_max=x_max).fit("mle")
    except (RuntimeError, ValueError):
        return {"mle_mu": np.nan, "mle_sigma": np.nan, "mle_sigma_dist": np.nan,
                "d_aic": np.nan, "mle_success": False}
    return {"mle_mu": r["mle_mu"], "mle_sigma": r["mle_sigma"],
            "mle_sigma_dist": r["mle_sigma_dist"], "d_aic": r["d_aic"],
            "mle_success": r["fit_success"]}


@lru_cache(maxsize=64)
def calibrate_null(n: int, n_sims: int = NULL_SIMS, seed: int = NULL_SEED,
                   q: float = NULL_Q) -> dict:
    """
    Thresholds for ``d_aic`` and ``skew`` at sample size ``n``, measured by
    simulating from a true Gaussian.

    This is the load-bearing piece of the whole module. Because the ceiling is
    pinned to the observed maximum, the truncated model gets a free advantage
    that grows as n shrinks; the only honest way to know how big it is, is to
    run the identical comparison on data that is Gaussian by construction.

    N(0,1) is sufficient for any gene: ``d_aic`` is a likelihood ratio between
    two location-scale families and ``skew`` is standardised, so both are
    invariant to the gene's own mean and variance. Cached because the cost is a
    few thousand MLE fits and every gene in a tissue shares one n.

    Returns ``{n, n_sims, d_aic_threshold, skew_lo, skew_hi, d_aic_median}``,
    where ``skew_lo``/``skew_hi`` are the central 90% of skew under the null --
    a gene inside that band is not measurably skewed at this sample size.
    """
    rng = np.random.default_rng(seed + n)
    d_aics = np.empty(n_sims, dtype=float)
    skews = np.empty(n_sims, dtype=float)
    for i in range(n_sims):
        x = rng.normal(0.0, 1.0, n)
        skews[i] = _skew(x)
        mu, sigma, nll_t, ok = _fit_mle_truncated(x, float(x.max()))
        if not ok:
            d_aics[i] = np.nan
            continue
        s_n = float(np.std(x, ddof=0))
        nll_n = -float(np.sum(
            -0.5 * ((x - x.mean()) / s_n) ** 2 - np.log(s_n) - 0.5 * np.log(2 * np.pi)))
        d_aics[i] = 2.0 * (nll_n - nll_t)

    good = d_aics[np.isfinite(d_aics)]
    return {
        "n": int(n),
        "n_sims": int(n_sims),
        "d_aic_threshold": float(np.quantile(good, q)) if good.size else float("nan"),
        "d_aic_median": float(np.median(good)) if good.size else float("nan"),
        "skew_lo": float(np.quantile(skews, 0.05)),
        "skew_hi": float(np.quantile(skews, 0.95)),
    }


# -- fit quality --------------------------------------------------------------

def r_squared(values, y0, A, x0, w, bins: int = HIST_BINS) -> float:
    """
    Coefficient of determination of the 4-parameter fit against its own bin
    counts, on exactly the binning ``BhuvanFitter`` used.

    ``sumsquarevalue`` in the tables is an unnormalised residual sum of squares:
    it scales with n and with peak height, so p50=58 for kidney cortex says
    nothing about whether any individual fit is good. R^2 divides that by the
    total sum of squares and is comparable between genes and between tissues.

    Can be negative -- that means the fit is worse than a flat line at the mean
    count, which is a real and useful thing to be able to report.
    """
    arr = _finite(values)
    if arr.size == 0 or not all(np.isfinite([y0, A, x0, w])) or w == 0:
        return float("nan")
    counts, edges = np.histogram(arr, bins=bins)
    centres = 0.5 * (edges[:-1] + edges[1:])
    counts = counts.astype(float)
    tss = float(np.sum((counts - counts.mean()) ** 2))
    if tss <= 0:
        return float("nan")
    pred = _fourparam_gaussian(centres, y0, A, x0, w)
    if not np.all(np.isfinite(pred)):
        return float("nan")
    return float(1.0 - np.sum((counts - pred) ** 2) / tss)


# Absolute, not percentile-based. A percentile rule flags a fixed fraction no
# matter how good the fits are; "sigma is wider than the entire observed range"
# is a statement about validity that a percentile cannot make.
SIGMA_SPAN_LO = 0.05     # below: the fit collapsed onto one bin
SIGMA_SPAN_HI = 1.00     # above: sigma exceeds the data it was fit to (the APP mode)
MIN_QC_OBS = 30


def qc_gates(x0, w, lo, hi, n_obs, frac_at_floor,
             sigma_span_lo=SIGMA_SPAN_LO, sigma_span_hi=SIGMA_SPAN_HI,
             min_obs=MIN_QC_OBS, floor_frac=FLOOR_FRAC) -> dict:
    """
    Fit-validity flags, evaluated before any biological question is asked.

    Two distinct degeneracies hide behind one ``fit_success = True``:

      * ``sigma << span`` -- the Gaussian collapsed onto a single bin. Around a
        quarter of kidney-cortex fits have sigma under 0.3% of their data span.
        These produce |sigma_dist| in the thousands.
      * ``sigma >> span`` -- the fit is a shallow slice off the top of an
        enormous Gaussian (APP: w=22.7, sigma=16.0, against a 3.95-unit span).
        These produce sigma_dist near zero, which reads as maximal truncation
        and is an artefact.

    ``x0`` outside ``[min, max]`` means the reported peak was never observed;
    that is 23.5% of converged kidney-cortex fits.

    Returns the individual flags plus ``qc_pass``, true only if all hold.
    """
    span = (hi - lo) if all(np.isfinite([lo, hi])) else np.nan
    ratio = (abs(w) / np.sqrt(2.0)) / span if (np.isfinite(w) and np.isfinite(span)
                                               and span > 0) else np.nan

    g = {
        "sigma_span_ratio": float(ratio),
        "sigma_span_ok": bool(np.isfinite(ratio) and sigma_span_lo <= ratio <= sigma_span_hi),
        "x0_in_range": bool(np.isfinite(x0) and np.isfinite(lo) and np.isfinite(hi)
                            and lo <= x0 <= hi),
        "n_obs_ok": bool(np.isfinite(n_obs) and n_obs >= min_obs),
        "support_ok": bool(np.isfinite(frac_at_floor) and frac_at_floor < floor_frac),
    }
    g["qc_pass"] = bool(g["sigma_span_ok"] and g["x0_in_range"]
                        and g["n_obs_ok"] and g["support_ok"])
    return g


# -- stage 5: the verdict -----------------------------------------------------

def classify(support: dict, n_modes: int, tests: dict, trunc: dict,
             null: dict, alpha: float = 0.05) -> str:
    """
    Collapse the cascade into one label.

    Order matters and is not arbitrary -- each check is only meaningful once the
    ones above it have been ruled out. A bimodal gene is also "non-normal" and
    also "skewed", but reporting it as either would be misleading.

    ``zero_inflated``   >= 20% of donors at the detection floor
    ``undetermined``    too few points, or no spread, to say anything
    ``multimodal``      more than one KDE mode
    ``right_truncated`` skew below the null band AND d_aic past the simulated
                        threshold -- both, because either alone is common noise
    ``right_skewed``    skewed above the null band; a floor effect, and the
                        opposite of a ceiling
    ``non_normal``      rejects Shapiro-Wilk without a characterised direction
    ``normal``          survives every stage
    """
    if support.get("is_zero_inflated"):
        return "zero_inflated"
    if support.get("n_obs", 0) < MIN_TEST_OBS or not np.isfinite(tests.get("skew", np.nan)):
        return "undetermined"
    if n_modes > 1:
        return "multimodal"

    sk = tests["skew"]
    d_aic = trunc.get("d_aic", np.nan)
    thr = null.get("d_aic_threshold", np.nan)

    left_skewed = np.isfinite(sk) and sk < null.get("skew_lo", -np.inf)
    trunc_favoured = np.isfinite(d_aic) and np.isfinite(thr) and d_aic > thr
    if left_skewed and trunc_favoured:
        return "right_truncated"
    if np.isfinite(sk) and sk > null.get("skew_hi", np.inf):
        return "right_skewed"

    p = tests.get("shapiro_p", np.nan)
    if np.isfinite(p) and p < alpha:
        return "non_normal"
    return "normal"


# -- multiple testing and robust outliers -------------------------------------

def bh_fdr(pvals, q: float = 0.05):
    """
    Benjamini-Hochberg rejections at false-discovery rate ``q``.

    Returns a boolean array aligned with ``pvals``; NaN entries are never
    rejected and are excluded from the ranking. At 74,628 genes per tissue an
    uncorrected p < 0.05 flags ~3,700 by chance, so this is not optional.
    """
    p = np.asarray(pvals, dtype=float)
    out = np.zeros(p.shape, dtype=bool)
    idx = np.flatnonzero(np.isfinite(p))
    if idx.size == 0:
        return out
    order = idx[np.argsort(p[idx])]
    m = order.size
    passing = p[order] <= (np.arange(1, m + 1) / m) * q
    if not passing.any():
        return out
    out[order[: np.flatnonzero(passing).max() + 1]] = True
    return out


def robust_z(values):
    """
    Iglewicz-Hoaglin modified z-score, ``0.6745 * (x - median) / MAD``.

    Median and MAD rather than mean and SD because these parameter distributions
    are exactly the ones that break moment-based rules: ``ti_fourparam_sigma_dist``
    reaches 1e5 from degenerate fits, and a single such value moves the SD enough
    to hide every genuine outlier. |z| > 3.5 is the conventional cut.

    Returns all-NaN when the MAD is zero -- which happens when over half the
    values are identical, as for ``truncationindex`` where 96% are exactly 0.
    That is a signal to threshold within the nonzero part instead, not to fall
    back to a mean.
    """
    x = np.asarray(values, dtype=float)
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    if not np.isfinite(mad) or mad <= 0:
        return np.full(x.shape, np.nan)
    return 0.6745 * (x - med) / mad
