# -*- coding: utf-8 -*-
"""
tailmodel.py

A replacement candidate for RTI / LTI that models each tail the way it was
actually lost, and reports truncation as **the estimated fraction of the
population that is missing**.

The two edges of an expression distribution are not the same kind of edge:

- The **left** edge is usually the detection floor. ``log2(TPM+1)-1`` maps
  TPM = 0 to exactly -1, and a donor sitting there still exists -- we know how
  many there are, only not their value. That is **censoring**. The excluded
  tables throw those donors away and then read the gap as truncation, which is
  how LTI ended up tracking expression level.
- The **right** edge, if a dosage ceiling exists, removes donors outright. We do
  not know how many are gone. That is **truncation**.

So the likelihood for one gene, with latent Y ~ N(mu, sigma), ceiling
T = x_max, and floor L = -1, is

    n_cens * log Phi((L-mu)/sigma)              # censored donors, counted
  + sum_i  log phi((y_i-mu)/sigma) - log sigma  # detected donors
  - N      * log Phi((T-mu)/sigma)              # everyone survived the ceiling

and the right truncation index is ``rti_missing = 1 - Phi((T-mu)/sigma)``.

When **no** donor is censored the left edge is not the assay, so it is treated
as a second truncation point at x_min and ``lti_missing = Phi((x_min-mu)/sigma)``.
When donors **are** censored, left truncation cannot be told apart from the
detection limit, and ``lti_missing`` is NaN -- unknown, not zero.

Both cut points are pinned at the observed extremes, which is the MLE of a
truncation point but biased toward "truncated": under a true Gaussian of N
draws, the mass beyond the max is Beta(1, N), mean 1/(N+1). Judge the numbers
against a simulated null at the same N (see ``simulate_truncation.py``), never
against zero.
"""

import numpy as np
from scipy.optimize import minimize
from scipy.special import log_ndtr

FLOOR = -1.0
MIN_DETECTED = 10

# Same finite stand-in for an impossible likelihood as bhuvanfitter: inf
# collapses Nelder-Mead's simplex instead of steering it away.
_NLL_PENALTY = 1e12
_LOG_SQRT_2PI = 0.5 * np.log(2.0 * np.pi)


def _log_mass_between(a, b):
    """log(Phi(b) - Phi(a)) for a < b, accurate in either tail."""
    if a > 0.0:                          # both in the upper tail: mirror it
        a, b = -b, -a
    lb, la = log_ndtr(b), log_ndtr(a)
    d = la - lb
    if not d < 0.0:
        return -np.inf
    return lb + np.log1p(-np.exp(d))


def _nll(theta, y, n_cens, floor, lo, hi):
    """
    Negative log-likelihood. ``lo`` is the left truncation point, or None when
    the left edge is censoring at ``floor``; ``hi`` is the right truncation
    point, or None for an untruncated right tail.
    """
    mu, log_sigma = theta
    sigma = np.exp(log_sigma)
    if not np.isfinite(sigma) or sigma <= 0.0:
        return _NLL_PENALTY

    z = (y - mu) / sigma
    ll = float(np.sum(-0.5 * z * z)) - y.size * (log_sigma + _LOG_SQRT_2PI)
    n_total = y.size + n_cens
    if lo is None:
        if n_cens:
            ll += n_cens * log_ndtr((floor - mu) / sigma)
        log_z = log_ndtr((hi - mu) / sigma) if hi is not None else 0.0
    else:
        b = (hi - mu) / sigma if hi is not None else np.inf
        log_z = (_log_mass_between((lo - mu) / sigma, b) if np.isfinite(b)
                 else log_ndtr((mu - lo) / sigma))
    ll -= n_total * log_z
    return -ll if np.isfinite(ll) else _NLL_PENALTY


def _minimize(y, n_cens, floor, lo, hi, theta0):
    args = (y, n_cens, floor, lo, hi)
    best = minimize(_nll, theta0, args=args, method="Nelder-Mead",
                    options=dict(maxiter=4000, fatol=1e-10, xatol=1e-10))
    try:
        alt = minimize(_nll, best.x, args=args, method="BFGS",
                       options=dict(maxiter=500))
        if np.isfinite(alt.fun) and alt.fun < best.fun:
            best = alt
    except (ValueError, FloatingPointError):
        pass
    return best


def _start(y, n_cens, floor):
    """
    Starting point. With little censoring the detected moments are close; past
    that, the median and IQR of the full sample (censored donors placed at the
    floor) stay unbiased as long as the floor is below the median.
    """
    n_total = y.size + n_cens
    if n_cens / n_total < 0.25:
        return np.array([y.mean(), np.log(max(y.std(ddof=1), 1e-3))])
    full = np.concatenate([np.full(n_cens, floor), y])
    q25, q50, q75 = np.quantile(full, [0.25, 0.5, 0.75])
    if n_cens / n_total >= 0.5:          # median is censored too: lean on the top
        q75 = np.quantile(full, 0.9)
        sigma = max((q75 - floor) / 2.0, 0.05)
        return np.array([floor, np.log(sigma)])
    return np.array([q50, np.log(max((q75 - q25) / 1.349, 0.05))])


def _empty(n_total, n_cens, left_mode):
    nan = float("nan")
    return {"mu": nan, "sigma": nan, "x_min": nan, "x_max": nan,
            "rti_missing": nan, "lti_missing": nan,
            "rti_lr": nan, "lti_lr": nan,
            "n_total": int(n_total), "n_censored": int(n_cens),
            "left_mode": left_mode, "success": False}


def fit_tail_model(values, floor: float = FLOOR) -> dict:
    """
    Fit one gene. ``values`` are **all** of its finite expression values,
    including the ones at or below ``floor`` -- pass the raw row, not the
    excluded one, since the censored count is the information.

    Returns ``mu``, ``sigma`` (of the latent, untruncated distribution),
    ``rti_missing`` / ``lti_missing`` (fraction of the population estimated to
    be missing past each edge), ``rti_lr`` / ``lti_lr`` (2 x log-likelihood
    gain of truncating that side over leaving it untruncated; same caveat about
    a simulated null), ``n_total``, ``n_censored``, ``left_mode``
    (``"censored"`` or ``"truncated"``) and ``success``.
    """
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    cens = arr <= floor
    n_cens = int(cens.sum())
    y = arr[~cens]
    left_mode = "censored" if n_cens else "truncated"
    if y.size < MIN_DETECTED or float(np.ptp(y)) <= 0.0:
        return _empty(arr.size, n_cens, left_mode)

    x_min, x_max = float(y.min()), float(y.max())
    lo = None if n_cens else x_min
    theta0 = _start(y, n_cens, floor)

    full = _minimize(y, n_cens, floor, lo, x_max, theta0)
    mu, sigma = float(full.x[0]), float(np.exp(full.x[1]))
    ok = bool(np.isfinite(full.fun) and full.fun < _NLL_PENALTY
              and np.isfinite(mu) and np.isfinite(sigma) and sigma > 0.0)
    if not ok:
        return _empty(arr.size, n_cens, left_mode)

    # Each side's likelihood-ratio: refit with that side left untruncated.
    no_right = _minimize(y, n_cens, floor, lo, None, full.x)
    rti_lr = 2.0 * (no_right.fun - full.fun)
    if lo is None:
        lti_missing, lti_lr = float("nan"), float("nan")
    else:
        no_left = _minimize(y, 0, floor, None, x_max, full.x)
        lti_lr = 2.0 * (no_left.fun - full.fun)
        lti_missing = float(np.exp(log_ndtr((x_min - mu) / sigma)))

    return {"mu": mu, "sigma": sigma, "x_min": x_min, "x_max": x_max,
            "rti_missing": float(np.exp(log_ndtr((mu - x_max) / sigma))),
            "lti_missing": lti_missing,
            "rti_lr": float(rti_lr), "lti_lr": float(lti_lr),
            "n_total": int(arr.size), "n_censored": n_cens,
            "left_mode": left_mode, "success": True}
