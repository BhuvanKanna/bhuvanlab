# -*- coding: utf-8 -*-
"""
drti_core.py

Donor-level objects and right-truncation statistics for dRTI (../drti/PLAN.md).

Input everywhere is one gene's residual matrix `Z`: donors x tissues (NaN where
the donor has no sample or the tissue does not express the gene), as written by
drti_prep.py. Every statistic is oriented so that **higher = more truncated**.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize
from scipy.special import log_ndtr
from scipy.stats import genpareto

MIN_DONOR_TISSUES = 5
LAMBDA_CLIP = (0.05, 0.95)
GPD_Q = 90.0
XI_CLIP = (-1.5, 1.5)


# ---------------------------------------------------------------- variance
def anova_components(Z: np.ndarray) -> tuple[float, float]:
    """(sigma_L^2, sigma_w^2) by one-way random-effects ANOVA over donors with
    >= 2 values. sigma_L^2 is clipped at a tiny positive number."""
    k = np.sum(np.isfinite(Z), axis=1)
    Z = Z[k >= 2]
    k = k[k >= 2].astype(float)
    if len(k) < 20:
        return np.nan, np.nan
    means = np.nanmean(Z, axis=1)
    N, D = k.sum(), len(k)
    grand = np.nansum(Z) / N
    msb = np.sum(k * (means - grand) ** 2) / (D - 1)
    msw = np.nansum((Z - means[:, None]) ** 2) / (N - D)
    k0 = (N - np.sum(k ** 2) / N) / (D - 1)
    return max((msb - msw) / k0, 1e-6), msw


def icc_of(Z: np.ndarray) -> float:
    sL, sw = anova_components(Z)
    return sL / (sL + sw) if np.isfinite(sL) else np.nan


# ----------------------------------------------------------------- objects
def objects(Z: np.ndarray, min_k: int = MIN_DONOR_TISSUES) -> dict:
    """The three donor-level estimates of L_d, over donors with >= min_k values.
    Also returns k_d and the noise variance of the plain mean, s_d^2."""
    k = np.sum(np.isfinite(Z), axis=1)
    keep = k >= min_k
    Zk, kk = Z[keep], k[keep].astype(float)
    sL, sw = anova_components(Z)
    mean = np.nanmean(Zk, axis=1)
    blup = mean * sL / (sL + sw / kk)

    # Factor weights: each tissue's loading on the leave-one-out donor mean.
    tot = np.nansum(Z, axis=1)
    w = np.zeros(Z.shape[1])
    for t in range(Z.shape[1]):
        has = np.isfinite(Z[:, t]) & (k >= 2)
        if has.sum() < 20:
            continue
        loo = (tot[has] - Z[has, t]) / (k[has] - 1)
        lam = np.corrcoef(Z[has, t], loo)[0, 1]
        lam = np.clip(lam if np.isfinite(lam) else 0.0, *LAMBDA_CLIP)
        w[t] = lam / (1 - lam ** 2)
    W = np.where(np.isfinite(Zk), w[None, :], 0.0)
    factor = np.nansum(np.nan_to_num(Zk) * W, axis=1) / np.where(W.sum(1) > 0, W.sum(1), np.nan)
    return {"mean": mean, "blup": blup, "factor": factor,
            "k": kk, "s2": sw / kk, "sigma_L2": sL, "sigma_w2": sw, "keep": keep}


# -------------------------------------------------------------- statistics
def _q(v, ps):
    return np.percentile(v, ps)


def utc_p(v: np.ndarray, p: float) -> float:
    q50, q75, qp = _q(v, [50, 75, p])
    if q75 - q50 <= 0 or qp - q50 <= 0:
        return np.nan
    return -float(np.log((qp - q50) / (q75 - q50)))


def asym(v: np.ndarray) -> float:
    q025, q25, q50, q75, q975 = _q(v, [2.5, 25, 50, 75, 97.5])
    if min(q75 - q50, q975 - q50, q50 - q25, q50 - q025) <= 0:
        return np.nan
    return -float(np.log((q975 - q50) / (q75 - q50)) - np.log((q50 - q025) / (q50 - q25)))


def _exceedances(v: np.ndarray) -> np.ndarray:
    u = np.percentile(v, GPD_Q)
    y = np.sort(v[v > u] - u)
    return y


def xi_pwm(v: np.ndarray) -> float:
    """-xi by probability-weighted moments (Hosking & Wallis 1987)."""
    y = _exceedances(v)
    m = len(y)
    if m < 15:
        return np.nan
    p = (np.arange(1, m + 1) - 0.35) / m
    a0 = y.mean()
    a1 = np.mean((1 - p) * y)
    if a0 - 2 * a1 <= 0:
        return np.nan
    return float(np.clip(a0 / (a0 - 2 * a1) - 2, -XI_CLIP[1], -XI_CLIP[0]))


def xi_gpd(v: np.ndarray) -> float:
    """-xi by GPD maximum likelihood on exceedances over q90."""
    y = _exceedances(v)
    if len(y) < 15 or y.max() <= 0:
        return np.nan
    try:
        c, _, _ = genpareto.fit(y, floc=0)
    except Exception:
        return np.nan
    return float(np.clip(-c, -XI_CLIP[1], -XI_CLIP[0]))


def edge_lr(m: np.ndarray, s2: np.ndarray) -> tuple[float, float]:
    """Gaussian latent level, optionally truncated above at c, observed with
    known noise variance s2 per donor. Returns (sqrt(2*dlogL), missing fraction).

    H0: m ~ N(mu, sigma^2 + s2)
    H1: L ~ N(mu, sigma^2) truncated above at c;  m = L + e,  e ~ N(0, s2)
        f(m) = N(m; mu, v) * Phi((c - mu_post)/tau) / Phi((c - mu)/sigma)
        v = sigma^2 + s2, mu_post = mu + sigma^2/v (m - mu), tau^2 = sigma^2 s2 / v

    Fit by L-BFGS-B from two starts for c (the data max and mu + 3 sd).
    """
    ok = np.isfinite(m) & np.isfinite(s2)
    m, s2 = m[ok], s2[ok]
    if len(m) < 50:
        return np.nan, np.nan

    def nll0(th):
        mu, ls = th
        ls = np.clip(ls, -10.0, 10.0)
        v = np.exp(2 * ls) + s2
        return 0.5 * np.sum(np.log(v) + (m - mu) ** 2 / v)

    def nll1(th):
        mu, ls, c = th
        ls = np.clip(ls, -10.0, 10.0)
        sig2 = np.exp(2 * ls)
        v = sig2 + s2
        mu_post = mu + sig2 / v * (m - mu)
        tau = np.sqrt(sig2 * s2 / v)
        return (0.5 * np.sum(np.log(v) + (m - mu) ** 2 / v)
                - np.sum(log_ndtr((c - mu_post) / tau))
                + len(m) * log_ndtr((c - mu) / np.sqrt(sig2)))

    opt = {"method": "L-BFGS-B", "options": {"maxiter": 500, "ftol": 1e-10, "gtol": 1e-6}}
    r0 = minimize(nll0, [np.mean(m), np.log(max(np.std(m), 1e-3))], **opt)
    mu0, ls0 = r0.x
    sd = np.exp(ls0)
    best = None
    for c0 in (np.max(m), mu0 + 3 * sd):
        r = minimize(nll1, [mu0, ls0, c0], **opt)
        if np.isfinite(r.fun) and (best is None or r.fun < best.fun):
            best = r
    if best is None:
        return np.nan, np.nan
    d = max(r0.fun - best.fun, 0.0)
    mu, ls, c = best.x
    miss = float(np.exp(log_ndtr(-(c - mu) / np.exp(ls))))
    return float(np.sqrt(2 * d)), miss


def _edge_lr_nelder_mead(m: np.ndarray, s2: np.ndarray) -> tuple[float, float]:
    """Gaussian latent level, optionally truncated above at c, observed with
    known noise variance s2 per donor. Returns (sqrt(2*dlogL), missing fraction).

    H0: m ~ N(mu, sigma^2 + s2)
    H1: L ~ N(mu, sigma^2) truncated above at c;  m = L + e,  e ~ N(0, s2)
        f(m) = N(m; mu, v) * Phi((c - mu_post)/tau) / Phi((c - mu)/sigma)
        v = sigma^2 + s2, mu_post = mu + sigma^2/v (m - mu), tau^2 = sigma^2 s2 / v
    """
    ok = np.isfinite(m) & np.isfinite(s2)
    m, s2 = m[ok], s2[ok]
    if len(m) < 50:
        return np.nan, np.nan

    def nll0(th):
        mu, ls = th
        ls = np.clip(ls, -10.0, 10.0)
        v = np.exp(2 * ls) + s2
        return 0.5 * np.sum(np.log(v) + (m - mu) ** 2 / v)

    r0 = minimize(nll0, [np.mean(m), np.log(max(np.std(m), 1e-3))], method="Nelder-Mead",
                  options={"xatol": 1e-6, "fatol": 1e-8, "maxiter": 2000})

    def nll1(th):
        mu, ls, c = th
        ls = np.clip(ls, -10.0, 10.0)
        sig2 = np.exp(2 * ls)
        v = sig2 + s2
        mu_post = mu + sig2 / v * (m - mu)
        tau = np.sqrt(sig2 * s2 / v)
        return (0.5 * np.sum(np.log(v) + (m - mu) ** 2 / v)
                - np.sum(log_ndtr((c - mu_post) / tau))
                + len(m) * log_ndtr((c - mu) / np.sqrt(sig2)))

    mu0, ls0 = r0.x
    best = None
    sd = np.exp(ls0)
    for c0 in (np.percentile(m, 99), np.max(m), np.max(m) + sd, mu0 + 3 * sd):
        r = minimize(nll1, [mu0, ls0, c0], method="Nelder-Mead",
                     options={"xatol": 1e-6, "fatol": 1e-8, "maxiter": 4000})
        if best is None or r.fun < best.fun:
            best = r
    d = max(r0.fun - best.fun, 0.0)
    mu, ls, c = best.x
    miss = float(np.exp(log_ndtr(-(c - mu) / np.exp(ls))))
    return float(np.sqrt(2 * d)), miss


STATS = {
    "utc95": lambda v: utc_p(v, 95.0),
    "utc975": lambda v: utc_p(v, 97.5),
    "utc99": lambda v: utc_p(v, 99.0),
    "asym": asym,
    "xi_pwm": xi_pwm,
    "xi_gpd": xi_gpd,
}
OBJECTS = ("mean", "blup", "factor")


def all_stats(Z: np.ndarray, min_k: int = MIN_DONOR_TISSUES, with_lr: bool = True) -> dict:
    """{(object, stat): value} plus ('mean','edge_lr'), ('mean','missing_frac')
    and 'n' (donors used)."""
    o = objects(Z, min_k)
    out = {"n": int(o["keep"].sum())}
    for obj in OBJECTS:
        v = o[obj]
        v = v[np.isfinite(v)]
        for name, f in STATS.items():
            out[(obj, name)] = f(v) if len(v) >= 50 else np.nan
    if with_lr:
        lr, miss = edge_lr(o["mean"], o["s2"])
        out[("mean", "edge_lr")] = lr
        out[("mean", "missing_frac")] = miss
    return out
