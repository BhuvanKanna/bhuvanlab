# -*- coding: utf-8 -*-
"""
bhuvanfitter.py

Core fitting library for the BhuvanFitter project: the module-level
4-parameter Gaussian and the ``BhuvanFitter`` class. This is the single source
of truth for the fitting logic — both ``newbhuvanfitter.ipynb`` and
``generate_fourparam_stats.py`` import from here.

Public API
----------
_fourparam_gaussian(x, y0, A, x0, w)   -- the model curve (module level so
                                          scipy.optimize.curve_fit can use it)
BhuvanFitter(data, gene_name, x_max)   -- fit + metrics + plotting for one gene
gene_peaks(values, ...)                -- KDE peak detection for one gene
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde, skew, kurtosis
from scipy.signal import find_peaks


def _fourparam_gaussian(x, y0, A, x0, w):
    """
    4-parameter Gaussian model used by curve_fit.

        y = y0 + A * exp(-((x - x0) / w)^2)

    Defined at module level (not as a method) because scipy.optimize.curve_fit
    requires a plain picklable callable.

    Parameters
    ----------
    x  : array-like  Input x values.
    y0 : float       Baseline offset.
    A  : float       Amplitude (peak height above baseline).
    x0 : float       Centre of the peak.
    w  : float       Width parameter (w = sigma * sqrt(2)).
    """
    return y0 + A * np.exp(-((x - x0) / w) ** 2)


def gene_peaks(values, min_prominence_frac=0.05, bw_method="silverman",
               grid_size=1000, pad_frac=0.05, round_to=6):
    """
    Detect the peaks (modes) of one gene's expression distribution via a
    Gaussian KDE — bin-independent, unlike a histogram.

    A KDE is evaluated on a fine grid spanning the data range (padded), then
    ``scipy.signal.find_peaks`` keeps local maxima whose prominence is at least
    ``min_prominence_frac`` of the maximum density. This is the same detection
    approach as the project's former ``find_density_peaks``.

    Parameters
    ----------
    values : array-like
        Expression values across strains for one gene.
    min_prominence_frac : float
        Minimum peak prominence as a fraction of the maximum density. Scale-free,
        so the same threshold works across genes regardless of magnitude.
    bw_method : str or float
        Bandwidth passed to ``scipy.stats.gaussian_kde`` ('silverman', 'scott',
        a scalar, or a callable). Larger = smoother = fewer peaks.
    grid_size : int
        Number of grid points the density is evaluated on.
    pad_frac : float
        Fraction of the data range to pad the grid on each side.
    round_to : int
        Decimal places to round each peak's expression-value (the dict key) to.

    Returns
    -------
    dict
        ``{peak_expression_value: {"height": <kde density at peak>,
                                   "prominence": <peak prominence>}}``,
        sorted by ascending expression value. The number of peaks is ``len()``
        of this dict. Returns an empty dict for degenerate input (fewer than 5
        finite points, no spread, a singular KDE, or no interior mode).
    """
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size < 5:
        return {}

    lo, hi = float(arr.min()), float(arr.max())
    if hi - lo < 1e-12 or float(np.std(arr)) < 1e-8:
        return {}

    pad = (hi - lo) * pad_frac
    grid = np.linspace(lo - pad, hi + pad, grid_size)

    try:
        density = gaussian_kde(arr, bw_method=bw_method)(grid)
    except np.linalg.LinAlgError:
        return {}

    dmax = float(density.max())
    if dmax <= 0:
        return {}

    peak_idx, props = find_peaks(density, prominence=min_prominence_frac * dmax)

    peaks = {}
    for i, prom in zip(peak_idx, props["prominences"]):
        peak_value = round(float(grid[i]), round_to)
        peaks[peak_value] = {
            "height": float(density[i]),
            "prominence": float(prom),
        }
    return peaks


class BhuvanFitter:
    """
    Fit a 4-parameter Gaussian to one gene's expression histogram and report
    the fit parameters together with two truncation-index metrics.

    The histogram is always computed with 40 bins. Two fit models are
    registered:

    - ``fit("fourparam")`` returns a dict with everything
      ``compute_fourparam_table`` needs for one gene (minus
      ``has_minus_one_peak``); ``hist(lines=["fourparam"])`` overlays the
      fitted curve.
    - ``fit("kde")`` runs a bin-independent Gaussian KDE and detects its modes
      (reusing the module-level ``gene_peaks``); ``hist(lines=["kde"])``
      overlays the density curve (scaled to bin counts) with markers at each
      detected peak. The two overlays compose: ``hist(lines=["fourparam", "kde"])``.

    Parameters
    ----------
    data      : array-like   Expression values across strains for one gene.
    gene_name : str          Name of the gene (used in plot titles and repr).
    x_max     : float, optional
        Truncation point used by the truncation-index metrics. Defaults to the
        observed maximum of the data.
    """

    BINS = 40
    _FIT_REGISTRY = {"fourparam": "_fit_fourparam", "kde": "_fit_kde"}

    def __init__(self, data, gene_name: str, x_max=None):
        self.gene_name = gene_name

        arr = np.asarray(data, dtype=float)
        arr = arr[np.isfinite(arr)]
        if arr.size == 0:
            raise ValueError(f"No finite data values for gene '{gene_name}'.")
        self._data = arr

        self._x_max = float(x_max) if x_max is not None else float(arr.max())

        counts, edges = np.histogram(arr, bins=self.BINS)
        self.hist_counts = counts.astype(float)
        self.hist_edges = edges
        self.hist_centers = 0.5 * (edges[:-1] + edges[1:])

        self.active_fits = {name: False for name in self._FIT_REGISTRY}

        # fourparam result storage
        self.fourparam_y0 = None
        self.fourparam_A = None
        self.fourparam_x0 = None
        self.fourparam_w = None
        self.fourparam_sumsquare = None

        # kde result storage
        self.kde_object = None
        self.kde_grid = None
        self.kde_density = None
        self.kde_peaks = None
        self.kde_bw_method = None

    # -- Simple statistics -----------------------------------------------------

    def min(self):
        return float(self._data.min())

    def max(self):
        return float(self._data.max())

    def mean(self):
        return float(self._data.mean())

    def std(self):
        """Sample standard deviation (ddof=1) of the gene's finite values."""
        return float(np.std(self._data, ddof=1)) if self._data.size > 1 else float("nan")

    def skew(self):
        """Skewness of the gene's finite values (scipy default: bias=True).
        Matches the notebook's build_shape_features `skew(...)`."""
        return float(skew(self._data)) if self._data.size > 1 else float("nan")

    def kurt(self):
        """Fisher excess kurtosis of the gene's finite values (scipy default:
        Fisher=True so a normal distribution is 0, bias=True). Matches the
        notebook's build_shape_features `kurtosis(...)`."""
        return float(kurtosis(self._data)) if self._data.size > 1 else float("nan")

    # -- Public dispatch -------------------------------------------------------

    def fit(self, model: str, **kwargs) -> dict:
        """
        Run the requested fit and return its results as a dict.

        Parameters
        ----------
        model : str
            'fourparam' or 'kde'.
        **kwargs
            Extra keyword arguments forwarded to the chosen fit. 'fourparam'
            accepts ``max_nfev`` (curve_fit evaluation cap, default 10000);
            'kde' accepts ``bw_method``, ``min_prominence_frac``,
            ``grid_size`` and ``pad_frac`` (see ``_fit_kde`` / ``gene_peaks``).

        Returns
        -------
        dict
            For ``"fourparam"`` — everything compute_fourparam_table records for
            one gene, except has_minus_one_peak, plus the geometry behind the
            ratio:

                gene, y0, A, x0, w, sumsquarevalue,
                ti_fourparam_sigma_dist, truncationindex,
                min, max, mean, std, skew, kurt, right, maxheight, rightheight,
                n_obs, fit_success

            (truncationindex == rightheight / maxheight).

            For ``"kde"`` — the detected modes of a Gaussian KDE:

                gene, n_peaks, peaks, bw_method, n_obs, fit_success

            where ``peaks`` is ``{expression_value: {"height", "prominence"}}``
            (identical to ``gene_peaks`` / ``peaks.json``) and
            ``n_peaks == len(peaks)``.
        """
        if model not in self._FIT_REGISTRY:
            raise ValueError(
                f"Unknown model '{model}'. Supported: {list(self._FIT_REGISTRY)}"
            )
        return getattr(self, self._FIT_REGISTRY[model])(**kwargs)

    # -- fourparam fit ---------------------------------------------------------

    def _fit_fourparam(self, max_nfev: int = 10_000) -> dict:
        """
        Fit the 4-parameter Gaussian to the histogram bin counts by ordinary
        least squares (Trust Region Reflective), minimising the residual sum
        of squares, then compute the truncation-index metrics.

        ``max_nfev`` caps the curve_fit function evaluations. Genuine fits
        converge well under the default; lowering it (e.g. 2000) mainly cuts
        wasted effort on non-converging genes (which otherwise burn the whole
        budget before raising), at the small risk of flipping a slow-converging
        gene to ``fit_success=False``.
        """
        x = self.hist_centers
        y = self.hist_counts

        # Initial guesses from the histogram shape.
        y0_0 = float(y.min())
        A_0 = float(y.max() - y0_0)
        x0_0 = float(x[np.argmax(y)])

        half_max = y0_0 + A_0 / 2.0
        above = x[y >= half_max]
        if len(above) > 1:
            fwhm = float(above[-1] - above[0])
            w_0 = fwhm / (2.0 * np.sqrt(np.log(2.0)))
        else:
            w_0 = float(np.std(self._data)) * np.sqrt(2.0)
        w_0 = max(w_0, 1e-6)

        p0 = [y0_0, A_0, x0_0, w_0]

        try:
            popt, _ = curve_fit(
                _fourparam_gaussian,
                x, y,
                p0=p0,
                bounds=([-np.inf, 0.0, -np.inf, 1e-6],
                        [np.inf, np.inf, np.inf, np.inf]),
                method="trf",          # ordinary least squares (linear loss)
                max_nfev=max_nfev,
            )
        except (RuntimeError, ValueError) as exc:
            raise RuntimeError(
                f"fourparam fit failed for gene '{self.gene_name}': {exc}"
            ) from exc

        y0_fit, A_fit, x0_fit, w_fit = (float(p) for p in popt)
        residuals = y - _fourparam_gaussian(x, *popt)
        sumsquare = float(np.sum(residuals ** 2))

        self.fourparam_y0 = y0_fit
        self.fourparam_A = A_fit
        self.fourparam_x0 = x0_fit
        self.fourparam_w = w_fit
        self.fourparam_sumsquare = sumsquare
        self.active_fits["fourparam"] = True   # enables the metric properties below

        return {
            "gene": self.gene_name,
            "y0": y0_fit,
            "A": A_fit,
            "x0": x0_fit,
            "w": w_fit,
            "sumsquarevalue": sumsquare,
            "ti_fourparam_sigma_dist": self.ti_fourparam_sigma_dist,
            "truncationindex": self.truncationindex,
            "min": self.min(),
            "max": self.max(),
            "mean": self.mean(),
            "std": self.std(),
            "skew": self.skew(),
            "kurt": self.kurt(),
            "right": self._x_max,
            "maxheight": self.maxheight,
            "rightheight": self.rightheight,
            "n_obs": int(self._data.size),
            "fit_success": True,
        }

    # -- kde fit ---------------------------------------------------------------

    def _fit_kde(self, bw_method="silverman", min_prominence_frac=0.05,
                 grid_size=1000, pad_frac=0.05):
        """
        Fit a Gaussian KDE to the gene's expression values (bin-independent) and
        detect its modes by reusing the module-level ``gene_peaks`` — so the
        peaks reported here match ``generate_peaks.py`` / ``peaks.json`` exactly.

        The KDE density is also evaluated on the same padded grid and cached so
        ``hist(lines=["kde"])`` can overlay the curve. If the KDE is singular
        (e.g. zero-spread data) the density is left as None and ``fit_success``
        is False, but peak detection still returns ``{}`` rather than raising.

        Parameters mirror ``gene_peaks``: ``bw_method`` (KDE bandwidth),
        ``min_prominence_frac`` (peak prominence threshold as a fraction of the
        max density), ``grid_size`` and ``pad_frac`` (grid resolution / padding).
        """
        arr = self._data
        lo, hi = float(arr.min()), float(arr.max())
        pad = (hi - lo) * pad_frac
        grid = np.linspace(lo - pad, hi + pad, grid_size)

        try:
            kde = gaussian_kde(arr, bw_method=bw_method)
            density = kde(grid)
            fit_ok = True
        except np.linalg.LinAlgError:
            kde, density, fit_ok = None, None, False

        # Single source of truth for peak detection.
        peaks = gene_peaks(
            arr,
            min_prominence_frac=min_prominence_frac,
            bw_method=bw_method,
            grid_size=grid_size,
            pad_frac=pad_frac,
        )

        self.kde_object = kde
        self.kde_grid = grid if fit_ok else None
        self.kde_density = density
        self.kde_peaks = peaks
        self.kde_bw_method = bw_method
        self.active_fits["kde"] = True

        return {
            "gene": self.gene_name,
            "n_peaks": len(peaks),
            "peaks": peaks,
            "bw_method": bw_method,
            "n_obs": int(arr.size),
            "fit_success": fit_ok,
        }

    # -- Truncation-index metrics ----------------------------------------------

    @property
    def ti_fourparam_sigma_dist(self):
        """
        sigma-distance truncation index:  (x_max - x0) / (w / sqrt(2)).

        How many fitted sigmas x_max lies above the fitted peak x0.
        Lower = ceiling closer to the peak = stronger truncation.
        """
        if not self.active_fits.get("fourparam"):
            raise RuntimeError("fourparam fit has not been run. Call fit('fourparam') first.")
        sigma_fp = self.fourparam_w / np.sqrt(2.0)
        return float((self._x_max - self.fourparam_x0) / sigma_fp)

    # Grid resolution used when scanning the fitted curve over the histogram
    # interval for its min / max.
    _CURVE_GRID = 600

    def _curve_baseline(self):
        """
        Baseline subtracted from both maxheight and rightheight: the **minimum
        value the fitted curve attains over the histogram interval**
        ``[hist_edges[0], hist_edges[-1]]``.

        Using the curve's own interval-minimum (rather than f evaluated at the
        data minimum) guarantees ``0 <= rightheight <= maxheight``, so the
        truncationindex ratio is bounded to [0, 1].
        """
        x_range = np.linspace(self.hist_edges[0], self.hist_edges[-1], self._CURVE_GRID)
        return float(self.fourparam_function(x_range).min())

    @property
    def maxheight(self):
        """
        Full height of the fitted curve over the histogram interval, above the
        curve's interval-minimum baseline:  max(f) - min(f).
        Denominator of the truncationindex ratio.
        """
        if not self.active_fits.get("fourparam"):
            raise RuntimeError("fourparam fit has not been run. Call fit('fourparam') first.")
        x_range = np.linspace(self.hist_edges[0], self.hist_edges[-1], self._CURVE_GRID)
        return float(self.fourparam_function(x_range).max() - self._curve_baseline())

    @property
    def rightheight(self):
        """
        Height of the fitted curve at the right ceiling x_max, above the
        curve's interval-minimum baseline:  f(x_max) - min(f).
        Numerator of the truncationindex ratio.
        """
        if not self.active_fits.get("fourparam"):
            raise RuntimeError("fourparam fit has not been run. Call fit('fourparam') first.")
        return float(self.fourparam_function(self._x_max) - self._curve_baseline())

    @property
    def truncationindex(self):
        """
        Height-ratio truncation index (formerly ti_fourparam_height_ratio):
        rightheight / maxheight  ==  f(x_max)/f(peak) with the curve's
        interval-minimum baseline min(f) subtracted from both. Because that
        baseline is the curve's true minimum over the interval, the ratio is
        **bounded to [0, 1]**: 0 = ceiling sits at the curve minimum, 1 =
        ceiling sits at the peak. Higher = stronger truncation.

        Returns NaN for the degenerate case where maxheight == 0 (the fitted
        curve is flat over the interval, so there is no height to form a ratio).
        """
        if not self.active_fits.get("fourparam"):
            raise RuntimeError("fourparam fit has not been run. Call fit('fourparam') first.")
        mh = self.maxheight
        if mh == 0:
            return float("nan")
        ratio = self.rightheight / mh
        # The ratio is mathematically in [0, 1]; clamp away floating-point noise
        # (e.g. ~-1e-17 when the ceiling coincides with the curve minimum) so the
        # documented bound holds exactly.
        return float(min(1.0, max(0.0, ratio)))

    # -- Evaluate fitted curve -------------------------------------------------

    def fourparam_function(self, x):
        """Evaluate the fitted 4-parameter Gaussian at x."""
        if not self.active_fits.get("fourparam"):
            raise RuntimeError("fourparam fit has not been run.")
        return _fourparam_gaussian(
            np.asarray(x, dtype=float),
            self.fourparam_y0, self.fourparam_A,
            self.fourparam_x0, self.fourparam_w,
        )

    def kde_function(self, x):
        """Evaluate the fitted Gaussian KDE density at x."""
        if not self.active_fits.get("kde"):
            raise RuntimeError("kde fit has not been run.")
        if self.kde_object is None:
            raise RuntimeError(
                f"kde fit did not converge for gene '{self.gene_name}'."
            )
        return self.kde_object(np.asarray(x, dtype=float))

    # -- Visualisation ---------------------------------------------------------

    def hist(self, lines=None):
        """
        Plot the 40-bin histogram and optionally overlay fitted curves.

        Parameters
        ----------
        lines : list of str, optional
            Fits to overlay, e.g. ['fourparam']. A fit is only drawn if it is
            recognised and has already been run via ``fit``.
        """
        fig, ax = plt.subplots(figsize=(9, 5))

        ax.bar(
            self.hist_centers,
            self.hist_counts,
            width=np.diff(self.hist_edges).mean(),
            color="steelblue", alpha=0.6,
            label="Observed", zorder=2,
        )

        x_smooth = np.linspace(self.hist_edges[0], self.hist_edges[-1], 600)

        if lines:
            for fit_name in lines:
                if fit_name not in self._FIT_REGISTRY:
                    print(f"Warning: '{fit_name}' not recognised -- skipping.")
                    continue
                if not self.active_fits.get(fit_name):
                    print(f"Warning: '{fit_name}' not fitted yet -- skipping.")
                    continue

                if fit_name == "fourparam":
                    label = (
                        f"4-param Gaussian (histogram fit)\n"
                        f"A={self.fourparam_A:.3g}, x0={self.fourparam_x0:.3g}, "
                        f"w={self.fourparam_w:.3g}, y0={self.fourparam_y0:.3g}\n"
                        f"sumsquare={self.fourparam_sumsquare:.4g}\n"
                        f"sigma_dist={self.ti_fourparam_sigma_dist:.4f}, "
                        f"truncationindex={self.truncationindex:.4f}"
                    )
                    ax.plot(x_smooth, self.fourparam_function(x_smooth),
                            color="crimson", linewidth=2, label=label, zorder=3)

                elif fit_name == "kde":
                    if self.kde_object is None:
                        print(f"Warning: 'kde' did not converge for "
                              f"'{self.gene_name}' -- skipping.")
                        continue
                    # Scale the density (integrates to 1) onto the bin-count axis.
                    scale = self._data.size * np.diff(self.hist_edges).mean()
                    label = (
                        f"Gaussian KDE (bw={self.kde_bw_method})\n"
                        f"peaks={len(self.kde_peaks)}"
                    )
                    ax.plot(x_smooth, self.kde_function(x_smooth) * scale,
                            color="darkgreen", linewidth=2, label=label, zorder=3)
                    for peak_value, info in self.kde_peaks.items():
                        ax.plot(peak_value, info["height"] * scale, "v",
                                color="darkgreen", markersize=10, zorder=4)

        ax.set_title(self.gene_name, fontsize=14, fontweight="bold")
        ax.set_xlabel("Gene Expression Level", fontsize=12)
        ax.set_ylabel("Frequency (bin count)", fontsize=12)
        ax.legend(fontsize=9)
        ax.grid(True, linestyle="--", alpha=0.4, zorder=1)
        plt.tight_layout()
        plt.show()
        return ax

    # -- Dunder ----------------------------------------------------------------

    def __repr__(self):
        active = [k for k, v in self.active_fits.items() if v]
        parts = [
            f"BhuvanFitter(gene='{self.gene_name}', n={len(self._data)}, "
            f"active_fits={active})"
        ]
        if self.active_fits.get("fourparam"):
            parts.append(
                f"  y0={self.fourparam_y0:.3g}, A={self.fourparam_A:.3g}, "
                f"x0={self.fourparam_x0:.3g}, w={self.fourparam_w:.3g}, "
                f"sumsquare={self.fourparam_sumsquare:.4g}\n"
                f"  ti_fourparam_sigma_dist={self.ti_fourparam_sigma_dist:.4f}, "
                f"truncationindex={self.truncationindex:.4f}"
            )
        if self.active_fits.get("kde"):
            parts.append(
                f"  kde(bw={self.kde_bw_method}): n_peaks={len(self.kde_peaks)}, "
                f"peaks={sorted(self.kde_peaks)}"
            )
        return "\n".join(parts)
