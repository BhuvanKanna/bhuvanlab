"""dRTI building blocks (../drti/PLAN.md): orientation, invariance, matching."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import drti_core as core  # noqa: E402
from drti_evaluate import match2  # noqa: E402
from drti_final import knn_med_mad  # noqa: E402


def _panel(n=1500, t=20, icc=0.3, cut=None, seed=0):
    rng = np.random.default_rng(seed)
    L = rng.normal(size=n)
    if cut is not None:
        L = L[L < np.quantile(L, 1 - cut)]
    e = rng.normal(scale=np.sqrt((1 - icc) / icc), size=(len(L), t))
    Z = L[:, None] + e
    Z[rng.random(Z.shape) < 0.2] = np.nan
    return Z


def test_truncation_raises_every_statistic():
    full, cut = core.all_stats(_panel()), core.all_stats(_panel(cut=0.10))
    for key in [("mean", s) for s in core.STATS] + [("mean", "edge_lr")]:
        assert cut[key] > full[key], key


def test_quantile_stats_are_location_scale_invariant():
    rng = np.random.default_rng(1)
    v = rng.normal(size=2000)
    for name, f in core.STATS.items():
        if name == "xi_gpd":          # MLE: invariant only up to optimizer tolerance
            continue
        assert np.isclose(f(v), f(3.0 * v - 7.0)), name


def test_edge_lr_missing_fraction_tracks_cut():
    Z = _panel(n=4000, t=30, icc=0.5, cut=0.15)
    o = core.objects(Z)
    lr, miss = core.edge_lr(o["mean"], o["s2"])
    assert lr > 3 and 0.05 < miss < 0.35


def test_icc_recovers_truth():
    assert abs(core.icc_of(_panel(n=3000, icc=0.3)) - 0.3) < 0.05


def test_match2_respects_both_calipers():
    rng = np.random.default_rng(2)
    pos = pd.DataFrame({"gene": [f"P{i}" for i in range(50)],
                        "expr": rng.uniform(0, 5, 50), "icc": rng.uniform(0, 0.4, 50)})
    neg = pd.DataFrame({"gene": [f"N{i}" for i in range(300)],
                        "expr": rng.uniform(0, 5, 300), "icc": rng.uniform(0, 0.4, 300)})
    pairs = match2(pos, neg)
    p, n = pos.set_index("gene"), neg.set_index("gene")
    assert pairs and len({b for _, b in pairs}) == len(pairs)
    for a, b in pairs:
        assert abs(p.expr[a] - n.expr[b]) <= 0.1 and abs(p.icc[a] - n.icc[b]) <= 0.03


def test_knn_calibration_removes_a_smooth_trend():
    rng = np.random.default_rng(3)
    k = 13000                         # about the real number of scored genes
    tab = pd.DataFrame({"expr": rng.uniform(0, 8, k), "icc": rng.uniform(0, 0.5, k),
                        "n": rng.integers(300, 1000, k)})
    tab["s"] = 0.3 * tab.expr + rng.normal(size=k)
    med, mad = knn_med_mad(tab, "s")
    z = (tab.s - med) / mad
    before = np.corrcoef(tab.s, tab.expr)[0, 1]
    after = np.corrcoef(z, tab.expr)[0, 1]
    assert before > 0.5 and abs(after) < 0.1
    assert 0.7 < np.std(z) < 1.4
