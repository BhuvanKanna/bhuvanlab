# -*- coding: utf-8 -*-
"""
simulate_utc.py

Validation V1-V3 from ../utc/PLAN.md, run before the evaluation.

    cd fourparam
    python simulate_utc.py          # -> ../utc/simulation.csv and a printed verdict

UTC is invariant to shift and scale, so the scenarios use unit-scale draws.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import truncnorm

from utc import tail_metrics

OUT = Path(__file__).resolve().parent.parent / "utc"
N, REPS, SEED = 300, 2000, 0


def scenarios(rng: np.random.Generator) -> dict[str, np.ndarray]:
    s = {"V1 gaussian": rng.standard_normal((REPS, N))}
    for k in (1.0, 1.5, 2.0):
        s[f"V2 truncated +{k}sd"] = truncnorm.rvs(-np.inf, k, size=(REPS, N), random_state=rng)
    low = rng.random((REPS, N)) < 0.15
    s["V3 contamination 15% at -2sd"] = rng.standard_normal((REPS, N)) - 2.0 * low
    return s


def main() -> None:
    rng = np.random.default_rng(SEED)
    rows = []
    for name, x in scenarios(rng).items():
        m = tail_metrics(x)
        for k in ("utc", "ltc"):
            v = m[k][np.isfinite(m[k])]
            rows.append({"scenario": name, "metric": k, "mean": v.mean(), "sd": v.std(ddof=1),
                         "p05": np.percentile(v, 5), "p95": np.percentile(v, 95)})
    r = pd.DataFrame(rows)
    OUT.mkdir(exist_ok=True)
    r.to_csv(OUT / "simulation.csv", index=False, float_format="%.4f")
    print(r.round(3).to_string(index=False))

    utc = r[r.metric == "utc"].set_index("scenario")["mean"]
    ref = utc["V1 gaussian"]
    v2 = all(utc[s] < ref for s in utc.index if s.startswith("V2"))
    v3 = abs(utc["V3 contamination 15% at -2sd"] - ref) <= 0.05
    print(f"\nV2 (truncation lowers UTC): {'PASS' if v2 else 'FAIL'}")
    print(f"V3 (contamination leaves UTC within 0.05): {'PASS' if v3 else 'FAIL'} "
          f"(diff {utc['V3 contamination 15% at -2sd'] - ref:+.3f})")


if __name__ == "__main__":
    main()
