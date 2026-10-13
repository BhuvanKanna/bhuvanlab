# -*- coding: utf-8 -*-
"""
utc.py

Upper-tail compression (UTC) and its mirror (LTC), as fixed in ../utc/PLAN.md.

    cd fourparam
    python utc.py --jobs 4          # every tissue with >= 150 donors -> ../utc/

Per gene x tissue, on x = log2(TPM + 1) over all donors:

    UTC = ln[(q97.5 - q50) / (q75 - q50)]   upper half only
    LTC = ln[(q50 - q2.5)  / (q50 - q25)]   lower half only
    VAR = ln(q75 - q25)                     overall spread

Each is adjusted for expression within its tissue (residual from a rolling
median over the 201 genes nearest in q50, scaled by MAD) and the gene score is
the median z across eligible tissues. Do not change a threshold here without
logging it under "Deviations" in PLAN.md.

Writes ../utc/gene_tissue.csv.gz (one row per eligible gene x tissue) and
../utc/gene_scores.csv (one row per scored gene).
"""
from __future__ import annotations

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
OUT = HERE.parent / "utc"
DOCS = HERE.parent.parent / "docs"

MIN_DONORS = 150
MIN_MEDIAN_TPM = 5.0
MIN_TISSUES = 5
WINDOW = 201
METRICS = ("utc", "ltc", "var")


def strip_version(gid: str) -> str:
    return gid.split(".", 1)[0]


def protein_coding_ids() -> set[str]:
    man = json.loads((DOCS / "manifest.json").read_text(encoding="utf-8"))["gene_type"]
    text = (DOCS / man["base_url"] / man["file"]).read_text(encoding="ascii")
    ids = [strip_version(l.split("\t", 1)[0])
           for l in (DOCS / "genes.tsv").open(encoding="utf-8") if l.strip()]
    return {g for g, c in zip(ids, text) if man["classes"].get(c) == "protein_coding"}


def tail_metrics(x: np.ndarray) -> dict[str, np.ndarray]:
    """UTC, LTC, VAR and q50 for each row of x (genes x donors, log2(TPM+1)).
    A metric is NaN where its denominator or spread is <= 0."""
    q025, q25, q50, q75, q975 = np.percentile(x, [2.5, 25, 50, 75, 97.5], axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        up_in, up_out = q75 - q50, q975 - q50
        dn_in, dn_out = q50 - q25, q50 - q025
        iqr = q75 - q25
        utc = np.where((up_in > 0) & (up_out > 0), np.log(up_out / up_in), np.nan)
        ltc = np.where((dn_in > 0) & (dn_out > 0), np.log(dn_out / dn_in), np.nan)
        var = np.where(iqr > 0, np.log(iqr), np.nan)
    return {"q50": q50, "utc": utc, "ltc": ltc, "var": var}


def residual_z(q50: np.ndarray, values: np.ndarray, window: int = WINDOW) -> np.ndarray:
    """values minus a rolling median over the `window` genes nearest in q50,
    divided by 1.4826 x MAD of the residuals. Returned in the input order."""
    order = np.argsort(q50, kind="stable")
    s = pd.Series(values[order])
    base = s.rolling(window, center=True, min_periods=1).median().to_numpy()
    resid = np.empty_like(values, dtype=float)
    resid[order] = values[order] - base
    finite = resid[np.isfinite(resid)]
    mad = 1.4826 * np.median(np.abs(finite - np.median(finite))) if finite.size else np.nan
    return resid / mad if mad and np.isfinite(mad) else np.full_like(resid, np.nan)


def score_tissue(path: Path, pc: set[str]) -> pd.DataFrame | None:
    d = pd.read_csv(path, index_col=0)
    d = d.drop(columns=["Description"])
    n_donors = d.shape[1]
    if n_donors < MIN_DONORS:
        return None
    d.index = d.index.map(strip_version)
    d = d[d.index.isin(pc)]
    x = d.to_numpy(float) + 1.0                 # log2(TPM+1) - 1  ->  log2(TPM+1)
    m = tail_metrics(x)
    keep = (2.0 ** m["q50"] - 1.0) >= MIN_MEDIAN_TPM
    out = pd.DataFrame({"gene": d.index[keep]})
    for k in ("q50",) + METRICS:
        out[k] = m[k][keep]
    for k in METRICS:
        out["z_" + k] = residual_z(out["q50"].to_numpy(), out[k].to_numpy())
    out.insert(1, "tissue", path.name[len("v11_log2_"):-len(".csv.gz")])
    out.insert(2, "n_donors", n_donors)
    return out


def aggregate(gt: pd.DataFrame) -> pd.DataFrame:
    g = gt.groupby("gene")
    s = pd.DataFrame({"n_tissues": g.size(), "expr": g["q50"].mean()})
    for k in METRICS:
        s[k] = g["z_" + k].median()
    return s[s.n_tissues >= MIN_TISSUES].reset_index()


def _score(path: str) -> pd.DataFrame | None:
    return score_tissue(Path(path), protein_coding_ids())


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--jobs", type=int, default=4)
    args = p.parse_args()

    paths = sorted(str(f) for f in DATA.glob("v11_log2_*.csv.gz"))
    with ProcessPoolExecutor(args.jobs) as ex:
        parts = [r for r in ex.map(_score, paths) if r is not None]
    gt = pd.concat(parts, ignore_index=True)
    OUT.mkdir(exist_ok=True)
    gt.to_csv(OUT / "gene_tissue.csv.gz", index=False, float_format="%.6g")
    scores = aggregate(gt)
    scores.to_csv(OUT / "gene_scores.csv", index=False, float_format="%.6g")
    print(f"{gt.tissue.nunique()} tissues with >= {MIN_DONORS} donors, "
          f"{len(gt):,} eligible gene x tissue values, "
          f"{len(scores):,} genes scored in >= {MIN_TISSUES} tissues")


if __name__ == "__main__":
    main()
