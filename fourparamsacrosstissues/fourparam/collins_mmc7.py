# -*- coding: utf-8 -*-
"""
collins_mmc7.py

Join our cross-tissue truncation scores onto Collins et al. 2022 Table S7
(``mmc7.xlsx``: pHaplo, pTriplo and 13 conventional dosage / constraint measures)
and correlate them.

Per gene, across the 54 excluded (<= -1) tables:

    rti_top10     mean of the top 10% of the gene's per-tissue RTI values
    idklti_top10  same for idkLTI (old LTI, 0 where the left edge is the floor)
    expr_mean     mean of the gene's per-tissue ``mean`` (expression covariate)

"Top 10%" = the k = ceil(0.1 * n) largest values, where n is the number of tissues
with a successful fit for that gene (so k = 6 when all 54 tissues fit). Only
tissues with ``fit_success`` and a finite score count.

Ensembl ids are mapped to symbols with the tables' ``genename``. ``_PAR_Y`` copies
are dropped; where one symbol still maps to several ids, the id that fits in the
most tissues is kept. Only genes present in both sets are written.

Correlations are Spearman, plus a partial Spearman controlling for ``expr_mean``
(rank-residualized), because both scores track expression (see CLAUDE.md, "Read
LTI beside mean").

    cd fourparam
    python collins_mmc7.py --mmc7 ../../../mmc7.xlsx
"""
import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
OUTPUTS = HERE.parent / "outputs"
IDKLTI = HERE.parent / "idklti"
OUT_DIR = HERE.parent / "results" / "collins_mmc7"
SUFFIX = "_fourparam_excluded_at_or_below_-1.csv"
SCORES = ["rti_top10", "idklti_top10"]


def top10_mean(values: pd.Series) -> float:
    v = np.sort(values.dropna().to_numpy())[::-1]
    if v.size == 0:
        return np.nan
    return float(v[:math.ceil(0.1 * v.size)].mean())


def load_long() -> pd.DataFrame:
    frames = []
    for table in sorted(OUTPUTS.glob("v11_log2_*" + SUFFIX)):
        tissue = table.name[len("v11_log2_"):-len(SUFFIX)]
        df = pd.read_csv(table, usecols=["gene", "genename", "rti", "mean",
                                         "fit_success"])
        idk = pd.read_csv(IDKLTI / f"v11_log2_{tissue}_idklti_excluded_at_or_below_-1.csv")
        df = df.merge(idk, on="gene", how="left")
        df = df[df["fit_success"].astype(str) == "True"]
        df["tissue"] = tissue
        frames.append(df.drop(columns="fit_success"))
    if len(frames) != 54:
        raise SystemExit(f"expected 54 excluded tables, found {len(frames)}")
    return pd.concat(frames, ignore_index=True)


def per_gene(long: pd.DataFrame) -> pd.DataFrame:
    g = long.groupby("gene")
    out = pd.DataFrame({
        "genename": g["genename"].first(),
        "rti_top10": g["rti"].agg(top10_mean),
        "idklti_top10": g["idklti"].agg(top10_mean),
        "n_tissues": g["rti"].count(),
        "expr_mean": g["mean"].mean(),
    }).reset_index()
    out = out[~out["gene"].str.endswith("_PAR_Y")]
    out = out.sort_values(["genename", "n_tissues"], ascending=[True, False])
    n_dup = int(out.duplicated("genename").sum())
    out = out.drop_duplicates("genename", keep="first")
    print(f"{len(out)} symbols ({n_dup} duplicate-symbol ids dropped)")
    return out


def partial_spearman(x, y, z):
    rx, ry, rz = (stats.rankdata(a) for a in (x, y, z))
    res = [a - np.polyval(np.polyfit(rz, a, 1), rz) for a in (rx, ry)]
    return stats.pearsonr(*res)


def correlate(merged: pd.DataFrame, measures: list) -> pd.DataFrame:
    rows = []
    for score in SCORES:
        for m in measures:
            d = merged[[score, m, "expr_mean"]].apply(pd.to_numeric, errors="coerce").dropna()
            rho, p = stats.spearmanr(d[score], d[m])
            prho, pp = partial_spearman(d[score], d[m], d["expr_mean"])
            erho, _ = stats.spearmanr(d["expr_mean"], d[m])
            rows.append({"score": score, "measure": m, "n": len(d),
                         "spearman_rho": rho, "p": p,
                         "partial_rho_given_expr": prho, "partial_p": pp,
                         "expr_mean_rho": erho})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mmc7", required=True, type=Path)
    args = ap.parse_args()

    mmc7 = pd.read_excel(args.mmc7)
    measures = [c for c in mmc7.columns if c != "Gene"]
    genes = per_gene(load_long())

    merged = mmc7.merge(genes.rename(columns={"gene": "ensembl_id"}),
                        left_on="Gene", right_on="genename", how="inner")
    merged = merged.drop(columns="genename")
    print(f"mmc7 {len(mmc7)} genes, ours {len(genes)}, common {len(merged)}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    merged.to_csv(OUT_DIR / "mmc7_with_rti_idklti_top10.csv", index=False)
    corr = correlate(merged, measures)
    corr.to_csv(OUT_DIR / "mmc7_correlations.csv", index=False)
    with pd.option_context("display.width", 200, "display.precision", 3):
        print(corr.to_string(index=False))


if __name__ == "__main__":
    main()
