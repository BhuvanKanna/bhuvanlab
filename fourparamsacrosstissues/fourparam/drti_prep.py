# -*- coding: utf-8 -*-
"""
drti_prep.py

Per-tissue residual expression for the donor-level right-truncation work
(see ../drti/PLAN.md).

For each tissue with >= MIN_DONORS samples, on x = log2(TPM + 1):
  1. keep protein-coding genes with median TPM >= MIN_MEDIAN_TPM
  2. regress out donor/sample covariates (sex, age, Hardy scale, RIN,
     ischemic time, exonic rate, center)
  3. regress out the top N_PCS principal components of the step-2 residuals
     (composition, batch and other broad technical axes)
  4. robust-standardize each gene: (r - median) / (1.4826 * MAD)

Also keeps the robust-standardized raw values (no step 2-3) for comparison.

    python drti_prep.py --meta <dir with GTEx v11 SampleAttributesDS/SubjectPhenotypesDS>
                        [--genes N] [--out DIR] [--jobs 4]

`--genes N` keeps a fixed random sample of N protein-coding genes (seed 0),
for viability runs; omit it for every gene. Writes one .npz per tissue.
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
DOCS = HERE.parent.parent / "docs"

MIN_DONORS = 100
MIN_MEDIAN_TPM = 1.0
N_PCS = 10


def strip_version(g: str) -> str:
    return g.split(".", 1)[0]


def protein_coding_ids() -> list[str]:
    man = json.loads((DOCS / "manifest.json").read_text(encoding="utf-8"))["gene_type"]
    text = (DOCS / man["base_url"] / man["file"]).read_text(encoding="ascii")
    ids = [strip_version(l.split("\t", 1)[0])
           for l in (DOCS / "genes.tsv").open(encoding="utf-8") if l.strip()]
    return sorted(g for g, c in zip(ids, text) if man["classes"].get(c) == "protein_coding")


def donor_of(sample: str) -> str:
    return "-".join(sample.split("-")[:2])


def load_covariates(meta: Path) -> pd.DataFrame:
    sa = pd.read_csv(meta / "v11_SampleAttributesDS.txt", sep="\t", low_memory=False,
                     usecols=["SAMPID", "SMRIN", "SMTSISCH", "SMEXNCRT", "SMCENTER"])
    sp = pd.read_csv(meta / "v11_SubjectPhenotypesDS.txt", sep="\t")
    sa["SUBJID"] = sa.SAMPID.map(donor_of)
    d = sa.merge(sp, on="SUBJID", how="left").set_index("SAMPID")
    d["age"] = d.AGE.str[:2].astype(float) + 4.5
    return d


def design(samples: list[str], cov: pd.DataFrame) -> np.ndarray:
    c = cov.reindex(samples)
    cols = [np.ones(len(samples))]
    for k in ("SMRIN", "SMTSISCH", "SMEXNCRT", "age"):
        v = c[k].astype(float).to_numpy()
        v = np.where(np.isfinite(v), v, np.nanmedian(v) if np.isfinite(v).any() else 0.0)
        if v.std() > 0:
            cols.append((v - v.mean()) / v.std())
    for k in ("SEX", "DTHHRDY", "SMCENTER"):
        lev = c[k].astype(str).fillna("NA")
        for level in sorted(lev.unique())[1:]:          # drop one level as reference
            ind = (lev == level).to_numpy(float)
            if 0 < ind.sum() < len(ind):
                cols.append(ind)
    return np.column_stack(cols)


def residualize(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """y: samples x genes. OLS residuals on X."""
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def robust_z(r: np.ndarray) -> np.ndarray:
    med = np.median(r, axis=0)
    mad = 1.4826 * np.median(np.abs(r - med), axis=0)
    mad[mad <= 0] = np.nan
    return (r - med) / mad


def prep_tissue(path: Path, cov: pd.DataFrame, keep_genes: set[str] | None,
                pc: set[str], keep_raw: bool = True) -> dict | None:
    d = pd.read_csv(path, index_col=0).drop(columns=["Description"])
    if d.shape[1] < MIN_DONORS:
        return None
    d.index = d.index.map(strip_version)
    d = d[d.index.isin(pc)]
    x = d.to_numpy(float).T + 1.0                     # samples x genes, log2(TPM+1)
    expressed = np.median(x, axis=0) >= np.log2(MIN_MEDIAN_TPM + 1)
    x, genes = x[:, expressed], d.index[expressed].to_numpy()
    samples = list(d.columns)
    X = design(samples, cov)
    r1 = residualize(x, X)
    s = (r1 - r1.mean(0)) / np.where(r1.std(0) > 0, r1.std(0), 1.0)
    u, sv, _ = np.linalg.svd(s, full_matrices=False)
    pcs = u[:, :N_PCS]
    r2 = residualize(r1, np.column_stack([X, pcs]))
    sel = np.isin(genes, list(keep_genes)) if keep_genes is not None else np.ones(len(genes), bool)
    out = {
        "tissue": path.name[len("v11_log2_"):-len(".csv.gz")],
        "genes": genes[sel],
        "donors": np.array([donor_of(s_) for s_ in samples]),
        "z": robust_z(r2[:, sel]).astype(np.float32),        # donors x genes
        "median_log": np.median(x[:, sel], axis=0).astype(np.float32),
        "pc_var": (sv[:N_PCS] ** 2 / (sv ** 2).sum()).astype(np.float32),
    }
    if keep_raw:
        out["z_raw"] = robust_z(x[:, sel]).astype(np.float32)
    return out


def _job(args):
    path, meta, keep, out, keep_raw = args
    pc = set(protein_coding_ids())
    res = prep_tissue(Path(path), load_covariates(Path(meta)), keep, pc, keep_raw)
    if res is None:
        return None
    np.savez_compressed(Path(out) / f"{res['tissue']}.npz", **res)
    return res["tissue"], len(res["donors"]), len(res["genes"])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--meta", type=Path, required=True)
    p.add_argument("--out", type=Path, default=HERE.parent / "drti" / "prep")
    p.add_argument("--genes", type=int, default=None)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--no-raw", action="store_true",
                   help="Skip the unadjusted z_raw copy (only the viability check reads it).")
    args = p.parse_args()

    keep = None
    if args.genes:
        rng = np.random.default_rng(0)
        allpc = protein_coding_ids()
        keep = set(rng.choice(allpc, size=args.genes, replace=False))
    args.out.mkdir(parents=True, exist_ok=True)
    paths = sorted(str(f) for f in DATA.glob("v11_log2_*.csv.gz"))
    jobs = [(p_, str(args.meta), keep, str(args.out), not args.no_raw) for p_ in paths]
    with ProcessPoolExecutor(args.jobs) as ex:
        for r in ex.map(_job, jobs):
            if r:
                print(f"{r[0]:<45} donors {r[1]:>4}  genes {r[2]:>6}", flush=True)


if __name__ == "__main__":
    main()
