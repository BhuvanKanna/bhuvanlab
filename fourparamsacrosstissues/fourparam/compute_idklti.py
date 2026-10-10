# -*- coding: utf-8 -*-
"""
compute_idklti.py

Add ``idklti`` to a tissue's tail-model side-car: the existing 4-parameter LTI
(from the excluded table), **forced to 0 for genes whose left edge sits on the
detection floor**.

The rule: if the gene's smallest detected value (``min`` in the excluded table)
is <= -0.75, ``idklti = 0``; otherwise ``idklti = lti``. A failed fit stays NaN.

-0.75 is TPM <= 0.19 under ``log2(TPM+1)-1`` -- 0.25 log2 units of slack above
the -1 floor, and the same trace-expression ``FLOOR`` the QC cascade
(``normality.py``) uses. For those genes the observed floor is the assay, not
biology, so their LTI says "not measurable" and is recorded as 0. Muscle
skeletal: the floor rule covers 80.6% of genes with an LTI (72.7% at -0.9,
84.1% at -0.5, so the cut sits on a plateau).

**Zeroing is itself expression-dependent.** The zeroed genes are the
low-expression ones, so a marginal comparison on ``idklti`` will reward the rule,
not the biology. Score it on expression-matched controls only.

    cd fourparam
    python compute_idklti.py --tissue muscle_skeletal
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUTPUTS = HERE.parent / "outputs"
TAIL_DIR = HERE.parent / "tailmodel"
FLOOR_CUT = -0.75


def idklti(lti, xmin, cut=FLOOR_CUT):
    """LTI with floor-bordering genes set to 0; NaN where LTI is NaN."""
    lti = np.asarray(lti, dtype=float)
    xmin = np.asarray(xmin, dtype=float)
    return np.where(np.isfinite(lti) & (xmin <= cut), 0.0, lti)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    p.add_argument("--tissue", required=True)
    a = p.parse_args()

    src = OUTPUTS / f"v11_log2_{a.tissue}_fourparam_excluded_at_or_below_-1.csv"
    side = TAIL_DIR / f"v11_log2_{a.tissue}_tailmodel.csv"
    fp = pd.read_csv(src, usecols=["gene", "min", "lti"], dtype={"gene": str})
    fp["idklti"] = idklti(fp["lti"], fp["min"])
    lti = pd.to_numeric(fp["lti"], errors="coerce")
    forced = int((np.isfinite(lti) & (fp["min"] <= FLOOR_CUT) & (lti != 0)).sum())

    tm = pd.read_csv(side, dtype={"gene": str}).drop(columns=["idklti"],
                                                      errors="ignore")
    out = tm.merge(fp[["gene", "idklti"]], on="gene", how="left")
    if len(out) != len(tm):
        raise SystemExit("gene join changed the row count")
    out.to_csv(side, index=False, lineterminator="\n")

    v = out["idklti"]
    # "forced" counts only what the rule changed; the old LTI is already exactly
    # 0 for some genes (the shared RTI/LTI baseline), and those are not the rule.
    print(f"{side.name}: idklti defined {int(v.notna().sum()):,}, "
          f"forced to 0 by the floor rule {forced:,}, "
          f"0 in total {int((v == 0).sum()):,}, non-zero {int((v > 0).sum()):,}")


if __name__ == "__main__":
    main()
