# -*- coding: utf-8 -*-
"""
compute_idklti.py

``idklti`` for every gene of every excluded table: the existing 4-parameter LTI,
**forced to 0 for genes whose left edge sits on the detection floor**. Writes a
side-car joined on ``gene``, like ``qc/`` and ``r2/``:

    ../idklti/v11_log2_<tissue>_idklti_excluded_at_or_below_-1.csv   (gene, idklti)

The rule: if the gene's smallest detected value (``min`` in the excluded table)
is <= -0.75, ``idklti = 0``; otherwise ``idklti = lti``. A failed fit stays blank.

-0.75 is TPM <= 0.19 under ``log2(TPM+1)-1`` -- 0.25 log2 units of slack above
the -1 floor, and the same trace-expression ``FLOOR`` the QC cascade
(``normality.py``) uses. For those genes the observed floor is the assay, not
biology, so their LTI says "not measurable" and is recorded as 0. Muscle
skeletal: the floor rule covers 80.6% of genes with an LTI (72.7% at -0.9,
84.1% at -0.5, so the cut sits on a plateau).

**Zeroing is itself expression-dependent.** The zeroed genes are the
low-expression ones, so a marginal comparison on ``idklti`` will reward the rule,
not the biology. Score it on expression-matched controls only.

Nothing is refit, and a kept ``lti`` is copied as **text**, so it is
byte-identical to the table it came from (parsing and re-writing a float
shortens it). The browser computes the same value from the loaded row
(``IDK_FLOOR`` in docs/index.html); ``tests/test_idklti.py`` pins the two
cut-offs together.

    cd fourparam
    python compute_idklti.py --all
    python compute_idklti.py --tissues muscle_skeletal,liver
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUTPUTS = HERE.parent / "outputs"
OUT_DIR = HERE.parent / "idklti"
SUFFIX = "_fourparam_excluded_at_or_below_-1.csv"
FLOOR_CUT = -0.75


def idklti(lti, xmin, cut=FLOOR_CUT):
    """LTI with floor-bordering genes set to 0; NaN where LTI is NaN."""
    lti = np.asarray(lti, dtype=float)
    xmin = np.asarray(xmin, dtype=float)
    return np.where(np.isfinite(lti) & (xmin <= cut), 0.0, lti)


def idklti_text(lti: pd.Series, xmin: pd.Series, cut=FLOOR_CUT) -> pd.Series:
    """Same rule on the tables' own text: "0", the lti string verbatim, or ""."""
    lti_v = pd.to_numeric(lti, errors="coerce")
    min_v = pd.to_numeric(xmin, errors="coerce")
    floor = np.isfinite(lti_v) & (min_v <= cut)
    out = lti.where(np.isfinite(lti_v), "")
    return out.mask(floor, "0")


def run_one(table: Path) -> dict:
    tissue = table.name[len("v11_log2_"):-len(SUFFIX)]
    df = pd.read_csv(table, usecols=["gene", "min", "lti"], dtype=str,
                     keep_default_na=False)
    df["idklti"] = idklti_text(df["lti"], df["min"])
    out = OUT_DIR / f"v11_log2_{tissue}_idklti_excluded_at_or_below_-1.csv"
    df[["gene", "idklti"]].to_csv(out, index=False, lineterminator="\n")

    lti_v = pd.to_numeric(df["lti"], errors="coerce")
    defined = np.isfinite(lti_v)
    forced = defined & (pd.to_numeric(df["min"], errors="coerce") <= FLOOR_CUT)
    return {"tissue": tissue, "genes": len(df), "defined": int(defined.sum()),
            "forced": int((forced & (lti_v != 0)).sum()),
            "zero": int((pd.to_numeric(df["idklti"], errors="coerce") == 0).sum()),
            "nonzero": int((defined & ~forced & (lti_v != 0)).sum())}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--all", action="store_true")
    g.add_argument("--tissues", type=str)
    a = p.parse_args()

    tables = sorted(OUTPUTS.glob(f"v11_log2_*{SUFFIX}"))
    if a.tissues:
        want = set(a.tissues.split(","))
        tables = [t for t in tables
                  if t.name[len("v11_log2_"):-len(SUFFIX)] in want]
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    stats = [run_one(t) for t in tables]
    s = pd.DataFrame(stats)
    print(s.to_string(index=False))
    tot = s[["defined", "forced", "nonzero"]].sum()
    print(f"\n{len(s)} tissues: {int(tot.defined):,} defined, "
          f"{int(tot.forced):,} forced to 0 by the floor rule "
          f"({tot.forced / tot.defined:.1%}), {int(tot.nonzero):,} non-zero")


if __name__ == "__main__":
    main()
