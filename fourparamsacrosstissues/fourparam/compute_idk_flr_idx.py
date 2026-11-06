# -*- coding: utf-8 -*-
"""
compute_idk_flr_idx.py

``idk_flr_idx`` for every gene of every excluded table: the existing 4-parameter flr_idx,
**forced to 0 for genes whose left edge sits on the detection floor**. Writes a
side-car joined on ``gene``, like ``qc/`` and ``r2/``:

    ../idk_flr_idx/v11_log2_<tissue>_idk_flr_idx_excluded_at_or_below_-1.csv   (gene, idk_flr_idx)

The rule: if the gene's smallest detected value (``min`` in the excluded table)
is <= -0.75, ``idk_flr_idx = 0``; otherwise ``idk_flr_idx = flr_idx``. A failed fit stays blank.

-0.75 is TPM <= 0.19 under ``log2(TPM+1)-1`` -- 0.25 log2 units of slack above
the -1 floor, and the same trace-expression ``FLOOR`` the QC cascade
(``normality.py``) uses. For those genes the observed floor is the assay, not
biology, so their flr_idx says "not measurable" and is recorded as 0. Muscle
skeletal: the floor rule covers 80.6% of genes with an flr_idx (72.7% at -0.9,
84.1% at -0.5, so the cut sits on a plateau).

**Zeroing is itself expression-dependent.** The zeroed genes are the
low-expression ones, so a marginal comparison on ``idk_flr_idx`` will reward the rule,
not the biology. Score it on expression-matched controls only.

Nothing is refit, and a kept ``flr_idx`` is copied as **text**, so it is
byte-identical to the table it came from (parsing and re-writing a float
shortens it). The browser computes the same value from the loaded row
(``IDK_FLOOR`` in docs/index.html); ``tests/test_idk_flr_idx.py`` pins the two
cut-offs together.

    cd fourparam
    python compute_idk_flr_idx.py --all
    python compute_idk_flr_idx.py --tissues muscle_skeletal,liver
    python compute_idk_flr_idx.py --worm     # -> ../idk_flr_idx/worm_idk_flr_idx_excluded_at_or_below_-1.csv

The C. elegans table in ``../worm/`` (same ``log2(TPM+1)-1`` scale, same -1
floor) takes the identical rule. Its ``gene`` is a transcript id, so the
side-car is per transcript. The cut does not sit on a plateau there: 31.8% of
worm transcripts have ``min`` <= -0.9, 41.4% <= -0.75, 48.9% <= -0.5.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUTPUTS = HERE.parent / "outputs"
OUT_DIR = HERE.parent / "idk_flr_idx"
SUFFIX = "_fourparam_excluded_at_or_below_-1.csv"
WORM_TABLE = HERE.parent / "worm" / "worm_fourparam_excluded_at_or_below_-1.csv"
FLOOR_CUT = -0.75


def idk_flr_idx(flr_idx, xmin, cut=FLOOR_CUT):
    """flr_idx with floor-bordering genes set to 0; NaN where flr_idx is NaN."""
    flr_idx = np.asarray(flr_idx, dtype=float)
    xmin = np.asarray(xmin, dtype=float)
    return np.where(np.isfinite(flr_idx) & (xmin <= cut), 0.0, flr_idx)


def idk_flr_idx_text(flr_idx: pd.Series, xmin: pd.Series, cut=FLOOR_CUT) -> pd.Series:
    """Same rule on the tables' own text: "0", the flr_idx string verbatim, or ""."""
    flr_idx_v = pd.to_numeric(flr_idx, errors="coerce")
    min_v = pd.to_numeric(xmin, errors="coerce")
    floor = np.isfinite(flr_idx_v) & (min_v <= cut)
    out = flr_idx.where(np.isfinite(flr_idx_v), "")
    return out.mask(floor, "0")


def run_one(table: Path, tissue: str = None, out: Path = None) -> dict:
    tissue = tissue or table.name[len("v11_log2_"):-len(SUFFIX)]
    df = pd.read_csv(table, usecols=["gene", "min", "flr_idx"], dtype=str,
                     keep_default_na=False)
    df["idk_flr_idx"] = idk_flr_idx_text(df["flr_idx"], df["min"])
    out = out or OUT_DIR / f"v11_log2_{tissue}_idk_flr_idx_excluded_at_or_below_-1.csv"
    df[["gene", "idk_flr_idx"]].to_csv(out, index=False, lineterminator="\n")

    flr_idx_v = pd.to_numeric(df["flr_idx"], errors="coerce")
    defined = np.isfinite(flr_idx_v)
    forced = defined & (pd.to_numeric(df["min"], errors="coerce") <= FLOOR_CUT)
    return {"tissue": tissue, "genes": len(df), "defined": int(defined.sum()),
            "forced": int((forced & (flr_idx_v != 0)).sum()),
            "zero": int((pd.to_numeric(df["idk_flr_idx"], errors="coerce") == 0).sum()),
            "nonzero": int((defined & ~forced & (flr_idx_v != 0)).sum())}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--all", action="store_true")
    g.add_argument("--tissues", type=str)
    g.add_argument("--worm", action="store_true",
                   help="the C. elegans table in ../worm/ instead of a tissue")
    a = p.parse_args()

    if a.worm:
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        s = run_one(WORM_TABLE, tissue="worm",
                    out=OUT_DIR / "worm_idk_flr_idx_excluded_at_or_below_-1.csv")
        print(pd.DataFrame([s]).to_string(index=False))
        return

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
