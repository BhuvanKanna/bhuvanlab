# -*- coding: utf-8 -*-
"""
compute_tailmodel.py

Fit ``tailmodel.fit_tail_model`` to every gene of ONE tissue matrix and write a
side-car table joined on ``gene``, like ``qc/`` and ``r2/``:

    ../tailmodel/v11_log2_<tissue>_tailmodel.csv

**Reads the raw rows, deliberately.** The model treats donors at TPM = 0 (exactly
-1) as censored, and their count is the information it uses. The excluded
tables throw that count away, so there is no excluded variant of this table.

    cd fourparam
    python compute_tailmodel.py --input ../data/v11_log2_muscle_skeletal.csv.gz
    python compute_tailmodel.py --input ... --force --jobs 8

Columns: gene, genename, mu, sigma (latent, untruncated), x_min, x_max
(detected extremes), rti_missing, lti_missing (estimated fraction of donors
missing past each edge; lti is NaN when any donor is censored), rti_lr, lti_lr,
n_total, n_censored, left_mode, success, usable. Read the metrics only where
``usable`` is True (see tailmodel.fit_tail_model). Under a true Gaussian of N
donors ``rti_missing`` is ~1/(N+1), not 0.
"""
import argparse
import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from tailmodel import fit_tail_model

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE.parent / "tailmodel"
COLUMNS = ["gene", "genename", "mu", "sigma", "x_min", "x_max", "rti_missing",
           "lti_missing", "rti_lr", "lti_lr", "n_total", "n_censored",
           "left_mode", "success", "usable"]
CHUNK_ROWS = 2000


def out_csv_for(input_path: Path) -> Path:
    stem = input_path.name.split(".csv")[0]
    return OUT_DIR / f"{stem}_tailmodel.csv"


def _fit(item):
    gene, name, values = item
    r = fit_tail_model(values)
    r.update(gene=gene, genename=name)
    return r


def _items(input_path: Path, id_col: str, name_col: str):
    for chunk in pd.read_csv(input_path, chunksize=CHUNK_ROWS):
        samples = [c for c in chunk.columns if c not in (id_col, name_col)]
        mat = chunk[samples].to_numpy(dtype=float)
        for gene, name, row in zip(chunk[id_col], chunk[name_col], mat):
            yield gene, name, row


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--id-col", default="Name")
    p.add_argument("--name-col", default="Description")
    p.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    p.add_argument("--force", action="store_true")
    a = p.parse_args()

    out = out_csv_for(a.input)
    if out.exists() and not a.force:
        raise SystemExit(f"{out} exists; pass --force to overwrite")

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        rows = list(ex.map(_fit, _items(a.input, a.id_col, a.name_col),
                           chunksize=64))
    df = pd.DataFrame(rows)[COLUMNS]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False, lineterminator="\n")
    ok = df[df.usable]
    print(f"{out.name}: {len(df):,} genes, {int(df.success.sum()):,} fit, "
          f"{len(ok):,} usable, "
          f"{int((ok.left_mode == 'truncated').sum()):,} with a defined LTI, "
          f"{time.time() - t0:.0f}s")
    print(f"median rti_missing {ok.rti_missing.median():.4g}  "
          f"(null ~1/(N+1) = {np.median(1 / (ok.n_total + 1)):.4g})")


if __name__ == "__main__":
    main()
