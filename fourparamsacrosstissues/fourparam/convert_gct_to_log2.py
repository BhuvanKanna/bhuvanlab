#!/usr/bin/env python
"""
convert_gct_to_log2.py

Convert a raw GTEx ``gene_tpm`` GCT matrix into a ``log2(TPM + 1) - 1`` CSV.

The transform is applied elementwise to every per-sample TPM value:

    log2(TPM + 1) - 1

so TPM = 0 maps to exactly -1 (the floor, no clamping) and every value is
>= -1 by construction. Integer TPMs 0,1,2,3,... -> -1, 0, 0.585, 1, ...

GCT layout (tab-separated): the first two lines are the version tag
(``#1.2`` / ``#1.3``) and the ``<nrows> <ncols> ...`` dims line, so both are
skipped. The header row is ``Name  Description  <samples>`` (GCT 1.2) or
``id  Name  Description  <samples>`` (GCT 1.3); the ``id`` column, if present,
is dropped. ``Name`` (Ensembl id) and ``Description`` (gene symbol) lead.

Usage
-----
    python convert_gct_to_log2.py --input raw.gct.gz --output ../data/v11_log2_liver.csv.gz
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def convert(input_gct: Path, output_csv: Path) -> None:
    print(f"Reading {input_gct.name} ...", flush=True)
    df = pd.read_csv(input_gct, sep="\t", skiprows=2)  # skip #1.x and dims lines
    if "id" in df.columns:
        df = df.drop(columns=["id"])

    meta = ["Name", "Description"]
    sample_cols = [c for c in df.columns if c not in meta]
    print(f"  {df.shape[0]} genes x {len(sample_cols)} samples "
          f"-> log2(TPM+1)-1", flush=True)

    vals = df[sample_cols].to_numpy(dtype=float)
    df[sample_cols] = np.log2(vals + 1.0) - 1.0

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    # .csv.gz is written gzipped automatically by pandas from the extension.
    df.to_csv(output_csv, index=False)
    print(f"Wrote {output_csv} ({df.shape[0]} rows, {df.shape[1]} columns)",
          flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    convert(args.input, args.output)


if __name__ == "__main__":
    main()
