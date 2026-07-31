#!/usr/bin/env python
"""
download_data.py

Download all 54 GTEx **v11 Gene-TPMs-by-tissue** matrices and convert each to a
``log2(TPM + 1) - 1`` CSV in ``../data/``. Resumable: a tissue whose converted
file already exists is skipped, and the raw ``.gct.gz`` is deleted after a
successful conversion, so ``../data/`` ends up holding **only** the 54 converted
matrices, named ``v11_log2_<tissue>.csv.gz``.

Source: https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/tpms-by-tissue/

Usage
-----
    python download_data.py                                  # all 54
    python download_data.py --tissues thyroid,whole_blood    # just these
"""
from __future__ import annotations

import argparse
import gzip
import subprocess
import sys
from pathlib import Path

from convert_gct_to_log2 import convert

LFS_POINTER_MAGIC = b"version https://git-lfs.github.com/spec/v1"


def _is_lfs_pointer(path: Path) -> bool:
    """True if `path` is a Git LFS pointer rather than the payload it stands for.

    ``data/`` is LFS-tracked and this repo is configured not to download the
    payload, so these ~134-byte text stubs are the normal on-disk state. They are
    not gzip, so `_valid_gz` would call them corrupt and this script would delete
    them and re-fetch gigabytes. Treat them as present instead."""
    try:
        with path.open("rb") as fh:
            return fh.read(len(LFS_POINTER_MAGIC)) == LFS_POINTER_MAGIC
    except OSError:
        return False


def _valid_gz(path: Path) -> bool:
    """True if `path` is a complete, readable gzip (a killed-mid-write partial
    fails here because its end-of-stream marker is missing)."""
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1 << 20):   # read to EOF; truncated files raise
                pass
        return True
    except Exception:  # noqa: BLE001
        return False

HERE = Path(__file__).resolve().parent
DATA_DIR = HERE.parent / "data"
BASE = ("https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/"
        "tpms-by-tissue/")

TISSUES = [
    "adipose_subcutaneous", "adipose_visceral_omentum", "adrenal_gland",
    "artery_aorta", "artery_coronary", "artery_tibial", "bladder",
    "brain_amygdala", "brain_anterior_cingulate_cortex_ba24",
    "brain_caudate_basal_ganglia", "brain_cerebellar_hemisphere",
    "brain_cerebellum", "brain_cortex", "brain_frontal_cortex_ba9",
    "brain_hippocampus", "brain_hypothalamus",
    "brain_nucleus_accumbens_basal_ganglia", "brain_putamen_basal_ganglia",
    "brain_spinal_cord_cervical_c-1", "brain_substantia_nigra",
    "breast_mammary_tissue", "cells_cultured_fibroblasts",
    "cells_ebv-transformed_lymphocytes", "cervix_ectocervix",
    "cervix_endocervix", "colon_sigmoid", "colon_transverse",
    "esophagus_gastroesophageal_junction", "esophagus_mucosa",
    "esophagus_muscularis", "fallopian_tube", "heart_atrial_appendage",
    "heart_left_ventricle", "kidney_cortex", "kidney_medulla", "liver", "lung",
    "minor_salivary_gland", "muscle_skeletal", "nerve_tibial", "ovary",
    "pancreas", "pituitary", "prostate", "skin_not_sun_exposed_suprapubic",
    "skin_sun_exposed_lower_leg", "small_intestine_terminal_ileum", "spleen",
    "stomach", "testis", "thyroid", "uterus", "vagina", "whole_blood",
]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tissues", default=None,
                   help="Comma-separated subset of TISSUES to fetch "
                        "(default: all 54).")
    args = p.parse_args()

    tissues = TISSUES
    if args.tissues:
        wanted = [t.strip() for t in args.tissues.split(",") if t.strip()]
        unknown = [t for t in wanted if t not in TISSUES]
        if unknown:
            raise SystemExit(f"Unknown tissue(s): {', '.join(unknown)}")
        tissues = wanted

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    done = 0
    for i, tissue in enumerate(tissues, 1):
        out = DATA_DIR / f"v11_log2_{tissue}.csv.gz"
        if out.exists() and _is_lfs_pointer(out):
            print(f"[{i}/{len(tissues)}] skip (LFS pointer): {out.name}",
                  flush=True)
            done += 1
            continue
        if out.exists() and _valid_gz(out):
            print(f"[{i}/{len(tissues)}] skip (valid): {out.name}", flush=True)
            done += 1
            continue
        if out.exists():   # present but truncated/corrupt -> discard and redo
            print(f"[{i}/{len(tissues)}] corrupt, re-fetching: {out.name}",
                  flush=True)
            out.unlink()

        fname = f"gene_tpm_adult_gtex_v11_{tissue}.gct.gz"
        url = BASE + fname
        raw = DATA_DIR / fname
        print(f"[{i}/{len(tissues)}] downloading {fname} ...", flush=True)
        rc = subprocess.call(
            ["curl", "-fsSL", "-o", str(raw), url])
        if rc != 0 or not raw.exists():
            print(f"  ERROR: download failed for {fname} (rc={rc})", flush=True)
            if raw.exists():
                raw.unlink()
            continue

        try:
            convert(raw, out)
        except Exception as exc:  # noqa: BLE001
            print(f"  ERROR: convert failed for {fname}: {exc}", flush=True)
            if out.exists():
                out.unlink()
            raw.unlink(missing_ok=True)
            continue

        raw.unlink(missing_ok=True)   # keep only the converted file
        done += 1

    print(f"Done. {done}/{len(tissues)} tissues present in {DATA_DIR}.",
          flush=True)
    if done != len(tissues):
        sys.exit(1)


if __name__ == "__main__":
    main()
