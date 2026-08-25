"""The R^2 side-car is read by fixed offset, so its geometry is load-bearing.

`docs/r2/*.txt` carries one `%6.3f` field per gene in genes.tsv order, and the
browser reads gene `gi` as `text.substr(gi * width, width)`. Nothing in that
lookup can detect a field that is one character wide or one gene out of step --
it just returns a number belonging to a different gene. These tests pin the
geometry, the round-trip, and the manifest contract the page relies on.

The published files are only checked when they exist, so a fresh clone that has
not run build_r2.py still passes.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import build_r2 as br

REPO = Path(__file__).resolve().parents[3]
DOCS = REPO / "docs"
DOCS_R2 = DOCS / "r2"
MANIFEST = DOCS / "manifest.json"


def _published():
    return sorted(DOCS_R2.glob("v11_log2_*_r2*.txt")) if DOCS_R2.exists() else []


def test_field_width_holds_the_clamp():
    """Every value the encoder can emit must fit the field exactly."""
    for v in (br.CLAMP, -br.CLAMP, 0.0, -0.0005, 0.9995):
        assert len(f"{v:{br.WIDTH}.{br.DECIMALS}f}") == br.WIDTH
    assert len(br.MISSING) == br.WIDTH


def test_encode_round_trips(tmp_path):
    """Decoding by fixed offset returns each gene's own value, and a gene with
    no R^2 decodes to blank rather than to some neighbouring number."""
    order = ["ENSG001", "ENSG002", "ENSG003", "ENSG004"]
    csv = tmp_path / "v11_log2_fake_r2_excluded_at_or_below_-1.csv"
    pd.DataFrame({
        # ENSG003 is absent entirely; ENSG002 is present but NaN.
        "gene": ["ENSG001.5", "ENSG002.1", "ENSG004.2"],
        "r_squared": [0.5987, np.nan, -0.0123],
    }).to_csv(csv, index=False)

    text, counts = br.encode(csv, order, strict=False)

    assert len(text) == br.WIDTH * len(order)
    assert counts == {"values": 2, "missing": 2, "clamped": 0}

    fields = [text[i * br.WIDTH:(i + 1) * br.WIDTH] for i in range(len(order))]
    assert float(fields[0]) == pytest.approx(0.599, abs=5e-4)
    assert fields[1].strip() == ""
    assert fields[2].strip() == ""
    assert float(fields[3]) == pytest.approx(-0.012, abs=5e-4)


def test_encode_clamps_or_refuses(tmp_path):
    """An R^2 too large for the field must never be written narrow -- one wide
    field would shift every gene after it."""
    order = ["ENSG001"]
    csv = tmp_path / "v11_log2_fake_r2_excluded_at_or_below_-1.csv"
    pd.DataFrame({"gene": ["ENSG001.1"], "r_squared": [-1e6]}).to_csv(csv, index=False)

    text, counts = br.encode(csv, order, strict=False)
    assert len(text) == br.WIDTH
    assert counts["clamped"] == 1
    assert float(text) == pytest.approx(-br.CLAMP)

    with pytest.raises(SystemExit):
        br.encode(csv, order, strict=True)


def test_version_mismatch_does_not_silently_blank(tmp_path):
    """genes.tsv ids are unversioned; matching on the versioned string is the
    easiest way to produce an all-blank file, so the stripping must hold."""
    order = ["ENSG001"]
    csv = tmp_path / "v11_log2_fake_r2_excluded_at_or_below_-1.csv"
    pd.DataFrame({"gene": ["ENSG001.99"], "r_squared": [0.42]}).to_csv(csv, index=False)
    text, counts = br.encode(csv, order, strict=False)
    assert counts["values"] == 1
    assert float(text) == pytest.approx(0.42, abs=5e-4)


@pytest.mark.skipif(not (DOCS / "genes.tsv").exists(), reason="no genes.tsv")
def test_published_files_are_gene_aligned():
    files = _published()
    if not files:
        pytest.skip("build_r2.py has not been run")
    n_genes = len(br.load_gene_order())
    for path in files:
        text = path.read_text(encoding="ascii")
        assert len(text) == br.WIDTH * n_genes, (
            f"{path.name} is {len(text)} chars, expected "
            f"{br.WIDTH * n_genes} for {n_genes:,} genes")


@pytest.mark.skipif(not MANIFEST.exists(), reason="no manifest.json")
def test_manifest_r2_block_matches_the_files():
    files = _published()
    if not files:
        pytest.skip("build_r2.py has not been run")
    man = json.loads(MANIFEST.read_text(encoding="utf-8"))
    r2 = man.get("r2")
    assert r2 and r2.get("available") is True
    # The page slices by manifest width rather than a literal, so a drift here
    # silently misreads every field.
    assert r2["width"] == br.WIDTH
    assert r2["decimals"] == br.DECIMALS
    assert r2["kind_ids"]["excluded"] == "excluded_at_or_below_-1"

    named = {p.name for p in files}
    for tissue, kinds in r2["files"].items():
        for kind_id, name in kinds.items():
            assert name in named, f"{tissue}/{kind_id} -> missing {name}"
