"""idk_flr_idx: the old flr_idx, forced to 0 where the gene's left edge is the floor.

Computed twice -- by ``compute_idk_flr_idx.py`` into ``idk_flr_idx/``, and by the browser
from the loaded row -- so the two cut-offs are pinned together here.
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import compute_idk_flr_idx as C

INDEX_HTML = Path(__file__).resolve().parents[3] / "docs" / "index.html"


def test_rule_on_text_keeps_kept_values_verbatim():
    flr_idx = pd.Series(["0.012596832467784065", "0.5", "", "0.3", "0.0"])
    mn = pd.Series(["-0.9", "-0.75", "-0.95", "-0.7", "2.0"])
    got = C.idk_flr_idx_text(flr_idx, mn).tolist()
    assert got == ["0", "0", "", "0.3", "0.0"]     # -0.75 itself is at the floor


def test_numeric_and_text_rules_agree():
    rng = np.random.default_rng(0)
    flr_idx = rng.uniform(0, 1, 500)
    flr_idx[::7] = np.nan
    mn = rng.uniform(-1, 1, 500)
    num = C.idk_flr_idx(flr_idx, mn)
    txt = pd.to_numeric(C.idk_flr_idx_text(pd.Series(flr_idx).map(
        lambda v: "" if np.isnan(v) else repr(float(v))), pd.Series(mn).map(repr)),
        errors="coerce").to_numpy()
    assert np.allclose(num, txt, equal_nan=True)


def test_browser_uses_the_same_cutoff():
    html = INDEX_HTML.read_text(encoding="utf-8")
    m = re.search(r"const IDK_FLOOR = (-?[0-9.]+);", html)
    assert m, "IDK_FLOOR not found in docs/index.html"
    assert float(m.group(1)) == C.FLOOR_CUT


def test_worm_sidecar_matches_rule():
    out = C.OUT_DIR / "worm_idk_flr_idx_excluded_at_or_below_-1.csv"
    src = pd.read_csv(C.WORM_TABLE, usecols=["gene", "min", "flr_idx"], dtype=str,
                      keep_default_na=False)
    got = pd.read_csv(out, dtype=str, keep_default_na=False)
    assert got["gene"].tolist() == src["gene"].tolist()
    assert got["idk_flr_idx"].tolist() == C.idk_flr_idx_text(src["flr_idx"], src["min"]).tolist()
