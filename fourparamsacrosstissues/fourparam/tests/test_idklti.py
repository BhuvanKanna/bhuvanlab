"""idklti: the old LTI, forced to 0 where the gene's left edge is the floor.

Computed twice -- by ``compute_idklti.py`` into ``idklti/``, and by the browser
from the loaded row -- so the two cut-offs are pinned together here.
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import compute_idklti as C

INDEX_HTML = Path(__file__).resolve().parents[3] / "docs" / "index.html"


def test_rule_on_text_keeps_kept_values_verbatim():
    lti = pd.Series(["0.012596832467784065", "0.5", "", "0.3", "0.0"])
    mn = pd.Series(["-0.9", "-0.75", "-0.95", "-0.7", "2.0"])
    got = C.idklti_text(lti, mn).tolist()
    assert got == ["0", "0", "", "0.3", "0.0"]     # -0.75 itself is at the floor


def test_numeric_and_text_rules_agree():
    rng = np.random.default_rng(0)
    lti = rng.uniform(0, 1, 500)
    lti[::7] = np.nan
    mn = rng.uniform(-1, 1, 500)
    num = C.idklti(lti, mn)
    txt = pd.to_numeric(C.idklti_text(pd.Series(lti).map(
        lambda v: "" if np.isnan(v) else repr(float(v))), pd.Series(mn).map(repr)),
        errors="coerce").to_numpy()
    assert np.allclose(num, txt, equal_nan=True)


def test_browser_uses_the_same_cutoff():
    html = INDEX_HTML.read_text(encoding="utf-8")
    m = re.search(r"const IDK_FLOOR = (-?[0-9.]+);", html)
    assert m, "IDK_FLOOR not found in docs/index.html"
    assert float(m.group(1)) == C.FLOOR_CUT
