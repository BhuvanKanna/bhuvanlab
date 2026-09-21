"""Migrating an existing table to the RTI/LTI schema.

The migration renames two columns, reorders to the new schema and computes the
four new left-side columns from the stored fit parameters. Nothing is refit.

The load-bearing property is that **no existing value is re-serialised**: the
byte-identity guarantee between outputs/, gene_major/ and extract_genes.py only
holds because every carried-over field is copied as text. Parsing a float and
writing it back silently shortens it (0.012596832467784065 -> 0.012596832467784),
which is pinned below.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import migrate_rti_lti as M
from bhuvanfitter import BhuvanFitter

OLD_HEADER = ("gene,genename,y0,A,x0,w,sumsquarevalue,ti_fourparam_sigma_dist,"
              "truncationindex,min,max,mean,std,skew,kurt,right,maxheight,"
              "rightheight,n_obs,fit_success,hist,hist_max,r_squared")


def fitted_gene(seed=0):
    rng = np.random.default_rng(seed)
    bf = BhuvanFitter(rng.lognormal(0.0, 0.6, size=3000), gene_name="G1")
    bf.fit("fourparam")
    return bf


def old_row_from(bf, r2="0.912345678901234"):
    """An old-schema row carrying this fit, with a deliberately long float."""
    return ",".join(str(v) for v in [
        "ENSG00000000001.1", "G1", bf.fourparam_y0, bf.fourparam_A,
        bf.fourparam_x0, bf.fourparam_w, bf.fourparam_sumsquare,
        bf.rti_sigma_dist, bf.rti, bf.min(), bf.max(), bf.mean(), bf.std(),
        bf.skew(), bf.kurt(), bf.max(), bf.maxheight, bf.rightheight,
        3000, True, "A" * 40, 210, r2])


def write_old(tmp_path, rows, name="v11_log2_test_fourparam.csv"):
    p = tmp_path / name
    p.write_text(OLD_HEADER + "\n" + "\n".join(rows) + "\n", encoding="utf-8",
                 newline="")
    return p


def test_header_becomes_the_new_schema_with_r_squared_still_last(tmp_path):
    src = write_old(tmp_path, [old_row_from(fitted_gene())])
    out = tmp_path / "out.csv"
    M.migrate_table(src, out)
    cols = out.read_text(encoding="utf-8").splitlines()[0].split(",")
    assert cols == M.NEW_COLUMNS + ["r_squared"]
    assert "truncationindex" not in cols and "ti_fourparam_sigma_dist" not in cols


def test_carried_over_values_are_copied_as_text_not_reparsed(tmp_path):
    """A float that loses digits through float64 round-tripping must survive."""
    long_r2 = "0.012596832467784065"
    src = write_old(tmp_path, [old_row_from(fitted_gene(), r2=long_r2)])
    out = tmp_path / "out.csv"
    M.migrate_table(src, out)
    head, row = out.read_text(encoding="utf-8").splitlines()
    got = dict(zip(head.split(","), row.split(",")))
    assert got["r_squared"] == long_r2
    old = dict(zip(OLD_HEADER.split(","), src.read_text(encoding="utf-8")
                   .splitlines()[1].split(",")))
    for new_name, old_name in [("rti", "truncationindex"),
                               ("rti_sigma_dist", "ti_fourparam_sigma_dist"),
                               ("y0", "y0"), ("A", "A"), ("x0", "x0"),
                               ("w", "w"), ("maxheight", "maxheight"),
                               ("rightheight", "rightheight"), ("hist", "hist")]:
        assert got[new_name] == old[old_name], new_name


def test_new_left_columns_match_the_fitter(tmp_path):
    bf = fitted_gene()
    src = write_old(tmp_path, [old_row_from(bf)])
    out = tmp_path / "out.csv"
    M.migrate_table(src, out)
    head, row = out.read_text(encoding="utf-8").splitlines()
    got = dict(zip(head.split(","), row.split(",")))
    assert float(got["left"]) == pytest.approx(bf.min())
    assert float(got["lti"]) == pytest.approx(bf.lti, abs=1e-9)
    assert float(got["lti_sigma_dist"]) == pytest.approx(bf.lti_sigma_dist, rel=1e-9)
    assert float(got["leftheight"]) == pytest.approx(bf.leftheight, rel=1e-9)


def test_failed_fit_rows_get_blank_left_columns(tmp_path):
    """fit_success=False rows carry no metrics; the new columns must stay empty
    rather than inventing a number from NaN parameters."""
    blank = ["ENSG00000000002.1", "G2"] + [""] * 6 + [""] * 11 + ["0", "False", "", "0", ""]
    row = ",".join(blank[:8] + [""] * (len(OLD_HEADER.split(",")) - len(blank)) + blank[8:])
    row = ",".join(["ENSG00000000002.1", "G2"] + [""] * 18 + ["", "0", ""])
    src = write_old(tmp_path, [row])
    out = tmp_path / "out.csv"
    M.migrate_table(src, out)
    head, got_row = out.read_text(encoding="utf-8").splitlines()
    got = dict(zip(head.split(","), got_row.split(",")))
    for c in ("left", "lti", "lti_sigma_dist", "leftheight"):
        assert got[c] == "", f"{c} should be blank on a failed-fit row, got {got[c]!r}"


def test_table_without_hist_columns_migrates(tmp_path):
    """Only 4 of the 108 tables carry hist/hist_max. The rest must still
    migrate, keeping schema order over whatever columns they do have."""
    cols = [c for c in OLD_HEADER.split(",") if c not in ("hist", "hist_max")]
    bf = fitted_gene()
    full = dict(zip(OLD_HEADER.split(","), old_row_from(bf).split(",")))
    src = write_old(tmp_path, [",".join(full[c] for c in cols)])
    src.write_text(",".join(cols) + "\n" + ",".join(full[c] for c in cols) + "\n",
                   encoding="utf-8", newline="")
    out = tmp_path / "out.csv"
    M.migrate_table(src, out)
    got_cols = out.read_text(encoding="utf-8").splitlines()[0].split(",")
    assert "hist" not in got_cols
    assert got_cols == [c for c in M.NEW_COLUMNS if c != "hist" and c != "hist_max"] \
        + ["r_squared"]


def test_recheck_passes_on_a_migrated_table(tmp_path):
    """recheck_table re-derives every metric from the stored parameters and
    compares against what the table says, so an already-migrated table can be
    audited without the original in hand."""
    src = write_old(tmp_path, [old_row_from(fitted_gene(s)) for s in range(5)])
    out = tmp_path / "out.csv"
    M.migrate_table(src, out)
    st = M.recheck_table(out)
    assert st["rows_checked"] == 5
    assert st["max_err"] < 1e-9
    assert st["mismatches"] == 0


def test_recheck_tolerates_a_numerically_flat_curve(tmp_path):
    """A fit whose curve is constant to floating point.

    ENSG00000172188 in whole blood fits A = 0.605 but with the peak so far
    outside [min, max] that the curve varies by 8.5e-12 across it. Stored and
    recomputed maxheight agree to five significant figures, yet scaling the
    error by maxheight itself -- a floating-point zero -- reports 3.3e-06 and
    condemns the row. The scale has to be the amplitude the curve is *evaluated*
    at (|A| + |y0|), because that is where the cancellation comes from.
    """
    import numpy as np
    from bhuvanfitter import _fourparam_gaussian
    y0, A, x0, w, lo, hi = 1.0, 0.605435, 1.27, 0.146277, 2.0, 2.0575
    grid = np.linspace(lo, hi, M.CURVE_GRID)
    curve = _fourparam_gaussian(grid, y0, A, x0, w)
    base = curve.min()
    maxheight = curve.max() - base
    # tiny but not zero, exactly like the row this reproduces
    assert 0 < maxheight < 1e-9 * (abs(A) + abs(y0)), maxheight

    cols = OLD_HEADER.split(",")
    vals = dict(zip(cols, old_row_from(fitted_gene()).split(",")))
    rh = float(_fourparam_gaussian(hi, y0, A, x0, w) - base)
    vals.update(y0=repr(y0), A=repr(A), x0=repr(x0), w=repr(w),
                min=repr(lo), max=repr(hi), right=repr(hi),
                # off by one ulp of the evaluation scale, as on disk
                maxheight=repr(float(maxheight) * (1 + 3e-6)),
                rightheight=repr(rh),
                # every other metric has to describe THIS fit too, or the test
                # fails on a stale value rather than on the scale rule
                ti_fourparam_sigma_dist=repr((hi - x0) / (w / np.sqrt(2.0))),
                truncationindex=repr(min(1.0, max(0.0, rh / float(maxheight)))))
    src = write_old(tmp_path, [",".join(vals[c] for c in cols)])
    out = tmp_path / "out.csv"
    M.migrate_table(src, out)
    st = M.recheck_table(out)
    assert st["mismatches"] == 0, f"max_err={st['max_err']:.3e}"


def test_recheck_catches_a_corrupted_lti(tmp_path):
    """The check has to be able to fail, or it proves nothing."""
    src = write_old(tmp_path, [old_row_from(fitted_gene())])
    out = tmp_path / "out.csv"
    M.migrate_table(src, out)
    head, row = out.read_text(encoding="utf-8").splitlines()
    fields = row.split(",")
    fields[head.split(",").index("lti")] = "0.5"
    out.write_text(head + "\n" + ",".join(fields) + "\n", encoding="utf-8",
                   newline="")
    st = M.recheck_table(out)
    assert st["mismatches"] == 1


def test_extra_columns_keep_their_original_position(tmp_path):
    """The worm table carries `wormbasegeneid` between gene and genename. An
    extra column is not part of the schema, so it must stay where it was rather
    than being swept to the end with r_squared."""
    bf = fitted_gene()
    old_cols = OLD_HEADER.split(",")
    vals = dict(zip(old_cols, old_row_from(bf).split(",")))
    vals["wormbasegeneid"] = "WBGene00000156"
    cols = ["gene", "wormbasegeneid"] + old_cols[1:]
    src = tmp_path / "worm.csv"
    src.write_text(",".join(cols) + "\n" + ",".join(vals[c] for c in cols) + "\n",
                   encoding="utf-8", newline="")
    out = tmp_path / "out.csv"
    M.migrate_table(src, out)
    got = out.read_text(encoding="utf-8").splitlines()[0].split(",")
    assert got[:3] == ["gene", "wormbasegeneid", "genename"]
    assert got[-1] == "r_squared"


def test_verify_recomputes_maxheight_and_agrees(tmp_path):
    """--verify proves the reconstructed curve is the one that was fit."""
    src = write_old(tmp_path, [old_row_from(fitted_gene(s)) for s in range(5)])
    out = tmp_path / "out.csv"
    stats = M.migrate_table(src, out, verify=True)
    assert stats["verified"] == 5
    assert stats["max_rel_err"] < 1e-9


def test_verify_measures_error_against_the_curve_height_not_the_value(tmp_path):
    """An edge height stored as 4.7e-16 on a curve ~1e3 tall is a floating-point
    zero, not a disagreement. Scaling the error by the stored value turns that
    noise into a 100% relative error and fails a migration that is exactly right,
    so the error is measured as a fraction of the curve's own height.
    """
    bf = fitted_gene()
    fields = old_row_from(bf).split(",")
    cols = OLD_HEADER.split(",")
    # A ceiling far out in the tail: rightheight underflows to ~0 while
    # maxheight stays large.
    fields[cols.index("rightheight")] = "4.718e-16"
    src = write_old(tmp_path, [",".join(fields)])
    stats = M.migrate_table(src, tmp_path / "out.csv", verify=True)
    assert stats["max_rel_err"] < 1e-9
