"""The thumbnail histogram encoding is lossy on purpose; these tests pin how
lossy, and pin the degenerate cases that would otherwise divide by zero."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bhuvanfitter import (HIST_ALPHABET, HIST_BINS, HIST_LEVELS,
                          encode_histogram)


def decode(s, hist_max):
    """What the browser does: char -> level -> count."""
    return np.array([HIST_ALPHABET.index(c) for c in s]) * hist_max / HIST_LEVELS


def test_alphabet_is_64_unique_csv_safe_chars():
    assert len(HIST_ALPHABET) == 64
    assert len(set(HIST_ALPHABET)) == 64
    assert "," not in HIST_ALPHABET
    assert '"' not in HIST_ALPHABET
    assert "\n" not in HIST_ALPHABET


def test_encoding_is_always_exactly_40_chars_from_the_alphabet():
    rng = np.random.default_rng(0)
    for _ in range(50):
        data = rng.normal(size=rng.integers(10, 500))
        s, hmax = encode_histogram(data)
        assert len(s) == HIST_BINS
        assert set(s) <= set(HIST_ALPHABET)
        assert hmax > 0


def test_tallest_bin_encodes_to_the_top_level():
    rng = np.random.default_rng(1)
    s, _ = encode_histogram(rng.normal(size=400))
    assert max(HIST_ALPHABET.index(c) for c in s) == HIST_LEVELS


def test_hist_max_is_the_true_peak_count():
    data = np.concatenate([np.zeros(37), np.linspace(1, 2, 9)])
    counts, _ = np.histogram(data, bins=HIST_BINS)
    _, hmax = encode_histogram(data)
    assert hmax == counts.max()


def test_round_trip_error_is_within_one_level_of_peak():
    rng = np.random.default_rng(2)
    for _ in range(50):
        data = rng.normal(size=rng.integers(20, 1000))
        counts, _ = np.histogram(data, bins=HIST_BINS)
        s, hmax = encode_histogram(data)
        assert np.max(np.abs(decode(s, hmax) - counts)) <= hmax / HIST_LEVELS


def test_all_values_identical_gives_one_spike():
    s, hmax = encode_histogram(np.full(64, -1.0))
    assert hmax == 64
    assert len(s) == HIST_BINS
    assert sum(1 for c in s if HIST_ALPHABET.index(c) > 0) == 1


def test_single_observation():
    s, hmax = encode_histogram([3.5])
    assert hmax == 1
    assert sum(1 for c in s if HIST_ALPHABET.index(c) > 0) == 1


def test_no_finite_data_returns_empty():
    assert encode_histogram([]) == ("", 0)
    assert encode_histogram([np.nan, np.inf, -np.inf]) == ("", 0)


def test_non_finite_values_are_dropped_not_counted():
    s_clean, max_clean = encode_histogram([1.0, 2.0, 3.0, 2.0])
    s_dirty, max_dirty = encode_histogram([1.0, 2.0, np.nan, 3.0, 2.0, np.inf])
    assert (s_clean, max_clean) == (s_dirty, max_dirty)


def test_matches_the_fitter_own_histogram():
    """encode_histogram must bin identically to BhuvanFitter, or the browser
    draws bars that do not line up with the curve fitted to them."""
    from bhuvanfitter import BhuvanFitter
    rng = np.random.default_rng(3)
    data = rng.normal(size=300)
    bf = BhuvanFitter(data, gene_name="X")
    s, hmax = encode_histogram(data)
    assert hmax == int(bf.hist_counts.max())
    assert np.allclose(decode(s, hmax), bf.hist_counts, atol=hmax / HIST_LEVELS)
