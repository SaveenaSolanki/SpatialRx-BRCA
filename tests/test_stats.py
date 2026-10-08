"""Regression tests for the exact permutation test and BH FDR (Gate 2).

The original ``exact_perm_paired`` doubled the two-sided p-value.  These tests
pin the correct behaviour against independent brute-force enumeration and against
statsmodels.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from spatialrx.stats import bh_fdr, exact_perm_paired, pairwise_permutation


def _brute_two_sided(d):
    d = np.asarray(d, float)
    obs = abs(d.mean())
    vals = [
        abs(np.mean(np.array(signs) * d))
        for signs in itertools.product([-1, 1], repeat=len(d))
    ]
    return float(np.mean(np.array(vals) >= obs - 1e-12))


def test_all_same_sign_matches_exact_two_sided():
    obs, p = exact_perm_paired([1, 1, 1, 1, 1, 1])
    assert obs == 1.0
    # only the all-+ and all- sign vectors reach |mean| = 1
    assert np.isclose(p, 2 / 64)


def test_monotone_values_analytical():
    obs, p = exact_perm_paired([1, 2, 3, 4, 5, 6])
    assert np.isclose(obs, 3.5)
    assert np.isclose(p, 2 / 64)


def test_not_doubled_against_brute_force():
    d = [0.09, 0.02, -0.01, 0.15, 0.05, 0.10]
    _, p = exact_perm_paired(d)
    expected = _brute_two_sided(d)
    assert np.isclose(p, expected), f"got {p}, brute force {expected}"
    # guard the specific bug: the doubled value would be 2x expected
    if expected < 0.5:
        assert not np.isclose(p, 2 * expected)


def test_symmetric_zero_mean_p_is_one():
    _, p = exact_perm_paired([1, -1, 1, -1, 1, -1])
    assert np.isclose(p, 1.0)


def test_one_sided_alternatives_are_ordered():
    d = [1, 1, 1, 1, 1, -1]
    _, p_greater = exact_perm_paired(d, alternative="greater")
    _, p_less = exact_perm_paired(d, alternative="less")
    assert p_greater < p_less
    assert 0 < p_greater < 1 and 0 < p_less < 1


def test_doubling_bug_is_absent():
    """Regression: the old code returned 2*s/N. For an extreme all-positive
    vector the correct exact two-sided p is 2/64; the buggy value was 4/64."""
    _, p = exact_perm_paired([1, 1, 1, 1, 1, 1])
    assert np.isclose(p, 2 / 64)
    assert not np.isclose(p, 4 / 64)


def test_too_few_observations():
    _, p = exact_perm_paired([1.0])
    assert np.isnan(p)
    _, p2 = exact_perm_paired([])
    assert np.isnan(p2)


def test_pairwise_permutation_ignores_nan():
    obs, p = pairwise_permutation([0, 1, 2, np.nan], [1, 2, 3, 4])
    assert np.isfinite(obs) and np.isfinite(p)


# ---------------------------------------------------------------------------
# BH FDR
# ---------------------------------------------------------------------------
def test_bh_matches_statsmodels_classic_example():
    sm = pytest.importorskip("statsmodels")
    from statsmodels.stats.multitest import multipletests

    p = np.array([
        0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216,
        0.222, 0.251, 0.269, 0.275, 0.34, 0.341, 0.384, 0.569, 0.594, 0.696,
        0.762, 0.94, 0.942, 0.975,
    ])
    assert np.allclose(bh_fdr(p), multipletests(p, method="fdr_bh")[1])


def test_bh_small_example():
    adj = bh_fdr([0.01, 0.02, 0.03, 0.04, 0.05])
    assert np.allclose(adj, 0.05)


def test_bh_monotone_non_decreasing():
    p = np.array([0.04, 0.01, 0.03, 0.02, 0.5])
    adj = bh_fdr(p)
    assert np.all(np.diff(adj[np.argsort(p)]) >= -1e-12)


def test_bh_bounds_and_nan():
    adj = bh_fdr([0.01, np.nan, 0.04])
    assert np.isnan(adj[1])
    assert 0 <= adj[0] <= 1 and 0 <= adj[2] <= 1
    assert len(bh_fdr([])) == 0
