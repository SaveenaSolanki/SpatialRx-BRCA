"""Statistical helpers with exact, unit-tested behaviour.

Centralises the paired sign-flip permutation test and BH FDR control so that
scripts and tests share one implementation.

The original analysis contained a bug: the two-sided sign-flip p-value was
multiplied by 2 even though the |statistic| comparison already counts both
tails.  ``exact_perm_paired`` here returns the correct value (``s / n_perm``).
"""

from __future__ import annotations

import numpy as np

__all__ = ["exact_perm_paired", "bh_fdr", "pairwise_permutation"]


def exact_perm_paired(values, alternative: str = "two-sided"):
    """Exact sign-flip permutation test for a paired mean.

    Parameters
    ----------
    values : array-like
        Per-patient paired differences (post - pre).
    alternative : {"two-sided", "greater", "less"}

    Returns
    -------
    (observed_mean, p_value)

    Notes
    -----
    For ``n <= 20`` the full 2**n sign-flip distribution is enumerated, so the
    test is exact.  For larger n a deterministic Monte-Carlo approximation with
    10,000 sign flips (seed 0) is used.  The two-sided p-value is the fraction of
    sign-flip means whose absolute value is at least the observed absolute mean
    (no additional factor of 2).
    """
    d = np.asarray(values, dtype=float)
    d = d[np.isfinite(d)]
    n = d.size
    if n < 3:
        obs = float(d.mean()) if n else float("nan")
        return obs, float("nan")
    obs = float(d.mean())

    if n <= 20:
        bits = np.arange(1 << n, dtype=np.uint64)
        idx = np.arange(n, dtype=np.uint64)
        signs = (((bits[:, None] >> idx) & 1).astype(np.float64) * 2.0) - 1.0
        means = signs @ d / n
        total = 1 << n
    else:
        rng = np.random.default_rng(0)
        signs = rng.choice(np.array([-1.0, 1.0]), size=(10_000, n))
        means = signs @ d / n
        total = signs.shape[0]

    tol = 1e-12
    if alternative == "two-sided":
        s = int(np.count_nonzero(np.abs(means) >= abs(obs) - tol))
    elif alternative == "greater":
        s = int(np.count_nonzero(means >= obs - tol))
    elif alternative == "less":
        s = int(np.count_nonzero(means <= obs + tol))
    else:
        raise ValueError(f"unknown alternative: {alternative!r}")

    return obs, float(min(1.0, s / total))


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR-adjusted p-values (matches statsmodels' fdr_bh).

    Handles NaN by excluding those entries from the correction and returning NaN
    for them.
    """
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    finite = np.isfinite(p)
    pf = p[finite]
    m = pf.size
    if m == 0:
        return out
    order = np.argsort(pf, kind="mergesort")
    ranked = pf[order] * m / (np.arange(m) + 1)
    # enforce monotonicity (step-up)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.clip(ranked, 0.0, 1.0)
    res = np.empty(m)
    res[order] = adj
    out[finite] = res
    return out


def pairwise_permutation(pre, post, alternative: str = "two-sided"):
    """Paired differences post - pre -> exact sign-flip test."""
    pre = np.asarray(pre, dtype=float)
    post = np.asarray(post, dtype=float)
    valid = np.isfinite(pre) & np.isfinite(post)
    return exact_perm_paired(post[valid] - pre[valid], alternative=alternative)
