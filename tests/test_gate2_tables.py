"""End-to-end regression tests tying the corrected statistics to committed
metrics tables (Gate 2).

These tests fail if the per-patient deltas or the corrected exact p-values drift,
and they encode the specific before/after values recorded in the gate report.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from spatialrx.config import PATHS
from spatialrx.stats import exact_perm_paired

METRICS = PATHS.tables / "07_paired_metrics.csv"
LR_SCORES = PATHS.tables / "10_lr_scores.csv"


@pytest.mark.skipif(not METRICS.exists(), reason="07 metrics not regenerated")
def test_bcell_paired_p_is_corrected_value():
    """B-cell activity: old stored p_perm = 0.1875 (=2*6/64); corrected = 0.09375.

    The original report prose already quoted 0.094, so the fix reconciles the
    table with the text.
    """
    t = pd.read_csv(METRICS).set_index("sample")
    piv = t.pivot_table(index="patient", columns="treatment",
                        values="score_B_cell_activity")
    d = (piv["post"] - piv["pre"]).dropna().to_numpy()
    obs, p = exact_perm_paired(d)
    assert np.isclose(obs, 0.0856055498084939)
    assert np.isclose(p, 0.09375)


@pytest.mark.skipif(not LR_SCORES.exists(), reason="10 LR scores not available")
def test_csf1_lr_paired_p_is_corrected_value():
    """CSF1-CSF1R (tumour->myeloid): committed table had p=0.0, which is
    impossible for an exact two-sided sign-flip test; committed script yields
    0.25 (doubled); corrected value is 0.125."""
    scores = pd.read_csv(LR_SCORES)
    sub = scores[
        (scores.ligand == "CSF1")
        & (scores.receptor == "CSF1R")
        & (scores.category == "Tumour->Immune(myeloid)")
    ]
    piv = sub.pivot_table(index="patient", columns="treatment",
                          values="interface_score").dropna()
    d = (piv["post"] - piv["pre"]).to_numpy()
    obs, p = exact_perm_paired(d)
    assert np.isclose(obs, 0.3325868, atol=1e-6)
    assert np.isclose(p, 0.125)
    # guard against the impossible zero
    assert p > 0


@pytest.mark.skipif(not METRICS.exists(), reason="07 metrics not regenerated")
def test_all_paired_p_values_are_multiples_of_reciprocal():
    """For n paired observations the exact p must be k / 2**n."""
    t = pd.read_csv(METRICS).set_index("sample")
    for col in t.columns:
        if col in ("patient", "treatment"):
            continue
        piv = t.pivot_table(index="patient", columns="treatment", values=col)
        if "pre" not in piv or "post" not in piv:
            continue
        d = (piv["post"] - piv["pre"]).dropna().to_numpy()
        if d.size < 3:
            continue
        _, p = exact_perm_paired(d)
        if np.isnan(p):
            continue
        assert np.isclose(p * (2 ** d.size), round(p * (2 ** d.size)))
