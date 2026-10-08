"""Normalisation regression tests (Gate 2).

Verifies the CLR and NNLS helpers against analytic expectations and, when the
cohort is available, against the stored ``protein_clr`` values produced by the
original pipeline.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from spatialrx.config import PATHS
from spatialrx.norm import clr_protein, log_cpm_gene, nnls_proportions


def _make_prot(n=40, p=5, seed=0):
    return np.random.default_rng(seed).poisson(5, size=(n, p)).astype(float)


def test_clr_rows_have_zero_mean():
    prot = _make_prot()
    out = clr_protein(prot)
    # CLR subtracts the per-row geometric mean, so each row sums to ~0
    assert np.allclose(out.mean(axis=1), 0.0, atol=1e-9)


def test_clr_is_deterministic():
    prot = _make_prot()
    assert np.allclose(clr_protein(prot), clr_protein(prot))


def test_clr_handles_zero_column():
    prot = _make_prot()
    prot[:, 0] = 0.0
    out = clr_protein(prot)
    assert np.isfinite(out).all()


def test_clr_empty():
    assert clr_protein(np.zeros((0, 3))).shape == (0, 3)


def test_nnls_recovers_known_mixture():
    # nnls_proportions column-normalises S, so build S with unit-sum columns;
    # then X = p_true @ S.T is recovered exactly by the solver.
    rng = np.random.default_rng(0)
    S = rng.random((8, 4))
    S = S / S.sum(axis=0, keepdims=True)
    p_true = rng.dirichlet(np.ones(4), size=30)
    X = p_true @ S.T
    P = nnls_proportions(X, S, [f"L{i}" for i in range(4)])
    assert np.allclose(P.to_numpy(), p_true, atol=1e-6)
    assert np.allclose(P.sum(axis=1), 1.0, atol=1e-9)


def test_nnls_proportions_sum_to_one():
    rng = np.random.default_rng(1)
    S = rng.random((6, 3))
    X = rng.random((10, 6))
    P = nnls_proportions(X, S, ["a", "b", "c"])
    assert np.allclose(P.sum(axis=1), 1.0, atol=1e-9)
    assert (P.to_numpy() >= 0).all()


def test_log_cpm_gene_known_value():
    import anndata as ad
    import pandas as pd
    from scipy import sparse

    counts = np.array([[10], [0], [100]], dtype=float)
    a = ad.AnnData(X=sparse.csr_matrix(counts), var=pd.DataFrame(index=["G"]))
    a.obs["n_counts"] = counts.sum(axis=1)
    got = log_cpm_gene(a, "G")
    expected = np.log1p(np.array([10, 0, 100]) / np.array([10, 1, 100]) * 1e6)
    assert np.allclose(got, expected)


@pytest.mark.skipif(
    not (PATHS.ecosystem_dir() / "HCC22-088-P4-S2_ecosystem.h5ad").exists(),
    reason="cohort processed objects not available",
)
def test_clr_matches_stored_pipeline_values():
    """Regression: recomputing CLR from raw protein counts must reproduce the
    ``obsm['protein_clr']`` stored by module 05."""
    from spatialrx.io import read_h5ad

    adata = read_h5ad(PATHS.ecosystem_dir() / "HCC22-088-P4-S2_ecosystem.h5ad")
    raw = adata.obsm["protein_counts"]
    raw = raw.toarray() if hasattr(raw, "toarray") else np.asarray(raw)
    recomputed = clr_protein(raw)
    stored = np.asarray(adata.obsm["protein_clr"])
    assert recomputed.shape == stored.shape
    assert np.allclose(recomputed, stored, atol=1e-8)
