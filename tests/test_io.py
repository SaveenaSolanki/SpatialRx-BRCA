"""Regression tests for deterministic AnnData IO (Gate 1)."""

from __future__ import annotations

import numpy as np
from scipy import sparse

from spatialrx.io import configure_anndata, read_h5ad, roundtrip_check, write_h5ad


def test_configure_anndata_returns_bool():
    assert isinstance(configure_anndata(), bool)


def test_write_read_roundtrip_dense(tmp_path, make_adata):
    adata = make_adata(n_obs=8, n_genes=5, n_ab=2)
    adata.X = np.asarray(adata.X.todense())
    p = tmp_path / "dense.h5ad"
    write_h5ad(adata, p)
    back = read_h5ad(p)
    assert back.shape == adata.shape
    assert np.allclose(back.X, adata.X)
    assert list(back.obs_names) == list(adata.obs_names)


def test_roundtrip_check_sparse_and_obsm(tmp_path, make_adata):
    adata = make_adata(n_obs=12, n_genes=6, n_ab=3)
    adata.obsm["protein_counts"] = sparse.csr_matrix(
        np.random.default_rng(1).poisson(1, size=(12, 3))
    )
    adata.obsm["protein_clr"] = np.random.default_rng(2).normal(size=(12, 3))
    checks = roundtrip_check(adata)
    assert checks, "expected at least one check"
    assert all(checks.values()), checks


def test_write_creates_parent_dirs(tmp_path, make_adata):
    adata = make_adata()
    p = tmp_path / "a" / "b" / "c.h5ad"
    write_h5ad(adata, p)
    assert p.exists()
