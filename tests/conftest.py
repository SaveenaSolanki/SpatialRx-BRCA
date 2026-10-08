"""Shared test fixtures."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
# make experiments/srx.py importable for consistency regression tests
sys.path.insert(0, str(ROOT / "experiments"))


@pytest.fixture()
def make_adata():
    """Factory returning a small AnnData with GEX + Antibody Capture features."""

    def _make(n_obs: int = 20, n_genes: int = 6, n_ab: int = 3, seed: int = 0):
        import anndata as ad
        import numpy as np
        import pandas as pd
        from scipy import sparse

        rng = np.random.default_rng(seed)
        total = n_genes + n_ab
        X = rng.poisson(2.0, size=(n_obs, total)).astype("float32")
        var = pd.DataFrame(
            {
                "feature_types": ["Gene Expression"] * n_genes
                + ["Antibody Capture"] * n_ab,
            },
            index=[f"G{i}" for i in range(n_genes)] + [f"AB{i}" for i in range(n_ab)],
        )
        obs = pd.DataFrame(index=[f"bc{i}" for i in range(n_obs)])
        adata = ad.AnnData(X=sparse.csr_matrix(X), obs=obs, var=var)
        adata.obsm["spatial"] = rng.random((n_obs, 2))
        adata.obs["sample"] = "SYNTH-S1"
        adata.obs["patient"] = "P1"
        adata.obs["treatment"] = "pre"
        return adata

    return _make
