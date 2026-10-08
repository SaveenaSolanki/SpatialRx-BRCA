"""Normalisation helpers (tested) shared by the analysis modules.

These implementations are copied verbatim from the original module 05 so that
existing results are unchanged; centralising them lets us regression-test the
normalisation actually used in the pipeline.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import nnls

__all__ = ["clr_protein", "nnls_proportions", "log_cpm_gene"]


def clr_protein(prot: np.ndarray) -> np.ndarray:
    """Winsorize the top/bottom 5% per protein, then per-spot CLR.

    Parameters
    ----------
    prot : array (n_spots, n_proteins), non-negative counts

    Returns
    -------
    array of the same shape, log(x / geometric_mean(x)) per row.
    """
    prot = np.asarray(prot, dtype=np.float64).copy()
    n = prot.shape[0]
    if n == 0:
        return prot
    k = max(1, int(n * 0.05))
    for j in range(prot.shape[1]):
        x = prot[:, j]
        if x.max() == 0:
            continue
        order = np.sort(x)
        lo, hi = order[min(k, n - 1)], order[max(0, n - 1 - k)]
        prot[:, j] = np.clip(x, lo, hi)
    prot = prot + 1.0
    gmean = np.exp(np.log(prot).mean(axis=1, keepdims=True))
    return np.log(prot / gmean)


def nnls_proportions(X: np.ndarray, S: np.ndarray, lineage_names) -> "object":
    """Solve ``argmin ||X_i - S p_i||, p >= 0`` per spot; return proportions.

    ``S`` columns are scaled to comparable magnitude before solving, matching the
    original module-05 behaviour.
    """
    import pandas as pd

    S = np.asarray(S, dtype=np.float64)
    X = np.asarray(X, dtype=np.float64)
    S = S / (S.sum(axis=0, keepdims=True) + 1e-12)
    P = np.zeros((X.shape[0], S.shape[1]))
    for i in range(X.shape[0]):
        sol, _ = nnls(S, X[i])
        tot = sol.sum()
        P[i] = sol / tot if tot > 0 else 0
    return pd.DataFrame(P, columns=lineage_names)


def log_cpm_gene(counts_adata, gene: str):
    """log1p(CPM) for a single gene from a raw-count AnnData."""
    if gene not in counts_adata.var_names:
        return None
    X = counts_adata[:, gene].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    n = np.asarray(counts_adata.obs["n_counts"]).ravel().astype(np.float64)
    n[n == 0] = 1
    return np.log1p(X.ravel() / n * 1e6)
