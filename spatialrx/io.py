"""Deterministic AnnData IO.

The original pipeline failed to write h5ad under the drifted environment because
pandas 3.x produces ``StringArray`` indexes that anndata refuses to serialise by
default.  The canonical environment avoids that, but we still defensively enable
nullable-string writing when the installed anndata exposes the switch, so IO is
robust across both environments.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    import anndata as ad

__all__ = ["configure_anndata", "read_h5ad", "write_h5ad", "roundtrip_check"]


def configure_anndata() -> bool:
    """Enable nullable-string serialisation if available. Returns True if set."""
    try:
        import anndata

        if hasattr(anndata, "settings") and hasattr(
            anndata.settings, "allow_write_nullable_strings"
        ):
            anndata.settings.allow_write_nullable_strings = True
            return True
    except Exception:  # pragma: no cover - settings API varies
        pass
    return False


def read_h5ad(path):
    import anndata as ad

    configure_anndata()
    return ad.read_h5ad(str(path))


def write_h5ad(adata, path, compression: str | None = None):
    configure_anndata()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write(str(path), compression=compression)
    return path


def _allclose(a, b) -> bool:
    from scipy import sparse

    if sparse.issparse(a) or sparse.issparse(b):
        a = sparse.csr_matrix(a)
        b = sparse.csr_matrix(b)
        if a.shape != b.shape:
            return False
        d = (a - b)
        if d.nnz == 0:
            return True
        return bool(np.abs(d.data).max() <= 1e-9)
    return bool(np.allclose(np.asarray(a), np.asarray(b), equal_nan=True))


def roundtrip_check(adata) -> dict:
    """Write ``adata`` to a temp file, read it back and compare core fields."""
    configure_anndata()
    checks: dict[str, bool] = {}
    fd, tmp = tempfile.mkstemp(suffix=".h5ad")
    os.close(fd)
    try:
        adata.write(tmp)
        import anndata as ad

        back = ad.read_h5ad(tmp)
        checks["shape"] = adata.shape == back.shape
        checks["X"] = _allclose(adata.X, back.X)
        checks["obs_names"] = list(adata.obs_names) == list(back.obs_names)
        checks["var_names"] = list(adata.var_names) == list(back.var_names)
        for key in ("protein_clr", "protein_counts", "spatial"):
            if key in adata.obsm:
                checks[f"obsm[{key}]"] = key in back.obsm and _allclose(
                    adata.obsm[key], back.obsm[key]
                )
    finally:
        try:
            os.remove(tmp)
        except OSError:
            pass
    return checks
