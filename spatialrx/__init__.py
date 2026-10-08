"""SpatialRx-BRCA reproducibility & analysis toolkit.

Importing this package is deliberately cheap (no scanpy/squidpy/torch imports).
Heavy scientific dependencies are imported lazily inside the individual modules
that need them.
"""

from __future__ import annotations

__version__ = "0.2.0"

from .config import PATHS, Paths  # noqa: F401

__all__ = ["PATHS", "Paths", "__version__"]
