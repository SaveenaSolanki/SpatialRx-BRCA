"""Reproducibility smoke test (``spatialrx reproduce --smoke``).

Exercises the full scientific stack on a *representative* pair of processed
objects without regenerating the cohort:

  1. load a per-sample domains object and an ecosystem object
  2. build the 6-neighbour grid graph and compute Moran's I for a fixed gene
  3. score a fixed gene program
  4. run a deterministic PyTorch forward pass (bitwise repeatable)
  5. h5ad write/read round-trip

Returns a structured dict with per-check status so it can be used in CI.
"""

from __future__ import annotations

import numpy as np

from .config import PATHS, Paths

SPECIMEN = "HCC22-088-P4-S2"
FIXED_GENE = "COL1A1"
FIXED_PROGRAM = ["ESR1", "PGR", "TFF1"]
SEED = 0


def _check(status: str, detail: str = "", **extra) -> dict:
    d = {"status": status, "detail": detail}
    d.update(extra)
    return d


def run_smoke(paths: Paths | None = None, specimen: str = SPECIMEN) -> dict:
    paths = paths or PATHS
    checks: dict[str, dict] = {}

    domains = paths.domains_dir() / f"{specimen}_domains.h5ad"
    ecosystem = paths.ecosystem_dir() / f"{specimen}_ecosystem.h5ad"
    if not domains.exists() or not ecosystem.exists():
        return {
            "ok": False,
            "checks": {
                "inputs": _check(
                    "SKIP",
                    f"missing {domains.name if not domains.exists() else ecosystem.name}",
                )
            },
        }

    import scanpy as sc
    import squidpy as sq
    import torch

    from .io import read_h5ad, roundtrip_check

    # 1. load
    try:
        dom = read_h5ad(domains)
        eco = read_h5ad(ecosystem)
        checks["load"] = _check(
            "PASS", f"domains {dom.shape}, ecosystem {eco.shape}"
        )
    except Exception as e:  # noqa: BLE001
        return {"ok": False, "checks": {"load": _check("FAIL", f"{type(e).__name__}: {e}")}}

    # 2. Moran's I on a fixed gene
    try:
        orig = np.asarray(eco.obsm["spatial"]).copy()
        eco.obsm["spatial"] = np.column_stack([eco.obs["array_row"], eco.obs["array_col"]])
        sq.gr.spatial_neighbors(eco, coord_type="grid", n_neighs=6)
        eco.obsm["spatial"] = orig
        mor = sq.gr.spatial_autocorr(
            eco, mode="moran", use_raw=True, genes=[FIXED_GENE],
            n_perms=None, copy=True,
        )
        val = float(mor["I"].iloc[0])
        checks["morans_I"] = _check(
            "PASS" if np.isfinite(val) else "FAIL",
            f"I({FIXED_GENE})={val:.4f}", value=val,
        )
    except Exception as e:  # noqa: BLE001
        checks["morans_I"] = _check("FAIL", f"{type(e).__name__}: {e}")

    # 3. gene-program score
    try:
        sc.tl.score_genes(eco, FIXED_PROGRAM, score_name="_smoke_score", use_raw=True)
        mean = float(np.nanmean(np.asarray(eco.obs["_smoke_score"])))
        checks["score_genes"] = _check(
            "PASS" if np.isfinite(mean) else "FAIL", f"mean={mean:.4f}", value=mean
        )
    except Exception as e:  # noqa: BLE001
        checks["score_genes"] = _check("FAIL", f"{type(e).__name__}: {e}")

    # 4. deterministic torch forward pass
    try:
        torch.manual_seed(SEED)
        net = torch.nn.Sequential(
            torch.nn.Linear(16, 8), torch.nn.ReLU(), torch.nn.Linear(8, 16)
        )
        x = torch.randn(32, 16)
        y1 = net(x).detach().numpy()
        y2 = net(x).detach().numpy()
        deterministic = bool(np.array_equal(y1, y2))
        checks["torch"] = _check(
            "PASS" if deterministic else "FAIL",
            f"torch {torch.__version__}, deterministic={deterministic}",
        )
    except Exception as e:  # noqa: BLE001
        checks["torch"] = _check("FAIL", f"{type(e).__name__}: {e}")

    # 5. h5ad round-trip
    try:
        rt = roundtrip_check(dom)
        checks["h5ad_roundtrip"] = _check(
            "PASS" if all(rt.values()) else "FAIL", str(rt)
        )
    except Exception as e:  # noqa: BLE001
        checks["h5ad_roundtrip"] = _check("FAIL", f"{type(e).__name__}: {e}")

    return {"ok": all(c["status"] == "PASS" for c in checks.values()), "checks": checks}
