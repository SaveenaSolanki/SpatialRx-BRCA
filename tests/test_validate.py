"""Tests for ``spatialrx validate-data`` (Gate 1)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from spatialrx import data as dm
from spatialrx.config import Paths
from spatialrx.validate import run_validate_data, sha256_file


def make_paths(tmp: Path) -> Paths:
    raw = tmp / "data" / "raw"
    processed = tmp / "data" / "processed"
    outputs = tmp / "outputs"
    return Paths(
        root=tmp,
        raw=raw,
        extracted=raw / "extracted",
        processed=processed,
        outputs=outputs,
        figures=outputs / "figures",
        tables=outputs / "tables",
        logs=tmp / "logs",
    )


def build_specimen_h5ad(path: Path, spec: str, patient: str, treatment: str, n_obs: int = 5):
    import anndata as ad

    n_genes, n_ab = dm.N_GENES_EXPECTED, dm.N_PROTEINS_EXPECTED
    var = pd.DataFrame(
        {"feature_types": ["Gene Expression"] * n_genes},
        index=[f"G{i}" for i in range(n_genes)],
    )
    obs = pd.DataFrame(
        {
            "sample": [spec] * n_obs,
            "patient": [patient] * n_obs,
            "treatment": [treatment] * n_obs,
        },
        index=[f"bc{i}" for i in range(n_obs)],
    )
    a = ad.AnnData(X=sparse.csr_matrix((n_obs, n_genes), dtype="float32"),
                   obs=obs, var=var)
    a.obsm["protein_counts"] = sparse.csr_matrix((n_obs, n_ab))
    a.obsm["spatial"] = np.zeros((n_obs, 2))
    path.parent.mkdir(parents=True, exist_ok=True)
    a.write(str(path))


def test_validate_non_strict_ok_with_partial_data(tmp_path):
    paths = make_paths(tmp_path)
    build_specimen_h5ad(
        paths.ingest_dir() / "HCC22-088-P1-S1.h5ad", "HCC22-088-P1-S1", "P1", "pre"
    )
    report = run_validate_data(paths, strict=False)
    assert report["ok"] is True
    by_name = {c["name"]: c for c in report["checks"]}
    assert by_name["sample_table"]["status"] == "PASS"
    assert by_name["processed:HCC22-088-P1-S1"]["status"] == "PASS"
    assert by_name["raw_tar"]["status"] == "SKIP"
    assert by_name["sample_qc"]["status"] == "SKIP"


def test_validate_strict_fails_without_raw_data(tmp_path):
    paths = make_paths(tmp_path)
    build_specimen_h5ad(
        paths.ingest_dir() / "HCC22-088-P1-S1.h5ad", "HCC22-088-P1-S1", "P1", "pre"
    )
    report = run_validate_data(paths, strict=True)
    assert report["ok"] is False
    assert report["n_fail"] > 0


def test_validate_detects_wrong_dimensions(tmp_path):
    paths = make_paths(tmp_path)
    import anndata as ad

    a = ad.AnnData(X=sparse.csr_matrix((3, 10), dtype="float32"))
    a.var["feature_types"] = ["Gene Expression"] * 10
    a.obs["sample"] = "HCC22-088-P1-S1"
    a.obs["patient"] = "P1"
    a.obs["treatment"] = "pre"
    a.obsm["protein_counts"] = sparse.csr_matrix((3, 1))
    a.obsm["spatial"] = np.zeros((3, 2))
    p = paths.ingest_dir() / "HCC22-088-P1-S1.h5ad"
    p.parent.mkdir(parents=True, exist_ok=True)
    a.write(str(p))
    report = run_validate_data(paths, strict=False)
    by_name = {c["name"]: c for c in report["checks"]}
    assert by_name["processed:HCC22-088-P1-S1"]["status"] == "FAIL"


def test_sha256_file(tmp_path):
    p = tmp_path / "x.bin"
    p.write_bytes(b"hello world")
    assert sha256_file(p) == (
        "b94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9"
    )
