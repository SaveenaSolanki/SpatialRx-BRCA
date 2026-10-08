"""Regression tests for merged-atlas integrity (Gate 3)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from spatialrx.atlas import atlas_compartment_counts, validate_atlas
from spatialrx.config import PATHS


def _synthetic_atlas(sample_ids, patients=None, treatments=None, compartments=None, n=6):
    import anndata as ad

    obs = pd.DataFrame({"sample": sample_ids[:n]})
    obs["patient"] = patients[:n] if patients is not None else "P1"
    obs["treatment"] = treatments[:n] if treatments is not None else "pre"
    obs["compartment"] = compartments[:n] if compartments is not None else "Stroma"
    return ad.AnnData(X=np.zeros((n, 3)), obs=obs)


def test_validate_atlas_flags_integer_sample_ids():
    a = _synthetic_atlas([str(i % 14) for i in range(6)])
    problems = validate_atlas(a)
    assert any("integer batch labels" in p for p in problems)
    assert any("sample IDs mismatch" in p for p in problems)


def test_validate_atlas_flags_patient_mismatch():
    # real specimen name but wrong patient
    a = _synthetic_atlas(
        ["HCC22-088-P4-S2"] * 6, patients=["P1"] * 6, treatments=["post"] * 6
    )
    problems = validate_atlas(a)
    assert any("patient mismatch" in p for p in problems)


def test_validate_atlas_flags_unknown_compartment():
    a = _synthetic_atlas(
        ["HCC22-088-P4-S2"] * 6, patients=["P4"] * 6, treatments=["post"] * 6,
        compartments=["Weird"] * 6,
    )
    problems = validate_atlas(a)
    assert any("unknown compartments" in p for p in problems)


@pytest.mark.skipif(
    not (PATHS.domains_dir() / "atlas_merged.h5ad").exists(),
    reason="atlas not available",
)
def test_regenerated_atlas_is_valid():
    from spatialrx.io import read_h5ad

    adata = read_h5ad(PATHS.domains_dir() / "atlas_merged.h5ad")
    assert validate_atlas(adata) == []
    assert adata.obs["sample"].nunique() == 14
    counts = atlas_compartment_counts(adata)
    assert counts, "atlas has no compartment annotation"
    # the corrected pipeline recovers a Tumour compartment (the broken atlas had none)
    assert "Tumour" in counts
    assert set(counts) <= {"Tumour", "Immune", "Stroma", "Mixed"}


def test_domain_composition_table_is_valid():
    path = PATHS.tables / "03_domain_composition_fractions.csv"
    if not path.exists():
        pytest.skip("composition table not regenerated")
    t = pd.read_csv(path).set_index("sample")
    assert t.index.nunique() == 14
    frac = [c for c in ["Tumour", "Immune", "Stroma", "Mixed"] if c in t.columns]
    assert frac, "no compartment columns"
    assert np.allclose(t[frac].sum(axis=1), 1.0, atol=1e-9)
