"""Regression tests for cohort metadata and the RNA/protein feature split."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

from spatialrx import data as dm

ROOT = Path(__file__).resolve().parent.parent


def _load_legacy_srx():
    spec = importlib.util.spec_from_file_location(
        "srx_legacy", ROOT / "experiments" / "srx.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_sample_table_consistent():
    assert dm.validate_sample_table(dm.sample_meta()) == []


def test_metadata_matches_legacy_srx():
    """The package table must not silently diverge from experiments/srx.py."""
    legacy = _load_legacy_srx()
    assert dm.sample_meta().values.tolist() == legacy.SAMPLE_META.values.tolist()
    assert list(dm.sample_meta().columns) == list(legacy.SAMPLE_META.columns)
    assert dm.PAIRED == legacy.PAIRED
    assert dm.REPLICATE_SECTIONS == legacy.REPLICATE_SECTIONS


def test_metadata_shape():
    df = dm.sample_meta()
    assert len(df) == 14
    assert df["specimen"].nunique() == 14
    assert df["patient"].nunique() == 6
    assert set(df["treatment"]) == {"pre", "post"}


def test_replicates_excluded_from_paired():
    for rep in dm.REPLICATE_SECTIONS:
        assert rep not in dm.PAIRED


def test_patient_treatment_helpers():
    assert dm.patient_of("HCC22-088-P4-S2_D") == "P4"
    assert dm.treatment_of("HCC22-088-P3-S1_A") == "pre"
    assert dm.treatment_of("HCC22-088-P6-S2_D") == "post"


def test_patient_of_correct_for_all_specimens_and_matches_legacy():
    """Regression: patient_of previously used index 3 and returned the section
    token (e.g. 'S2_D') instead of the patient ('P4')."""
    legacy = _load_legacy_srx()
    for spec in dm.sample_meta()["specimen"]:
        expected = spec.split("-")[2]
        assert dm.patient_of(spec) == expected
        assert legacy.patient_of(spec) == expected
        assert dm.patient_of(spec).startswith("P")
    assert dm.patient_of("HCC22-088-P4-S2_1i_rep") == "P4"
    assert dm.patient_of("HCC22-088-P3-S1_A") == "P3"


def test_validate_sample_table_detects_duplicates():
    df = dm.sample_meta()
    bad = pd.concat([df, df.iloc[[0]]], ignore_index=True)
    problems = dm.validate_sample_table(bad)
    assert any("duplicate" in p for p in problems)


def test_validate_sample_table_detects_patient_mismatch():
    df = dm.sample_meta().copy()
    df.loc[0, "patient"] = "PX"
    problems = dm.validate_sample_table(df)
    assert any("patient mismatch" in p for p in problems)


def test_split_rna_protein(make_adata):
    adata = make_adata(n_obs=10, n_genes=6, n_ab=3)
    rna, prot = dm.split_rna_protein(adata)
    assert rna.shape == (10, 6)
    assert prot.shape == (10, 3)
    assert list(rna.var_names) == ["G0", "G1", "G2", "G3", "G4", "G5"]


def test_split_rna_protein_does_not_treat_dashed_gene_as_protein(make_adata):
    adata = make_adata(n_obs=4, n_genes=3, n_ab=1)
    # rename a gene to a dashed symbol and a duplicate-symbol suffix
    adata.var_names = ["NKX2-1", "TBCE-1", "KRTAP5-1", "CD3-1"]
    rna, prot = dm.split_rna_protein(adata)
    assert list(rna.var_names) == ["NKX2-1", "TBCE-1", "KRTAP5-1"]
    assert prot.shape == (4, 1)


def test_cohort_spot_totals(tmp_path):
    csv = tmp_path / "qc.csv"
    pd.DataFrame(
        {
            "gsm": ["a", "b"],
            "specimen": ["s1", "s2"],
            "patient": ["P1", "P1"],
            "treatment": ["pre", "post"],
            "spots": [100, 300],
            "genes": [18085, 18085],
            "proteins": [35, 35],
        }
    ).to_csv(csv, index=False)
    totals = dm.cohort_spot_totals(csv)
    assert totals["n_sections"] == 2
    assert totals["total_spots"] == 400
    assert totals["total_proteins"] == 35
