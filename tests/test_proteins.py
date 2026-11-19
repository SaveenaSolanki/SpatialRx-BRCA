"""Regression tests for antibody-name handling (ISSUE-001)."""

from __future__ import annotations

import pytest

from spatialrx.proteins import antibody_index, isotype, protein_symbols, protein_to_symbol


@pytest.mark.parametrize(
    "feature,symbol",
    [
        ("CD68-1", "CD68"),
        ("PTPRC-2", "PTPRC"),
        ("CD3E-1", "CD3E"),
        ("HLA-DRA", "HLA-DRA"),   # dashed symbol, no panel suffix -> must not become "HLA"
        ("mouse_IgG2a", "mouse_IgG2a"),
    ],
)
def test_protein_to_symbol(feature, symbol):
    assert protein_to_symbol(feature) == symbol


def test_ambiguous_dashed_symbol_resolved_by_known_symbols():
    """NKX2-1 is a real symbol that looks like a -1 suffix.

    Default (antibody-only) behaviour strips the suffix; supplying the gene-symbol
    universe preserves the symbol. Documented, not silent.
    """
    assert protein_to_symbol("NKX2-1") == "NKX2"
    assert protein_to_symbol("NKX2-1", known_symbols={"NKX2-1", "CD68"}) == "NKX2-1"
    # a suffixed panel name is still stripped even when the universe is supplied
    assert protein_to_symbol("CD68-1", known_symbols={"NKX2-1", "CD68"}) == "CD68"


def test_hla_dra_is_not_truncated():
    """The ISSUE-001 regression: split('-')[0] gave 'HLA', dropping the pair."""
    assert protein_to_symbol("HLA-DRA") != "HLA"


def test_antibody_index_prefers_first_isoform():
    names = ["CD68-1", "PTPRC-1", "PTPRC-2", "HLA-DRA", "mouse_IgG1k", "CD3E-1"]
    idx = antibody_index(names)
    assert idx["PTPRC"] == 1          # first isoform
    assert idx["HLA-DRA"] == 3
    assert idx["CD3E"] == 5
    assert "mouse_IgG1k" not in idx   # isotypes excluded
    assert "HLA" not in idx


def test_antibody_index_covers_real_panel():
    """All 31 non-isotype panel features must resolve to a symbol."""
    panel = [
        "ACTA2-1", "BCL2-1", "CCR7-1", "CD14-1", "CD163-1", "CD19-1", "CD27-1",
        "CD274-1", "CD3E-1", "CD4-1", "CD40-1", "CD68-1", "CD8A-1", "CEACAM8-1",
        "CR2-1", "CXCR5-1", "EPCAM-1", "FCGR3A-1", "HLA-DRA", "ITGAM-1", "ITGAX-1",
        "KRT5-1", "MS4A1-1", "PAX5-1", "PCNA-1", "PDCD1-1", "PECAM1-1", "PTPRC-1",
        "PTPRC-2", "SDC1-1", "VIM-1",
    ]
    idx = antibody_index(panel)
    assert len(idx) == 30            # 31 features - 1 duplicate isoform
    for expected in ("HLA-DRA", "EPCAM", "PTPRC", "VIM", "ACTA2"):
        assert expected in idx


def test_isotype():
    assert isotype("rat_IgG2a")
    assert not isotype("CD68-1")


def test_protein_symbols_vectorised():
    assert protein_symbols(["CD68-1", "HLA-DRA"]) == ["CD68", "HLA-DRA"]
