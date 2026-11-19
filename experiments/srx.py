"""
Shared constants/utilities for SpatialRx-BRCA.

Run everything with:  conda run -n spatialrx python experiments/XX_*.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# pandas>=3 writes Arrow-backed string columns, which anndata<0.11 refuses by
# default; opt in once here because every pipeline script imports this module.
# (R1 mitigation: shared entry point rather than per-script boilerplate.)
# ---------------------------------------------------------------------------
try:  # pragma: no cover - environment shim
    import anndata as _ad

    _ad.settings.allow_write_nullable_strings = True
except Exception:  # pragma: no cover
    pass

RAW_DIR = "data/raw"
EXTRACT_DIR = "data/raw/extracted"
PROC_DIR = "data/processed"
FIG_DIR = "outputs/figures"
TABLE_DIR = "outputs/tables"

# ---------------------------------------------------------------------------
# Sample metadata from GSE289326 series matrix + CITEgeist paper
# (bioRxiv 2025.02.15.638331; trial NCT05914792; P#-S1 = pre-treatment biopsy,
#  P#-S2 = post primary endocrine therapy surgery; 3-50 months of pET)
# ---------------------------------------------------------------------------
SAMPLE_META = pd.DataFrame(
    [
        # gsm, specimen, patient, treatment, section_note
        ("GSM8789203", "HCC22-088-P1-S1",   "P1", "pre",  ""),
        ("GSM8789204", "HCC22-088-P1-S2",   "P1", "post", ""),
        ("GSM8789205", "HCC22-088-P2-S1",   "P2", "pre",  ""),
        ("GSM8789206", "HCC22-088-P2-S2",   "P2", "post", ""),
        ("GSM8789207", "HCC22-088-P3-S1_A", "P3", "pre",  "S1 section A"),
        ("GSM8789208", "HCC22-088-P3-S2",   "P3", "post", ""),
        ("GSM8789209", "HCC22-088-P4-S1",   "P4", "pre",  ""),
        ("GSM8789210", "HCC22-088-P4-S2",   "P4", "post", "ESR1 D538G case (paper)"),
        ("GSM8789211", "HCC22-088-P4-S2_1i_rep", "P4", "post", "replicate section of S2"),
        ("GSM8789212", "HCC22-088-P5-S1",   "P5", "pre",  ""),
        ("GSM8789213", "HCC22-088-P5-S2",   "P5", "post", ""),
        ("GSM8789214", "HCC22-088-P5-S2_F_rep", "P5", "post", "replicate section of S2"),
        ("GSM8789215", "HCC22-088-P6-S1",   "P6", "pre",  ""),
        ("GSM8789216", "HCC22-088-P6-S2_D", "P6", "post", ""),
    ],
    columns=["gsm", "specimen", "patient", "treatment", "note"],
)
# Unique tissue sections that enter the paired analysis (replicates kept in QC,
# excluded from paired pre/post stats to avoid patient double-counting).
PAIRED = ["HCC22-088-P1-S1", "HCC22-088-P1-S2",
          "HCC22-088-P2-S1", "HCC22-088-P2-S2",
          "HCC22-088-P3-S1_A", "HCC22-088-P3-S2",
          "HCC22-088-P4-S1", "HCC22-088-P4-S2",
          "HCC22-088-P5-S1", "HCC22-088-P5-S2",
          "HCC22-088-P6-S1", "HCC22-088-P6-S2_D"]
REPLICATE_SECTIONS = ["HCC22-088-P4-S2_1i_rep", "HCC22-088-P5-S2_F_rep"]

# Response status per patient is NOT public per-patient in the paper main text;
# do not invent it. (Paper states 4 responding / 2 progressing overall.)
PATIENT_RESPONSE = {}  # filled only if found in supplement

# ---------------------------------------------------------------------------
# Canonical lineage markers (RNA) used for scoring / interpretation
# ---------------------------------------------------------------------------
MARKERS = {
    "epithelial_tumour": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7"],
    "luminal_ER": ["ESR1", "PGR", "AR", "FOXA1", "GATA3", "TFF3"],
    "basal": ["KRT5", "KRT14", "KRT17", "TP63"],
    "myoepithelial": ["ACTA2", "MYH11", "KRT5", "KRT14"],
    "proliferation": ["MKI67", "TOP2A", "PCNA", "CCNB1", "BIRC5", "UBE2C"],
    "T_cells": ["CD3D", "CD3E", "CD3G", "TRBC1", "TRBC2"],
    "CD8_T": ["CD8A", "CD8B", "GZMB", "GZMA", "PRF1"],
    "CD4_T": ["CD4", "IL7R", "CCR7", "FOXP3"],
    "B_cells": ["CD79A", "CD79B", "MS4A1", "JCHAIN", "IGHG1"],
    "myeloid": ["CD68", "LST1", "C1QA", "C1QB", "AIF1", "LYZ", "CD14"],
    "macrophage_M1": ["TNF", "IL1B", "CXCL9", "CXCL10", "NOS2"],
    "macrophage_M2": ["CD163", "MRC1", "MSR1", "TGFBI"],
    "fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB", "FAP"],
    "CAF": ["FAP", "ACTA2", "POSTN", "COL11A1"],
    "endothelial": ["PECAM1", "VWF", "EMCN", "FLT1", "CDH5"],
    "mast": ["TPSAB1", "TPSB2", "CPA3", "KIT"],
    "adipocyte": ["ADIPOQ", "LEP", "PLIN1", "FABP4"],
    "hypoxia": ["VEGFA", "SLC2A1", "LDHA", "PGK1", "CA9", "ADM"],
    "interferon": ["ISG15", "MX1", "OAS1", "IFI6", "STAT1", "IRF7"],
    "EMT": ["VIM", "FN1", "TWIST1", "SNAI2", "ZEB1", "CDH2"],
    "angiogenesis": ["VEGFA", "ANGPT1", "ANGPT2", "PDGFB", "NRP1"],
    "estrogen_response": ["ESR1", "PGR", "TFF1", "TFF3", "GREB1", "STC2", "XBP1", "CA12"],
    "cytotoxicity": ["GZMB", "GZMA", "PRF1", "NKG7", "GNLY"],
}

# Protein panel actually present in GSE289326 (10x Visium CytAssist
# Immuno-Oncology): 31 antibody features = 30 distinct targets (PTPRC has two
# isoforms, PTPRC-1/-2) + 4 isotype controls = 35 features total.
# Feature names carry a "-<n>" suffix except HLA-DRA. NOTE: reference metadata
# only; the pipeline derives symbols via spatialrx.proteins.antibody_index.
PROTEIN_LINEAGE = {
    "ACTA2-1": "Stroma/myofibroblast", "VIM-1": "Mesenchymal",
    "PECAM1-1": "Endothelial",
    "EPCAM-1": "Epithelial/tumour", "KRT5-1": "Basal epithelial/tumour",
    "CD3E-1": "T cells", "CD4-1": "CD4 T", "CD8A-1": "CD8 T",
    "CD19-1": "B cells", "MS4A1-1": "B cells", "PAX5-1": "B cells",
    "PTPRC-1": "Leukocytes (CD45)", "PTPRC-2": "Leukocytes (CD45 isoform)",
    "CR2-1": "B cells/follicular dendritic", "SDC1-1": "Plasma cells",
    "CD68-1": "Macrophage", "CD163-1": "Macrophage (M2-like)", "CD14-1": "Monocyte",
    "ITGAM-1": "Myeloid", "ITGAX-1": "Dendritic/myeloid", "HLA-DRA": "APC (MHC-II)",
    "FCGR3A-1": "NK/monocyte", "CEACAM8-1": "Neutrophil",
    "CCR7-1": "Naive T/DC", "CD27-1": "T/NK", "CXCR5-1": "B/Tfh",
    "PDCD1-1": "Exhausted T", "CD274-1": "Checkpoint (PD-L1)",
    "CD40-1": "APC/T", "BCL2-1": "Survival", "PCNA-1": "Proliferation",
}
PROTEIN_ISOTYPES = ("mouse_IgG1k", "mouse_IgG2a", "mouse_IgG2bk", "rat_IgG2a")

# Gene programs used by Module 5 (name -> genes) — curated, MSigDB-inspired
PROGRAMS = {
    "Estrogen_response": ["ESR1", "PGR", "TFF1", "TFF3", "GREB1", "STC2", "XBP1", "CA12"],
    "Proliferation": ["MKI67", "TOP2A", "PCNA", "CCNB1", "BIRC5", "UBE2C"],
    "EMT": ["VIM", "FN1", "TWIST1", "SNAI2", "ZEB1", "CDH2"],
    "Hypoxia": ["VEGFA", "SLC2A1", "LDHA", "PGK1", "CA9", "ADM"],
    "Interferon_response": ["ISG15", "MX1", "OAS1", "IFI6", "STAT1", "IRF7"],
    "T_cell_activity": ["CD3D", "CD3E", "GZMB", "GZMA", "PRF1", "NKG7"],
    "Cytotoxicity": ["GZMB", "GZMA", "PRF1", "NKG7", "GNLY"],
    "Myeloid_inflammation": ["CD68", "LST1", "C1QA", "C1QB", "LYZ", "IL1B", "TNF"],
    "Fibroblast_activation": ["FAP", "ACTA2", "POSTN", "COL11A1", "COL1A1"],
    "Angiogenesis": ["VEGFA", "ANGPT1", "ANGPT2", "PDGFB", "NRP1", "PECAM1"],
    "B_cell_activity": ["CD79A", "CD79B", "MS4A1", "JCHAIN"],
}

# Domain annotation palette (module 03) — consistent colors across figures
DOMAIN_COLORS = {
    "Tumour": "#d62728", "Tumour/immune mix": "#ff9896",
    "Immune": "#1f77b4", "Fibroblast-rich": "#2ca02c",
    "Endothelial/stroma": "#9467bd", "Adipose": "#ffd92f",
    "Low-quality/edge": "#7f7f7f",
}

TREATMENT_COLORS = {"pre": "#1f77b4", "post": "#d62728"}


def patient_of(specimen: str) -> str:
    # NOTE: previously used index 3, which returned the section token
    # (e.g. "S2_D"). The patient is the 3rd dash-separated field.
    return specimen.split("-")[2]


def treatment_of(specimen: str) -> str:
    return "pre" if "-S1" in specimen else "post"