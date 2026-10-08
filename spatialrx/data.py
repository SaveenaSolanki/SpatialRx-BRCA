"""Cohort metadata and data-integrity helpers for GSE289326.

This module is the single source of truth for:
  * GSM -> specimen -> patient -> treatment mapping
  * which sections are paired vs replicate
  * the raw-tar checksum
  * the RNA/protein feature split
  * lightweight validation of the sample table

It must stay importable without scanpy (feature splitting imports numpy only).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    import anndata as ad
    import pandas as pd

# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------
RAW_TAR = "GSE289326_RAW.tar"
RAW_TAR_SHA256 = "0292ba6990606e52fb8153abfd35ad662c1c266dbe7641a6c51f206b95e75843"
RAW_TAR_BYTES = 1_824_737_280

N_GENES_EXPECTED = 18_085
N_PROTEINS_EXPECTED = 35

# ---------------------------------------------------------------------------
# Sample table (verified against the GSE289326 series metadata)
# ---------------------------------------------------------------------------
_SAMPLE_ROWS = [
    ("GSM8789203", "HCC22-088-P1-S1", "P1", "pre", ""),
    ("GSM8789204", "HCC22-088-P1-S2", "P1", "post", ""),
    ("GSM8789205", "HCC22-088-P2-S1", "P2", "pre", ""),
    ("GSM8789206", "HCC22-088-P2-S2", "P2", "post", ""),
    ("GSM8789207", "HCC22-088-P3-S1_A", "P3", "pre", "S1 section A"),
    ("GSM8789208", "HCC22-088-P3-S2", "P3", "post", ""),
    ("GSM8789209", "HCC22-088-P4-S1", "P4", "pre", ""),
    ("GSM8789210", "HCC22-088-P4-S2", "P4", "post", "ESR1 D538G case (paper)"),
    ("GSM8789211", "HCC22-088-P4-S2_1i_rep", "P4", "post", "replicate section of S2"),
    ("GSM8789212", "HCC22-088-P5-S1", "P5", "pre", ""),
    ("GSM8789213", "HCC22-088-P5-S2", "P5", "post", ""),
    ("GSM8789214", "HCC22-088-P5-S2_F_rep", "P5", "post", "replicate section of S2"),
    ("GSM8789215", "HCC22-088-P6-S1", "P6", "pre", ""),
    ("GSM8789216", "HCC22-088-P6-S2_D", "P6", "post", ""),
]

_SAMPLE_COLUMNS = ["gsm", "specimen", "patient", "treatment", "note"]


def sample_meta() -> "pd.DataFrame":
    """Return the cohort table (a fresh copy each call)."""
    import pandas as pd

    return pd.DataFrame(_SAMPLE_ROWS, columns=_SAMPLE_COLUMNS)


# Unique tissue sections entering the paired analysis. Replicate sections are
# kept in QC but excluded here so a patient is never counted twice.
PAIRED = [
    "HCC22-088-P1-S1", "HCC22-088-P1-S2",
    "HCC22-088-P2-S1", "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A", "HCC22-088-P3-S2",
    "HCC22-088-P4-S1", "HCC22-088-P4-S2",
    "HCC22-088-P5-S1", "HCC22-088-P5-S2",
    "HCC22-088-P6-S1", "HCC22-088-P6-S2_D",
]
REPLICATE_SECTIONS = ["HCC22-088-P4-S2_1i_rep", "HCC22-088-P5-S2_F_rep"]

# Per-patient responder status is not published per patient; never infer it.
PATIENT_RESPONSE: dict[str, str] = {}


def patient_of(specimen: str) -> str:
    # "HCC22-088-P4-S2_D" -> "P4"  (patient is the 3rd dash-separated field)
    return specimen.split("-")[2]


def treatment_of(specimen: str) -> str:
    return "pre" if "-S1" in specimen else "post"


def paired_patients() -> list[str]:
    return sorted({patient_of(s) for s in PAIRED})


# ---------------------------------------------------------------------------
# Feature split
# ---------------------------------------------------------------------------
def split_rna_protein(adata: "ad.AnnData"):
    """Split an ingested AnnData into (rna, protein_csr).

    Protein features are selected **strictly** by ``var['feature_types']`` being
    "Antibody Capture" (case-insensitive). We deliberately do *not* use the
    "-1" name suffix: many human gene symbols contain dashes or have duplicate
    disambiguation suffixes.

    Returns
    -------
    (rna, protein_matrix) where ``protein_matrix`` is a scipy CSR matrix with
    shape (n_obs, n_proteins).
    """
    import scipy.sparse as sp

    ft = adata.var["feature_types"].astype(str)
    is_ab = ft.str.fullmatch("Antibody Capture|Protein", case=False)
    mask = is_ab.to_numpy()
    rna = adata[:, ~mask].copy()
    x = adata.X
    if sp.issparse(x):
        prot = x.tocsc()[:, mask].tocsr()
    else:
        prot = sp.csr_matrix(np.asarray(x)[:, mask])
    return rna, prot


# ---------------------------------------------------------------------------
# Sample-table validation
# ---------------------------------------------------------------------------
def validate_sample_table(df: "pd.DataFrame") -> list[str]:
    """Return a list of problems with a cohort table (empty = valid)."""
    problems: list[str] = []
    required = set(_SAMPLE_COLUMNS)
    if not required.issubset(df.columns):
        return [f"missing columns: {sorted(required - set(df.columns))}"]
    if df["gsm"].duplicated().any():
        problems.append("duplicate GSM ids")
    if df["specimen"].duplicated().any():
        problems.append("duplicate specimen ids")
    for _, r in df.iterrows():
        if patient_of(r["specimen"]) != r["patient"]:
            problems.append(f"{r['specimen']}: patient mismatch ({r['patient']})")
        if treatment_of(r["specimen"]) != r["treatment"]:
            problems.append(f"{r['specimen']}: treatment mismatch ({r['treatment']})")
    # each patient must have exactly one pre and one non-replicate post section
    for pat, sub in df.groupby("patient"):
        pre = sub[sub["treatment"] == "pre"]
        post = sub[(sub["treatment"] == "post")
                   & (~sub["specimen"].isin(REPLICATE_SECTIONS))]
        if len(pre) != 1:
            problems.append(f"{pat}: expected 1 pre section, found {len(pre)}")
        if len(post) != 1:
            problems.append(f"{pat}: expected 1 paired post section, found {len(post)}")
    return problems


def cohort_spot_totals(qc_csv) -> dict:
    """Read a module-00 QC table and return totals used by the report.

    Kept dependency-light (pandas only) so it can run in tests.
    """
    import pandas as pd

    qc = pd.read_csv(qc_csv)
    return {
        "n_sections": int(len(qc)),
        "n_patients": int(qc["patient"].nunique()),
        "total_spots": int(qc["spots"].sum()),
        "total_proteins": int(qc["proteins"].iloc[0]) if len(qc) else 0,
        "total_genes": int(qc["genes"].iloc[0]) if len(qc) else 0,
    }
