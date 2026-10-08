"""Integrity validation for the merged atlas.

The original atlas stored integer batch indices as sample IDs and contained no
Tumour compartment because the annotation ran on scaled data.  ``validate_atlas``
makes those failure modes explicit and testable.
"""

from __future__ import annotations

from .data import sample_meta

ALLOWED_COMPARTMENTS = {"Tumour", "Immune", "Stroma", "Mixed"}


def validate_atlas(adata, meta=None) -> list[str]:
    """Return a list of integrity problems with a merged atlas (empty = valid)."""
    problems: list[str] = []
    meta = sample_meta() if meta is None else meta
    expected = set(meta["specimen"].astype(str))

    if "sample" not in adata.obs:
        return ["missing obs['sample']"]
    got = set(adata.obs["sample"].astype(str).unique())
    if got != expected:
        problems.append(
            "sample IDs mismatch: "
            f"missing={sorted(expected - got)[:5]} extra={sorted(got - expected)[:5]}"
        )
        if got and all(s.isdigit() for s in got):
            problems.append("sample IDs are integer batch labels (concat keys missing)")

    if "patient" not in adata.obs:
        problems.append("missing obs['patient']")
    else:
        pmap = dict(zip(meta["specimen"].astype(str), meta["patient"].astype(str)))
        tmap = dict(zip(meta["specimen"].astype(str), meta["treatment"].astype(str)))
        sample = adata.obs["sample"].astype(str)
        bad_p = sample.map(pmap) != adata.obs["patient"].astype(str)
        if bad_p.any():
            problems.append(f"{int(bad_p.sum())} rows with patient mismatch")
        if "treatment" in adata.obs:
            bad_t = sample.map(tmap) != adata.obs["treatment"].astype(str)
            if bad_t.any():
                problems.append(f"{int(bad_t.sum())} rows with treatment mismatch")

    if "compartment" not in adata.obs:
        problems.append("missing obs['compartment']")
    else:
        vals = set(adata.obs["compartment"].astype(str).unique())
        extra = vals - ALLOWED_COMPARTMENTS
        if extra:
            problems.append(f"unknown compartments: {sorted(extra)}")
    return problems


def atlas_compartment_counts(adata) -> dict:
    if "compartment" not in adata.obs:
        return {}
    return {
        str(k): int(v)
        for k, v in adata.obs["compartment"].value_counts().items()
    }
