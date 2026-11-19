"""
Module 09 — CITEgeist Step-2 (spot-level QP deconvolution) for one section.

Reproducible replacement for the ad-hoc run that produced
`data/processed/09_citegeist/` (ISSUE-005). It
  1. stages a SpaceRanger-style directory from the GEO files,
  2. runs CITEgeist's protein-first QP solver with a curated panel profile dict,
  3. writes proportions + a comparison against our transparent NNLS output,
  4. writes a JSON provenance record (versions, inputs, spot counts).

Requirements (isolated env — see environment/requirements-citegeist.txt):
    pip install citegeist cuopt-cu12 --extra-index-url https://pypi.nvidia.com
The QP solver needs an NVIDIA GPU + cuOPT; without them the script exits with
status 2 and an explanatory message (it does not silently degrade).

Usage:
    python experiments/09_citegeist_step2.py [--specimen HCC22-088-P4-S2] [--min-counts 100]
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from srx import PROC_DIR, TABLE_DIR, SAMPLE_META  # noqa: E402

RAW = Path("data/raw/extracted")
STAGE = Path("data/tmp_cg")
OUT = Path(PROC_DIR) / "09_citegeist"

# Curated protein-panel profiles (markers carry the "-1" antibody suffix, as CITEgeist
# requires for real data). Panel = 10x Visium CytAssist Immuno-Oncology (31 markers).
CELL_PROFILE_DICT = {
    "Endothelial": {"Major": ["PECAM1-1"]},
    "Fibroblasts": {"Major": ["ACTA2-1"]},
    "B_Cells": {"Major": ["CD19-1"]},
    "Macrophages": {"Major": ["CD68-1", "CD163-1"]},
    "Monocytes": {"Major": ["CD14-1"]},
    "CD8_T_Cells": {"Major": ["CD3E-1", "CD8A-1"]},
    "CD4_T_Cells": {"Major": ["CD4-1", "CD3E-1"]},
    "Cancer_Luminal": {"Major": ["EPCAM-1"]},
    "Cancer_Basal": {"Major": ["KRT5-1", "SDC1-1"]},
    "Dendritic_Cells": {"Major": ["ITGAX-1"]},
    "NK": {"Major": ["FCGR3A-1"]},
    "Plasma": {"Major": ["SDC1-1"]},
}

# CITEgeist cell types -> our NNLS lineages (for the comparison table)
AGGREGATE = {
    "epithelial_tumour": ["Cancer_Luminal", "Cancer_Basal"],
    "T_cells": ["CD8_T_Cells", "CD4_T_Cells"],
    "myeloid": ["Macrophages", "Monocytes", "Dendritic_Cells"],
    "fibroblast": ["Fibroblasts"],
    "endothelial": ["Endothelial"],
    "B_cells": ["B_Cells"],
}


def stage_spaceranger(specimen: str) -> Path:
    """Build a SpaceRanger-style directory from the GEO tar contents."""
    gsm = SAMPLE_META.loc[SAMPLE_META.specimen == specimen, "gsm"].iloc[0]
    base = STAGE / f"{specimen.split(chr(45))[-2]}-{specimen.split(chr(45))[-1]}"
    spatial = base / "spatial"
    spatial.mkdir(parents=True, exist_ok=True)
    shutil.copy(RAW / f"{gsm}_{specimen}_filtered_feature_bc_matrix.h5",
                base / "filtered_feature_bc_matrix.h5")
    pos = pd.read_csv(RAW / f"{gsm}_{specimen}_tissue_positions.csv.gz", compression="gzip")
    pos.to_csv(spatial / "tissue_positions_list.csv", index=False)
    with gzip.open(RAW / f"{gsm}_{specimen}_scalefactors_json.json.gz", "rt") as fh:
        json.dump(json.load(fh), open(spatial / "scalefactors_json.json", "w"))
    for img in ("tissue_hires_image.png", "tissue_lowres_image.png"):
        with gzip.open(RAW / f"{gsm}_{specimen}_{img}.gz", "rb") as fin, \
                open(spatial / img, "wb") as fout:
            fout.write(fin.read())
    with gzip.open(RAW / f"{gsm}_{specimen}_cytassist_image.tiff.gz", "rb") as fin, \
            open(spatial / "tissue_fullres_image.tif", "wb") as fout:
        fout.write(fin.read())
    return base


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--specimen", default="HCC22-088-P4-S2")
    ap.add_argument("--min-counts", type=int, default=100,
                    help="CITEgeist filter_gex min_counts (100 biopsy, 25 surgical)")
    args = ap.parse_args()

    try:
        import squidpy as sq
        from CITEgeist.model.citegeist_model import CitegeistModel
    except ImportError as exc:  # pragma: no cover - environment guard
        print(f"[09] CITEgeist stack unavailable ({exc}).\n"
              "     Install the isolated env: environment/requirements-citegeist.txt",
              file=sys.stderr)
        return 2
    try:
        import cuopt  # noqa: F401
    except ImportError:
        print("[09] cuOPT not installed; the QP solver requires it "
              "(pip install cuopt-cu12 --extra-index-url https://pypi.nvidia.com)",
              file=sys.stderr)
        return 2

    OUT.mkdir(parents=True, exist_ok=True)
    stage = stage_spaceranger(args.specimen)
    print(f"[09] staged SpaceRanger dir: {stage}")

    adata = sq.read.visium(str(stage), counts_file="filtered_feature_bc_matrix.h5",
                           load_images=True, gex_only=False)
    print(f"[09] loaded {adata.shape[0]} spots x {adata.shape[1]} features")

    model = CitegeistModel(sample_name=args.specimen.replace("-", "_"), adata=adata,
                           output_folder=str(OUT), simulation=False)
    model.split_adata()
    model.filter_gex(nonzero_percentage=0.01, mean_expression_threshold=1.1,
                     min_counts=args.min_counts)
    model.copy_gex_to_protein_adata()
    model.preprocess_gex(target_sum=10000)
    model.preprocess_antibody()
    model.load_cell_profile_dict(CELL_PROFILE_DICT)
    model.run_cell_proportion_model(method="qp", use_detection_gating=True,
                                    use_gex_detection=True, gex_detection_k=10,
                                    validation_warn_only=True)
    model.append_proportions_to_adata(key="finetuned")

    # CITEgeist auto-saves <sample>_cell_prop_finetuned_results.csv; normalise its name
    produced = sorted(OUT.glob("*cell_prop_finetuned_results.csv"))
    if not produced:  # pragma: no cover
        print("[09] no proportions CSV produced", file=sys.stderr)
        return 3

    # ---- comparison against our NNLS deconvolution (module 05) ----
    eco_path = Path(PROC_DIR) / "05_ecosystem" / f"{args.specimen}_ecosystem.h5ad"
    comparison = None
    if eco_path.exists():
        import scanpy as sc
        from scipy import stats
        cg = pd.read_csv(produced[-1], index_col=0)
        ad = sc.read_h5ad(eco_path)
        common = cg.index.intersection(ad.obs_names)
        cg, ad = cg.loc[common], ad[common]
        rows = []
        for lineage, cols in AGGREGATE.items():
            cols = [c for c in cols if c in cg.columns]
            if not cols or f"rna_{lineage}" not in ad.obs:
                continue
            x = cg[cols].sum(axis=1)
            y = ad.obs[f"rna_{lineage}"].to_numpy()
            rows.append({"lineage": lineage, "n_spots": len(common),
                         "citegeist_mean": float(x.mean()),
                         "nnls_mean": float(np.mean(y)),
                         "pearson_r": float(stats.pearsonr(x, y).statistic),
                         "spearman_r": float(stats.spearmanr(x, y).statistic)})
        comparison = pd.DataFrame(rows)
        comparison.to_csv(Path(TABLE_DIR) / "09_citegeist_vs_nnls.csv", index=False)
        print(comparison.round(3).to_string(index=False))

    record = {
        "specimen": args.specimen,
        "min_counts": args.min_counts,
        "spots_input": int(adata.shape[0]),
        "spots_with_proportions": int(len(pd.read_csv(produced[-1], index_col=0))),
        "cell_types": list(CELL_PROFILE_DICT),
        "proportions_csv": str(produced[-1]),
        "comparison_csv": str(Path(TABLE_DIR) / "09_citegeist_vs_nnls.csv")
        if comparison is not None else None,
        "note": "QP (cuOPT barrier solver); cell-level assignment (Step 3) not run - "
                "requires StarDist nuclei segmentation weights not shipped with the package.",
    }
    (OUT / "09_run_record.json").write_text(json.dumps(record, indent=2))
    print(f"[09] wrote {produced[-1]} and run record")
    return 0


if __name__ == "__main__":
    sys.exit(main())
