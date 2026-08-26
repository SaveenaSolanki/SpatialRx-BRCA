"""
Module 02 — Spatial foundations.
Plot expression back onto the tissue section (X/Y coordinates + H&E image).
Also prints numeric summaries (fraction of spots expressing markers) so results
are verifiable without eyeballing images.

Usage: conda run -n spatialrx python experiments/02_spatial_foundations.py
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import squidpy as sq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from srx import PROC_DIR, FIG_DIR

sc.settings.verbosity = 0
os.makedirs(FIG_DIR, exist_ok=True)

SAMPLE = os.path.join(PROC_DIR, "01_rna_foundations", "P4-S2_mod1.h5ad")
MARKERS = ["EPCAM", "ESR1", "MKI67", "CD3D", "CD68", "COL1A1", "PECAM1"]


def main():
    ad = sc.read_h5ad(SAMPLE)
    lib = ad.uns.get("library_id", "visium")

    # numeric summary
    X = ad.raw[:, MARKERS].X.toarray() if hasattr(ad.raw.X, "toarray") else \
        np.asarray(ad.raw[:, MARKERS].X)
    exp = X > 0
    print("[m2] % spots with >0 counts (P4-S2):")
    for i, g in enumerate(MARKERS):
        print(f"      {g:10s}: {exp[:, i].mean():6.1%}   mean counts {X[:, i].mean():.2f}")

    # spatial overlay panels
    fig, axes = plt.subplots(2, 4, figsize=(16, 9))
    cells = [ad, ad, ad, ad, ad, ad, ad, ad]
    colors = ["n_counts", "EPCAM", "ESR1", "MKI67", "CD3D", "CD68", "COL1A1", "PECAM1"]
    for ax, c in zip(axes.ravel(), colors):
        sq.pl.spatial_scatter(ad, color=c, ax=ax, img=True, library_id=lib,
                              frameon=False, colorbar=False,
                              title=c, size=4)
    fig.suptitle("HCC22-088-P4-S2 (post endocrine therapy) — modules 1&2: RNA on tissue",
                 y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "02_spatial_markers_P4S2.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # crop zoom on a central region to see spot-grid + tissue texture
    coords = ad.obsm["spatial"]  # pxl_row, pxl_col
    cx, cy = coords.mean(axis=0)
    r = 500
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    try:
        for ax, c in zip(axes, ["EPCAM", "CD3D", "leiden"]):
            sq.pl.spatial_scatter(ad, color=c, ax=ax, img=True,
                                  library_id=lib, frameon=False,
                                  crop_coord=[(cx - r, cx + r, cy - r, cy + r)],
                                  colorbar=False, size=6, title=f"{c} (zoomed)")
        fig.suptitle("Zoomed view — molecular domains vs H&E", y=1.0)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, "02_spatial_zoom_P4S2.png"), dpi=150,
                    bbox_inches="tight")
    except Exception as e:
        print(f"[m2] zoom figure skipped: {e}")
    plt.close(fig)

    # cluster map on tissue
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(["#d62728", "#1f77b4", "#2ca02c", "#ff7f0e",
                           "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"])
    fig, ax = plt.subplots(figsize=(7, 7))
    sq.pl.spatial_scatter(ad, color="leiden", ax=ax, img=True, library_id=lib,
                          frameon=False, size=4, palette=cmap)
    ax.set_title("Leiden clusters on tissue")
    fig.savefig(os.path.join(FIG_DIR, "02_leiden_on_tissue_P4S2.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)
    print("[m2] figures saved")


if __name__ == "__main__":
    main()