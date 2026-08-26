"""
Module 01 — RNA foundations on a single sample (conceptual scRNA-seq workflow).
Learn: a spot is an observation, a gene is a feature.
Pipeline: QC -> normalize -> log1p -> HVG -> scale -> PCA -> neighbours -> UMAP -> Leiden.

Sample choice: HCC22-088-P4-S2 (post-pET surgical section, ESR1 D538G case; ~5k spots).
Also renders the smallest pre-biopsy sample (P6-S1) for contrast in the QC plot.

Usage: conda run -n spatialrx python experiments/01_rna_foundations.py
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from srx import PROC_DIR, FIG_DIR

sc.settings.verbosity = 1
sc.settings.n_jobs = 8
os.makedirs(FIG_DIR, exist_ok=True)

SAMPLE = os.path.join(PROC_DIR, "00_ingest", "HCC22-088-P4-S2.h5ad")


def qc_panel(pth: str, tag: str):
    ad = sc.read_h5ad(pth)
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.5))
    axes[0].hist(np.asarray(ad.obs["n_counts"]), bins=60, color="#1f77b4")
    axes[0].set_title(f"{tag}: UMI counts per spot")
    axes[0].set_xlabel("n_counts")
    axes[1].hist(np.asarray(ad.obs["n_genes"]), bins=60, color="#2ca02c")
    axes[1].set_title("Genes detected per spot")
    axes[1].set_xlabel("n_genes")
    axes[2].scatter(ad.obs["n_counts"], ad.obs["n_genes"], s=2, alpha=0.3)
    axes[2].set_xlabel("n_counts"); axes[2].set_ylabel("n_genes")
    axes[2].set_title("Counts vs genes")
    mt = ad.obsm["protein_counts"].sum(axis=1)
    axes[3].hist(np.asarray(mt).ravel(), bins=60, color="#9467bd")
    axes[3].set_title("Protein (antibody) UMIs per spot")
    axes[3].set_xlabel("protein UMIs")
    fig.suptitle(f"Module 1 QC — {tag} ({ad.n_obs} spots, {ad.n_vars} genes)", y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, f"01_qc_{tag}.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    return ad


def main():
    ad = qc_panel(SAMPLE, "P4-S2_post")

    # minimal spot QC (keep tissue spots; remove extreme low-depth)
    ad.obs["n_genes"] = np.asarray((ad.X > 0).sum(axis=1)).ravel()
    low = ad.obs["n_genes"] < 200
    print(f"[m1] spots below 200 genes: {low.sum()} ({low.mean():.1%})")
    ad = ad[~low].copy()
    # gene QC: expressed in >= 10 spots
    sc.pp.filter_genes(ad, min_cells=10)
    print(f"[m1] after filtering: {ad.n_obs} spots x {ad.n_vars} genes")

    # normalisation + log (HVG runs on raw counts in X before log)
    sc.pp.highly_variable_genes(ad, n_top_genes=2000, flavor="seurat_v3", n_bins=20)
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    ad.layers["lognorm"] = ad.X.copy()
    n_hvg = int(ad.var["highly_variable"].sum())
    print(f"[m1] {n_hvg} HVGs")

    fig, ax = plt.subplots(figsize=(6, 4))
    sc.pl.highly_variable_genes(ad, show=False)
    plt.title(f"HVGs (n={n_hvg})")
    fig.savefig(os.path.join(FIG_DIR, "01_hvg.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    ad.raw = ad
    ad = ad[:, ad.var["highly_variable"]].copy()
    sc.pp.scale(ad, max_value=10)
    sc.tl.pca(ad, n_comps=30, svd_solver="arpack")
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(np.arange(1, 31), ad.uns["pca"]["variance_ratio"], "o-")
    ax.set_xlabel("PC"); ax.set_ylabel("variance ratio")
    ax.set_title("PCA variance explained")
    fig.savefig(os.path.join(FIG_DIR, "01_pca_variance.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    sc.pp.neighbors(ad, n_neighbors=10, n_pcs=20)
    sc.tl.umap(ad)
    sc.tl.leiden(ad, resolution=0.8, flavor="igraph", n_iterations=2,
                 directed=False, key_added="leiden")
    print("[m1] Leiden clusters:", sorted(ad.obs["leiden"].unique().astype(int).tolist()))

    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.umap(ad, color="leiden", show=False, ax=ax, palette="tab20", size=6)
    fig.savefig(os.path.join(FIG_DIR, "01_umap_leiden.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # top markers per cluster (on raw counts via .raw)
    sc.tl.rank_genes_groups(ad, "leiden", method="t-test", use_raw=True, n_genes=25)
    top = pd_markers(ad)
    fig, ax = plt.subplots(figsize=(5, 4))
    sc.pl.rank_genes_groups(ad, n_genes=8, sharey=False, show=False)
    plt.savefig(os.path.join(FIG_DIR, "01_ranked_genes.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    outdir = os.path.join(PROC_DIR, "01_rna_foundations")
    os.makedirs(outdir, exist_ok=True)
    ad.write(os.path.join(outdir, "P4-S2_mod1.h5ad"))
    print("[m1] saved P4-S2_mod1.h5ad")


def pd_markers(ad):
    import pandas as pd
    rows = []
    for grp in ad.uns["rank_genes_groups"]["names"].dtype.names:
        names = ad.uns["rank_genes_groups"]["names"][grp][:10]
        scores = ad.uns["rank_genes_groups"]["scores"][grp][:10]
        for n, s in zip(names, scores):
            rows.append({"cluster": grp, "gene": n, "score": s})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    main()