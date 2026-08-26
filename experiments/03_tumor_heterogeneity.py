"""
Module 03 — Tumour heterogeneity across the cohort.
Per-sample: normalized pipeline -> PCA -> UMAP -> Leiden -> domain annotation.
Merged atlas: all 14 sections, harmony batch correction, shared UMAP + Leiden.

Domain annotation per cluster is *rule-based and auditable*: mean z-score of curated
lineage marker sets -> argmax domain, then broad compartment collapse.
Outputs tables + figures; each cluster's marker table is saved for verification.

Usage: conda run -n spatialrx python experiments/03_tumor_heterogeneity.py
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from srx import PROC_DIR, FIG_DIR, TABLE_DIR, SAMPLE_META, MARKERS, DOMAIN_COLORS

sc.settings.verbosity = 0
sc.settings.n_jobs = 8
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)

INGEST = os.path.join(PROC_DIR, "00_ingest")
OUT = os.path.join(PROC_DIR, "03_domains")
os.makedirs(OUT, exist_ok=True)

DOMAIN_RULES = {
    "epithelial_tumour": ["EPCAM", "KRT8", "KRT18", "KRT19", "ESR1", "PGR"],
    "proliferation": ["MKI67", "TOP2A", "BIRC5", "UBE2C"],
    "T_cells": ["CD3D", "CD3E", "CD8A", "CD4", "TRBC1", "TRBC2"],
    "myeloid": ["CD68", "LST1", "C1QA", "C1QB", "LYZ"],
    "B_cells": ["CD79A", "CD79B", "MS4A1"],
    "fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "FAP"],
    "endothelial": ["PECAM1", "VWF", "EMCN", "CDH5"],
}
COMPARTMENT = {
    "epithelial_tumour": "Tumour", "proliferation": "Tumour",
    "T_cells": "Immune", "myeloid": "Immune", "B_cells": "Immune",
    "fibroblast": "Stroma", "endothelial": "Stroma",
}
DOMAIN_PALETTE = {
    "Tumour": "#d62728", "Immune": "#1f77b4", "Stroma": "#2ca02c",
    "Mixed": "#bcbd22",
}


def cluster_zscore_table(ad):
    """Mean z-scored log-normalized expression of lineage markers per cluster."""
    genes = sorted({g for gs in DOMAIN_RULES.values() for g in gs})
    genes = [g for g in genes if g in ad.raw.var_names]
    X = ad.raw[:, genes].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.log1p(X / X.sum(axis=1, keepdims=True) * 1e4)
    df = pd.DataFrame(X, index=ad.obs_names, columns=genes)
    cl = ad.obs["leiden"].astype(str)
    mean = df.groupby(cl).mean()
    z = (mean - mean.mean()) / mean.std().replace(0, np.nan)
    return z


def annotate_domains(ad, sample):
    z = cluster_zscore_table(ad)
    rule_scores = {}
    for dom, genes in DOMAIN_RULES.items():
        gs = [g for g in genes if g in z.columns]
        rule_scores[dom] = z[gs].mean(axis=1)
    rs = pd.DataFrame(rule_scores).fillna(0)
    # compartment-level decision FIRST (Tumour/Immune/Stroma), then pick the
    # winning fine rule inside the chosen compartment. This avoids collapsing
    # proliferative-tumour clusters (epithelial ~ proliferation) to "Mixed".
    comp_score = pd.DataFrame(
        {c: rs[[d for d, cc in COMPARTMENT.items() if cc == c]].max(axis=1)
         for c in ["Tumour", "Immune", "Stroma"]})
    gap = comp_score.apply(
        lambda r: r.nlargest(2).iloc[0] - r.nlargest(2).iloc[1], axis=1)
    weak = comp_score.max(axis=1) < 0.2
    comp_choice = comp_score.idxmax(axis=1)
    comp_choice = comp_choice.mask((gap < 0.25) | weak, "Mixed")
    dom = pd.Series("mixed", index=rs.index)
    for c in ["Tumour", "Immune", "Stroma"]:
        sel = comp_choice == c
        cands = [d for d, cc in COMPARTMENT.items() if cc == c]
        dom.loc[sel] = rs.loc[sel, cands].idxmax(axis=1)
    ad.obs["domain"] = dom.loc[ad.obs["leiden"].astype(str)].to_numpy()
    ad.obs["compartment"] = ad.obs["domain"].map(COMPARTMENT).fillna("Mixed")
    out = pd.DataFrame({"sample": sample, "leiden": dom.index, **rs.to_dict("list")})
    return out


def per_sample_pipeline():
    tab = []
    for _, m in SAMPLE_META.iterrows():
        p = os.path.join(INGEST, f"{m.specimen}.h5ad")
        ad = sc.read_h5ad(p)
        ad.obs["n_genes"] = np.asarray((ad.X > 0).sum(axis=1)).ravel()
        keep = ad.obs["n_genes"] >= 100
        if keep.sum() < 50:  # sanity: keep at least 50 spots
            keep = ad.obs["n_genes"] >= 50
        ad = ad[keep].copy()
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)
        ad.layers["lognorm"] = ad.X.copy()
        ad.raw = ad.copy()
        hvg = sc.pp.highly_variable_genes(ad, n_top_genes=2000, flavor="seurat",
                                          n_bins=20, inplace=False)
        ad = ad[:, hvg["highly_variable"]].copy()
        sc.pp.scale(ad, max_value=10)
        sc.tl.pca(ad, n_comps=30, svd_solver="arpack")
        sc.pp.neighbors(ad, n_neighbors=12, n_pcs=20)
        sc.tl.umap(ad)
        sc.tl.leiden(ad, resolution=0.7, flavor="igraph", n_iterations=2,
                     directed=False, key_added="leiden")
        ztab = annotate_domains(ad, m.specimen)
        tab.append(ztab)
        ad.write(os.path.join(OUT, f"{m.specimen}_domains.h5ad"))
        print(f"[m3] {m.specimen}: {ad.n_obs} spots, {ad.obs['leiden'].nunique()} clusters, "
              f"compartments: {ad.obs['compartment'].value_counts().to_dict()}")
    ztabs = pd.concat(tab)
    ztabs.to_csv(os.path.join(TABLE_DIR, "03_cluster_domain_scores.csv"), index=False)
    return None


def merged_atlas():
    import anndata
    ads = []
    for _, m in SAMPLE_META.iterrows():
        ad = sc.read_h5ad(os.path.join(OUT, f"{m.specimen}_domains.h5ad"))
        ads.append(ad)
    ad = anndata.concat(ads, join="outer", index_unique="_", label="section")
    ad.obs = ad.obs.drop(columns=["batch"], errors="ignore")
    ad.obs["sample"] = ad.obs["section"]
    # outer join fills missing genes with NaN; convert to explicit zeros + sparse
    import scipy.sparse as sp
    ad.X = sp.csr_matrix(np.nan_to_num(ad.X))
    ad.layers.clear()
    # recompute var flags
    ad.var["hvg_seurat"] = ad.var.get("highly_variable", False)
    ad.layers["lognorm"] = ad.X.copy()
    ad.raw = ad
    hvg = sc.pp.highly_variable_genes(ad, n_top_genes=3000, flavor="seurat",
                                      n_bins=20, inplace=False)
    ad = ad[:, hvg["highly_variable"]].copy()
    sc.pp.scale(ad, max_value=10)
    sc.tl.pca(ad, n_comps=30, svd_solver="arpack")
    try:
        import harmonypy as hm
        Z = np.asarray(ad.obsm["X_pca"][:, :20], dtype=np.float64)
        meta = pd.DataFrame({"sample": ad.obs["sample"].values})
        ho = hm.run_harmony(Z, meta, "sample", max_iter_harmony=10)
        ad.obsm["X_pca_harmony"] = np.asarray(ho.Z_corr)  # already (n_obs, n_pcs)
        emb = "X_pca_harmony"
        print("[m3] harmony integration used")
    except Exception as e:
        print(f"[m3] harmony failed ({e}); using plain PCA")
        emb = "X_pca"
    sc.pp.neighbors(ad, n_neighbors=12, n_pcs=20, use_rep=emb)
    sc.tl.umap(ad)
    sc.tl.leiden(ad, resolution=0.9, flavor="igraph", n_iterations=2,
                 directed=False, key_added="leiden_atlas")
    # atlas-level annotation via same rules (on shared clusters)
    ztab = annotate_domains(ad, "atlas")
    ad.obs["domain_atlas"] = ad.obs["domain"]
    ad.obs["compartment_atlas"] = ad.obs["compartment"]
    ad.uns["domain_colors"] = DOMAIN_PALETTE
    ad.write(os.path.join(OUT, "atlas_merged.h5ad"))
    print("[m3] atlas saved:", ad.shape)

    # figure: atlas UMAP by sample/treatment/compartment
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    sc.pl.umap(ad, color="sample", ax=axes[0], show=False, frameon=False, size=5)
    axes[0].set_title("Sections (14)")
    sc.pl.umap(ad, color="treatment", ax=axes[1], show=False, frameon=False, size=5,
               palette=["#1f77b4", "#d62728"])
    axes[1].set_title("Pre vs post endocrine therapy")
    comp_colors = [DOMAIN_PALETTE.get(c, "#bcbd22") for c in
                   sorted(ad.obs["compartment"].unique())]
    sc.pl.umap(ad, color="compartment", ax=axes[2], show=False, frameon=False,
               size=5, palette=comp_colors)
    axes[2].set_title("Compartments")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "03_atlas_umap.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    return ad


def domain_composition_table(ad):
    tab = ad.obs.groupby(["sample", "patient", "treatment", "compartment"]).size()
    tab = tab.unstack(fill_value=0)
    tab = tab.div(tab.sum(axis=1), axis=0)
    tab.to_csv(os.path.join(TABLE_DIR, "03_domain_composition_fractions.csv"))
    print("\n[m3] mean compartment fraction by treatment (of analysed spots):")
    print(tab.groupby("treatment").mean().round(3).to_string())
    return tab


def main():
    per_sample_pipeline()
    ad = merged_atlas()
    domain_composition_table(ad)


if __name__ == "__main__":
    main()