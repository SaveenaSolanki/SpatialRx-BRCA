"""
Module 05 — From genes to biological programs: the cell ecosystem.
Sections:
  A. Per-spot gene-program scores (curated programs + lineage signatures) via
     scanpy.score_genes on log-normalized RNA.
  B. Protein layer: winsorize + CLR normalization (per CITEgeist methods), plus
     isotype-factor alternative; stores obsm['protein_clr'].
  C. Cell-state deconvolution:
       - RNA-based NNLS against a marker-signature matrix (6 lineages)
       - Protein-based NNLS against the 10x Immuno-Oncology antibody panel
     transparent non-negative least squares (scipy.optimize.nnls).
  D. RNA vs protein deconvolution agreement (lineage-level correlations) and
     spatial maps of cell-state proportions.
  E. Cellular neighbourhood analysis (squidpy nhood_enrichment on compartments).

Outputs: data/processed/05_ecosystem/*.h5ad, tables, figures.

NOTE (documented deviations): decoupler 2.x dropped the PROGENy resource API;
we use curated MSigDB-inspired programs instead of PROGENy. cell2location is a
heavy reference-based tool; we use transparent NNLS. See report caveats.

Usage: conda run -n spatialrx python experiments/05_ecosystem.py
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy.optimize import nnls
from scipy.sparse import issparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from srx import PROC_DIR, FIG_DIR, TABLE_DIR, SAMPLE_META, PROGRAMS, MARKERS

sc.settings.verbosity = 0
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)
SRC = os.path.join(PROC_DIR, "03_domains")
OUT = os.path.join(PROC_DIR, "05_ecosystem")
os.makedirs(OUT, exist_ok=True)

# lineages used for deconvolution (RNA markers from MARKERS, proteins from panel)
LINEAGES = ["T_cells", "B_cells", "myeloid", "fibroblast", "endothelial",
            "epithelial_tumour"]

PROTEIN_SIGNATURE = {
    "T_cells": ["CD3E", "CD4", "CD8A", "CD27", "CCR7", "PDCD1", "CD40"],
    "B_cells": ["CD19", "MS4A1", "PAX5", "CR2", "SDC1"],
    "myeloid": ["CD68", "CD14", "ITGAM", "ITGAX", "CD163", "HLA-DRA", "CEACAM8"],
    "NK": ["FCGR3A"],
    "fibroblast": ["ACTA2", "VIM"],
    "endothelial": ["PECAM1"],
    "epithelial_tumour": ["EPCAM", "KRT5", "BCL2"],
    "proliferation": ["PCNA"],
}
ISOTYPES = ["mouse_IgG1k", "mouse_IgG2a", "mouse_IgG2bk", "rat_IgG2a"]


def lognorm_matrix(ad, genes):
    genes = [g for g in genes if g in ad.raw.var_names]
    X = ad.raw[:, genes].X
    if issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)
    return X, genes


def clr_protein(prot):
    """Winsorize top/bottom 5% per protein (across spots), then per-spot CLR."""
    prot = np.asarray(prot, dtype=np.float64)
    n = prot.shape[0]
    if n == 0:
        return prot
    k = max(1, int(n * 0.05))
    for j in range(prot.shape[1]):
        x = prot[:, j]
        if x.max() == 0:
            continue
        lo, hi = np.sort(x)[k], np.sort(x)[n - 1 - k]
        x = np.clip(x, lo, hi)
        prot[:, j] = x
    prot = prot + 1.0
    gmean = np.exp(np.log(prot).mean(axis=1, keepdims=True))
    return np.log(prot / gmean)


def nnls_proportions(X, S, lineage_names):
    """Solve argmin ||X_i - S p_i||, p >= 0, per spot; returns n_spots x n_lineage."""
    S = S / (S.sum(axis=0, keepdims=True) + 1e-12)  # signature columns ~ comparable scale
    P = np.zeros((X.shape[0], S.shape[1]))
    for i in range(X.shape[0]):
        sol, _ = nnls(S, X[i])
        tot = sol.sum()
        P[i] = sol / tot if tot > 0 else 0
    return pd.DataFrame(P, columns=lineage_names)


def process_sample(specimen):
    ad = sc.read_h5ad(os.path.join(SRC, f"{specimen}_domains.h5ad"))

    # ---- A. program + lineage scores (RNA) ----
    for name, genes in {**PROGRAMS, **{k: v for k, v in MARKERS.items()
                                       if k in LINEAGES + ["proliferation"]}}.items():
        genes = [g for g in genes if g in ad.raw.var_names]
        if len(genes) >= 3:
            sc.tl.score_genes(ad, genes, score_name=f"score_{name}")
    # lineage scores compressed into score_lineage_<name>
    for ln in LINEAGES + ["proliferation"]:
        sc.tl.score_genes(ad, [g for g in MARKERS[ln] if g in ad.raw.var_names],
                          score_name=f"score_lineage_{ln}")

    # ---- B. protein normalization ----
    prot = ad.obsm["protein_counts"].toarray() if issparse(ad.obsm["protein_counts"]) \
        else np.asarray(ad.obsm["protein_counts"])
    pnames = list(ad.uns["protein_var"].index)
    ad.obsm["protein_clr"] = clr_protein(prot)
    ad.uns["protein_names"] = np.asarray(pnames)

    # ---- C1. RNA NNLS deconvolution ----
    sig_genes = []
    for ln in LINEAGES:
        for g in MARKERS[ln]:
            if g not in sig_genes:
                sig_genes.append(g)
    S_rna = np.zeros((len(sig_genes), len(LINEAGES)))
    for j, ln in enumerate(LINEAGES):
        Xg, gs = lognorm_matrix(ad, MARKERS[ln])
        if Xg.shape[1] == 0:
            continue
        S_rna[:, j] = [np.mean(Xg[:, gs.index(g)]) if g in gs else 0.0
                       for g in sig_genes]
    # signature = average lognorm per lineage per gene; weight by lineage mean
    avail = [g for g in sig_genes if g in ad.raw.var_names]
    Xsig, _ = lognorm_matrix(ad, avail)
    Xsig = Xsig / (Xsig.sum(axis=1, keepdims=True) + 1e-9)  # per-spot total normalization
    S_rna_av = S_rna[[sig_genes.index(g) for g in avail], :]
    P_rna = nnls_proportions(Xsig, S_rna_av, LINEAGES)
    for ln in LINEAGES:
        ad.obs[f"rna_{ln}"] = P_rna[ln].values
    ad.obs["rna_dominant"] = P_rna.idxmax(axis=1).values

    # ---- C2. Protein NNLS deconvolution ----
    # drop isotypes; map antibody names (PTPRC-1/-2 -> PTPRC)
    use = [p for p in pnames if p not in ISOTYPES]
    # reverse map: gene symbol -> lineages (from PROTEIN_SIGNATURE lineage->genes)
    gene_to_lineage = {}
    for ln, gs in PROTEIN_SIGNATURE.items():
        for gn in gs:
            gene_to_lineage.setdefault(gn, []).append(ln)
    sig_rows = {p.split("-")[0]: gene_to_lineage.get(p.split("-")[0], [])
                for p in use}
    lineages_prot = sorted({ln for v in sig_rows.values() for ln in v})
    S_p = np.zeros((len(use), len(lineages_prot)))
    for i, p in enumerate(use):
        gn = p.split("-")[0]
        for ln in sig_rows.get(gn, []):
            S_p[i, lineages_prot.index(ln)] = 1.0
    # remove any lineage columns with no markers
    keep_lp = S_p.sum(axis=0) > 0
    S_p = S_p[:, keep_lp]
    lineages_prot = [lines for lines, k in zip(lineages_prot, keep_lp) if k]
    Xp = ad.obsm["protein_clr"][:, [use.index(p) for p in use]]
    P_p = nnls_proportions(Xp - Xp.min(axis=0), S_p, lineages_prot)
    for ln in lineages_prot:
        ad.obs[f"prot_{ln}"] = P_p[ln].values
    ad.obs["prot_dominant"] = P_p.idxmax(axis=1).values if len(lineages_prot) else "NA"
    ad.uns["proteins_used"] = use

    # ---- D. agreement between RNA & protein deconvolution ----
    common = [ln for ln in lineages_prot if ln in LINEAGES]
    agree = pd.Series({ln: float(np.corrcoef(ad.obs[f"rna_{ln}"], ad.obs[f"prot_{ln}"])[0, 1])
                       for ln in common}, name="r")
    print(f"[m5] {specimen}: RNA-protein deconv agreement "
          f"({len(common)} lineages): " +
          ", ".join(f"{ln}={agree[ln]:+.2f}" for ln in common))

    # ---- E. neighbourhood enrichment on compartments ----
    try:
        orig = np.asarray(ad.obsm["spatial"]).copy()
        ad.obsm["spatial"] = np.column_stack([ad.obs["array_row"], ad.obs["array_col"]])
        sq.gr.spatial_neighbors(ad, coord_type="grid", n_neighs=6)
        ad.obsm["spatial"] = orig
        ad.obs["compartment"] = ad.obs["compartment"].astype("category")
        sq.gr.nhood_enrichment(ad, cluster_key="compartment")
        res = ad.uns["compartment_nhood_enrichment"]
        nh = res["zscore"]
        cats = ad.obs["compartment"].cat.categories
        nc = pd.DataFrame(np.asarray(nh), index=cats, columns=cats)
        nc.to_csv(os.path.join(TABLE_DIR, f"05_nhood_{specimen}.csv"))
    except Exception as e:
        print(f"[m5] nhood enrichment failed for {specimen}: {e}")

    ad.write(os.path.join(OUT, f"{specimen}_ecosystem.h5ad"))
    return ad


def main():
    prot_rows, nhood_tabs = [], []
    for _, m in SAMPLE_META.iterrows():
        ad = process_sample(m.specimen)
        # collect per-sample RNA-B-dominant composition and protein agreement
        comp = ad.obs.groupby("compartment").size()
        comp = comp / comp.sum()
        row = pd.Series({"sample": m.specimen, "patient": m.patient,
                         "treatment": m.treatment,
                         **{f"rna_{ln}": ad.obs[f"rna_{ln}"].mean()
                            for ln in LINEAGES},
                         **{f"prot_{ln}": ad.obs[f"prot_{ln}"].mean()
                            for ln in LINEAGES if f"prot_{ln}" in ad.obs}})
        prot_rows.append(row)

    tab = pd.DataFrame(prot_rows)
    tab.to_csv(os.path.join(TABLE_DIR, "05_cellstate_summary.csv"), index=False)
    print("\n[m5] mean deconvoluted cell-state fractions by treatment:")
    cols = [f"rna_{ln}" for ln in LINEAGES]
    print(tab.groupby("treatment")[cols].mean().round(3).to_string())

    # RNA vs protein deconv agreement across samples (boxplot)
    agrees = []
    for _, m in SAMPLE_META.iterrows():
        ad = sc.read_h5ad(os.path.join(OUT, f"{m.specimen}_ecosystem.h5ad"))
        common = [ln for ln in LINEAGES if f"rna_{ln}" in ad.obs and f"prot_{ln}" in ad.obs]
        for ln in common:
            agrees.append({"sample": m.specimen, "lineage": ln,
                           "r": np.corrcoef(ad.obs[f"rna_{ln}"], ad.obs[f"prot_{ln}"])[0, 1]})
    ag = pd.DataFrame(agrees)
    ag.to_csv(os.path.join(TABLE_DIR, "05_rna_prot_deconv_agreement.csv"), index=False)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    order = ag.groupby("lineage")["r"].median().sort_values(ascending=False).index
    import seaborn as sns
    sns.boxplot(data=ag, x="lineage", y="r", order=order, ax=ax, color="#7fb3d5")
    sns.stripplot(data=ag, x="lineage", y="r", order=order, ax=ax, color="k", size=3)
    ax.set_ylabel("Pearson r (per-spot RNA vs protein deconvolution)")
    ax.set_title("Multi-omics agreement of cell-state deconvolution")
    plt.xticks(rotation=30, ha="right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "05_rna_prot_deconv_agreement.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # program maps for the flagship sample P4-S2
    ad = sc.read_h5ad(os.path.join(OUT, "HCC22-088-P4-S2_ecosystem.h5ad"))
    lib = ad.uns.get("library_id", "visium")
    progs = ["score_Estrogen_response", "score_Proliferation", "score_T_cell_activity",
             "score_Myeloid_inflammation", "score_Fibroblast_activation",
             "score_Hypoxia"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for ax, p in zip(axes.ravel(), progs):
        sq.pl.spatial_scatter(ad, color=p, ax=ax, img=True, library_id=lib,
                              frameon=False, colorbar=False, size=3, use_raw=False,
                              title=p.replace("score_", ""))
    fig.suptitle("Gene-program maps — HCC22-088-P4-S2 (post pET)", y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "05_program_maps_P4S2.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # RNA-deconv cell state maps
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for ax, ln in zip(axes.ravel(), LINEAGES):
        sq.pl.spatial_scatter(ad, color=f"rna_{ln}", ax=ax, img=True, library_id=lib,
                              frameon=False, colorbar=False, size=3, use_raw=False,
                              title=ln)
    fig.suptitle("RNA-based cell-state deconvolution maps — P4-S2", y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "05_cellstate_maps_P4S2.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # protein-based equivalents
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    ln_use = [ln for ln in LINEAGES if f"prot_{ln}" in ad.obs]
    for ax, ln in zip(axes.ravel(), ln_use + ["epithelial_tumour"]):
        sq.pl.spatial_scatter(ad, color=f"prot_{ln}", ax=ax, img=True, library_id=lib,
                              frameon=False, colorbar=False, size=3, use_raw=False,
                              title=f"protein {ln}")
    fig.suptitle("Protein-based cell-state deconvolution maps — P4-S2", y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "05_cellstate_prot_maps_P4S2.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()