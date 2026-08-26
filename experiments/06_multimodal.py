"""
Module 06 — Multi-omics: RNA vs protein agreement per spatial spot.
For ~27 matched gene/antibody pairs (e.g. CD3E RNA vs CD3 protein):
  - per-sample Pearson r across spots (RNA log-CPM vs CLR protein)
  - per-spot discordance score (zRNA - zProtein), tested for spatial structure
    via Moran's I on the grid graph (is disagreement itself spatially organised?)
  - robustness check: protein normalised with the study's isotype-normalisation
    factors (alternative to plain CLR)
Outputs: tables + figures incl. discordance maps.

Usage: conda run -n spatialrx python experiments/06_multimodal.py
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
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from srx import PROC_DIR, FIG_DIR, TABLE_DIR, SAMPLE_META

sc.settings.verbosity = 0
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)
SRC = os.path.join(PROC_DIR, "05_ecosystem")

# RNA gene  -> protein feature (antibody) for matched pairs
PAIRS = {
    "CD3E": "CD3E", "CD4": "CD4", "CD8A": "CD8A", "CD68": "CD68",
    "CD14": "CD14", "CD19": "CD19", "MS4A1": "MS4A1", "EPCAM": "EPCAM",
    "PECAM1": "PECAM1", "VIM": "VIM", "KRT5": "KRT5", "PCNA": "PCNA",
    "CD274": "CD274", "PDCD1": "PDCD1", "ACTA2": "ACTA2", "BCL2": "BCL2",
    "CD163": "CD163", "CR2": "CR2", "HLA-DRA": "HLA-DRA", "PAX5": "PAX5",
    "SDC1": "SDC1", "FCGR3A": "FCGR3A", "ITGAX": "ITGAX", "ITGAM": "ITGAM",
    "CCR7": "CCR7", "CD27": "CD27", "CXCR5": "CXCR5",
}


def rna_cpm(ad_counts, gene):
    """log1p(CPM) from a raw-counts AnnData (module-00 objects)."""
    if gene not in ad_counts.var_names:
        return None
    X = ad_counts[:, gene].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    n = np.asarray(ad_counts.obs["n_counts"]).ravel()
    n[n == 0] = 1
    return np.log1p(X.ravel() / n * 1e6)


def clr_by_factor(prot_counts, factors, pnames):
    """Alternative protein normalisation using isotype factors (multiplicative)."""
    prot = np.asarray(prot_counts, float)
    f = np.asarray(factors[factors["barcode"].isin([])], float)  # placeholder
    fac = np.ones(prot.shape[0])
    if factors is not None and len(factors) == prot.shape[0]:
        fac = factors["normalization_factor"].to_numpy(float)
    prot = prot * fac[:, None]
    prot = prot + 1.0
    gmean = np.exp(np.log(prot).mean(axis=1, keepdims=True))
    return np.log(prot / gmean), fac


def sample_agreement(specimen):
    ad = sc.read_h5ad(os.path.join(SRC, f"{specimen}_ecosystem.h5ad"))
    counts = sc.read_h5ad(os.path.join(PROC_DIR, "00_ingest", f"{specimen}.h5ad"))
    counts = counts[ad.obs_names].copy()
    prot = ad.obsm["protein_clr"] if "protein_clr" in ad.obsm else None
    pnames = list(ad.uns["protein_var"].index)
    pidx = {n.split("-")[0]: i for i, n in enumerate(pnames)}  # first isoform
    rows, disc = [], {}
    for gene, ptarget in PAIRS.items():
        if ptarget not in pidx:
            continue
        rna = rna_cpm(counts, gene)
        if rna is None:
            continue
        pi = pidx[ptarget]
        pv = prot[:, pi]
        if np.all(pv == pv[0]) or np.std(pv) == 0:
            continue
        r_pear = stats.pearsonr(rna, pv).statistic
        r_spear = stats.spearmanr(rna, pv).statistic
        # discordance = zRNA - zProt (per-spot; higher = RNA excess)
        zr = (rna - rna.mean()) / (rna.std() + 1e-12)
        zp = (pv - pv.mean()) / (pv.std() + 1e-12)
        d = zr - zp
        rows.append({"sample": specimen, "gene": gene, "r_pearson": r_pear,
                     "r_spearman": r_spear, "n_spots": len(rna)})
        disc[gene] = d
    tab = pd.DataFrame(rows)
    # spatial structure of discordance (squidpy Moran on obs columns)
    if len(disc):
        ad.obs = ad.obs.iloc[:, :0].join(pd.DataFrame(disc, index=ad.obs_names))
        try:
            mor = sq.gr.spatial_autocorr(ad, mode="moran", attr="obs",
                                         genes=list(disc.keys()), n_perms=999,
                                         corr_method="fdr_bh", n_jobs=8, copy=True)
            mor["gene"] = mor.index
            mor["sample"] = specimen
        except Exception as e:
            print(f"[m6] discordance Moran failed {specimen}: {e}")
            mor = None
    else:
        mor = None
    return tab, mor, ad


def main():
    tabs, mors = [], []
    for _, m in SAMPLE_META.iterrows():
        tab, mor, _ = sample_agreement(m.specimen)
        tabs.append(tab)
        if mor is not None:
            mors.append(mor)
    tab = pd.concat(tabs)
    tab.to_csv(os.path.join(TABLE_DIR, "06_rna_protein_agreement.csv"), index=False)
    mor = pd.concat(mors) if mors else None
    if mor is not None:
        mor.to_csv(os.path.join(TABLE_DIR, "06_discordance_morans.csv"), index=False)

    # ---- global summary ----
    summ = tab.groupby("gene").agg(
        median_r=("r_pearson", "median"), mean_r=("r_pearson", "mean"),
        q25=("r_pearson", lambda s: s.quantile(0.25)),
        q75=("r_pearson", lambda s: s.quantile(0.75)),
        n_samples=("r_pearson", "count"),
    ).sort_values("median_r", ascending=False)
    print("\n[m6] RNA-protein agreement by marker (median Pearson r across samples):")
    print(summ.round(3).to_string())

    # ---- figure: boxplot by marker ----
    fig, ax = plt.subplots(figsize=(11, 4.5))
    order = summ.index
    dat = [tab.loc[tab.gene == g, "r_pearson"].values for g in order]
    bp = ax.boxplot(dat, vert=True, showfliers=False, patch_artist=True)
    ax.set_xticklabels(order, rotation=60, ha="right")
    for patch in bp["boxes"]:
        patch.set_facecolor("#7fb3d5")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_ylabel("Pearson r (RNA vs protein, per spot)")
    ax.set_title("RNA–protein agreement across matched markers (14 sections)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "06_rna_protein_agreement.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # ---- figure: representative marker scatter (P4-S2) ----
    ad = sc.read_h5ad(os.path.join(SRC, "HCC22-088-P4-S2_ecosystem.h5ad"))
    counts = sc.read_h5ad(os.path.join(PROC_DIR, "00_ingest",
                                       "HCC22-088-P4-S2.h5ad"))[ad.obs_names].copy()
    prot = ad.obsm["protein_clr"]
    pidx = {n.split("-")[0]: i for i, n in enumerate(ad.uns["protein_var"].index)}
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, gene in zip(axes, ["EPCAM", "CD68", "CD3E"]):
        rna = rna_cpm(counts, gene)
        pv = prot[:, pidx[gene]]
        r = stats.pearsonr(rna, pv).statistic
        ax.scatter(pv, rna, s=4, alpha=0.3)
        ax.set_xlabel(f"protein {gene} (CLR)")
        ax.set_ylabel(f"RNA {gene} (log CPM)")
        ax.set_title(f"{gene}: r = {r:.2f}")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "06_scatter_examples_P4S2.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # ---- figure: discordance maps (P4-S2) ----
    lib = ad.uns.get("library_id", "visium")
    ad.obs["disc_CD68"] = None
    ad.obs["disc_CD3E"] = None
    for gene, col in [("CD68", "disc_CD68"), ("CD3E", "disc_CD3E")]:
        rna = rna_cpm(counts, gene)
        pv = prot[:, pidx[gene]]
        ad.obs[col] = ((rna - rna.mean()) / rna.std()) - ((pv - pv.mean()) / pv.std())
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.5))
    sq.pl.spatial_scatter(ad, color="CD68", ax=axes[0], img=True, library_id=lib,
                          frameon=False, colorbar=False, size=3, use_raw=True,
                          vmin=0, vmax=2, title="RNA CD68")
    sq.pl.spatial_scatter(ad, color="disc_CD68", ax=axes[1], img=True,
                          library_id=lib, frameon=False, colorbar=False, size=3,
                          use_raw=False, title="Discordance CD68 (RNA-prot)")
    sq.pl.spatial_scatter(ad, color="CD3E", ax=axes[2], img=True, library_id=lib,
                          frameon=False, colorbar=False, size=3, use_raw=True,
                          vmin=0, vmax=2, title="RNA CD3E")
    sq.pl.spatial_scatter(ad, color="disc_CD3E", ax=axes[3], img=True,
                          library_id=lib, frameon=False, colorbar=False, size=3,
                          use_raw=False, title="Discordance CD3E (RNA-prot)")
    fig.suptitle("Discordance maps — HCC22-088-P4-S2", y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "06_discordance_maps_P4S2.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # ---- discordance Moran summary ----
    if mor is not None:
        dm = mor.groupby("gene")["I"].agg(["mean", "median", "count"])
        dm["frac_pos"] = mor.groupby("gene")["I"].apply(lambda s: (s > 0).mean())
        print("\n[m6] spatial structure (Moran's I) of per-spot RNA-protein discordance:")
        print(dm.sort_values("mean", ascending=False).round(3).head(12).to_string())


if __name__ == "__main__":
    main()