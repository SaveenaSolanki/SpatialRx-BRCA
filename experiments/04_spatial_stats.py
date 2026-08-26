"""
Module 04 — Spatial statistics: Moran's I for spatially structured expression.
1. Builds spatial (grid) neighbour graphs per sample.
2. Runs Squidpy Moran's I on log-normalized RNA for all genes.
3. Cross-validates against the authors' uploaded spatial_enrichment.csv (their own
   per-feature Moran's I table), reporting the correlation between the two.
4. Ranks spatially variable genes (SVGs) & writes per-sample + cohort tables.

Usage: conda run -n spatialrx python experiments/04_spatial_stats.py
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from srx import PROC_DIR, FIG_DIR, TABLE_DIR, SAMPLE_META, MARKERS

sc.settings.verbosity = 0
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)
OUT = os.path.join(PROC_DIR, "03_domains")


def morans_per_sample(specimen, gsm, perm_genes=None):
    ad = sc.read_h5ad(os.path.join(OUT, f"{specimen}_domains.h5ad"))
    # grid neighbours from array coords (Visium hexagonal lattice): squidpy reads
    # obsm[spatial_key], so temporarily swap in lattice coordinates
    orig = np.asarray(ad.obsm["spatial"]).copy()
    ad.obsm["spatial"] = np.column_stack([ad.obs["array_row"], ad.obs["array_col"]])
    sq.gr.spatial_neighbors(ad, coord_type="grid", n_neighs=6, key_added="spatial_neighbors")
    ad.obsm["spatial"] = orig
    genes_all = list(ad.raw.var_names)
    # stage 1: exact Moran's I for ALL genes (no permutations = fast)
    mor = sq.gr.spatial_autocorr(ad, mode="moran", use_raw=True,
                                 genes=genes_all, n_perms=None,
                                 connectivity_key="spatial_neighbors_connectivities",
                                 corr_method="fdr_bh", n_jobs=8, copy=True)
    # stage 2: permutations only for candidate genes (top by |I|)
    if perm_genes is not None:
        cand = [g for g in perm_genes if g in ad.raw.var_names]
        morp = sq.gr.spatial_autocorr(ad, mode="moran", use_raw=True,
                                      genes=cand, n_perms=999,
                                      connectivity_key="spatial_neighbors_connectivities",
                                      corr_method="fdr_bh", n_jobs=8, copy=True)
        mor = mor.drop(columns=[c for c in ["pval_norm", "pval_z_sim", "pval_sim",
                                            "pval_norm_fdr_bh", "pval_z_sim_fdr_bh",
                                            "pval_sim_fdr_bh"] if c in mor.columns])
        mor = mor.merge(morp[["pval_norm_fdr_bh", "pval_sim_fdr_bh", "pval_z_sim_fdr_bh"]],
                        left_index=True, right_index=True, how="left")
    mor["gene"] = mor.index
    mor["sample"] = specimen
    # cross-check with authors' table
    auth = pd.read_csv(os.path.join(PROC_DIR, "00_ingest",
                                    f"{specimen}__spatial_enrichment.csv"))
    auth = auth[auth["Feature Type"] == "Gene Expression"]
    m = auth[["Feature Name", "I", "P value"]].rename(
        columns={"Feature Name": "gene", "I": "I_auth", "P value": "p_auth"})
    cmp = mor.merge(m, on="gene")
    r = (cmp[["I", "I_auth"]].corr(method="pearson").iloc[0, 1]
         if len(cmp) > 10 else np.nan)
    print(f"[m4] {specimen}: n={len(mor)} genes; Moran's I vs authors r={r:.3f} "
          f"({len(cmp)} matched genes)")
    return mor, cmp, ad


def main():
    rows, authall = [], []
    # first pass: exact I for all genes (no perms)
    for _, m in SAMPLE_META.iterrows():
        mor, cmp, _ = morans_per_sample(m.specimen, m.gsm, perm_genes=None)
        rows.append(mor)
        authall.append(cmp)
    mor = pd.concat(rows)
    # top candidate genes by mean |I| -> second pass with permutations
    top_cand = (mor.groupby("gene")["I"]
                   .agg(lambda s: np.mean(np.abs(s)))
                   .nlargest(500).index)
    rows2 = []
    for _, m in SAMPLE_META.iterrows():
        mor2, _cmp, _ = morans_per_sample(m.specimen, m.gsm, perm_genes=top_cand)
        rows2.append(mor2)
    mor = pd.concat(rows2)
    mor.to_csv(os.path.join(TABLE_DIR, "04_morans_I_per_sample.csv"), index=False)

    # ---- SVG ranking across samples ----
    svg = (mor.groupby("gene")["I"]
              .agg(["mean", "median", lambda s: np.percentile(s, 25), "count"])
              .rename(columns={"<lambda_0>": "I_Q25"}))
    svg = svg[svg["count"] >= 8].sort_values("mean", ascending=False)
    svg.to_csv(os.path.join(TABLE_DIR, "04_svg_ranking.csv"))
    print("\n[m4] top 25 spatially variable genes (mean Moran's I across samples):")
    print(svg.head(25).round(3).to_string())
    print("\n[m4] bottom 10 (least structured):")
    print(svg.tail(10).round(3).to_string())

    # canonical markers: Moran's I summary
    canon = {g: k for k, gs in MARKERS.items() for g in gs}
    mm = mor[mor["gene"].isin(canon)].copy()
    mm["program"] = mm["gene"].map(canon)
    summ = mm.groupby("gene").agg(mean_I=("I", "mean"), sd_I=("I", "std"),
                                  n=("I", "count")).sort_values("mean_I", ascending=False)
    print("\n[m4] canonical lineage markers by mean Moran's I:")
    print(summ.round(3).to_string())

    # ---- author cross-validation figure ----
    cmpall = pd.concat(authall)
    r_all = cmpall[["I", "I_auth"]].corr(method="pearson").iloc[0, 1]
    r_sp = cmpall.groupby("sample")[["I", "I_auth"]].corr().iloc[0::2, 1].reset_index()
    r_sp = r_sp.rename(columns={r_sp.columns[-1]: "I"})
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].scatter(cmpall["I_auth"], cmpall["I"], s=3, alpha=0.25, color="#1f77b4")
    axes[0].set_xlabel("Moran's I (authors' spatial_enrichment.csv)")
    axes[0].set_ylabel("Moran's I (our Squidpy run)")
    axes[0].set_title(f"All samples pooled: Pearson r = {r_all:.3f}")
    axes[1].bar(range(len(r_sp)), r_sp["I"], color="#2ca02c")
    axes[1].set_xticks(range(len(r_sp)))
    axes[1].set_xticklabels([s.split("-")[-2] + "-" + s.split("-")[-1]
                             for s in r_sp["sample"]], rotation=45, ha="right")
    axes[1].set_ylabel("per-sample r")
    axes[1].set_ylim(0, 1)
    axes[1].set_title("Per-sample agreement with authors' Moran's I")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "04_morans_validation.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # ---- SVG heatmap (top 30 genes x samples) ----
    top = svg.head(30).index
    pivot = mor[mor.gene.isin(top)].pivot_table(index="gene", columns="sample",
                                                values="I")
    fig, ax = plt.subplots(figsize=(10, 8))
    sns_ = __import__("seaborn")
    sns_.heatmap(pivot, cmap="RdBu_r", center=0, ax=ax)
    ax.set_title("Top 30 SVGs — Moran's I per sample")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "04_svg_heatmap.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # ---- spatial map of a few markers with high Moran's I ----
    ad = sc.read_h5ad(os.path.join(OUT, "HCC22-088-P4-S2_domains.h5ad"))
    sq.gr.spatial_neighbors(ad, coord_type="grid", n_neighs=6)
    picks = ["COL1A1", "EPCAM", "ESR1", "CD68"]
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.5))
    for ax, g in zip(axes, picks):
        sq.pl.spatial_scatter(ad, color=g, ax=ax, img=True,
                              library_id=ad.uns["library_id"],
                              frameon=False, colorbar=False, size=3,
                              use_raw=True, title=f"{g} (I={svg.loc[g,'mean']:.2f} "
                              f"if in top-SVGs else local)")
    fig.suptitle("Highly structured genes on tissue (HCC22-088-P4-S2)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "04_svg_spatial_map.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()