"""
Module 06b — Spatial structure of the protein (antibody) layer.
Computes Moran's I of CLR-normalised proteins per section and validates against the
authors' spatial_enrichment.csv entries for the same antibody feature names.
Also reports how many proteins are spatially structured (FDR<0.05) per section.

Usage: conda run -n spatialrx python experiments/06b_protein_spatial.py
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
from srx import PROC_DIR, FIG_DIR, TABLE_DIR, SAMPLE_META

sc.settings.verbosity = 0
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)
ECO = os.path.join(PROC_DIR, "05_ecosystem")


def main():
    rows, cmp_rows = [], []
    for _, m in SAMPLE_META.iterrows():
        ad = sc.read_h5ad(os.path.join(ECO, f"{m.specimen}_ecosystem.h5ad"))
        pnames = list(ad.uns["protein_var"].index)
        X = ad.obsm["protein_clr"]
        # store protein matrix as an obs field for squidpy
        tmp = pd.DataFrame(X, index=ad.obs_names, columns=pnames)
        ad.obs = ad.obs.join(tmp)
        orig = np.asarray(ad.obsm["spatial"]).copy()
        ad.obsm["spatial"] = np.column_stack([ad.obs["array_row"], ad.obs["array_col"]])
        if "spatial_connectivities" not in ad.obsp:
            sq.gr.spatial_neighbors(ad, coord_type="grid", n_neighs=6)
        ad.obsm["spatial"] = orig
        mor = sq.gr.spatial_autocorr(ad, mode="moran", attr="obs", genes=pnames,
                                     n_perms=199, corr_method="fdr_bh",
                                     copy=True, n_jobs=8)
        mor["protein"] = mor.index
        mor["sample"] = m.specimen
        mor["feature"] = [p.split("-")[0] for p in mor.index]
        rows.append(mor)
        # authors' table
        auth = pd.read_csv(os.path.join(PROC_DIR, "00_ingest",
                                        f"{m.specimen}__spatial_enrichment.csv"))
        auth_ab = auth[auth["Feature Type"] == "Antibody Capture"]
        aa = auth_ab[["Feature Name", "I"]].rename(
            columns={"Feature Name": "feature", "I": "I_auth"})
        cmp = mor.merge(aa, on="feature", how="inner")
        if len(cmp) > 5:
            r = cmp[["I", "I_auth"]].corr().iloc[0, 1]
        else:
            r = np.nan
        print(f"[m6b] {m.specimen}: {len(mor)} proteins; n_sig(FDR<0.05)="
              f"{int((mor.pval_norm_fdr_bh < 0.05).sum())}; "
              f"vs authors r={r:.3f} ({len(cmp)} matched)")
        cmp_rows.append(cmp)
    mor_all = pd.concat(rows)
    mor_all.to_csv(os.path.join(TABLE_DIR, "06b_protein_morans.csv"), index=False)
    cmp_all = pd.concat(cmp_rows)
    r_all = cmp_all[["I", "I_auth"]].corr().iloc[0, 1]
    # aggregate protein spatial structure
    summ = mor_all.groupby("protein")["I"].agg(["mean", "median", "count"])
    summ = summ.sort_values("mean", ascending=False)
    print("\n[m6b] top-15 spatially structured proteins (mean Moran's I):")
    print(summ.head(15).round(3).to_string())
    print(f"\n[m6b] protein Moran's I vs authors (all matched): r={r_all:.3f}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].scatter(cmp_all["I_auth"], cmp_all["I"], s=6, alpha=0.4, color="#9467bd")
    axes[0].set_xlabel("authors' Moran's I (antibody)")
    axes[0].set_ylabel("our Moran's I (CLR protein)")
    axes[0].set_title(f"protein validation r={r_all:.3f}")
    axes[1].bar(np.arange(len(summ.head(15))), summ.head(15)["mean"], color="#2ca02c")
    axes[1].set_xticks(np.arange(len(summ.head(15))))
    axes[1].set_xticklabels(summ.head(15).index, rotation=60, ha="right", fontsize=7)
    axes[1].set_ylabel("mean Moran's I")
    axes[1].set_title("Most spatially structured proteins")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "06b_protein_spatial.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()