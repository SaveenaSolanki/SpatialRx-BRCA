"""
Module 11 — Are RNA–protein discordance hotspots biological (cell-state-driven) or
technical (tissue-edge/antibody-penetration artifacts)?
Hypothesis test: per-spot discordance (zRNA - zProt) for key markers should correlate
with cell-state proportions and pathway programs IF discordance reflects regional
post-transcriptional regulation; it should correlate with distance-to-tissue-edge
IF it is an antibody-penetration artifact.

Controls/metrics per marker:
  1. Spearman rho between discordance and rna_<lineage> proportions + program scores
     (per sample, pooled with sign consistency).
  2. Spearman rho between discordance and distance to nearest non-tissue spot
     (edge-distance) — artifact check.
  3. Hotspot (top quartile) vs elsewhere: cell-state composition differences.
  4. Pre vs post comparison of hotspot cell states (paired, n<=6).

Usage: conda run -n spatialrx python experiments/11_discordance_cellstate.py
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
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from srx import PROC_DIR, FIG_DIR, TABLE_DIR, SAMPLE_META, PAIRED

sc.settings.verbosity = 0
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)
ECO = os.path.join(PROC_DIR, "05_ecosystem")
ING = os.path.join(PROC_DIR, "00_ingest")

MARKERS = ["CD68", "CD3E", "CD163", "HLA-DRA", "CD14", "EPCAM", "ACTA2", "VIM", "PCNA"]
LINEAGES = ["T_cells", "B_cells", "myeloid", "fibroblast", "endothelial",
            "epithelial_tumour"]
PROGRAMS = ["score_M2_macrophage", "score_exhaustion", "score_Proliferation",
            "score_Hypoxia", "score_Interferon_response", "score_Fibroblast_activation",
            "score_Estrogen_response"]

M2_GENES = ["CD163", "MRC1", "MSR1", "TGFBI"]
EXH_GENES = ["PDCD1", "CTLA4", "LAG3", "HAVCR2"]


def rna_cpm(counts_ad, gene):
    if gene not in counts_ad.var_names:
        return None
    X = counts_ad[:, gene].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    n = np.asarray(counts_ad.obs["n_counts"]).ravel().astype(np.float64)
    n[n == 0] = 1
    return np.log1p(X.ravel() / n * 1e6)


def edge_exposure(ad):
    """Continuous tissue-edge exposure: 6 - grid degree (0 = deep interior,
    6 = fully exposed). Valid for core biopsies whose spots are all peripheral."""
    if "spatial_connectivities" not in ad.obsp:
        import squidpy as sq
        orig = np.asarray(ad.obsm["spatial"]).copy()
        ad.obsm["spatial"] = np.column_stack([ad.obs["array_row"], ad.obs["array_col"]])
        sq.gr.spatial_neighbors(ad, coord_type="grid", n_neighs=6)
        ad.obsm["spatial"] = orig
    C = ad.obsp["spatial_connectivities"]
    deg = np.asarray(C.sum(axis=1)).ravel()
    return (6.0 - deg).astype(float)

def main():
    rows = []
    for spec in PAIRED:
        ad = sc.read_h5ad(os.path.join(ECO, f"{spec}_ecosystem.h5ad"))
        counts = sc.read_h5ad(os.path.join(ING, f"{spec}.h5ad"))[ad.obs_names].copy()
        pnames = list(ad.uns["protein_var"].index)
        pidx = {p.split("-")[0]: i for i, p in enumerate(pnames)}
        prot = ad.obsm["protein_clr"]
        edist = edge_exposure(ad)
        # ensure programs: M2 and exhaustion (may not exist in obs)
        sc.tl.score_genes(ad, [g for g in M2_GENES if g in ad.raw.var_names],
                          score_name="score_M2_macrophage", use_raw=True)
        sc.tl.score_genes(ad, [g for g in EXH_GENES if g in ad.raw.var_names],
                          score_name="score_exhaustion", use_raw=True)
        for met in MARKERS:
            if met not in pidx:
                continue
            rna = rna_cpm(counts, met)
            if rna is None:
                continue
            pv = prot[:, pidx[met]]
            zr = (rna - rna.mean()) / (rna.std() + 1e-12)
            zp = (pv - pv.mean()) / (pv.std() + 1e-12)
            disc = zr - zp
            entry = {"sample": spec, "patient": spec.split("-")[2],
                     "treatment": "pre" if "-S1" in spec else "post", "marker": met,
                     "disc_edge_rho": stats.spearmanr(disc, edist).statistic,
                     "hotspot_frac": float((disc > np.percentile(disc, 75)).mean())}
            for ln in LINEAGES:
                if f"rna_{ln}" in ad.obs:
                    entry[f"rho_{ln}"] = stats.spearmanr(
                        disc, ad.obs[f"rna_{ln}"]).statistic
            for pr in PROGRAMS:
                if pr in ad.obs:
                    entry[f"rho_{pr}"] = stats.spearmanr(
                        disc, ad.obs[pr]).statistic
            rows.append(entry)
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(TABLE_DIR, "11_discordance_cellstate.csv"), index=False)

    # pooled summary per marker
    rho_cols = [c for c in tab.columns if c.startswith("rho_")]
    summ = []
    for met in MARKERS:
        sub = tab[tab["marker"] == met]
        if sub.empty:
            continue
        r = {"marker": met, "n_samples": len(sub),
             "disc_edge_rho_median": sub["disc_edge_rho"].median()}
        for c in rho_cols:
            v = sub[c]
            r[f"{c}_median"] = v.median()
            r[f"{c}_signfrac"] = (v > 0).mean()
        summ.append(r)
    sm = pd.DataFrame(summ).set_index("marker")
    sm.to_csv(os.path.join(TABLE_DIR, "11_discordance_cellstate_summary.csv"))
    pd.set_option("display.width", 220)
    print("\n[m11] discordance correlates with edge-distance (artifact check) "
          "and cell states (median rho, frac samples with same sign):")
    print(sm.round(3).to_string())

    # key: which cell states track CD68 / CD3E / EPCAM discordance
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for ax, met in zip(axes.ravel(), ["CD68", "CD3E", "EPCAM", "ACTA2"]):
        sub = tab[tab["marker"] == met]
        cols = ["rho_T_cells", "rho_myeloid", "rho_epithelial_tumour",
                "rho_fibroblast", "rho_score_Hypoxia", "rho_score_M2_macrophage",
                "rho_score_exhaustion", "disc_edge_rho"]
        cols = [c for c in cols if c in sub.columns]
        med = [sub[c].median() for c in cols]
        ax.barh(range(len(cols)), med, color=["#d62728" if v > 0 else "#1f77b4"
                                              for v in med])
        ax.set_yticks(range(len(cols)))
        ax.set_yticklabels([c.replace("rho_score_", "").replace("rho_", "") for c in cols])
        ax.axvline(0, color="k", lw=0.7)
        ax.set_title(met)
    fig.suptitle("What drives RNA–protein discordance? (median Spearman rho)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "11_discordance_drivers.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # hotspots (top quartile) composition vs elsewhere, pre vs post
    hot_rows = []
    for met in MARKERS:
        sub = tab[tab["marker"] == met]
        for spec in PAIRED:
            a = sub[(sub["sample"] == spec) & (sub["marker"] == met)]
            if a.empty:
                continue
            ad = sc.read_h5ad(os.path.join(ECO, f"{spec}_ecosystem.h5ad"))
            counts = sc.read_h5ad(os.path.join(ING, f"{spec}.h5ad"))[ad.obs_names].copy()
            pnames = list(ad.uns["protein_var"].index)
            pidx = {p.split("-")[0]: i for i, p in enumerate(pnames)}
            if met not in pidx:
                continue
            rna = rna_cpm(counts, met)
            pv = ad.obsm["protein_clr"][:, pidx[met]]
            disc = ((rna - rna.mean()) / rna.std()) - ((pv - pv.mean()) / pv.std())
            hot = disc > np.percentile(disc, 75)
            for ln in LINEAGES:
                if f"rna_{ln}" in ad.obs:
                    hot_rows.append({"sample": spec, "marker": met, "state": ln,
                                     "hotspot_delta": ad.obs.loc[hot, f"rna_{ln}"].mean()
                                     - ad.obs.loc[~hot, f"rna_{ln}"].mean()})
    ht = pd.DataFrame(hot_rows)
    ht.to_csv(os.path.join(TABLE_DIR, "11_hotspot_state_deltas.csv"), index=False)
    hm = ht.groupby("marker")["state"].apply(
        lambda s: pd.Series({st: ht.loc[(ht["marker"] == s.name) & (ht["state"] == st),
                                        "hotspot_delta"].median() for st in LINEAGES}))
    pp = hm.T
    print("\n[m11] median cell-state delta in discordance hotspots vs elsewhere "
          "(per 100 proportion points):")
    print((pp * 100).round(2).to_string())


if __name__ == "__main__":
    main()