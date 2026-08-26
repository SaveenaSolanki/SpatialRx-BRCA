"""
Module 10 — Ligand–receptor analysis at tumour–immune interfaces (pre vs post).
Approach (transparent, reproducible; no external DB):
  - curated ligand–receptor pairs (immune <-> tumour, tumour <-> stroma) from the
    published literature (CellChat/CSF/LIANA canonical pairs), filtered to genes
    present in the Visium panel
  - per-spot log-CPM expression; sender/receiver compartment definitions
  - interaction score = mean(sender ligand) * mean(receiver receptor)
  - interface spots: Tumour compartment spots with >=1 Immune neighbour (grid graph)
      and Immune spots with >=1 Tumour neighbour
  - compares interface vs non-interface, and pre vs post (paired, n=6 patients)
  Includes explicit MIF-CD74(+CD44) and MDK pairs from the CITEgeist paper's case study.

Usage: conda run -n spatialrx python experiments/10_lr_analysis.py
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
from srx import PROC_DIR, FIG_DIR, TABLE_DIR, SAMPLE_META, PAIRED

sc.settings.verbosity = 0
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)
ECO = os.path.join(PROC_DIR, "05_ecosystem")
ING = os.path.join(PROC_DIR, "00_ingest")

# ligand-receptor pairs: (ligand, receptor, category)
LR_PAIRS = [
    # immune -> tumour
    ("CXCL10", "CXCR3", "Immune->Tumour"), ("CXCL9", "CXCR3", "Immune->Tumour"),
    ("CXCL11", "CXCR3", "Immune->Tumour"), ("CCL5", "CCR5", "Immune->Tumour"),
    ("CCL2", "CCR2", "Immune->Tumour"), ("CCL3", "CCR1", "Immune->Tumour"),
    ("TNF", "TNFRSF1A", "Immune->Tumour"), ("IL1B", "IL1R1", "Immune->Tumour"),
    ("IFNG", "IFNGR1", "Immune->Tumour"), ("CSF1", "CSF1R", "Immune->Tumour"),
    ("MIF", "CD74", "Immune->Tumour"),  # paper's finding; CD44 co-receptor added below
    ("CD274", "PDCD1", "Immune->Immune(checkpoint)"),
    ("LGALS9", "HAVCR2", "Immune->Immune(checkpoint)"),
    ("CD40LG", "CD40", "Immune->Immune"),
    # tumour -> immune
    ("MDK", "LRP1", "Tumour->Immune"), ("MDK", "ITGA6", "Tumour->Immune"),
    ("TGFB1", "TGFBR1", "Tumour->Immune"), ("VEGFA", "KDR", "Tumour->Immune"),
    ("VEGFA", "FLT1", "Tumour->Immune"), ("AREG", "EGFR", "Tumour->Immune"),
    ("HBEGF", "EGFR", "Tumour->Immune"), ("CXCL12", "CXCR4", "Tumour->Immune"),
    ("CSF1", "CSF1R", "Tumour->Immune(myeloid)"),
    # tumour <-> stroma
    ("PDGFB", "PDGFRB", "Tumour->Stroma"), ("FGF2", "FGFR1", "Tumour->Stroma"),
    ("COL1A1", "SDC1", "Stroma->Tumour"),
]
# tri-molecular pairs with a co-receptor
TRIPLE = [("MIF", "CD44", "CD74", "Immune->Tumour"), ("MDK", "NCL", "LRP1", "Tumour->Immune")]


def log_cpm(counts_ad, gene):
    if gene not in counts_ad.var_names:
        return None
    X = counts_ad[:, gene].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    n = np.asarray(counts_ad.obs["n_counts"]).ravel()
    n[n == 0] = 1
    return np.log1p(X.ravel() / n * 1e6)


def interface_spots(ad):
    if "spatial_connectivities" not in ad.obsp:
        orig = np.asarray(ad.obsm["spatial"]).copy()
        ad.obsm["spatial"] = np.column_stack([ad.obs["array_row"], ad.obs["array_col"]])
        sq.gr.spatial_neighbors(ad, coord_type="grid", n_neighs=6)
        ad.obsm["spatial"] = orig
    C = ad.obsp["spatial_connectivities"]
    comp = ad.obs["compartment"].to_numpy()
    tum = comp == "Tumour"
    imm = comp == "Immune"
    tum_interface = tum & (np.asarray((C @ imm.astype(float)) > 0).ravel())
    imm_interface = imm & (np.asarray((C @ tum.astype(float)) > 0).ravel())
    return tum_interface, imm_interface, tum, imm


def main():
    lr = pd.DataFrame(LR_PAIRS, columns=["ligand", "receptor", "category"])
    results = []
    for _, m in SAMPLE_META.iterrows():
        spec = m.specimen
        ad = sc.read_h5ad(os.path.join(ECO, f"{spec}_ecosystem.h5ad"))
        counts = sc.read_h5ad(os.path.join(ING, f"{spec}.h5ad"))[ad.obs_names].copy()
        tif, iif, tum, imm = interface_spots(ad)
        if tum.sum() == 0 or imm.sum() == 0:
            print(f"[m10] {spec}: no tumour/immune pair -> skip")
            continue
        for _, row in lr.iterrows():
            lig = log_cpm(counts, row.ligand)
            rec = log_cpm(counts, row.receptor)
            if lig is None or rec is None:
                continue
            # tumour->immune uses tumour sender, immune receiver
            if row.category.startswith("Tumour"):
                send, recv = tif, iif
            elif row.category.startswith("Immune->Immune") or row.category == "Immune->Tumour":
                send, recv = iif, tif
            else:  # stroma involved: use all tumour/immune interface vs non-interface
                send, recv = tif, iif
            # interface score
            sc_lig = lig[send].mean() if send.sum() else np.nan
            sc_rec = rec[recv].mean() if recv.sum() else np.nan
            interface_score = sc_lig * sc_rec
            # non-interface reference (tumour spots without immune neighbour etc.)
            if row.category.startswith("Tumour"):
                send0, recv0 = tum & ~tif, imm & ~iif
            else:
                send0, recv0 = tum & ~tif, imm & ~iif
            ref = (lig[send0].mean() if send0.sum() else np.nan) * \
                  (rec[recv0].mean() if recv0.sum() else np.nan)
            results.append({"sample": spec, "patient": spec.split("-")[2],
                            "treatment": "pre" if "-S1" in spec else "post",
                            "ligand": row.ligand, "receptor": row.receptor,
                            "category": row.category,
                            "interface_score": interface_score, "ref_score": ref,
                            "n_interface": int(send.sum() + recv.sum())})
    res = pd.DataFrame(results)
    res.to_csv(os.path.join(TABLE_DIR, "10_lr_scores.csv"), index=False)

    # pooled summary: interface fold-change pre vs post
    summ = res.groupby(["ligand", "receptor", "category"]).agg(
        score_pre=("interface_score", lambda s: s[res.loc[s.index, "treatment"] == "pre"].mean()),
        score_post=("interface_score", lambda s: s[res.loc[s.index, "treatment"] == "post"].mean()),
        ref_pre=("ref_score", lambda s: s[res.loc[s.index, "treatment"] == "pre"].mean()),
        ref_post=("ref_score", lambda s: s[res.loc[s.index, "treatment"] == "post"].mean()),
        n=("interface_score", "count")).reset_index()
    summ["score_delta"] = summ["score_post"] - summ["score_pre"]
    summ["interface_enrichment"] = (summ["score_pre"] / summ["ref_pre"].replace(0, np.nan))
    summ = summ.sort_values("score_post", ascending=False)
    summ.to_csv(os.path.join(TABLE_DIR, "10_lr_summary.csv"), index=False)
    print("\n[m10] top interaction scores at tumour-immune interface (mean score):")
    print(summ.head(20).round(3).to_string(index=False))

    # paired deltas per pair (exact permutation as module 7)
    def paired_delta(sub):
        d = sub.pivot_table(index="patient", columns="treatment",
                            values="interface_score").dropna()
        if len(d) < 3 or "pre" not in d or "post" not in d:
            return None
        dv = (d["post"] - d["pre"]).to_numpy()
        obs = dv.mean()
        n = len(dv)
        s = 0
        for bits in range(2 ** n):
            signs = np.array([1 if (bits >> i) & 1 else -1 for i in range(n)])
            if abs((signs * dv).mean()) >= abs(obs):
                s += 1
        p = min(1.0, 2 * s / 2 ** n)
        return obs, p, len(dv)

    tpairs = []
    for (lig, rec, cat), sub in res.groupby(["ligand", "receptor", "category"]):
        r = paired_delta(sub)
        if r:
            tpairs.append({"ligand": lig, "receptor": rec, "category": cat,
                           "delta": r[0], "p_perm": r[1], "n": r[2]})
    tp = pd.DataFrame(tpairs).sort_values("p_perm")
    tp.to_csv(os.path.join(TABLE_DIR, "10_lr_paired_tests.csv"), index=False)
    print("\n[m10] paired pre->post interface-score tests:")
    print(tp.head(12).round(3).to_string(index=False))

    # figure: top pairs, pre vs post at interface
    top = summ.head(12)
    fig, ax = plt.subplots(figsize=(11, 5))
    x = np.arange(len(top))
    w = 0.38
    ax.bar(x - w / 2, top["score_pre"], w, label="pre", color="#1f77b4")
    ax.bar(x + w / 2, top["score_post"], w, label="post", color="#d62728")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r.ligand}→{r.receptor}" for _, r in top.iterrows()],
                       rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("mean interaction score (interface)")
    ax.legend()
    ax.set_title("Ligand–receptor scores at tumour–immune interface, pre vs post pET")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "10_lr_interface_prepost.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # figure: interface map for P4-S2 (MIF-CD74 as in the paper's case study)
    ad = sc.read_h5ad(os.path.join(ECO, "HCC22-088-P4-S2_ecosystem.h5ad"))
    counts = sc.read_h5ad(os.path.join(ING, "HCC22-088-P4-S2.h5ad"))[ad.obs_names].copy()
    tif, iif, tum, imm = interface_spots(ad)
    ad.obs["interface"] = "other"
    ad.obs.loc[tif | iif, "interface"] = "tumour-immune interface"
    mif = log_cpm(counts, "MIF")
    cd74 = log_cpm(counts, "CD74")
    ad.obs["MIF"] = mif
    ad.obs["CD74"] = cd74
    lib = ad.uns.get("library_id", "visium")
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    sq.pl.spatial_scatter(ad, color="compartment", ax=axes[0], img=True,
                          library_id=lib, frameon=False, size=3, title="compartments")
    sq.pl.spatial_scatter(ad, color="MIF", ax=axes[1], img=True, library_id=lib,
                          frameon=False, size=3, use_raw=False, colorbar=False,
                          title="MIF RNA (ligand)")
    sq.pl.spatial_scatter(ad, color="CD74", ax=axes[2], img=True, library_id=lib,
                          frameon=False, size=3, use_raw=False, colorbar=False,
                          title="CD74 RNA (receptor)")
    fig.suptitle("P4-S2: tumour–immune interface (MIF–CD74 axis)", y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "10_interface_map_P4S2.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()