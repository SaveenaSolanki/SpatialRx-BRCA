"""
Module 12 — Matched-region subsampling for treatment tests.
The pre sections are small core biopsies (294-2,256 spots), post sections are whole
surgical specimens (3,598-4,992). Any pre/post comparison of *absolute* quantities is
confounded by section size; we test whether conclusions hold when post sections are
subsampled to the size of the paired pre biopsy:
  Strategy A (random): uniformly subsample n_pre spots from the post section.
  Strategy B (contiguous): sample a contiguous "needle-core-like" region (BFS on the
    grid graph) of n_pre spots from the post section.
Metrics recomputed on each subsample; deltas (post-pre) evaluated per patient and pooled.
Report: full-section delta, median subsampled delta, 95% CI, sign stability.

Usage: conda run -n spatialrx python experiments/12_subsampling.py
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
from srx import PROC_DIR, FIG_DIR, TABLE_DIR, PAIRED

sc.settings.verbosity = 0
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)
ECO = os.path.join(PROC_DIR, "05_ecosystem")

METRICS = ["frac_Tumour", "frac_Immune", "score_Estrogen_response",
           "score_Proliferation", "score_T_cell_activity",
           "score_Myeloid_inflammation", "tumour_immune_contact"]
B_RANDOM = 200
B_BFS = 100


def graph(ad):
    if "spatial_connectivities" not in ad.obsp:
        orig = np.asarray(ad.obsm["spatial"]).copy()
        ad.obsm["spatial"] = np.column_stack([ad.obs["array_row"], ad.obs["array_col"]])
        sq.gr.spatial_neighbors(ad, coord_type="grid", n_neighs=6)
        ad.obsm["spatial"] = orig
    return ad.obsp["spatial_connectivities"].tocsr()


def metrics_subset(ad, idx):
    idx = np.asarray(idx)
    comp = ad.obs["compartment"].to_numpy()[idx]
    m = {"frac_Tumour": (comp == "Tumour").mean(),
         "frac_Immune": (comp == "Immune").mean()}
    for c in [c for c in METRICS if c.startswith("score_")]:
        m[c] = ad.obs[c].to_numpy()[idx].mean()
    # tumour-immune contact within subset (via subset adjacency)
    C = ad.obsp["spatial_connectivities"].tocsr()
    subC = C[idx][:, idx]
    tum = comp == "Tumour"
    imm = comp == "Immune"
    if tum.sum() == 0:
        m["tumour_immune_contact"] = np.nan
    else:
        neigh_imm = np.asarray((subC @ imm.astype(float)) > 0).ravel()
        m["tumour_immune_contact"] = neigh_imm.mean()
    return m


def bfs_core(C, seed, k):
    n = C.shape[0]
    seen = np.zeros(n, dtype=bool)
    order = [seed]
    seen[seed] = True
    frontier = [seed]
    Cc = C.tocsr()
    while len(order) < k and frontier:
        nxt = []
        for f in frontier:
            nbrs = Cc.indices[Cc.indptr[f]:Cc.indptr[f + 1]]
            for nb in nbrs:
                if not seen[nb]:
                    seen[nb] = True
                    order.append(nb)
                    nxt.append(nb)
                    if len(order) >= k:
                        break
            if len(order) >= k:
                break
        frontier = nxt
    return np.array(order[:k])


def main():
    rng = np.random.default_rng(0)
    pre_metrics = {}
    rows = []
    for spec in PAIRED:
        patient = spec.split("-")[2]
        trt = "pre" if "-S1" in spec else "post"
        ad = sc.read_h5ad(os.path.join(ECO, f"{spec}_ecosystem.h5ad"))
        m = metrics_subset(ad, np.arange(ad.n_obs))
        if trt == "pre":
            pre_metrics[patient] = m
            continue
        # post: subsample
        C = graph(ad)
        n_pre = sum(1 for s in PAIRED if s.split("-")[2] == patient and "-S1" in s)
        # n_pre spots on the paired pre section
        pre_spec = [s for s in PAIRED if s.split("-")[2] == patient and "-S1" in s][0]
        pre_ad = sc.read_h5ad(os.path.join(ECO, f"{pre_spec}_ecosystem.h5ad"))
        k = pre_ad.n_obs
        for met in METRICS:
            rows.append({"patient": patient, "metric": met, "strategy": "full",
                         "delta": m[met] - pre_metrics[patient][met]})
        for rep in range(B_RANDOM):
            idx = rng.choice(ad.n_obs, size=k, replace=False)
            sm = metrics_subset(ad, idx)
            for met in METRICS:
                rows.append({"patient": patient, "metric": met,
                             "strategy": "random", "rep": rep,
                             "delta": sm[met] - pre_metrics[patient][met]})
        for rep in range(B_BFS):
            seed = int(rng.integers(0, ad.n_obs))
            idx = bfs_core(C, seed, k)
            if len(idx) < max(10, k // 2):
                continue
            sm = metrics_subset(ad, idx)
            for met in METRICS:
                rows.append({"patient": patient, "metric": met,
                             "strategy": "contiguous", "rep": rep,
                             "delta": sm[met] - pre_metrics[patient][met]})
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(TABLE_DIR, "12_subsampling_deltas.csv"), index=False)

    # pooled summaries
    summ = []
    for met in METRICS:
        full = tab[(tab.metric == met) & (tab.strategy == "full")]["delta"]
        for strat in ["random", "contiguous"]:
            d = tab[(tab.metric == met) & (tab.strategy == strat)]["delta"].dropna()
            if len(d) == 0:
                continue
            summ.append({"metric": met, "strategy": strat,
                         "full_mean_delta": full.mean(),
                         "sub_median_delta": d.median(),
                         "sub_ci_lo": d.quantile(0.025), "sub_ci_hi": d.quantile(0.975),
                         "sign_stability_gt0": (d > 0).mean(),
                         "full_in_ci": float(full.mean() >= d.quantile(0.025)
                                             and full.mean() <= d.quantile(0.975))})
    sm = pd.DataFrame(summ)
    sm.to_csv(os.path.join(TABLE_DIR, "12_subsampling_summary.csv"), index=False)
    pd.set_option("display.width", 160)
    print("\n[m12] pre/post deltas: full section vs size-matched subsampling")
    print(sm.round(3).to_string(index=False))

    # figure: distribution of subsampled deltas vs full-section delta
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    for ax, met in zip(axes.ravel(), METRICS):
        full = tab[(tab.metric == met) & (tab.strategy == "full")]
        for strat, color in [("random", "#1f77b4"), ("contiguous", "#2ca02c")]:
            d = tab[(tab.metric == met) & (tab.strategy == strat)]["delta"].dropna()
            if len(d):
                ax.hist(d, bins=40, alpha=0.45, color=color,
                        label=strat, density=True)
        ax.axvline(full["delta"].mean(), color="#d62728", lw=2, label="full section")
        ax.axvline(0, color="k", lw=0.7)
        ax.set_title(met, fontsize=9)
        ax.legend(fontsize=7)
    fig.suptitle("Does the pre/post effect survive size-matched subsampling of post sections?")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "12_subsampling.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # per-patient medians (which patients carry each effect)
    per = tab[tab.strategy == "random"].groupby(["metric", "patient"])["delta"].median()
    per = per.unstack().T
    per.to_csv(os.path.join(TABLE_DIR, "12_subsampling_per_patient.csv"))
    print("\n[m12] per-patient median deltas (random subsampling):")
    print(per.round(3).to_string())


if __name__ == "__main__":
    main()