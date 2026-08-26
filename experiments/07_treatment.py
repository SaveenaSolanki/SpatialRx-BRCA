"""
Module 07 — Endocrine therapy remodelling: paired pre/post spatial analysis.
Paired metrics per patient (n=6; replicate sections excluded from paired stats):
  1. compartment fractions (Tumour / Immune / Stroma / Mixed)
  2. RNA-deconv cell-state proportions
  3. gene-program scores (estrogen, proliferation, immune, hypoxia, ...)
  4. ESR1 expression (mean CPM within Tumour-compartment spots)
  5. tumour–immune contact (fraction of Tumour spots adjacent to Immune spots)
  6. neighbourhood enrichment z-scores (compartment + dominant-cell-state)
  7. spatial structure of programs (Moran's I per sample)
Statistics: paired per-patient deltas; Wilcoxon signed-rank + exact permutations;
BH correction across metrics (transparent about n=6).

CAVEAT: pre = core biopsy, post = surgical specimen -> absolute counts/areas differ;
only proportions / per-spot-normalized metrics are compared.

Usage: conda run -n spatialrx python experiments/07_treatment.py
"""
from __future__ import annotations

import os
import sys
from itertools import combinations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from srx import PROC_DIR, FIG_DIR, TABLE_DIR, SAMPLE_META, PAIRED, MARKERS, PROGRAMS

sc.settings.verbosity = 0
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)
SRC = os.path.join(PROC_DIR, "05_ecosystem")

PROGRAM_COLS = [f"score_{p}" for p in PROGRAMS]
LINEAGE_COLS = [f"rna_{ln}" for ln in ["T_cells", "B_cells", "myeloid",
                                       "fibroblast", "endothelial", "epithelial_tumour"]]


def tumour_immune_contact(ad):
    """Fraction of Tumour-compartment spots with >=1 Immune neighbour (grid graph)."""
    if "spatial_connectivities" not in ad.obsp:
        orig = np.asarray(ad.obsm["spatial"]).copy()
        ad.obsm["spatial"] = np.column_stack([ad.obs["array_row"], ad.obs["array_col"]])
        sq.gr.spatial_neighbors(ad, coord_type="grid", n_neighs=6)
        ad.obsm["spatial"] = orig
    comp = ad.obs["compartment"].to_numpy()
    C = ad.obsp["spatial_connectivities"]
    tum = np.where(comp == "Tumour")[0]
    if len(tum) == 0:
        return np.nan, 0
    imm = np.zeros(ad.n_obs, dtype=bool)
    imm[comp == "Immune"] = True
    neigh_has_imm = np.asarray((C[tum] @ imm.astype(float)) > 0).ravel()
    return neigh_has_imm.mean(), len(tum)


def per_sample_metrics(ad, specimen):
    m = {}
    comp = ad.obs["compartment"].value_counts(normalize=True)
    for c in ["Tumour", "Immune", "Stroma", "Mixed"]:
        m[f"frac_{c}"] = comp.get(c, 0.0)
    for c in LINEAGE_COLS:
        m[c] = ad.obs[c].mean()
    for c in PROGRAM_COLS:
        if c in ad.obs:
            m[c] = ad.obs[c].mean()
    # ESR1 within tumour compartment (RNA CPM on raw counts)
    counts = sc.read_h5ad(os.path.join(PROC_DIR, "00_ingest", f"{specimen}.h5ad"))
    counts = counts[ad.obs_names].copy()
    tum_mask = (ad.obs["compartment"] == "Tumour").to_numpy()
    if "ESR1" in counts.var_names and tum_mask.sum() > 5:
        X = counts[:, "ESR1"].X
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        n = np.asarray(counts.obs["n_counts"]).ravel()
        n[n == 0] = 1
        cpm = np.log1p(X.ravel() / n * 1e6)
        m["ESR1_tumour_cpm"] = cpm[tum_mask].mean()
    else:
        m["ESR1_tumour_cpm"] = np.nan
    contact, n_tum = tumour_immune_contact(ad)
    m["tumour_immune_contact"] = contact
    m["n_tumour_spots"] = n_tum
    return m


def program_morans(ad):
    """Moran's I of program scores per sample (grid graph)."""
    orig = np.asarray(ad.obsm["spatial"]).copy()
    ad.obsm["spatial"] = np.column_stack([ad.obs["array_row"], ad.obs["array_col"]])
    if "spatial_connectivities" not in ad.obsp:
        sq.gr.spatial_neighbors(ad, coord_type="grid", n_neighs=6)
    ad.obsm["spatial"] = orig
    cols = PROGRAM_COLS + ["rna_epithelial_tumour", "rna_T_cells"]
    cols = [c for c in cols if c in ad.obs]
    mor = sq.gr.spatial_autocorr(ad, mode="moran", attr="obs", genes=cols,
                                 n_perms=None, copy=True)
    return mor["I"]


def lineage_nhood(ad):
    """nhood enrichment on dominant RNA lineage (cell-organisation view)."""
    lab = "rna_dominant"
    if lab not in ad.obs:
        return None
    ad.obs["_dom_tmp"] = ad.obs[lab].astype("category")
    orig = np.asarray(ad.obsm["spatial"]).copy()
    ad.obsm["spatial"] = np.column_stack([ad.obs["array_row"], ad.obs["array_col"]])
    if "spatial_connectivities" not in ad.obsp:
        sq.gr.spatial_neighbors(ad, coord_type="grid", n_neighs=6)
    ad.obsm["spatial"] = orig
    try:
        sq.gr.nhood_enrichment(ad, cluster_key="_dom_tmp")
        res = ad.uns["_dom_tmp_nhood_enrichment"]
        cats = ad.obs["_dom_tmp"].cat.categories
        z = pd.DataFrame(np.asarray(res["zscore"]), index=cats, columns=cats)
        return z
    except Exception:
        return None


def exact_perm_paired(d):
    """Exact permutation test for paired mean (n=6); two-sided."""
    d = np.asarray(d, float)
    d = d[~np.isnan(d)]
    if len(d) < 3:
        return np.nan, np.nan
    n = len(d)
    obs = d.mean()
    s = 0
    total = 0
    # enumerate sign flips (2^n); exact two-sided for n<=6
    for bits in range(2 ** n):
        signs = np.array([1 if (bits >> i) & 1 else -1 for i in range(n)])
        if abs(float((signs * d).mean())) >= abs(float(obs)):
            s += 1
        total += 1
    p = 2 * s / total
    p = min(p, 1.0)
    return obs, p


def main():
    # ---- metric table for paired sections ----
    rows = []
    for spec in PAIRED:
        ad = sc.read_h5ad(os.path.join(SRC, f"{spec}_ecosystem.h5ad"))
        m = per_sample_metrics(ad, spec)
        m.update({"sample": spec, "patient": spec.split("-")[2],
                  "treatment": "pre" if "-S1" in spec else "post"})
        rows.append(m)
    tab = pd.DataFrame(rows).set_index("sample")
    tab.to_csv(os.path.join(TABLE_DIR, "07_paired_metrics.csv"))

    # ---- paired tests ----
    metric_cols = [c for c in tab.columns if c not in ("patient", "treatment")]
    tests = []
    for c in metric_cols:
        pre = tab.loc[tab.treatment == "pre", c].to_numpy()
        post = tab.loc[tab.treatment == "post", c].to_numpy()
        valid = ~np.isnan(pre) & ~np.isnan(post)
        pre, post = pre[valid], post[valid]
        if len(pre) < 3:
            continue
        d = post - pre
        obs, p_perm = exact_perm_paired(d)
        w, p_wil = stats.wilcoxon(pre, post) if len(pre) >= 5 else (np.nan, np.nan)
        pooled = np.sqrt((np.var(pre) + np.var(post)) / 2)
        es = (post.mean() - pre.mean()) / pooled if pooled > 0 else np.nan
        tests.append({"metric": c, "delta_mean": obs, "p_perm": p_perm,
                      "p_wilcoxon": p_wil, "effect_size": es,
                      "n": len(pre), "pre_mean": pre.mean(), "post_mean": post.mean()})
    tests = pd.DataFrame(tests)
    # BH correction across metrics
    tests["p_perm_BH"] = tests["p_perm"].apply(
        lambda p: min(1.0, p * len(tests) / max(1, (tests["p_perm"] <= p).sum())))
    tests = tests.sort_values("p_perm")
    tests.to_csv(os.path.join(TABLE_DIR, "07_paired_tests.csv"), index=False)
    print("\n[m7] paired pre->post tests (n=6 patients; permutation p, BH-corrected):")
    print(tests.round(3).to_string(index=False))

    # ---- figures ----
    panels = [("frac_Tumour", "Tumour fraction"), ("frac_Immune", "Immune fraction"),
              ("frac_Stroma", "Stroma fraction"),
              ("score_Estrogen_response", "Estrogen-response program"),
              ("score_Proliferation", "Proliferation program"),
              ("score_Interferon_response", "Interferon program"),
              ("tumour_immune_contact", "Tumour–immune contact"),
              ("ESR1_tumour_cpm", "ESR1 in tumour (log CPM)")]
    panels = [(c, t) for c, t in panels if c in tab.columns]
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    for ax, (c, t) in zip(axes.ravel(), panels):
        pre_pts = tab.loc[tab.treatment == "pre", c]
        post_pts = tab.loc[tab.treatment == "post", c]
        for pt in set(tab["patient"]):
            vals = tab.loc[tab.patient == pt, c]
            if len(vals) == 2:
                ax.plot([0, 1], vals.values, "o-", color="grey", alpha=0.6)
        ax.scatter(np.zeros(len(pre_pts)), pre_pts.values, color="#1f77b4", s=30, zorder=3)
        ax.scatter(np.ones(len(post_pts)), post_pts.values, color="#d62728", s=30, zorder=3)
        r = tests[tests.metric == c]
        if len(r):
            star = "***" if r.p_perm_BH.iloc[0] < 0.001 else "**" if r.p_perm_BH.iloc[0] < 0.01 \
                else "*" if r.p_perm_BH.iloc[0] < 0.05 else "ns"
            ax.set_title(f"{t}\nperm p={r.p_perm.iloc[0]:.3f} ({star})", fontsize=9)
        else:
            ax.set_title(t, fontsize=9)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["pre", "post"])
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Paired pre vs post primary endocrine therapy (6 patients)", y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "07_paired_metrics.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # ---- spatial maps for an illustrative pair (P1 and P4) ----
    for pt in ["P1", "P4"]:
        fig, axes = plt.subplots(2, 4, figsize=(19, 8.5))
        specs = [s for s in PAIRED if s.split("-")[2] == pt]
        for axr, spec in zip(axes, specs):
            ad = sc.read_h5ad(os.path.join(SRC, f"{spec}_ecosystem.h5ad"))
            lib = ad.uns.get("library_id", "visium")
            sq.pl.spatial_scatter(ad, color="compartment", ax=axr[0], img=True,
                                  library_id=lib, frameon=False, size=3,
                                  title=f"{spec} compartments")
            for ax, col, t in zip(axr[1:], ["score_Estrogen_response",
                                            "score_Proliferation",
                                            "score_T_cell_activity"],
                                  ["Estrogen", "Proliferation", "T-cell activity"]):
                sq.pl.spatial_scatter(ad, color=col, ax=ax, img=True, library_id=lib,
                                      frameon=False, size=3, use_raw=False,
                                      colorbar=False, title=t)
        fig.suptitle(f"Patient {pt}: pre (top) vs post endocrine therapy (bottom)", y=1.0)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, f"07_spatial_prepost_{pt}.png"), dpi=150,
                    bbox_inches="tight")
        plt.close(fig)

    # ---- neighbourhood enrichment: post - pre ----
    nh_tabs = {}
    for spec in PAIRED:
        f = os.path.join(TABLE_DIR, f"05_nhood_{spec}.csv")
        if os.path.exists(f):
            nh_tabs[spec] = pd.read_csv(f, index_col=0)
    cats = ["Tumour", "Immune", "Stroma", "Mixed"]
    delta = pd.DataFrame(np.nan, index=cats, columns=cats)
    counts = pd.DataFrame(0, index=cats, columns=cats)
    for pt in set(tab["patient"]):
        pre_s = [s for s in PAIRED if s.split("-")[2] == pt and "-S1" in s
                 and s in nh_tabs]
        post_s = [s for s in PAIRED if s.split("-")[2] == pt and "-S2" in s
                  and s in nh_tabs]
        if not pre_s or not post_s:
            continue
        for a in cats:
            for b in cats:
                if a in nh_tabs[pre_s[0]].columns and b in nh_tabs[pre_s[0]].index:
                    v_pre = nh_tabs[pre_s[0]].loc[b, a] if (b in nh_tabs[pre_s[0]].index and a in nh_tabs[pre_s[0]].columns) else np.nan
                    v_post = nh_tabs[post_s[0]].loc[b, a] if (b in nh_tabs[post_s[0]].index and a in nh_tabs[post_s[0]].columns) else np.nan
                    if not np.isnan(v_pre) and not np.isnan(v_post):
                        delta.loc[a, b] = np.nansum([delta.loc[a, b], v_post - v_pre]) if not np.isnan(delta.loc[a, b]) else v_post - v_pre
                        counts.loc[a, b] += 1
    delta = delta / counts.replace(0, np.nan)
    delta.to_csv(os.path.join(TABLE_DIR, "07_nhood_delta_post_pre.csv"))
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(delta.values, cmap="RdBu_r", vmin=-np.nanmax(np.abs(delta.values)),
                   vmax=np.nanmax(np.abs(delta.values)))
    ax.set_xticks(range(len(cats))); ax.set_yticks(range(len(cats)))
    ax.set_xticklabels(cats); ax.set_yticklabels(cats)
    for i in range(len(cats)):
        for j in range(len(cats)):
            v = delta.values[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:+.1f}", ha="center", va="center",
                        color="white" if abs(v) > 5 else "black")
    ax.set_xlabel("neighbour compartment"); ax.set_ylabel("central compartment")
    ax.set_title("Δ neighbourhood enrichment z (post − pre)")
    fig.colorbar(im, ax=ax, shrink=0.8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "07_nhood_delta.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- Moran's I of programs: pre vs post ----
    mor_rows = []
    for spec in PAIRED:
        ad = sc.read_h5ad(os.path.join(SRC, f"{spec}_ecosystem.h5ad"))
        mi = program_morans(ad)
        for g, v in mi.items():
            mor_rows.append({"sample": spec, "patient": spec.split("-")[2],
                             "treatment": "pre" if "-S1" in spec else "post",
                             "program": g, "moranI": v})
    mors = pd.DataFrame(mor_rows)
    mors.to_csv(os.path.join(TABLE_DIR, "07_program_morans.csv"), index=False)
    piv = mors.pivot_table(index="program", columns="treatment", values="moranI",
                           aggfunc="mean")
    piv["delta"] = piv.get("post", 0) - piv.get("pre", 0)
    print("\n[m7] mean Moran's I of programs (pre vs post):")
    print(piv.round(3).to_string())

    # paired permutation tests on Moran's I deltas
    rows_t = []
    for prog in mors["program"].unique():
        pre = mors[(mors.program == prog) & (mors.treatment == "pre")].set_index("patient")["moranI"]
        post = mors[(mors.program == prog) & (mors.treatment == "post")].set_index("patient")["moranI"]
        d = (post - pre).dropna()
        if len(d) < 3:
            continue
        obs, p = exact_perm_paired(d.to_numpy())
        rows_t.append({"program": prog, "delta_moranI": obs, "p_perm": p,
                       "n": len(d)})
    mt = pd.DataFrame(rows_t).sort_values("p_perm")
    mt["p_BH"] = mt["p_perm"].apply(
        lambda p: min(1.0, p * len(mt) / max(1, (mt["p_perm"] <= p).sum())))
    mt.to_csv(os.path.join(TABLE_DIR, "07_program_morans_paired_tests.csv"), index=False)
    print("\n[m7] paired tests on spatial-structure (Moran's I) change post-pre:")
    print(mt.round(3).to_string(index=False))


if __name__ == "__main__":
    main()