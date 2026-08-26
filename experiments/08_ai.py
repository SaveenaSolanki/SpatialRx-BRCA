"""
Module 08 — Multimodal spatial representation (small, transparent models).
Question: does adding protein and spatial-graph context recover biologically
coherent spatial tumour states better than RNA alone?

Per sample (12 paired sections, replicate sections included as extra robustness):
  features:
    RNA     = log-normalized, PCA(64) of the top genes
    protein = CLR-normalized antibody panel, PCA(16)
    graph   = grid adjacency (row-normalised)
  models (self-supervised 32-d embeddings, MSE reconstruction):
    M1 RNA only (MLP AE)
    M2 RNA + protein (MLP AE)
    M3 RNA + protein + spatial graph (GCN-style: 1 graph-conv pre-layer + MLP AE)
  metrics (within-sample, averaged across samples):
    - silhouette vs annotated compartment
    - spatial coherence = mean Moran's I of embedding dimensions
    - biological alignment = mean |corr| between embedding dims and 11 gene programs
    - niche purity (AUROC-style): separation of program-high spots by embedding
Honest framing: NOT a response-prediction model (n=6). Compare configs, report CI.

Usage: conda run -n spatialrx python experiments/08_ai.py
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
import torch
import torch.nn as nn
from scipy.sparse import csr_matrix, issparse
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from srx import PROC_DIR, FIG_DIR, TABLE_DIR, SAMPLE_META, PROGRAMS

sc.settings.verbosity = 0
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)
SRC = os.path.join(PROC_DIR, "03_domains")
ECO = os.path.join(PROC_DIR, "05_ecosystem")
DEVICE = "cpu"
torch.manual_seed(0)
np.random.seed(0)

PROGRAM_COLS = [f"score_{p}" for p in PROGRAMS]


class AE(nn.Module):
    def __init__(self, d_in, d_hid=128, d_z=32):
        super().__init__()
        self.enc = nn.Sequential(nn.Linear(d_in, d_hid), nn.ReLU(),
                                 nn.Linear(d_hid, d_z))
        self.dec = nn.Sequential(nn.Linear(d_z, d_hid), nn.ReLU(),
                                 nn.Linear(d_hid, d_in))

    def forward(self, x):
        z = self.enc(x)
        return z, self.dec(z)


class GCN_AE(nn.Module):
    """Graph-conv prelayer (A_hat X W) then MLP AE on the smoothed features."""

    def __init__(self, d_in, d_hid=128, d_z=32):
        super().__init__()
        self.w = nn.Linear(d_in, d_in, bias=False)  # learned feature transform
        self.ae = AE(d_in, d_hid, d_z)

    def forward(self, x, ahat):
        h = torch.relu(self.w(x))
        h = ahat @ h
        return self.ae(h)


def embed(model, X, A=None):
    model.eval()
    with torch.no_grad():
        Xt = torch.tensor(np.asarray(X, dtype=np.float32))
        if A is not None:
            At = torch.tensor(np.asarray(A.todense(), dtype=np.float32))
            z, _ = model(Xt, At)
        else:
            z, _ = model(Xt)
    return z.numpy()


def fit(model, X, A=None, epochs=300):
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    Xt = torch.tensor(np.asarray(X, dtype=np.float32))
    Xt = Xt.to(DEVICE)
    A_t = torch.tensor(np.asarray(A.todense(), dtype=np.float32)) if A is not None else None
    model.to(DEVICE)
    for ep in range(epochs):
        opt.zero_grad()
        if A_t is not None:
            z, rec = model(Xt, A_t)
        else:
            z, rec = model(Xt)
        loss = ((rec - Xt) ** 2).mean()
        loss.backward()
        opt.step()
    return float(loss.item())


def graph_adj(ad):
    """Row-normalised grid adjacency from array coords."""
    orig = np.asarray(ad.obsm["spatial"]).copy()
    ad.obsm["spatial"] = np.column_stack([ad.obs["array_row"], ad.obs["array_col"]])
    if "spatial_connectivities" not in ad.obsp:
        import squidpy as sq
        sq.gr.spatial_neighbors(ad, coord_type="grid", n_neighs=6)
    ad.obsm["spatial"] = orig
    A = ad.obsp["spatial_connectivities"].tocsr().astype(np.float32)
    d = np.asarray(A.sum(axis=1)).ravel()
    d[d == 0] = 1
    Dinv = csr_matrix((1.0 / d, (np.arange(A.shape[0]), np.arange(A.shape[0]))))
    Ahat = (Dinv @ A).tocsr()
    return Ahat


def morans_obsm(z, Ahat):
    """Closed-form Moran's I per embedding dimension (row-normalised W)."""
    z = np.asarray(z, float)
    n = z.shape[0]
    W = Ahat.tocsr()
    out = []
    for j in range(z.shape[1]):
        col = z[:, j]
        zm = col - col.mean()
        num = (zm @ (W @ zm)) / W.sum()
        den = (zm ** 2).sum() / n
        out.append(num / den if den > 0 else np.nan)
    return np.nanmean(out)


def sample_features(ad):
    # RNA: lognorm PCA
    X = ad.X
    if issparse(X):
        X = X.toarray()
    rna_pca = np.asarray(ad.obsm.get("X_pca", np.eye(0)))[:, :64] if "X_pca" in ad.obsm else None
    if rna_pca is None or rna_pca.shape[0] != ad.n_obs:
        # quick PCA from scratch (module-03 objects already have X_pca usually)
        mu = X.mean(0)
        Xc = X - mu
        u, s, vt = np.linalg.svd(Xc, full_matrices=False)
        rna_pca = (Xc @ vt[:64].T) / (s[:64] + 1e-9)
    prot = ad.obsm["protein_clr"] if "protein_clr" in ad.obsm else None
    if prot is not None:
        pm = prot.mean(0)
        pc = prot - pm
        u2, s2, vt2 = np.linalg.svd(pc, full_matrices=False)
        prot_pca = (pc @ vt2[:16].T) / (s2[:16] + 1e-9)
    else:
        prot_pca = None
    return rna_pca.astype(np.float32), prot_pca


def program_matrix(ad):
    prog = np.column_stack([np.asarray(ad.obs[c]) for c in PROGRAM_COLS])
    return prog


def niche_separation(z, prog):
    """Program-high vs low spots separation: mean AUROC over programs
    (per-spot embedding -> linear score)."""
    from sklearn.linear_model import LogisticRegression
    aucs = []
    for j in range(prog.shape[1]):
        y = prog[:, j] > np.percentile(prog[:, j], 80)
        if y.sum() < 5 or (~y).sum() < 5:
            continue
        clf = LogisticRegression(max_iter=500)
        clf.fit(z, y)
        aucs.append(roc_auc_score(y, clf.decision_function(z)))
    return float(np.mean(aucs)) if aucs else np.nan


def run_sample(specimen):
    ad = sc.read_h5ad(os.path.join(SRC, f"{specimen}_domains.h5ad"))
    eco = sc.read_h5ad(os.path.join(ECO, f"{specimen}_ecosystem.h5ad"))
    rna, prot = sample_features(ad)
    Ahat = graph_adj(ad)
    comp = np.asarray(eco.obs["compartment"])
    prog = program_matrix(eco)

    X2 = np.hstack([rna, prot]) if prot is not None else rna
    out = {"sample": specimen}

    emb = {}
    # M1 RNA
    m = AE(rna.shape[1]).to(DEVICE)
    fit(m, rna, epochs=300)
    emb["RNA"] = embed(m, rna)
    # M2 RNA+protein
    m2 = AE(X2.shape[1]).to(DEVICE)
    fit(m2, X2, epochs=300)
    emb["RNA+protein"] = embed(m2, X2)
    # M3 RNA+protein+graph
    m3 = GCN_AE(X2.shape[1]).to(DEVICE)
    fit(m3, X2, A=Ahat, epochs=300)
    emb["RNA+protein+graph"] = embed(m3, X2, A=Ahat)

    for name, z in emb.items():
        sil = silhouette_score(z, comp) if len(set(comp)) > 1 else np.nan
        mor = morans_obsm(z, Ahat)
        mark = niche_separation(z, prog)
        out[f"silhouette_{name}"] = sil
        out[f"moran_{name}"] = mor
        out[f"niche_{name}"] = mark
    out["n_spots"] = ad.n_obs
    return out, emb, ad


def main():
    rows, all_emb = [], {}
    for _, m in SAMPLE_META.iterrows():
        r, emb, ad = run_sample(m.specimen)
        rows.append(r)
        all_emb[m.specimen] = (emb, ad)
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(TABLE_DIR, "08_embedding_metrics.csv"), index=False)

    print("\n[m8] embedding quality metrics (per-sample, then mean +/- SD):")
    for key in ["silhouette_RNA", "silhouette_RNA+protein", "silhouette_RNA+protein+graph",
                "moran_RNA", "moran_RNA+protein", "moran_RNA+protein+graph",
                "niche_RNA", "niche_RNA+protein", "niche_RNA+protein+graph"]:
        v = tab[key].dropna()
        print(f"  {key:38s}: {v.mean():.3f} +/- {v.std():.3f}  (n={len(v)})")

    # paired per-sample comparison
    print("\n[m8] per-sample deltas (multimodal - RNA-only):")
    for key in ["silhouette", "moran", "niche"]:
        d = tab[f"{key}_RNA+protein+graph"] - tab[f"{key}_RNA"]
        print(f"  {key:12s}: mean delta {d.mean():+.3f} +/- {d.std():.3f} "
              f"(+fixed win count {int((d > 0).sum())}/{len(d)})")

    # combine table to long form for the figure
    long = pd.melt(tab, id_vars=["sample", "n_spots"],
                   value_vars=[c for c in tab.columns if c.startswith((
                       "silhouette_", "moran_", "niche_"))],
                   var_name="metric_config", value_name="value")
    long["metric"] = long["metric_config"].str.split("_").str[0]
    long["config"] = long["metric_config"].str.split("_", n=1).str[1]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for ax, met in zip(axes, ["silhouette", "moran", "niche"]):
        sub = long[long.metric == met]
        order = ["RNA", "RNA+protein", "RNA+protein+graph"]
        dat = [sub.loc[sub.config == c, "value"].values for c in order]
        bp = ax.boxplot(dat, vert=True, showfliers=False, patch_artist=True)
        ax.set_xticklabels(order, rotation=15)
        for patch in bp["boxes"]:
            patch.set_facecolor("#a0c8e0")
        ax.set_title({"silhouette": "Silhouette vs compartments",
                      "moran": "Spatial coherence (mean Moran's I of dims)",
                      "niche": "Biological niche separation (AUROC)"}[met])
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "08_embedding_metrics.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # UMAP of the three embeddings for the flagship sample
    emb, ad = all_emb["HCC22-088-P4-S2"]
    eco = sc.read_h5ad(os.path.join(ECO, "HCC22-088-P4-S2_ecosystem.h5ad"))
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    for ax, (name, z) in zip(axes, emb.items()):
        u, s, vt = np.linalg.svd(z - z.mean(0), full_matrices=False)
        um = (z - z.mean(0)) @ vt[:2].T
        sc_ = ax.scatter(um[:, 0], um[:, 1], c=eco.obs["compartment"].map(
            {"Tumour": 0, "Immune": 1, "Stroma": 2, "Mixed": 3}).to_numpy(),
            s=3, cmap="tab10", vmax=4)
        ax.set_title(name)
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("Embedding space (P4-S2), colored by compartment", y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "08_embedding_umap_P4S2.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()