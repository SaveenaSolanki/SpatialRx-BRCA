"""Leakage-free benchmark for Gate 4 (Module 08 reconstruction).

Design is pre-registered in ``benchmark/BENCHMARK_CONTRACT.md``.  Key properties:

* patient-held-out (primary) and section-held-out (secondary) splits;
* HVG selection, scaling, PCA, protein scaling and every model are fit on the
  training folds only;
* all methods are evaluated on identical held-out partitions;
* negative controls: random projection, degree-preserving graph randomisation,
  and a label-shuffled probe.

The module exposes small, testable primitives; the runner lives in
``benchmark/run_benchmark.py``.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import sparse

from .config import PATHS, Paths

# 11 curated programs scored in module 05 (order matters for labels)
PROGRAMS = [
    "Estrogen_response", "Proliferation", "EMT", "Hypoxia",
    "Interferon_response", "T_cell_activity", "Cytotoxicity",
    "Myeloid_inflammation", "Fibroblast_activation", "Angiogenesis",
    "B_cell_activity",
]

METHODS = ["pca", "smoothing", "random", "rna_ae", "rna_protein_ae",
           "graph_ae", "graph_ae_randomized"]
STOCHASTIC = {"random", "rna_ae", "rna_protein_ae", "graph_ae", "graph_ae_randomized"}

# Training budget. Declared deviation from the contract's 128-hidden / 200-epoch
# specification to keep the benchmark tractable on the shared, heavily loaded
# host (see BENCHMARK_CONTRACT.md addendum). Identical for every learned method.
DEFAULT_EPOCHS = 50
DEFAULT_HIDDEN = 64


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
@dataclass
class SectionData:
    specimen: str
    patient: str
    treatment: str
    rna: sparse.csr_matrix          # spots x common genes, log-normalised
    genes: list[str]
    protein: np.ndarray             # spots x 35 (CLR)
    programs: np.ndarray            # spots x 11
    compartment: np.ndarray         # spots (str)
    rows: np.ndarray
    cols: np.ndarray
    adj: sparse.csr_matrix | None = None


def _align_genes(sections: dict[str, SectionData]) -> None:
    """Reindex every section to the common gene set (order of the first section)."""
    gene_sets = [set(s.genes) for s in sections.values()]
    common = set.intersection(*gene_sets)
    ref = [g for g in next(iter(sections.values())).genes if g in common]
    index = {g: i for i, g in enumerate(ref)}
    for s in sections.values():
        if s.genes == ref:
            s.genes = ref
            continue
        pos = np.array([s.genes.index(g) for g in ref])
        s.rna = s.rna[:, pos].tocsr()
        s.genes = ref
    _ = index


def build_graph(rows: np.ndarray, cols: np.ndarray, n_neighs: int = 6) -> sparse.csr_matrix:
    import anndata as ad
    import squidpy as sq

    tmp = ad.AnnData(X=np.zeros((len(rows), 1), dtype=np.float32))
    tmp.obsm["spatial"] = np.column_stack([rows, cols]).astype(float)
    sq.gr.spatial_neighbors(tmp, coord_type="grid", n_neighs=n_neighs)
    return tmp.obsp["spatial_connectivities"].tocsr()


def load_sections(paths: Paths | None = None, specimens: Iterable[str] | None = None) -> dict[str, SectionData]:
    from .data import sample_meta
    from .io import read_h5ad

    paths = paths or PATHS
    meta = sample_meta().set_index("specimen")
    specs = list(specimens) if specimens is not None else list(meta.index)
    out: dict[str, SectionData] = {}
    for spec in specs:
        dom = read_h5ad(paths.domains_dir() / f"{spec}_domains.h5ad")
        eco = read_h5ad(paths.ecosystem_dir() / f"{spec}_ecosystem.h5ad")[dom.obs_names]
        raw = dom.raw
        genes = list(raw.var_names)
        rna = sparse.csr_matrix(raw[:, :].X)
        prog = eco.obs[[f"score_{p}" for p in PROGRAMS]].to_numpy(dtype=np.float32)
        sec = SectionData(
            specimen=spec,
            patient=str(meta.loc[spec, "patient"]),
            treatment=str(meta.loc[spec, "treatment"]),
            rna=rna,
            genes=genes,
            protein=np.asarray(eco.obsm["protein_clr"], dtype=np.float32),
            programs=prog,
            compartment=eco.obs["compartment"].astype(str).to_numpy(),
            rows=dom.obs["array_row"].to_numpy(float),
            cols=dom.obs["array_col"].to_numpy(float),
        )
        sec.adj = build_graph(sec.rows, sec.cols)
        out[spec] = sec
    _align_genes(out)
    return out


# ---------------------------------------------------------------------------
# Splits
# ---------------------------------------------------------------------------
@dataclass
class Fold:
    name: str
    analysis: str
    train_sections: list[str]
    test_sections: list[str]
    train_patients: list[str]
    test_patients: list[str]


def make_folds(sections: dict[str, SectionData], analysis: str) -> list[Fold]:
    by_patient: dict[str, list[str]] = {}
    for spec, s in sections.items():
        by_patient.setdefault(s.patient, []).append(spec)
    folds: list[Fold] = []
    if analysis == "patient":
        for pat, specs in sorted(by_patient.items()):
            tr = [x for p, xs in by_patient.items() if p != pat for x in xs]
            folds.append(Fold(f"LOPO_{pat}", analysis, tr, sorted(specs), 
                              sorted(p for p in by_patient if p != pat), [pat]))
    elif analysis == "section":
        for spec in sorted(sections):
            tr = [x for x in sections if x != spec]
            folds.append(Fold(f"LOSO_{spec}", analysis, tr, [spec],
                              sorted({sections[x].patient for x in tr}),
                              [sections[spec].patient]))
    else:
        raise ValueError(f"unknown analysis {analysis!r}")
    return folds


def assert_disjoint(folds: list[Fold]) -> None:
    if folds and folds[0].analysis == "patient":
        for f in folds:
            assert not (set(f.train_patients) & set(f.test_patients)), f"patient leakage in {f.name}"
    for f in folds:
        assert not (set(f.train_sections) & set(f.test_sections)), f"section leakage in {f.name}"


# ---------------------------------------------------------------------------
# Preprocessing (fit on train only)
# ---------------------------------------------------------------------------
class Preprocessor:
    def __init__(self, n_hvg: int = 2000, n_pcs: int = 64):
        self.n_hvg = n_hvg
        self.n_pcs = n_pcs

    def fit(self, train: list[SectionData]):
        import anndata as ad
        import scanpy as sc
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler

        X = sparse.vstack([s.rna for s in train], format="csr")
        a = ad.AnnData(X=X)
        a.var_names = list(train[0].genes)
        sc.pp.highly_variable_genes(a, n_top_genes=self.n_hvg, flavor="seurat", n_bins=20)
        self.hvg_idx = np.where(a.var["highly_variable"].to_numpy())[0]
        Xd = np.asarray(a[:, self.hvg_idx].X.todense(), dtype=np.float32)
        self.rna_scaler = StandardScaler().fit(Xd)
        Xs = np.clip(self.rna_scaler.transform(Xd), -10, 10).astype(np.float32)
        self.pca = PCA(n_components=self.n_pcs, random_state=0).fit(Xs)
        prot = np.vstack([s.protein for s in train])
        self.prot_scaler = StandardScaler().fit(prot)
        return self

    def transform(self, s: SectionData):
        Xd = np.asarray(s.rna[:, self.hvg_idx].todense(), dtype=np.float32)
        Xs = np.clip(self.rna_scaler.transform(Xd), -10, 10).astype(np.float32)
        rna = self.pca.transform(Xs).astype(np.float32)
        prot = self.prot_scaler.transform(s.protein).astype(np.float32)
        return rna, prot


def featurize(prep: Preprocessor, sections: list[SectionData]):
    rna_list, prot_list = [], []
    for s in sections:
        r, p = prep.transform(s)
        rna_list.append(r)
        prot_list.append(p)
    rna = np.vstack(rna_list)
    prot = np.vstack(prot_list)
    return {"rna": rna, "prot": prot, "full": np.hstack([rna, prot])}


# ---------------------------------------------------------------------------
# Graphs
# ---------------------------------------------------------------------------
def normalize_adj(A: sparse.csr_matrix) -> sparse.csr_matrix:
    A = A.tocsr().astype(np.float32)
    d = np.asarray(A.sum(axis=1)).ravel()
    d[d == 0] = 1.0
    Dinv = sparse.diags(1.0 / d)
    return (Dinv @ A).tocsr()


def block_diag_graph(adjs: list[sparse.csr_matrix]) -> sparse.csr_matrix:
    return sparse.block_diag(adjs, format="csr")


def randomize_graph_degree_preserving(A: sparse.csr_matrix, seed: int) -> sparse.csr_matrix:
    """Degree-preserving double-edge swaps (igraph), returns a symmetric CSR."""
    import igraph as ig

    A = A.tocoo()
    mask = A.row < A.col
    edges = list(zip(A.row[mask].tolist(), A.col[mask].tolist()))
    n = A.shape[0]
    g = ig.Graph(n=n, edges=edges, directed=False)
    if g.ecount() > 0:
        g.rewire(n=max(1, 5 * g.ecount()), mode="simple")
    el = np.array(g.get_edgelist())
    if el.size == 0:
        return sparse.csr_matrix((n, n), dtype=np.float32)
    rows = np.concatenate([el[:, 0], el[:, 1]])
    cols = np.concatenate([el[:, 1], el[:, 0]])
    data = np.ones(rows.size, dtype=np.float32)
    return sparse.csr_matrix((data, (rows, cols)), shape=(n, n))


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------
def _torch():
    import torch
    return torch


class AE:
    """Small MLP autoencoder (numpy wrapper around torch)."""

    def __init__(self, d_in: int, d_hid: int = DEFAULT_HIDDEN, d_z: int = 32,
                 seed: int = 0, epochs: int | None = None, lr: float = 1e-3):
        import torch.nn as nn
        torch = _torch()
        torch.manual_seed(seed)
        self.net = nn.Sequential(nn.Linear(d_in, d_hid), nn.ReLU(),
                                 nn.Linear(d_hid, d_z))
        self.dec = nn.Sequential(nn.Linear(d_z, d_hid), nn.ReLU(),
                                 nn.Linear(d_hid, d_in))
        self.epochs = DEFAULT_EPOCHS if epochs is None else epochs
        self.lr = lr

    def fit(self, X: np.ndarray):
        torch = _torch()
        Xt = torch.tensor(X, dtype=torch.float32)  # convert once
        opt = torch.optim.Adam(list(self.net.parameters()) + list(self.dec.parameters()),
                               lr=self.lr)
        for _ in range(self.epochs):
            opt.zero_grad()
            z = self.net(Xt)
            loss = ((self.dec(z) - Xt) ** 2).mean()
            loss.backward()
            opt.step()
        self.loss_ = float(loss.item())
        return self

    def embed(self, X: np.ndarray) -> np.ndarray:
        torch = _torch()
        with torch.no_grad():
            z = self.net(torch.tensor(X, dtype=torch.float32))
        return z.numpy()


def _sparse_torch(Ahat: sparse.csr_matrix):
    torch = _torch()
    A = Ahat.tocoo()
    idx = torch.tensor(np.vstack([A.row, A.col]), dtype=torch.long)
    val = torch.tensor(A.data, dtype=torch.float32)
    return torch.sparse_coo_tensor(idx, val, A.shape).coalesce()


class GCNAE(AE):
    """Graph-conv pre-layer (A_hat X W) then MLP autoencoder."""

    def __init__(self, d_in: int, seed: int = 0, **kw):
        super().__init__(d_in, seed=seed, **kw)
        torch = _torch()
        import torch.nn as nn
        torch.manual_seed(seed)
        self.w = nn.Linear(d_in, d_in, bias=False)

    def _set_graph(self, Ahat):
        self._At = _sparse_torch(Ahat)

    def _smooth_t(self, Xt):
        torch = _torch()
        return torch.sparse.mm(self._At, torch.relu(self.w(Xt)))

    def fit(self, X, Ahat):
        torch = _torch()
        self._set_graph(Ahat)
        Xt = torch.tensor(X, dtype=torch.float32)  # convert once
        opt = torch.optim.Adam(list(self.w.parameters()) + list(self.net.parameters())
                               + list(self.dec.parameters()), lr=self.lr)
        for _ in range(self.epochs):
            opt.zero_grad()
            h = self._smooth_t(Xt)
            z = self.net(h)
            loss = ((self.dec(z) - h) ** 2).mean()
            loss.backward()
            opt.step()
        self.loss_ = float(loss.item())
        return self

    def embed(self, X, Ahat):
        torch = _torch()
        self._set_graph(Ahat)
        with torch.no_grad():
            h = self._smooth_t(torch.tensor(X, dtype=torch.float32))
        return super().embed(h.numpy())


def fit_method(method: str, feats: dict, adjs: list[sparse.csr_matrix], seed: int):
    if method == "pca":
        return None
    if method == "smoothing":
        return None
    if method == "random":
        rng = np.random.default_rng(seed)
        return rng.normal(size=(feats["rna"].shape[1], 32)).astype(np.float32)
    if method == "rna_ae":
        return AE(feats["rna"].shape[1], seed=seed).fit(feats["rna"])
    if method == "rna_protein_ae":
        return AE(feats["full"].shape[1], seed=seed).fit(feats["full"])
    if method == "graph_ae":
        A = normalize_adj(block_diag_graph(adjs))
        return GCNAE(feats["full"].shape[1], seed=seed).fit(feats["full"], A)
    if method == "graph_ae_randomized":
        Ar = normalize_adj(block_diag_graph(
            [randomize_graph_degree_preserving(a, seed) for a in adjs]))
        return {"gcn": GCNAE(feats["full"].shape[1], seed=seed).fit(feats["full"], Ar),
                "seed": seed}
    raise ValueError(method)


def embed_method(method: str, model, feats: dict, adjs: list[sparse.csr_matrix]) -> np.ndarray:
    if method == "pca":
        return feats["rna"][:, :32]
    if method == "smoothing":
        # non-learned operator: smooth each split with its own within-section graph
        return (normalize_adj(block_diag_graph(adjs)) @ feats["rna"])[:, :32]
    if method == "random":
        return feats["rna"] @ model
    if method == "rna_ae":
        return model.embed(feats["rna"])
    if method == "rna_protein_ae":
        return model.embed(feats["full"])
    if method == "graph_ae":
        A = normalize_adj(block_diag_graph(adjs))
        return model.embed(feats["full"], A)
    if method == "graph_ae_randomized":
        Ar = normalize_adj(block_diag_graph(
            [randomize_graph_degree_preserving(a, model["seed"]) for a in adjs]))
        return model["gcn"].embed(feats["full"], Ar)
    raise ValueError(method)


# ---------------------------------------------------------------------------
# Labels, probe, metrics
# ---------------------------------------------------------------------------
def binarize_stack(mats: list[np.ndarray], q: int = 80) -> np.ndarray:
    out = []
    for m in mats:
        y = np.zeros(m.shape, dtype=np.int64)
        for j in range(m.shape[1]):
            col = m[:, j]
            thr = np.percentile(col, q)
            y[:, j] = (col > thr).astype(np.int64)
        out.append(y)
    return np.vstack(out)


def probe_auroc(z_tr, y_tr, z_te, y_te, C: float = 1.0) -> np.ndarray:
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score

    aucs = []
    for j in range(y_tr.shape[1]):
        if y_tr[:, j].min() == y_tr[:, j].max() or y_te[:, j].min() == y_te[:, j].max():
            aucs.append(np.nan)
            continue
        clf = LogisticRegression(C=C, max_iter=1000)
        clf.fit(z_tr, y_tr[:, j])
        aucs.append(roc_auc_score(y_te[:, j], clf.decision_function(z_te)))
    return np.asarray(aucs, dtype=float)


def select_C(z_tr, y_tr, groups, Cs=(0.1, 1.0)) -> float:
    """Inner leave-one-group-out selection of the probe regularisation."""
    best_C, best_score = Cs[0], -np.inf
    for C in Cs:
        scores = []
        for g in np.unique(groups):
            te = groups == g
            tr = ~te
            if y_tr[tr].sum() == 0:
                continue
            scores.append(np.nanmean(probe_auroc(z_tr[tr], y_tr[tr], z_tr[te], y_tr[te], C)))
        score = np.nanmean(scores) if scores else -np.inf
        if score > best_score:
            best_C, best_score = C, score
    return best_C


def label_shuffle_auroc(z_tr, y_tr, z_te, y_te, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    y_perm = y_tr.copy()
    for j in range(y_perm.shape[1]):
        y_perm[:, j] = rng.permutation(y_perm[:, j])
    return probe_auroc(z_tr, y_perm, z_te, y_te)


def silhouette_per_section(z, comp, sections_index) -> float:
    from sklearn.metrics import silhouette_score

    vals = []
    for zz, cc in zip(z, comp):
        if len(set(cc)) > 1 and min(Counter(cc).values()) > 1:
            vals.append(silhouette_score(zz, cc))
    return float(np.nanmean(vals)) if vals else np.nan


def moran_per_section(z, adjs) -> float:
    vals = []
    for zz, A in zip(z, adjs):
        W = normalize_adj(A)
        wsum = W.sum()
        dims = []
        for j in range(zz.shape[1]):
            x = zz[:, j]
            xm = x - x.mean()
            den = (xm ** 2).sum() / len(x)
            if den <= 0:
                continue
            dims.append((xm @ (W @ xm)) / wsum / den)
        vals.append(np.nanmean(dims))
    return float(np.nanmean(vals)) if vals else np.nan


def split_sections(z: np.ndarray, sections: list[SectionData]):
    out, i = [], 0
    for s in sections:
        n = s.rna.shape[0]
        out.append((z[i:i + n], s.compartment, s.adj))
        i += n
    return out


def seed_correlation(z0: np.ndarray, z: np.ndarray) -> float:
    """Mean absolute Pearson correlation between matched latent dimensions."""
    if z0.shape != z.shape:
        return np.nan
    cors = []
    for j in range(z0.shape[1]):
        a, b = z0[:, j], z[:, j]
        if a.std() == 0 or b.std() == 0:
            continue
        cors.append(abs(np.corrcoef(a, b)[0, 1]))
    return float(np.nanmean(cors)) if cors else np.nan
