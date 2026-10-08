"""Tests for the Gate-4 benchmark primitives (leakage-safety and controls)."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import sparse

from spatialrx import benchmark as B


def _dummy_sections(patients):
    """sections dict keyed by name with only a .patient attribute (for splits)."""
    out = {}
    for i, p in enumerate(patients):
        out[f"S{i}"] = SimpleNamespace(patient=p)
    return out


# ---------------------------------------------------------------------------
# Splits
# ---------------------------------------------------------------------------
def test_patient_folds_are_patient_disjoint():
    secs = _dummy_sections(["P1", "P1", "P2", "P2", "P3", "P3"])
    folds = B.make_folds(secs, "patient")
    assert len(folds) == 3
    B.assert_disjoint(folds)
    for f in folds:
        assert not (set(f.train_patients) & set(f.test_patients))
        assert set(f.test_patients) == {secs[f.test_sections[0]].patient}


def test_section_folds_hold_out_one_section():
    secs = _dummy_sections(["P1", "P1", "P2", "P2"])
    folds = B.make_folds(secs, "section")
    assert len(folds) == 4
    B.assert_disjoint(folds)
    for f in folds:
        assert len(f.test_sections) == 1


# ---------------------------------------------------------------------------
# Graph controls
# ---------------------------------------------------------------------------
def _path_graph(n=6):
    rows, cols = [], []
    for i in range(n - 1):
        rows += [i, i + 1]
        cols += [i + 1, i]
    return sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))


def test_randomize_graph_preserves_degree_sequence():
    A = _path_graph(8)
    Ar = B.randomize_graph_degree_preserving(A, seed=0)
    d0 = np.sort(np.asarray(A.sum(1)).ravel())
    d1 = np.sort(np.asarray(Ar.sum(1)).ravel())
    assert np.array_equal(d0, d1)
    assert Ar.shape == A.shape
    assert (Ar != Ar.T).nnz == 0  # symmetric


def test_block_diag_has_no_cross_edges():
    A1 = _path_graph(3)
    A2 = _path_graph(3)
    A = B.block_diag_graph([A1, A2])
    assert A.shape == (6, 6)
    assert A[0, 3] == 0 and A[2, 5] == 0
    assert A[0, 1] == 1 and A[3, 4] == 1


def test_normalize_adj_rows_sum_to_one():
    A = _path_graph(5)
    W = B.normalize_adj(A)
    row_sums = np.asarray(W.sum(1)).ravel()
    assert np.allclose(row_sums[row_sums > 0], 1.0)


# ---------------------------------------------------------------------------
# Labels / probe / controls
# ---------------------------------------------------------------------------
def test_binarize_stack_expected():
    m = np.array([[1.0], [2.0], [3.0], [4.0], [5.0]])
    y = B.binarize_stack([m], q=80)
    assert y[:, 0].tolist() == [0, 0, 0, 0, 1]


def test_probe_auroc_chance_on_random():
    rng = np.random.default_rng(0)
    ztr = rng.normal(size=(400, 8))
    ytr = (rng.normal(size=(400, 1)) > 0).astype(int)
    zte = rng.normal(size=(200, 8))
    yte = (rng.normal(size=(200, 1)) > 0).astype(int)
    auc = B.probe_auroc(ztr, ytr, zte, yte)[0]
    assert 0.35 < auc < 0.65


def test_label_shuffle_removes_signal():
    rng = np.random.default_rng(1)
    ztr = rng.normal(size=(400, 8))
    ytr = (ztr[:, 0] > 0).astype(int).reshape(-1, 1)  # perfectly learnable
    zte = rng.normal(size=(200, 8))
    yte = (zte[:, 0] > 0).astype(int).reshape(-1, 1)
    real = B.probe_auroc(ztr, ytr, zte, yte)[0]
    # average over permutations for a stable null estimate
    shuffled = np.mean([B.label_shuffle_auroc(ztr, ytr, zte, yte, seed=s)[0]
                        for s in range(20)])
    assert real > 0.9
    assert 0.42 < shuffled < 0.58


def test_select_C_returns_grid_member():
    rng = np.random.default_rng(2)
    z = rng.normal(size=(300, 6))
    y = (rng.normal(size=(300, 2)) > 0).astype(int)
    groups = np.array(["P1"] * 100 + ["P2"] * 100 + ["P3"] * 100)
    assert B.select_C(z, y, groups, Cs=(0.1, 1.0)) in (0.1, 1.0)


# ---------------------------------------------------------------------------
# Preprocessing: training-only fitting
# ---------------------------------------------------------------------------
def _synth_section(name, patient, n=40, g=80, seed=0):
    rng = np.random.default_rng(seed)
    X = sparse.csr_matrix(rng.poisson(1.0, size=(n, g)).astype(np.float32))
    return B.SectionData(
        specimen=name, patient=patient, treatment="pre", rna=X,
        genes=[f"g{i}" for i in range(g)],
        protein=rng.normal(size=(n, 4)).astype(np.float32),
        programs=rng.normal(size=(n, 11)).astype(np.float32),
        compartment=np.array(["A", "B"] * (n // 2)),
        rows=np.arange(n).astype(float), cols=np.arange(n).astype(float),
    )


def test_preprocessor_fit_is_invariant_to_test_data():
    train = [_synth_section("S1", "P1", seed=1), _synth_section("S2", "P1", seed=2)]
    test_a = _synth_section("T1", "P2", seed=3)
    test_b = _synth_section("T1", "P2", seed=4)
    prep = B.Preprocessor(n_hvg=30, n_pcs=5).fit(train)
    snapshot = (prep.hvg_idx.copy(), prep.rna_scaler.mean_.copy(),
                prep.pca.components_.copy())
    prep.transform(test_a)
    prep.transform(test_b)
    assert np.array_equal(prep.hvg_idx, snapshot[0])
    assert np.array_equal(prep.rna_scaler.mean_, snapshot[1])
    assert np.array_equal(prep.pca.components_, snapshot[2])


def test_featurize_shapes():
    train = [_synth_section("S1", "P1", seed=1)]
    prep = B.Preprocessor(n_hvg=30, n_pcs=5).fit(train)
    feats = B.featurize(prep, train)
    assert feats["rna"].shape == (40, 5)
    assert feats["full"].shape == (40, 9)


def test_embed_method_shapes():
    n, g = 50, 12
    feats = {"rna": np.random.default_rng(0).normal(size=(n, 64)).astype(np.float32)}
    feats["prot"] = np.random.default_rng(1).normal(size=(n, 4)).astype(np.float32)
    feats["full"] = np.hstack([feats["rna"], feats["prot"]])
    adjs = [_path_graph(25), _path_graph(25)]
    for method in ["pca", "smoothing", "random"]:
        model = B.fit_method(method, feats, adjs, seed=0)
        z = B.embed_method(method, model, feats, adjs)
        assert z.shape == (n, 32), method
