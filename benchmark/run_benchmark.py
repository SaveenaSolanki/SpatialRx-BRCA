#!/usr/bin/env python
"""Run the Gate-4 leakage-free benchmark (see benchmark/BENCHMARK_CONTRACT.md).

Folds are independent, so they are executed in parallel worker processes.
Preprocessing and models are still fit strictly per fold on that fold's
training sections.

    python benchmark/run_benchmark.py --analyses patient section --seeds 5 --jobs 10
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

# Pin BLAS threads BEFORE importing numpy; each parallel worker uses 2.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "2")

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from spatialrx import benchmark as B  # noqa: E402
from spatialrx.config import PATHS  # noqa: E402
from spatialrx.stats import bh_fdr, exact_perm_paired  # noqa: E402

_SECTIONS = None  # inherited by forked workers


# ---------------------------------------------------------------------------
# Aggregation helpers
# ---------------------------------------------------------------------------
def mean_ci(values):
    v = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
    if len(v) == 0:
        return np.nan, np.nan, np.nan, 0
    m = v.mean()
    if len(v) == 1:
        return m, m, m, 1
    from scipy import stats
    half = stats.t.ppf(0.975, len(v) - 1) * v.std(ddof=1) / np.sqrt(len(v))
    return m, m - half, m + half, len(v)


def summarize(long: pd.DataFrame, value: str, metric: str) -> pd.DataFrame:
    if long.empty:
        return pd.DataFrame()
    per_ps = (long.groupby(["analysis", "method", "patient", "seed"])[value]
              .mean().reset_index())
    rows = []
    for (analysis, method), sub in per_ps.groupby(["analysis", "method"]):
        per_patient = sub.groupby("patient")[value].mean()
        m, lo, hi, n = mean_ci(per_patient.values)
        seed_sd = sub.groupby("patient")[value].std().mean()
        rows.append({"analysis": analysis, "method": method, "metric": metric,
                     "mean": m, "ci_lo": lo, "ci_hi": hi, "n_patients": n,
                     "seed_sd": seed_sd})
    return pd.DataFrame(rows)


def paired_vs_reference(long: pd.DataFrame, value: str, metric: str,
                        reference: str) -> pd.DataFrame:
    per_ps = (long.groupby(["analysis", "method", "patient", "seed"])[value]
              .mean().reset_index())
    per_p = per_ps.groupby(["analysis", "method", "patient"])[value].mean()
    rows = []
    for analysis in per_p.index.get_level_values("analysis").unique():
        ref = per_p.loc[analysis, reference]
        for method in per_p.loc[analysis].index.get_level_values("method").unique():
            if method == reference:
                continue
            comp = per_p.loc[analysis, method]
            joined = pd.concat([ref, comp], axis=1, keys=["ref", "method"]).dropna()
            if len(joined) < 3:
                continue
            d = (joined["method"] - joined["ref"]).to_numpy()
            obs, p = exact_perm_paired(d)
            pooled = np.sqrt((joined["ref"].var(ddof=1) + joined["method"].var(ddof=1)) / 2)
            es = d.mean() / pooled if pooled > 0 else np.nan
            rows.append({"analysis": analysis, "method": method, "metric": metric,
                         "reference": reference, "mean_diff": obs, "p_perm": p,
                         "n_patients": len(joined), "effect_size": es})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["p_BH"] = bh_fdr(out["p_perm"].to_numpy(dtype=float))
    return out


# ---------------------------------------------------------------------------
# Worker
# ---------------------------------------------------------------------------
def _init_worker(threads: int):
    global _SECTIONS
    import torch
    torch.set_num_threads(threads)
    if _SECTIONS is None:
        _SECTIONS = B.load_sections(PATHS)


def _run_fold(task):
    analysis, fold_name, test_patients, train_specs, test_specs, seeds, methods = task
    sections = _SECTIONS
    train = [sections[s] for s in train_specs]
    test = [sections[s] for s in test_specs]
    manifest = {"analysis": analysis, "fold": fold_name,
                "test_patients": "|".join(test_patients),
                "n_train_sections": len(train), "n_test_sections": len(test),
                "n_train_spots": int(sum(s.rna.shape[0] for s in train)),
                "n_test_spots": int(sum(s.rna.shape[0] for s in test))}

    prep = B.Preprocessor().fit(train)
    ftr, fte = B.featurize(prep, train), B.featurize(prep, test)
    ytr = B.binarize_stack([s.programs for s in train])
    yte = B.binarize_stack([s.programs for s in test])
    atr = [s.adj for s in train]
    ate = [s.adj for s in test]
    groups = np.concatenate([np.full(s.rna.shape[0], s.patient) for s in train])
    tp = test_patients[0]

    rows_auroc, rows_emb, rows_stab, rows_ctrl = [], [], [], []
    for method in methods:
        seed_list = list(range(seeds)) if method in B.STOCHASTIC else [0]
        model0 = B.fit_method(method, ftr, atr, seed_list[0])
        ztr0 = B.embed_method(method, model0, ftr, atr)
        C = B.select_C(ztr0, ytr, groups)
        z0 = None
        for seed in seed_list:
            if seed == seed_list[0]:
                model, ztr = model0, ztr0
            else:
                model = B.fit_method(method, ftr, atr, seed)
                ztr = B.embed_method(method, model, ftr, atr)
            zte = B.embed_method(method, model, fte, ate)
            auc = B.probe_auroc(ztr, ytr, zte, yte, C=C)
            for j, prog in enumerate(B.PROGRAMS):
                rows_auroc.append([analysis, method, seed, fold_name, tp, prog, auc[j], C])
            parts = B.split_sections(zte, test)
            sil = B.silhouette_per_section([p[0] for p in parts],
                                           [p[1] for p in parts], test)
            mor = B.moran_per_section([p[0] for p in parts], [p[2] for p in parts])
            rows_emb.append([analysis, method, seed, fold_name, tp, sil, mor])
            if seed == 0:
                z0 = zte
            else:
                rows_stab.append([analysis, method, seed, fold_name, tp,
                                  B.seed_correlation(z0, zte)])
            if method == "rna_ae":
                sh = B.label_shuffle_auroc(ztr, ytr, zte, yte, seed=seed)
                for j, prog in enumerate(B.PROGRAMS):
                    rows_ctrl.append([analysis, "label_shuffle", seed, fold_name, tp,
                                      prog, sh[j]])
    return manifest, rows_auroc, rows_emb, rows_stab, rows_ctrl


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
AUROC_COLS = ["analysis", "method", "seed", "fold", "patient", "program", "auroc", "C"]
EMB_COLS = ["analysis", "method", "seed", "fold", "patient", "silhouette", "moran"]
STAB_COLS = ["analysis", "method", "seed", "fold", "patient", "stability_vs_seed0"]
CTRL_COLS = ["analysis", "control", "seed", "fold", "patient", "program", "value"]


def build_tasks(sections, analyses, seeds, methods, max_folds):
    tasks = []
    for analysis in analyses:
        folds = B.make_folds(sections, analysis)
        B.assert_disjoint(folds)
        if max_folds:
            folds = folds[:max_folds]
        for f in folds:
            tasks.append((analysis, f.name, f.test_patients,
                          f.train_sections, f.test_sections, seeds, methods))
    return tasks


def run(analyses, seeds, methods, max_folds=None, out_dir: Path = None, jobs: int = 10,
        resume: bool = True):
    import pickle
    global _SECTIONS
    import torch
    torch.set_num_threads(1)  # keep parent quiet before forking
    out_dir = out_dir or (ROOT / "benchmark")
    out_dir.mkdir(parents=True, exist_ok=True)
    folds_dir = out_dir / "folds"
    folds_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    print("[bench] loading sections ...", flush=True)
    _SECTIONS = B.load_sections(PATHS)
    print(f"[bench] {len(_SECTIONS)} sections loaded in {time.time()-t0:.1f}s", flush=True)

    tasks = build_tasks(_SECTIONS, analyses, seeds, methods, max_folds)
    done = {}
    if resume:
        for p in folds_dir.glob("*.pkl"):
            try:
                done[p.stem] = pickle.loads(p.read_bytes())
            except Exception:
                pass
    pending = [t for t in tasks if f"{t[0]}__{t[1]}" not in done]
    print(f"[bench] {len(tasks)} folds ({len(done)} cached, {len(pending)} pending), "
          f"jobs={jobs}", flush=True)

    if pending:
        if jobs <= 1:
            _init_worker(2)
            iterator = iter(_run_fold(t) for t in pending)
            for i, res in enumerate(iterator, 1):
                key = f"{res[0]['analysis']}__{res[0]['fold']}"
                done[key] = res
                (folds_dir / f"{key}.pkl").write_bytes(pickle.dumps(res))
                print(f"[bench]   fold {i}/{len(pending)} done ({res[0]['fold']})", flush=True)
        else:
            ctx = mp.get_context("spawn")
            with ctx.Pool(jobs, initializer=_init_worker, initargs=(2,)) as pool:
                for i, res in enumerate(pool.imap_unordered(_run_fold, pending, chunksize=1), 1):
                    key = f"{res[0]['analysis']}__{res[0]['fold']}"
                    done[key] = res
                    (folds_dir / f"{key}.pkl").write_bytes(pickle.dumps(res))
                    print(f"[bench]   fold {i}/{len(pending)} done ({res[0]['fold']})", flush=True)

    manifests = [r[0] for r in done.values()]
    aurocs = [row for r in done.values() for row in r[1]]
    embs = [row for r in done.values() for row in r[2]]
    stabs = [row for r in done.values() for row in r[3]]
    ctrls = [row for r in done.values() for row in r[4]]

    auroc = pd.DataFrame(aurocs, columns=AUROC_COLS)
    emb = pd.DataFrame(embs, columns=EMB_COLS)
    stab = pd.DataFrame(stabs, columns=STAB_COLS)
    ctrl = pd.DataFrame(ctrls, columns=CTRL_COLS)
    auroc.to_csv(out_dir / "results_auroc.csv", index=False)
    emb.to_csv(out_dir / "results_embedding.csv", index=False)
    pd.DataFrame(manifests).to_csv(out_dir / "split_manifest.csv", index=False)
    stab.to_csv(out_dir / "stability.csv", index=False)
    ctrl.to_csv(out_dir / "negative_controls.csv", index=False)

    summary = pd.concat([summarize(auroc, "auroc", "auroc"),
                         summarize(emb, "silhouette", "silhouette"),
                         summarize(emb, "moran", "moran_spatial_coherence")],
                        ignore_index=True)
    summary.to_csv(out_dir / "summary.csv", index=False)

    ab = pd.concat([paired_vs_reference(auroc, "auroc", "auroc", "pca"),
                    paired_vs_reference(auroc, "auroc", "auroc", "rna_ae"),
                    paired_vs_reference(emb, "silhouette", "silhouette", "pca"),
                    paired_vs_reference(emb, "moran", "moran_spatial_coherence", "rna_ae")],
                   ignore_index=True)
    ab.to_csv(out_dir / "ablation_results.csv", index=False)

    ctrl_summary = []
    for (analysis, control), sub in ctrl.groupby(["analysis", "control"]):
        m, lo, hi, n = mean_ci(sub.groupby("patient")["value"].mean().values)
        ctrl_summary.append({"analysis": analysis, "control": control, "mean": m,
                             "ci_lo": lo, "ci_hi": hi, "n_patients": n})
    for analysis, sub in auroc[auroc.method == "random"].groupby("analysis"):
        m, lo, hi, n = mean_ci(sub.groupby("patient")["auroc"].mean().values)
        ctrl_summary.append({"analysis": analysis, "control": "random_projection",
                             "mean": m, "ci_lo": lo, "ci_hi": hi, "n_patients": n})
    pd.DataFrame(ctrl_summary).to_csv(out_dir / "negative_controls_summary.csv", index=False)

    runtime = round(time.time() - t0, 1)
    record = {"runtime_s": runtime, "analyses": analyses, "seeds": list(range(seeds)),
              "methods": methods, "jobs": jobs, "n_sections": len(_SECTIONS),
              "n_folds": len(tasks),
              "outputs": sorted([p.name for p in out_dir.glob("*.csv")]
                                + [p.name for p in out_dir.glob("*.md")])}
    (out_dir / "run_record.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"[bench] done in {runtime}s -> {out_dir}", flush=True)
    return summary, ab


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analyses", nargs="+", default=["patient", "section"])
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--methods", nargs="+", default=B.METHODS)
    ap.add_argument("--max-folds", type=int, default=None)
    ap.add_argument("--jobs", type=int, default=10)
    ap.add_argument("--no-resume", action="store_true")
    ap.add_argument("--out-dir", default=str(ROOT / "benchmark"))
    args = ap.parse_args()
    summary, ab = run(args.analyses, args.seeds, args.methods, args.max_folds,
                      Path(args.out_dir), args.jobs, resume=not args.no_resume)
    print("\n=== summary ===")
    print(summary.round(4).to_string(index=False))
    print("\n=== paired vs reference ===")
    print(ab.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
