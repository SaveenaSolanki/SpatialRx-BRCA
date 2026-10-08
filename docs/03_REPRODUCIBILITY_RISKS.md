# 03 — Reproducibility Risks

**Audited commit:** `fd3de57` · **Audit date:** 2026-10-08 · **Branch:** `audit/phase0`

Risk register for re-running and releasing SpatialRx-BRCA. Severity combines likelihood and
scientific impact. Each risk has a concrete mitigation and the phase in which it should be
resolved.

---

## 0. Live reproducibility test performed in this audit

A non-destructive smoke test was run in `/storage2/savi/conda/envs/spatialrx`
(Python 3.11.15, anndata 0.12.19, pandas 3.0.3) against the committed processed objects.
No repository outputs were overwritten.

| Test | Result |
|---|---|
| Load `03_domains/HCC22-088-P4-S2_domains.h5ad` + `05_ecosystem/..._ecosystem.h5ad` | ✅ pass |
| Build 6-neighbour grid graph + Squidpy Moran's I on an `obs` program score | ✅ pass (I=0.585) |
| `scanpy.tl.score_genes` | ✅ pass |
| PyTorch model forward pass | ✅ pass |
| AnnData `write()` round-trip | ❌ **fail** — `RuntimeError: anndata.settings.allow_write_nullable_strings is False` |
| AnnData `write()` with `anndata.settings.allow_write_nullable_strings = True` | ✅ pass |
| Recompute Moran validation r from `04_morans_I_per_sample.csv` + authors' CSV | ✅ pass, mean 0.833 |

**Interpretation:** read/compute paths still work, but **every module that writes an h5ad
(00, 01, 03, 05) fails by default in the current environment**. This is the single most
important reproducibility finding.

---

## 1. Risk register

### R1 — Unpinned, drifted environment (blocker)

| Field | Detail |
|---|---|
| Severity | **Critical** |
| Evidence | `README.md` pins `numpy<2.2`; installed numpy `2.4.6`, pandas `3.0.3`. `CHANGELOG.md` admits the drift. h5ad write fails (smoke test). |
| Impact | Modules 00/01/03/05 are not executable as committed; "reproduce from clean env" is impossible; review-level claim of reproducibility fails. |
| Mitigation | Add a pinned `environment.yml` + `pyproject.toml` (or `requirements.lock`) capturing the *working* versions; set `anndata.settings.allow_write_nullable_strings = True` in a shared helper/entry point; add `spatialrx validate-data` to assert versions. Optionally pin pandas `<3` for exact byte-reproduction of modules 00–08 and keep a second env for CITEgeist. |
| Phase | 2 / 7 |

### R2 — No tests and no CI

| Field | Detail |
|---|---|
| Severity | **Critical** (engineering) |
| Evidence | no `tests/`, no `.github/`; zero assertions in code |
| Impact | regressions undetectable; reviewers cannot verify metadata/normalisation/statistics helpers |
| Mitigation | unit tests for metadata parsing, feature split, CLR, NNLS, permutation test, BH; integration test on a tiny synthetic fixture; CI running those without the cohort. |
| Phase | 7 |

### R3 — Corrupted merged atlas (D1/D2)

| Field | Detail |
|---|---|
| Severity | **High** |
| Evidence | `atlas_merged.h5ad.obs["sample"]` = {0…13}; compartments = {Immune, Stroma, Mixed}, no Tumour; `03_domain_composition_fractions.csv` is indexed by Leiden id with almost all values NaN |
| Impact | `03_atlas_umap.png` "Sections (14)" panel is mislabeled; atlas composition table is invalid; any future use of `atlas_merged.h5ad` inherits the bug |
| Mitigation | pass `keys=` to `anndata.concat`; preserve specimen names in a dedicated obs column; re-annotate atlas using atlas-level clusters (`leiden_atlas`), not per-sample `leiden`; add a shape/label assertion before writing; regenerate `03_domain_composition_fractions.csv`. |
| Phase | 2 / 7 |

### R4 — Doubled exact permutation p-value (D3)

| Field | Detail |
|---|---|
| Severity | **High** (statistics correctness) |
| Evidence | `exact_perm_paired` in `07_treatment.py` returns `2*s/total`; stored B-cell `p_perm=0.1875` while report quotes `0.094` (the correct exact value / Wilcoxon). |
| Impact | all permutation p-values in `07_paired_tests.csv` and `10_lr_paired_tests.csv` are inflated ≈2×; uncorrected trend claims are wrong; BH conclusion unchanged. |
| Mitigation | return `s/total` (already two-sided via `|stat| ≥ |obs|`); add a unit test with a known small-n example; recompute tables and update report. |
| Phase | 3 / 7 |

### R5 — Leaky/circular ML benchmark metric (D4)

| Field | Detail |
|---|---|
| Severity | **High** (ML validity) |
| Evidence | `08_ai.py::niche_separation` calls `clf.fit(z, y)` then `roc_auc_score(y, clf.decision_function(z))` on the same spots; `y` = program-high derived from the same RNA used to build `z`. |
| Impact | headline "niche AUROC 0.765→0.752" is not a valid estimate; no train/test separation; no patient/section holdout; no baselines; no edge-randomisation control; no seed variation. |
| Mitigation | split spots train/test (with patient-level grouping for inductive tasks); fit transformations only on train; add PCA, spatial-smoothing, and graph-edge-randomisation baselines; report across seeds and patients; rename spatial-coherence metric as a manipulation check. Preserve the negative finding only after re-evaluation. |
| Phase | 4 |

### R6 — Claim/table inconsistencies (D9–D12)

| Field | Detail |
|---|---|
| Severity | **Medium** |
| Evidence | protein FDR range 30–35 vs reported 34–35; 26 vs 27 markers (HLA-DRA dropped); 22 vs 24 metrics; isotype robustness never run (`clr_by_factor` dead code, no-op). |
| Impact | erodes reviewer trust; some "robustness" claims are unsupported. |
| Mitigation | reconcile report to regenerated tables; either implement isotype normalisation correctly (or delete the claim); fix marker matching to handle `HLA-DRA`. |
| Phase | 3 / 6 |

### R7 — Non-determinism (D7)

| Field | Detail |
|---|---|
| Severity | **Medium** |
| Evidence | no seeds in modules 03/05/07/11; UMAP/Leiden/PCA randomised components; module 08 seeds once but does not vary seeds. |
| Impact | figures/tables not bit-reproducible; uncertainty across seeds not quantified. |
| Mitigation | central seed config; set numpy/torch/python/scanpy seeds; run ≥5 seeds for ML metrics; report seed variance. |
| Phase | 2 / 4 |

### R8 — CITEgeist step not reproducible (D13)

| Field | Detail |
|---|---|
| Severity | **Medium** |
| Evidence | `data/processed/09_citegeist/` artifacts exist; no script; `cg_step2.log` ends in an AnnData write error. |
| Impact | the NNLS-vs-CITEgeist validation (r=0.76) cannot be regenerated or reviewed. |
| Mitigation | add `experiments/09_citegeist_qp.py` (or document as an external step with exact commands/versions), mark optional/heavy, and keep it out of the default smoke path. |
| Phase | 2 / 6 |

### R9 — Heuristic deconvolution / LR scores without external resources

| Field | Detail |
|---|---|
| Severity | **Medium** |
| Evidence | NNLS signatures built from hand-curated markers; LR list hand-curated; no LIANA/CellChat resource; no shuffle controls. |
| Impact | biological interpretation is weaker than it appears; not directly comparable to published methods. |
| Mitigation | document as heuristic; add expression-shuffled and spatially-shuffled controls (Phase 5); optionally add a traceable resource; keep language as hypothesis. |
| Phase | 5 |

### R10 — Large intermediates, no data policy

| Field | Detail |
|---|---|
| Severity | **Low–Medium** |
| Evidence | `data/processed/` = 8.8 GB, 45 h5ad, untracked; no retention policy or manifest of hashes. |
| Impact | disk churn; reviewers cannot verify intermediate provenance; risk of accidental commit. |
| Mitigation | add `data/MANIFEST.md` with hashes/sizes; separate raw/intermediate/release dirs; add a `make clean`-style target; keep `.gitignore` strict. |
| Phase | 2 |

### R11 — No external validation

| Field | Detail |
|---|---|
| Severity | **Medium** (scientific completeness) |
| Evidence | `next_steps.md` documents a cohort scan but no candidate dataset was actually used; no compatibility table artifact. |
| Impact | Phase 6 requirement unmet; can be honestly declared unavailable, but must be recorded with a search manifest. |
| Mitigation | produce a machine-readable compatibility table + search log; either run one defensible cross-cohort analysis (e.g., GSE331245 signature transfer) or mark `EXTERNAL_VALIDATION_UNAVAILABLE` with reasons. |
| Phase | 6 |

### R12 — Logs are a noisy provenance record

| Field | Detail |
|---|---|
| Severity | **Low** |
| Evidence | `m4_run.log`, `m7_run.log`, `m8_run.log`, `m11_run*.log`, `cg_step2.log` contain tracebacks from superseded runs; run numbering is inconsistent. |
| Impact | cannot cite a log as the evidence for a final number. |
| Mitigation | structured execution logs (module, git SHA, env hash, start/end, status) written to `logs/runs/`; keep failed runs separate. |
| Phase | 2 |

### R13 — Missing publication/engineering scaffolding

| Field | Detail |
|---|---|
| Severity | **Low–Medium** |
| Evidence | no `environment.yml`, `CITATION.cff` (LICENSE present), no CLI, no figure-generation scripts separate from analysis, figures are dpi=150 raster only. |
| Impact | fails Phase 7/8/9 acceptance gates. |
| Mitigation | add packaging + CLI (`spatialrx …`), a `figures/` generation layer exporting PDF/SVG + 600-dpi PNG, `CITATION.cff`, README upgrades. |
| Phase | 7 / 8 / 9 |

---

## 2. Risk heat map

| | Low impact | Medium impact | High impact | Critical |
|---|---|---|---|---|
| **High likelihood** | R12 | R6, R7, R9, R10 | R3, R4 | R1, R2 |
| **Medium likelihood** | R13 | R11 | R5 | — |
| **Low likelihood** | — | R8 | — | — |

---

## 3. Reproducibility checklist (current state)

| Requirement | Status |
|---|---|
| Version-controlled code | ✅ |
| Documented data source + checksum | ✅ |
| Raw data present and verified | ✅ |
| Pinned environment | ❌ |
| One-command smoke test | ❌ |
| Configurable paths | ❌ (hard-coded relative paths) |
| Deterministic seeds | ❌ (partial) |
| Structured logs | ❌ |
| Restart/resume for long jobs | ❌ |
| Unit/integration tests | ❌ |
| CI | ❌ |
| Separate raw / intermediate / release | ⚠️ partial (`data/` only) |
| Data manifest with hashes | ⚠️ tar only |
| Figures regenerable by script | ⚠️ yes, but 150-dpi raster and script-coupled |

---

## 4. Highest-priority fixes (ordered)

1. **R1** pin environment + make h5ad writes work everywhere.
2. **R2** tests + CI for metadata/normalisation/statistics helpers.
3. **R3** fix merged atlas and its invalid composition table.
4. **R4** fix exact permutation p-value and recompute all paired tables.
5. **R5** replace the leaky niche metric with a split/patient-held-out benchmark and baselines.
6. **R6–R7** reconcile report to tables; add seeds.
7. **R8–R11** CITEgeist script, traceable LR resources + controls, data manifest, external-validation manifest.
8. **R12–R13** structured logs, packaging/CLI, publication figures.

These map directly onto the execution plan presented after this audit.
