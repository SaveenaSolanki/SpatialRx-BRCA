# 02 — Implementation Status

**Audited commit:** `fd3de57` · **Audit date:** 2026-10-08 · **Branch:** `audit/phase0`

Per-module implementation, executability, testing and reproduction status. Colours:

- State ∈ {`NOT_IMPLEMENTED`, `IMPLEMENTED`, `EXECUTABLE`, `TESTED`, `REPRODUCED`, `INDEPENDENTLY_VALIDATED`}.
- "EXECUTABLE (historical)" = the committed artifact set proves it ran once in the pre-drift
  environment; it is **not** demonstrated to run in the current environment, and it has no tests.

---

## Dependency graph

```
srx.py (constants)
  └─ 00_ingest ──▶ data/processed/00_ingest/*.h5ad  +  00_sample_qc.csv  +  companion CSVs
        ├─ 01_rna_foundations ─▶ 01_rna_foundations/P4-S2_mod1.h5ad ─▶ 02_spatial_foundations
        └─ 03_tumor_heterogeneity ─▶ 03_domains/*_domains.h5ad + atlas_merged.h5ad
              ├─ 04_spatial_stats
              ├─ 05_ecosystem ─▶ 05_ecosystem/*_ecosystem.h5ad
              │     ├─ 06_multimodal
              │     ├─ 06b_protein_spatial
              │     ├─ 07_treatment
              │     ├─ 10_lr_analysis
              │     ├─ 11_discordance_cellstate
              │     └─ 12_subsampling
              └─ 08_ai  (also reads 05_ecosystem)
```

`01`/`02` are a separate teaching track on a single section and are **not** prerequisites for
`03`–`12`. `09` (CITEgeist) exists only as artifacts.

---

## Module-by-module

### `experiments/srx.py` — shared constants (135 LOC)

| Field | Value |
|---|---|
| State | `IMPLEMENTED` (used, but partially stale) |
| Provides | paths, `SAMPLE_META`, `PAIRED`, `REPLICATE_SECTIONS`, `MARKERS`, `PROGRAMS`, `PROTEIN_LINEAGE`, colours, helpers |
| Issues | `PROTEIN_LINEAGE` (D12) is **unused** and describes a 44-entry panel (`Vimentin-1`, `SMA-1`, `ER-alpha-1`, …) that does not match the actual 35-feature panel (`VIM-1`, `ACTA2-1`, `HLA-DRA`, …). Dead and misleading. |

### `00_ingest.py` — ingers GSE289326 (226 LOC)

| Field | Value |
|---|---|
| State | `EXECUTABLE (historical)`; blocked in current env at write time (D5) |
| Inputs | `data/raw/GSE289326_RAW.tar` |
| Outputs | `data/processed/00_ingest/<specimen>.h5ad` (RNA counts, `obsm["protein_counts"]`, spatial coords, images), `00_sample_qc.csv`, companion `__isotype_factors.csv` / `__spatial_enrichment.csv` |
| Strengths | correct `feature_types=="Antibody Capture"` split (the dash-suffix trap is handled); idempotent extraction; per-sample QC |
| Issues | **no random seed needed**; writes h5ad without `allow_write_nullable_strings` → fails under pandas 3 / anndata 0.12; no checksum verification inside the script |

### `01_rna_foundations.py` — single-section RNA workflow (129 LOC)

| Field | Value |
|---|---|
| State | `EXECUTABLE (historical)` |
| Inputs | `00_ingest/HCC22-088-P4-S2.h5ad` |
| Outputs | `01_*` figures, `data/processed/01_rna_foundations/P4-S2_mod1.h5ad` |
| Method | QC → HVG (`seurat_v3`) → norm → log → scale → PCA → neighbours → UMAP → Leiden(0.8) |
| Issues | UMAP/Leiden unseeded (D7); writes h5ad (D5) |

### `02_spatial_foundations.py` — spatial overlay (91 LOC)

| Field | Value |
|---|---|
| State | `EXECUTABLE (historical)` |
| Inputs | `01_rna_foundations/P4-S2_mod1.h5ad` |
| Outputs | `02_spatial_markers_P4S2.png`, `02_leiden_on_tissue_P4S2.png` |
| Issues | `02_spatial_zoom_P4S2.png` is referenced but silently skipped by a bare `try/except` (D11-adjacent); no numeric assertion |

### `03_tumor_heterogeneity.py` — clustering + domain annotation (214 LOC)

| Field | Value |
|---|---|
| State | `EXECUTABLE (historical)`; per-sample part sound, **merged-atlas part defective** |
| Inputs | `00_ingest/*.h5ad` |
| Outputs | `03_domains/<specimen>_domains.h5ad`, `atlas_merged.h5ad`, `03_*` tables/figures |
| Method | per-sample PCA/UMAP/Leiden(0.7) + auditable rule-based domains; merged Harmony atlas + Leiden(0.9) |
| Strengths | compartment-first rule avoids the epithelial≈proliferation collapse; documented |
| Issues | **D1:** `anndata.concat` is called without `keys`, so `obs["section"]` = integers 0–13 and `obs["sample"]` is overwritten with them; the atlas UMAP panel titled "Sections (14)" is therefore coloured by integer id. **D2:** `annotate_domains` is re-applied to the atlas using the *per-sample* `leiden` column, producing compartments {Immune, Stroma, Mixed} with no Tumour and an invalid `03_domain_composition_fractions.csv`. Unseeded (D7). Writes h5ad (D5). |

### `04_spatial_stats.py` — Moran's I (167 LOC)

| Field | Value |
|---|---|
| State | `EXECUTABLE (historical)`; `INDEPENDENTLY VALIDATED` (metric-level) |
| Inputs | `03_domains/*_domains.h5ad`, `00_ingest/*__spatial_enrichment.csv` |
| Outputs | `04_morans_I_per_sample.csv`, `04_svg_ranking.csv`, `04_*` figures |
| Method | 6-neighbour grid graph on array coords; Squidpy Moran's I on all genes; permutation p for top-500 |
| Verification | per-sample r vs study table **recomputed = 0.695–0.938, mean 0.833** ✓ |
| Issues | the authors' graph/normalisation is not documented in the repo, so this is metric agreement not exact replication; permutation p-values exist only for top-500 genes while the report implies all genes |

### `05_ecosystem.py` — programs, CLR, NNLS deconvolution (291 LOC)

| Field | Value |
|---|---|
| State | `EXECUTABLE (historical)` |
| Inputs | `03_domains/*_domains.h5ad` |
| Outputs | `05_ecosystem/<specimen>_ecosystem.h5ad`, `05_*` tables/figures, `05_nhood_*.csv` |
| Method | `score_genes` programs; winsorise+CLR protein; per-spot NNLS (RNA 6 lineages, protein 8 lineages); squidpy neighbourhood enrichment |
| Issues | **D8:** the claimed isotype-factor alternative is absent (only plain CLR is computed); NNLS signatures are heuristic and not independently validated except for one lineage; unseeded (D7); writes h5ad (D5) |

### `06_multimodal.py` — RNA vs protein agreement/discordance (210 LOC)

| Field | Value |
|---|---|
| State | `EXECUTABLE (historical)` |
| Inputs | `05_ecosystem/*_ecosystem.h5ad`, `00_ingest/*.h5ad` |
| Outputs | `06_rna_protein_agreement.csv`, `06_discordance_morans.csv`, `06_*` figures |
| Method | per-spot Pearson/Spearman for matched pairs; discordance z-score difference; Moran's I of discordance |
| Issues | **D10:** 26/27 markers analysed — `HLA-DRA` is dropped because `split("-")[0]` yields `HLA`; **D8:** `clr_by_factor` never called (dead code, mathematically a no-op) |

### `06b_protein_spatial.py` — protein Moran's I (94 LOC)

| Field | Value |
|---|---|
| State | `EXECUTABLE (historical)` |
| Inputs | `05_ecosystem/*_ecosystem.h5ad`, authors' CSV |
| Outputs | `06b_protein_morans.csv`, `06b_protein_spatial.png` |
| Verification | pooled r vs authors **recomputed = 0.694** ✓; per-sample 0.22–0.94 |
| Issues | **D9:** report says 34–35/35 FDR<0.05; true range 30–35 |

### `07_treatment.py` — paired pre/post (326 LOC)

| Field | Value |
|---|---|
| State | `EXECUTABLE (historical)`; design `INDEPENDENTLY VALIDATED`-quality (patient-level) |
| Inputs | `05_ecosystem/*_ecosystem.h5ad`, `00_ingest/*.h5ad`, `05_nhood_*.csv` |
| Outputs | `07_paired_metrics.csv`, `07_paired_tests.csv`, `07_program_morans*.csv`, `07_*` figures |
| Method | per-patient paired deltas; exact sign-flip + Wilcoxon; BH; tumour–immune contact; program Moran's I |
| Strengths | correct independent unit (patient); replicates excluded; n reported honestly |
| Issues | **D3:** `exact_perm_paired` returns `2·s/N`; stored `p_perm` is doubled and inconsistent with the report (C18). **D11:** "22 metrics" vs 24 rows. Contact metric rests on n=3 pairs. Unseeded (D7) |

### `08_ai.py` — multimodal embeddings (289 LOC)

| Field | Value |
|---|---|
| State | `EXECUTABLE (historical)` |
| Inputs | `03_domains/*_domains.h5ad`, `05_ecosystem/*_ecosystem.h5ad` |
| Outputs | `08_embedding_metrics.csv`, `08_*` figures |
| Method | 3 autoencoders (RNA MLP; RNA+protein MLP; GCN-prelayer+MLP) 32-d, MSE, 300 epochs |
| Issues | **D4 (critical):** `niche_separation` fits `LogisticRegression` and evaluates `roc_auc_score` on the **same** spots; labels are program-high spots derived from the same RNA input. The AUROC is optimistic and circular. **No baselines** (PCA, smoothing, edge randomisation), **no patient/section holdout**, **no seeds across runs** (torch/np seeded once but not varied), **no inductive evaluation**. Spatial-coherence metric is largely a manipulation check (acknowledged in report). |

### `10_lr_analysis.py` — ligand–receptor at interfaces (219 LOC)

| Field | Value |
|---|---|
| State | `EXECUTABLE (historical)` |
| Inputs | `05_ecosystem/*_ecosystem.h5ad`, `00_ingest/*.h5ad` |
| Outputs | `10_lr_*.csv`, `10_*` figures |
| Method | hand-curated 26 LR pairs; score = mean(sender ligand)×mean(receiver receptor); exact sign-flip on n≤4 |
| Issues | no traceable LR resource (CellChat/LIANA) — Phase 5 requirement; no expression-shuffled / spatially-shuffled controls; no multiple-testing correction; n=4. Report already labels as hypothesis-generating |

### `11_discordance_cellstate.py` — discordance drivers (191 LOC)

| Field | Value |
|---|---|
| State | `EXECUTABLE (historical)` |
| Inputs | `05_ecosystem/*_ecosystem.h5ad`, `00_ingest/*.h5ad` |
| Outputs | `11_*.csv`, `11_discordance_drivers.png` |
| Method | Spearman vs edge-exposure and cell-state proportions; hotspot composition |
| Issues | edge metric is grid-degree based (a proxy, not physical distance); no multiple-testing correction; unseeded |

### `12_subsampling.py` — size-matched control (186 LOC)

| Field | Value |
|---|---|
| State | `EXECUTABLE (historical)` |
| Inputs | `05_ecosystem/*_ecosystem.h5ad` |
| Outputs | `12_*.csv`, `12_subsampling.png` |
| Method | random + contiguous BFS subsampling of post sections to paired pre size; B=200/100 |
| Strengths | deterministic (`default_rng(0)`); good methodological control |
| Issues | `n_pre` computed then unused (harmless); contact metric undefined where a patient has no tumour/immune pair (P3/P5/P6) |

### CITEgeist Step 2 (`data/processed/09_citegeist/`)

| Field | Value |
|---|---|
| State | `NOT_IMPLEMENTED` (as repository code) |
| Evidence | CSV + log artifacts exist; no script in `experiments/` |
| Issues | `logs/cg_step2.log` ends in an AnnData write error; result cannot be regenerated from the repo |

---

## Test coverage

| Test type | Present? | Count |
|---|---|---|
| Unit tests | No | 0 |
| Integration tests | No | 0 |
| Data-shape / feature-mapping assertions | No (only inline `print`s) | 0 |
| CI workflow | No | 0 |
| Smoke-test entry point | No | 0 |

---

## Status summary

| State | Modules |
|---|---|
| `INDEPENDENTLY_VALIDATED` | 04 (metric-level vs study table) |
| `REPRODUCED` (artifact-level, this audit) | none fully; 04's headline number recomputed |
| `EXECUTABLE (historical)` | 00, 01, 02, 03, 04, 05, 06, 06b, 07, 08, 10, 11, 12 |
| `IMPLEMENTED` only | `srx.py` (partly stale) |
| `NOT_IMPLEMENTED` | CITEgeist `09` script; Module-08 baselines/holdout (required by brief) |
| `TESTED` | none |

**No module can be marked `REPRODUCED` until the environment is pinned and at least a smoke
run is executed in it.** Remediation priorities are ordered in the execution plan.
