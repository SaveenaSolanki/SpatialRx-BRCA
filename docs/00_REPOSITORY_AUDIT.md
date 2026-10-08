# 00 — Repository Audit (Phase 0)

**Audit date:** 2026-10-08
**Repository:** `/storage/saveena/savi_pro/SpatialRx-BRCA`
**Remote:** https://github.com/SaveenaSolanki/SpatialRx-BRCA
**Audited commit:** `fd3de57b6991f7bb302e0baedca883b4bf87bc1c` (`fd3de57`, branch `main`)
**Audit branch:** `audit/phase0`
**Auditor role:** computational biology / ML / software / reproducibility review

This document records the state of the repository **before** remediation. Nothing in the
existing analysis code or outputs was modified during this audit. Two non-destructive
verification scripts were executed (read-only compute; no output overwrite).

---

## 1. How this audit was performed

| Step | Method |
|---|---|
| Structure & git | `find`, `git log/status/branch/remote`, `git ls-files` |
| Documentation | full read of `README.md`, `PLAN.md`, `CHANGELOG.md`, `data/README.md`, both report files |
| Code | full read of all 14 `experiments/*.py` files (2,768 LOC) |
| Environment | direct import + version probe in `/storage2/savi/conda/envs/spatialrx` |
| Data | presence + SHA-256 of raw tar; h5ad inventory; feature/spot counts |
| Claim tracing | matched every headline number in the report to a stored table, and recomputed it where cheap |
| Executability | non-destructive smoke test (load → graph → Moran → score → torch → h5ad write) |

### Audit status taxonomy

| State | Meaning |
|---|---|
| `NOT_IMPLEMENTED` | No code exists for the claimed analysis |
| `IMPLEMENTED` | Code exists but has not been shown to run |
| `EXECUTABLE` | Code runs without error in a defined environment **and** produces its outputs |
| `TESTED` | Automated tests (unit/integration) exercise the code |
| `REPRODUCED` | Recorded outputs can be regenerated from the committed code and match the report |
| `INDEPENDENTLY_VALIDATED` | Result checked against an external/independent source |

**Key rule applied:** an existing script is **not** evidence of a validated scientific result.

---

## 2. Repository inventory

### 2.1 Tracked content (89 files)

| Area | Files | Notes |
|---|---|---|
| Docs | `README.md`, `PLAN.md`, `CHANGELOG.md`, `LICENSE`, `data/README.md` | no `environment.yml`, no `CITATION.cff` |
| Code | `experiments/srx.py` + 13 numbered modules | all Python; no package, no tests |
| Reports | `outputs/SpatialRx-BRCA_report.md`, `outputs/SpatialRx-BRCA_next_steps.md` | canonical results |
| Figures | 29 PNGs (`outputs/figures/`) | dpi=150, raster only |
| Tables | 38 CSVs (`outputs/tables/`) | result tables |

### 2.2 Untracked / local-only content

| Path | Size | Status |
|---|---|---|
| `data/raw/GSE289326_RAW.tar` | 1,824,737,280 B (1.7 GiB) | present, SHA-256 matches `data/README.md` |
| `data/raw/extracted/` | 211 files (14×15 section files + `_layout.txt`) | present |
| `data/processed/` | 8.8 GB (45 h5ad) | present, gitignored |
| `logs/` | 3.2 MB (32 logs) | present, gitignored; **noisy historical runs** |
| `notes/` | empty | placeholder only |
| `experiments/__pycache__/` | — | bytecode cache |

### 2.3 Infrastructure that is absent

`environment.yml` / `pyproject.toml` / `requirements.txt`, `CITATION.cff`, `.github/` (CI),
`tests/`, `docs/` (this is the first), `Makefile`, config file, CLI entry points, data manifest
with per-file hashes.

---

## 3. Data inventory and integrity

| Item | Finding |
|---|---|
| GEO accession | `GSE289326` (CITEgeist cohort) |
| Raw tar SHA-256 | `0292ba69…e75843` — **matches** the value documented in `data/README.md` |
| GSM entries | 14 (`GSM8789203`–`GSM8789216`) |
| Unique tissue sections | 12 + 2 replicate sections (`P4-S2_1i_rep`, `P5-S2_F_rep`) |
| Patients | 6 (`P1`–`P6`); S1 = pre-treatment biopsy, S2 = post-pET surgery |
| Features/section | 18,085 genes + 35 antibody-capture features (`feature_type`-based) |
| Spots ingested | 45,326 (sum of `00_sample_qc.csv`) |
| Spots after module-03 QC | 39,428 (HVG atlas) |
| RNA/protein split | strict `feature_types=="Antibody Capture"` — correct; avoids the dash-suffix trap |
| Companion files | per-section `isotype_normalization_factors.csv` and `spatial_enrichment.csv` extracted |

Data availability is **good**: the primary dataset is present, checksum-verified, and the
intermediate h5ads needed to re-run modules 04–12 without re-ingesting are present.

---

## 4. Environment

Two `spatialrx` environments exist (`~/miniconda3/envs/spatialrx`, `/storage2/savi/conda/envs/spatialrx`).
Probed environment (Python 3.11.15):

```
numpy 2.4.6   pandas 3.0.3   scipy 1.15.3   scikit-learn 1.9.0   statsmodels 0.14.6
matplotlib 3.11.1   seaborn 0.13.2   scanpy 1.11.5   squidpy 1.8.2   decoupler 2.2.0
torch 2.5.1+cpu   harmonypy 2.0.0   leidenalg 0.12.0   igraph 1.0.0   anndata 0.12.19
```

**The environment does not match the README pins.** `README.md` requests `numpy<2.2`; the
installed numpy is `2.4.6` and pandas `3.0.3`. `CHANGELOG.md`/`next_steps.md` already disclose
that installing `citegeist`+`cuopt-cu12` bumped these. Consequence: the environment that
produced modules 00–08 no longer exists as such, and the current environment cannot write
AnnData by default (see defect D5 / risk R1).

---

## 5. Module inventory and status (summary)

Detailed per-module evidence is in `02_IMPLEMENTATION_STATUS.md`.

| # | Module | State | One-line assessment |
|---|---|---|---|
| — | `srx.py` | IMPLEMENTED | constants used by all modules; contains **unused, stale** `PROTEIN_LINEAGE` |
| 00 | `00_ingest.py` | EXECUTABLE (historical) / blocked now | correct feature split; writes h5ad → fails under current pandas |
| 01 | `01_rna_foundations.py` | EXECUTABLE (historical) | single-section RNA workflow; correct |
| 02 | `02_spatial_foundations.py` | EXECUTABLE (historical) | one zoom figure silently skipped |
| 03 | `03_tumor_heterogeneity.py` | EXECUTABLE (historical) | per-sample annotation sound; **merged atlas corrupted** (D1/D2) |
| 04 | `04_spatial_stats.py` | EXECUTABLE (historical) | Moran's I validated vs study table (recomputed ✓) |
| 05 | `05_ecosystem.py` | EXECUTABLE (historical) | programs + heuristic NNLS; isotype alternative not actually implemented |
| 06 | `06_multimodal.py` | EXECUTABLE (historical) | 26/27 pairs analysed; isotype robustness check never called |
| 06b | `06b_protein_spatial.py` | EXECUTABLE (historical) | protein Moran's I; report overstates FDR range (D9) |
| 07 | `07_treatment.py` | EXECUTABLE (historical) | paired design correct to patient; **permutation p doubled** (D3) |
| 08 | `08_ai.py` | EXECUTABLE (historical) | **leaky/circular niche AUROC; no baselines** (D4) |
| 10 | `10_lr_analysis.py` | EXECUTABLE (historical) | curated LR list, no traceable DB; heuristic score |
| 11 | `11_discordance_cellstate.py` | EXECUTABLE (historical) | edge/niche controls implemented |
| 12 | `12_subsampling.py` | EXECUTABLE (historical) | size-matched control implemented correctly |
| 09 | CITEgeist QP run | **NOT_IMPLEMENTED as code** | artifact exists; no script in repo regenerates it |

"No module is `TESTED` (there are no tests) and none is `REPRODUCED` end-to-end from a clean
environment during this audit."

---

## 6. Artifact inventory

- **Figures (29 PNG, dpi=150).** All figures referenced in the report exist except
  `02_spatial_zoom_P4S2.png` (code path wrapped in `try/except`, silently skipped).
- **Tables (38 CSV).** All exist. `03_domain_composition_fractions.csv` is **invalid**
  (index is Leiden cluster IDs, almost all values NaN — see D2).
- **Processed objects (45 h5ad, 8.8 GB).** Present and loadable. `atlas_merged.h5ad` loads
  but its `obs["sample"]` is corrupted (D1).

---

## 7. Defect register (summary)

Full details and severities are in `03_REPRODUCIBILITY_RISKS.md`.

| ID | Defect | Severity |
|---|---|---|
| D1 | `atlas_merged.h5ad` `obs["sample"]` = integer batch labels 0–13; atlas UMAP "Sections" mislabeled | High |
| D2 | Atlas domain annotation yields **no Tumour compartment**; `03_domain_composition_fractions.csv` invalid | High |
| D3 | `exact_perm_paired` returns `2·s/N`; report quotes the un-doubled value while the stored `p_perm` is doubled | High |
| D4 | Module 08 `niche_separation` trains and evaluates logistic regression on the same spots; labels derived from the model's own inputs | High |
| D5 | Current env cannot write h5ad by default (pandas 3.0.3 nullable strings vs anndata 0.12.19) | High (reproducibility) |
| D6 | No pinned environment / no tests / no CI / no config / no CLI | High (engineering) |
| D7 | No random seeds in modules 03/05/07/11 → non-deterministic outputs | Medium |
| D8 | `clr_by_factor` (isotype-normalisation robustness check) is defined but never called; mathematically a no-op as written | Medium |
| D9 | Report says "34–35/35 proteins FDR<0.05"; true range is 30–35 | Medium (claim accuracy) |
| D10 | Module 06 analyses 26 markers (HLA-DRA silently dropped), report says 27 | Low |
| D11 | Report "22 metrics" vs 24 rows in `07_paired_tests.csv` | Low |
| D12 | `PROTEIN_LINEAGE` in `srx.py` is unused and describes a different panel than the actual 35 features | Low |
| D13 | CITEgeist QP (`09_citegeist/`) has no generating script | Medium |
| D14 | Logs contain failed historical runs; not a clean provenance record | Low |
| D15 | `data/processed/` is 8.8 GB of intermediates with no retention/regeneration policy | Low |

---

## 8. What is strong (preserve)

1. **Real data, real execution.** 45,326 spots, 14 sections, honest cohort framing.
2. **Correct cohort hierarchy.** Patient is the unit for paired inference; replicates excluded
   from paired tests (`srx.PAIRED`).
3. **Moran's I external comparison** against the study's own per-feature table reproduces
   exactly (per-sample r 0.695–0.938, mean 0.833; recomputed this audit).
4. **Honest negative result preserved** for the multimodal AI question.
5. **Documented limitations** (n=6, biopsy vs surgical section, PROGENy/cell2location not run).
6. **Data checksum** documented and verified.
7. **Size-matched subsampling control** (module 12) is a genuinely good methodological move.

---

## 9. Overall readiness verdict

**NOT READY.** The scientific core is real and partially reproducible at the artifact level,
but the repository currently fails the Phase-12 acceptance gates on: environment pinning,
tests, CI, reproducibility of the CITEgeist step, several claim/table inconsistencies, a
corrupted merged atlas, a doubled permutation p-value, and a leaky ML benchmark metric.
None of these are fatal to the scientific story, but all must be addressed before the
repository can be presented as a professional, reproducible release.

The next phase (remediation) is described in the execution plan presented after this audit.
