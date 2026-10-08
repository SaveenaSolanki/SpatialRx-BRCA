# Benchmark contract — Module 08 reconstruction (Gate 4)

**Status:** agreed pre-registration. Written **before** the benchmark was implemented.
**Branch:** `fix/repro-correctness-v1`  ·  **Baseline:** `1d7d2f7e72f66b34db6eb5b3050ed3b409efbed9`
**Owner:** SpatialRx-BRCA remediation  ·  **Date:** 2026-10-08

This contract fixes the evaluation design in advance so that the final numbers cannot be
tuned after seeing held-out performance. Any deviation must be recorded in
`benchmark/statistical_analysis.md` under "deviations".

---

## 1. Question

Does adding protein and/or spatial-graph context to a learned representation improve
*biologically meaningful, transferable* spot-level structure compared with an RNA-only
representation — under a leakage-free, patient-held-out evaluation?

This replaces the original Module 08 question, whose metric was evaluated on the same spots
used to fit the model (issue ISSUE-004).

## 2. Data and objects

- Cohort: GSE289326, 14 sections (12 paired + 2 replicates) from 6 patients.
- Frozen inputs (never re-derived for the benchmark):
  - `data/processed/03_domains/<specimen>_domains.h5ad` → `.raw` (full-gene log-normalised
    expression), spatial `array_row`/`array_col`, `compartment` annotation.
  - `data/processed/05_ecosystem/<specimen>_ecosystem.h5ad` → `obsm["protein_clr"]`,
    `score_<program>` columns (11 programs).
- No new data are downloaded. Raw counts/ingest are untouched.

## 3. Biological replicate units

- **Primary unit = patient.** Sections are nested within patients; replicate sections
  (`P4-S2_1i_rep`, `P5-S2_F_rep`) always travel with their patient.
- Spots are **not** independent replicates and are never treated as such.
- Secondary unit = section (exploratory, because a patient may appear in both train and test).

## 4. Leakage boundaries (hard constraints)

1. A held-out **patient's** spots, labels, graph, protein and gene-percentile thresholds are
   never used to fit any preprocessing, representation model, or probe.
2. HVG selection, scaling, PCA, and protein scaling are fit on the training patients only.
3. For the GCN, spatial edges are **within-section only**; no edge crosses a train/test
   boundary. Using the held-out section's own (label-free) graph at inference is a declared
   **transductive** boundary, reported explicitly; the degree-randomised control removes it.
4. The supervised probe is trained only on training-patient embeddings and training-patient
   binarised labels; it is evaluated on the held-out patient.
5. Test-invariance check: mutating held-out expression must not change any fitted preprocessing
   parameter or model weight (enforced by a regression test).
6. Thresholds for "program-high" labels are computed per section (train thresholds from train
   sections, test thresholds from the test section) so no global statistic leaks across folds.

## 5. Splits

- **Primary — Leave-One-Patient-Out (LOPO):** 6 folds; each fold holds out all sections of one
  patient; the other 5 patients (including their replicates) are training data.
- **Secondary — Leave-One-Section-Out (LOSO):** 14 folds; one section held out at a time.
  Explicitly exploratory (patient leakage possible); never used for the primary conclusion.
- All methods are evaluated on **identical** partitions (one `split_manifest.csv`).
- Nested validation: within each training set, an inner leave-one-training-patient-out selects
  the probe's regularisation `C`; the outer held-out patient is touched once, for testing.

## 6. Preprocessing (training-only)

| Step | Detail |
|---|---|
| HVG | top 2000 `seurat` HVGs computed on the concatenated **training** log-normalised matrix |
| Scaling | per-gene z-score, train mean/SD |
| RNA features | PCA(64) fit on training spots |
| Protein features | 35 `protein_clr` dims, z-scored on train |
| Spatial graph | 6-neighbour Visium grid from `array_row`/`array_col`, built per section |

## 7. Methods and baselines (identical folds)

| ID | Description | Kind |
|---|---|---|
| `pca` | first 32 dims of the training-fitted PCA(64) | linear baseline |
| `smoothing` | graph-average the 64 RNA PCs, then take 32 dims | non-learned spatial baseline |
| `random` | Gaussian random projection of the 64 RNA PCs → 32 | negative control |
| `rna_ae` | MLP autoencoder on the 64 RNA PCs → 32 | learned, RNA-only |
| `rna_protein_ae` | MLP autoencoder on (64 RNA PCs + 35 proteins) → 32 | multimodal |
| `graph_ae` | graph-conv pre-layer + MLP AE on (RNA + protein), real graph | multimodal + spatial |
| `graph_ae_randomized` | as `graph_ae` with a degree-preserving randomised adjacency | spatial negative control |

Model hyperparameters are fixed a priori (32-d latent, 128 hidden, Adam lr 1e-3, 200 epochs,
full batch). No hyperparameter is selected on the held-out patient.

## 8. Negative controls and expected behaviour

| Control | Expectation if the benchmark is valid |
|---|---|
| `random` projection | primary AUROC ≈ chance (0.5) |
| label-shuffle probe (train labels permuted) | AUROC ≈ chance |
| `graph_ae_randomized` | loses any spatial-coherence advantage of `graph_ae` |
| `pca` vs `rna_ae` | comparable unless the non-linear model adds real structure |

## 9. Metrics

- **Primary:** held-out **AUROC** of a logistic probe for "program-high" spots (binarised at the
  per-section 80th percentile), averaged over the 11 programs, per held-out patient.
- **Secondary (reported, not the basis of the primary conclusion):**
  - silhouette score of held-out embeddings against `compartment` labels;
  - spatial coherence = mean Moran's I of embedding dimensions on the held-out section graph;
  - cross-seed embedding stability (mean absolute correlation between matched dimensions).
- Interpretive guard: spatial coherence is **not** evidence of biological superiority when the
  graph is included by construction; it is reported as a manipulation check.

## 10. Uncertainty and statistics

- 5 independent seeds (0–4) for stochastic methods; deterministic methods use seed 0.
- Per-patient metric values are averaged over seeds, then summarised as mean ± 95 % CI
  (t-based over the 6 patients) and per-patient distributions.
- Method comparisons are **paired by patient** against the `pca` and `rna_ae` references, using
  the corrected exact two-sided sign-flip test (`spatialrx.stats.exact_perm_paired`); effect
  size = paired mean difference / pooled SD.
- Multiple-testing: BH across method×metric comparisons within the primary analysis.
- Variance across seeds is reported separately from variance across patients.

## 11. Validation vs final testing

- **Validation** (model-free hyperparameters): fixed a priori; probe `C` chosen by inner LOPO on
  training patients only.
- **Final test:** the outer held-out patient, touched once per fold; aggregated exactly once
  after all methods complete.
- No metric is optimised against the held-out test data.

## 12. Acceptance criteria (benchmark validity)

The benchmark is **valid** only if all hold:

1. All primary folds are patient-disjoint (asserted).
2. Preprocessing test-invariance passes (asserted by test).
3. `random` and label-shuffle controls are near chance (95 % CI contains 0.5).
4. Cross-seed SD of the primary metric is reported and small relative to the patient spread.
5. Every method is evaluated on identical partitions.

The scientific conclusion will be stated **only** from the leakage-free evaluation. If
protein/spatial context does not improve the primary metric, that negative result is reported
transparently and the original Module-08 conclusion is not preserved.

## 13. Deliverables

`benchmark/results.csv` · `benchmark/summary.csv` · `benchmark/ablation_results.csv` ·
`benchmark/negative_controls.csv` · `benchmark/split_manifest.csv` ·
`benchmark/methods.md` · `benchmark/statistical_analysis.md` · this contract.
Historical Module-08 outputs are preserved unchanged in `outputs/tables/08_embedding_metrics.csv`.

## 14. Declared limitations (known before running)

- n = 6 patients → low power; CIs will be wide.
- Program scores are derived from RNA, so the primary probe is not fully independent of the
  RNA modality; it measures transferable representational structure, not external ground truth.
- Pre sections are core biopsies and post sections are whole sections; the benchmark is run
  within this heterogeneous collection and per-patient spread will reflect it.
- Protein panel is 35 features (31 markers + 4 isotypes); isotypes are included as features.

---

## Addendum A — compute-driven reductions (recorded before the final run)

The benchmark host is shared and was running at load ≈ 68 on 32 cores during preparation.
To keep the design executable without changing any leakage or evaluation rule, the following
*capacity* parameters were reduced equally for all methods:

| Parameter | Contract | Executed | Rationale |
|---|---|---|---|
| AE hidden width | 128 | **64** | halves per-fit cost |
| AE epochs | 200 | **50** | identical for every learned method |
| Primary seeds | 5 | **3** | repeated independent runs retained |
| Secondary seeds | 5 | **2** | exploratory arm |
| Parallel workers | — | 8 | avoid oversubscribing the shared host |

These changes affect the *capacity* of the learned representations and the precision of the
seed-variance estimate; they do **not** affect split construction, training-only fitting, the
primary metric, the label-shuffle control, or the paired statistical test. Every remaining
requirement of the contract (patient-held-out primary, section-held-out secondary, PCA /
spatial-smoothing baselines, degree-preserving graph randomisation, identical partitions,
separation of validation and test) is preserved.
