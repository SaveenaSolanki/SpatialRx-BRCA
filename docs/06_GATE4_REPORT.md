# 06 — Gate 4 report: leakage-free Module-08 benchmark

**Branch:** `fix/repro-correctness-v1`
**Starting point:** `1d7d2f7e72f66b34db6eb5b3050ed3b409efbed9` (Stage A reconciled)
**Contract:** `benchmark/BENCHMARK_CONTRACT.md` (+ Addendum A)
**Detailed analysis:** `benchmark/statistical_analysis.md` · **Methods:** `benchmark/methods.md`

## What changed

The original `experiments/08_ai.py` metric fitted and scored a logistic probe on the **same
spots**, so it could not support a generalization claim (issue ISSUE-004). It is **left
unchanged and preserved as historical evidence**. A new, independent benchmark now evaluates
representations under contamination-free splits.

- `spatialrx/benchmark.py` — tested primitives: splits, training-only preprocessing,
  baselines, AE/GCN models, degree-preserving graph randomisation, metrics, controls.
- `benchmark/run_benchmark.py` — resumable, checkpointed, parallel runner.
- `benchmark/results_*.csv`, `summary.csv`, `ablation_results.csv`, `negative_controls*.csv`,
  `stability.csv`, `split_manifest.csv`, `run_record.json`.
- `tests/test_benchmark.py` — 12 tests (splits, graph controls, probe null, preprocessing
  invariance, shapes).
- `docs/ISSUES.md` — tracked-issue register (HLA-DRA parsing defect recorded `OPEN`).

## Design (leakage boundaries)

Primary = **leave-one-patient-out** (6 folds); secondary = **leave-one-section-out** (14 folds).
HVG selection, scaling, PCA, protein scaling and every model are fit on training folds only;
the test section is only transformed. Probe `C` is chosen by inner leave-one-training-patient-out.
Every method uses identical partitions (`assert_disjoint`). The spatial graph is within-section
only; the randomised-graph control preserves degree and never adds cross-section edges.

## Results (mean over patients, 95 % CI)

Primary held-out AUROC (higher better):

| Method | Patient-held-out | Section-held-out |
|---|---|---|
| **pca** | **0.7035** | **0.7082** |
| smoothing | 0.6862 | 0.6915 |
| random (linear baseline) | 0.6821 | 0.6852 |
| rna_ae | 0.6763 | 0.6792 |
| rna_protein_ae | 0.6745 | 0.6745 |
| graph_ae | 0.6334 | 0.6325 |
| graph_ae_randomized | 0.5053 | 0.5008 |

- Every method is worse than `pca` in all 6 patients (exact p = 2/64 = 0.0312; BH-significant).
- Protein adds nothing: `rna_protein_ae − rna_ae` = −0.0018 (p = 0.25).
- Spatial graph hurts the biological metric: `graph_ae − rna_ae` = −0.043 (p = 0.031).
- Spatial coherence (Moran's I): graph_ae 0.525 ≫ pca 0.273, but the randomised-graph twin is
  0.001 — a construction effect, not biological superiority.
- Controls: label-shuffle AUROC 0.503 / 0.493 (CI contains 0.5); degree-preserving graph
  randomisation AUROC ≈ 0.50 and Moran ≈ 0.

## Interpretation

**The original negative finding is confirmed and strengthened.** Under a patient-held-out,
contamination-free evaluation, neither the protein layer, the spatial graph, nor a non-linear
autoencoder improves transferable niche separation over linear RNA PCA. The graph model's
advantage on spatial coherence is mechanical. This supersedes the original Module-08
conclusion; the original claim is not preserved.

## Known failures / caveats

- The contract's expectation that a random projection would be a chance control was wrong
  (AUROC ≈ 0.68): a Gaussian projection of PCA space preserves linear structure. The true null
  is the label shuffle; documented in `statistical_analysis.md` §6.
- Capacity reductions (hidden 128→64, epochs 200→50, seeds 5→3/2) due to shared-host load;
  declared in Addendum A and applied equally to all learned methods.
- n = 6 patients; wide CIs; the program probe is RNA-derived (transferable structure, not
  external ground truth).
- ISSUE-001 (HLA-DRA parser) and ISSUE-003 (`cluster_zscore_table` NaN rows) remain `OPEN`.

## Reproducibility

```
python benchmark/run_benchmark.py --analyses patient section --seeds 3 --jobs 6
# secondary arm declared with --seeds 2
```
Runtime: primary 6 m 44 s, secondary 7 m 46 s (shared host, load ≈ 68/32 cores). Per-fold
checkpoints (`benchmark/folds/`, gitignored) allow exact resume. 63 tests pass.
