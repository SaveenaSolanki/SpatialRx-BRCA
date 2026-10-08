# Benchmark methods

**Contract:** `BENCHMARK_CONTRACT.md` (+ Addendum A, compute-driven reductions)
**Runner:** `benchmark/run_benchmark.py` · **Library:** `spatialrx/benchmark.py`
**Results:** `benchmark/results_auroc.csv`, `results_embedding.csv`, `summary.csv`,
`ablation_results.csv`, `negative_controls*.csv`, `stability.csv`, `split_manifest.csv`

## Data

Frozen per-sample objects (never re-ingested): `data/processed/03_domains/<specimen>_domains.h5ad`
(full-gene log-normalised `.raw`, spatial coordinates, compartment labels) and
`data/processed/05_ecosystem/<specimen>_ecosystem.h5ad` (`protein_clr`, 11 program scores).
14 sections, 6 patients, 12 paired + 2 replicates.

## Splits

- **Primary — Leave-One-Patient-Out (LOPO), 6 folds.** Replicate sections travel with their patient.
- **Secondary — Leave-One-Section-Out (LOSO), 14 folds.** Exploratory; patient leakage possible.
- `assert_disjoint` enforces no patient/section overlap between train and test.
- Every method is evaluated on the same `split_manifest.csv`.

## Preprocessing (fit on training folds only)

`spatialrx.benchmark.Preprocessor`:
1. top-2000 `seurat` HVGs on the concatenated **training** log-normalised matrix;
2. `StandardScaler` + clip ±10 on training spots;
3. `PCA(64)`;
4. protein `protein_clr` standardised on training spots.

Test sections are only ever transformed. A regression test
(`test_preprocessor_fit_is_invariant_to_test_data`) asserts that fitted parameters do not
change when test data changes.

## Methods (identical folds, 32-d embeddings)

| ID | Definition |
|---|---|
| `pca` | first 32 dims of the training-fitted RNA PCA(64) |
| `smoothing` | row-normalised within-section graph applied to the RNA PCA(64), first 32 dims |
| `random` | Gaussian random linear projection of RNA PCA(64) → 32 |
| `rna_ae` | MLP autoencoder on RNA PCA(64) (64 hidden, 50 epochs, Adam 1e-3) |
| `rna_protein_ae` | MLP AE on RNA PCA(64) ⊕ protein(35) |
| `graph_ae` | graph-conv pre-layer (A_hat X W) + MLP AE on RNA ⊕ protein, real grid graph |
| `graph_ae_randomized` | as `graph_ae`, adjacency degree-preserving double-edge swapped per section |

No model or preprocessing parameter is fit on the held-out unit. The probe regularisation `C`
is chosen by an inner leave-one-training-patient-out on training data only, once per
method/fold, and reused across seeds.

## Metrics

- **Primary:** held-out AUROC of a logistic probe for program-high spots (per-section 80th
  percentile), averaged over 11 programs, per held-out patient.
- **Secondary:** silhouette against compartment labels; mean Moran's I of embedding dims on
  the held-out section graph.
- **Stability:** mean |Pearson r| between matched embedding dims of seed 0 and later seeds.

## Controls

- **label-shuffle probe** — the proper chance baseline for the primary metric.
- **`graph_ae_randomized`** — degree-preserving graph null.
- **`random`** — a linear baseline (see the note in `statistical_analysis.md`: because it is a
  linear projection of the PCA space it is **not** a chance control, despite the contract's
  expectation).

## Statistics

Per-patient values (mean over seeds/folds) → mean ± 95 % t-CI over the 6 patients; seed SD
reported separately. Method comparisons are paired by patient against `pca` and `rna_ae`,
using the corrected exact two-sided sign-flip test (`spatialrx.stats.exact_perm_paired`) with
BH across comparisons.
