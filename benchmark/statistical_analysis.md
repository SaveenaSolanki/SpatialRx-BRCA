# Statistical analysis — leaky-free Module-08 benchmark (Gate 4)

**Branch:** `fix/repro-correctness-v1`
**Contract:** `benchmark/BENCHMARK_CONTRACT.md` (+ Addendum A)
**Runner:** `python benchmark/run_benchmark.py --analyses patient section --seeds 3 --jobs 6`
(primary seeds = 3; secondary seeds = 2, as declared in Addendum A)
**Runtime:** primary 6 m 44 s; secondary 7 m 46 s (shared host at load ≈ 68 / 32 cores)

All numbers below are produced by the runner; `benchmark/*.csv` are the machine-readable
artefacts. The historical, leaky Module-08 table
(`outputs/tables/08_embedding_metrics.csv`) is preserved unchanged as evidence.

---

## 1. Benchmark validity checklist

| Criterion (contract §12) | Result |
|---|---|
| Primary folds patient-disjoint | **PASS** — `assert_disjoint`; 6 LOPO folds, 14 LOSO folds in `split_manifest.csv` |
| Preprocessing test-invariance | **PASS** — `test_preprocessor_fit_is_invariant_to_test_data` |
| Label-shuffle control ≈ chance | **PASS** — patient 0.503 (95 % CI 0.487–0.519); section 0.493 (0.483–0.502) |
| Graph-randomisation control ≈ chance | **PASS** — AUROC 0.505 / 0.501; Moran's I 0.001 / −0.002 |
| Cross-seed SD reported | **PASS** — AUROC seed SD 0.0025–0.0107 |
| Identical partitions for all methods | **PASS** — one `split_manifest.csv` |
| `random` projection ≈ chance | **NOT MET (expected design flaw, not a validity failure)** — see §6 |

---

## 2. Primary metric — held-out AUROC (program-high detection)

Mean over patients with 95 % t-CI (higher is better).

| Method | Patient-held-out | Section-held-out |
|---|---|---|
| **pca** | **0.7035** (0.6807–0.7262) | **0.7082** (0.6725–0.7440) |
| smoothing | 0.6862 (0.6590–0.7134) | 0.6915 (0.6539–0.7290) |
| random (linear projection) | 0.6821 (0.6596–0.7046) | 0.6852 (0.6514–0.7190) |
| rna_ae | 0.6763 (0.6528–0.6999) | 0.6792 (0.6453–0.7131) |
| rna_protein_ae | 0.6745 (0.6494–0.6996) | 0.6745 (0.6420–0.7070) |
| graph_ae | 0.6334 (0.5883–0.6784) | 0.6325 (0.5905–0.6745) |
| graph_ae_randomized | 0.5053 (0.4999–0.5107) | 0.5008 (0.4960–0.5057) |

The ordering is identical in the primary and secondary analyses. **The linear RNA baseline
(`pca`) is the strongest method; no learned representation beats it.**

## 3. Paired comparisons against `pca` (primary, n = 6 patients)

| Method | Δ AUROC | exact p | BH p | effect size |
|---|---|---|---|---|
| graph_ae | −0.0701 | 0.0312 | 0.0312 | −2.06 |
| graph_ae_randomized | −0.1982 | 0.0312 | 0.0312 | −12.58 |
| rna_protein_ae | −0.0290 | 0.0312 | 0.0312 | −1.27 |
| rna_ae | −0.0272 | 0.0312 | 0.0312 | −1.23 |
| random | −0.0214 | 0.0312 | 0.0312 | −0.99 |
| smoothing | −0.0173 | 0.0312 | 0.0312 | −0.72 |

Every method is worse than `pca` in all 6 patients (exact two-sided p = 2/64 = 0.0312;
BH-significant). Section-held-out gives the same conclusion.

## 4. Multimodality / spatial-context effect sizes

**Does protein help?** `rna_protein_ae − rna_ae` = **−0.0018** (p = 0.25, 95 % CI includes 0):
no evidence of improvement from adding the 35-protein layer; if anything a negligible decrease.

**Does the spatial graph help the primary metric?** `graph_ae − rna_ae` = **−0.0430**
(p = 0.0312, BH 0.0417): adding graph message passing **significantly reduces** transferable
niche separation, even though it raises spatial coherence by construction.

## 5. Secondary metrics

Spatial coherence (mean Moran's I of embedding dims), patient-held-out:

| Method | Moran's I |
|---|---|
| graph_ae | **0.5250** (0.375–0.676) |
| rna_protein_ae | 0.4438 |
| smoothing | 0.4429 |
| rna_ae | 0.4241 |
| pca | 0.2730 |
| random | 0.2615 |
| graph_ae_randomized | **0.0013** |

`graph_ae` has the highest coherence — but this is a **manipulation check, not evidence of
biological superiority**: the graph is injected by construction, and the same model is the
worst on the primary metric. The degree-preserving randomisation drives coherence to ≈ 0,
confirming the effect is graph-driven.

Silhouette (vs compartment labels) is low for all methods (≈ 0.01–0.05); `pca` and `smoothing`
are nominally highest; graph models are near 0. No method separates compartments well.

**Seed stability** (mean |r| vs seed 0): rna_ae 0.36, rna_protein_ae 0.37, graph_ae 0.36,
random 0.18, graph_ae_randomized 0.15 — learned embeddings are only moderately stable across
seeds, reinforcing that cross-seed SD (≤ 0.011 AUROC) must accompany every number.

## 6. The `random` projection control did not behave as expected

The contract expected a Gaussian projection of the RNA PCA space to be near chance. It scores
0.682 (patient) / 0.685 (section) — close to the AE methods. This is **not** a validity failure:
a random linear projection approximately preserves the linear structure of the PCA space
(Johnson–Lindenstrauss), and the probe is itself linear, so it can still read the program
signal. The correct chance control for this probe is the **label shuffle**, which is ≈ 0.5.
We therefore reclassify `random` as a *linear baseline*, not a null model, and report the
mis-specification openly.

## 7. Scientific interpretation (leakage-free)

1. **The original Module-08 negative result is confirmed and strengthened.** Under a
   patient-held-out, contamination-free evaluation, neither the protein layer nor a spatial
   graph nor a non-linear autoencoder improves transferable niche separation over a plain
   linear RNA PCA baseline.
2. **More capacity/modalities did not help; the spatial graph actively hurt** the primary
   metric while mechanically inflating a spatial-coherence metric — exactly the confound the
   original brief warned about.
3. **Spatial coherence is not evidence of biological benefit** here: the graph model with the
   best coherence is the worst on the held-out biological metric, and its randomised-graph
   twin drops to chance on both.
4. At n = 6 patients (12 paired sections) these comparisons are powered to detect only large
   effects, and the consistent, BH-significant degradation of learned methods relative to PCA
   is the clearest available signal.

**Consequently the original claim that a multimodal/spatial embedding recovers biologically
coherent spatial states better than RNA alone is not supported and is superseded by this
analysis.**

## 8. Deviations from the pre-registered contract

| # | Deviation | Reason / impact |
|---|---|---|
| 1 | AE hidden 128 → 64; epochs 200 → 50; primary seeds 5 → 3; secondary seeds 5 → 2 | Addendum A — shared-host load. Applies equally to all learned methods; no split/metric/leakage rule changed. |
| 2 | `random` treated as a linear baseline, not a chance control (its AUROC was not ≈ 0.5) | §6 — the original expectation was mis-specified; the true null is the label shuffle. |
| 3 | Probe `C` selected once per method/fold on the first seed and reused across seeds | Avoids repeating identical inner selection; still training-only. |
| 4 | `graph_ae_randomized` uses the same randomisation for train and test sections (independent per section, degree-preserving) | Keeps the control well defined; no cross-section edges are ever created. |

## 9. Limitations

- n = 6 patients; CIs are wide and the design is underpowered for small effects.
- Program scores are derived from RNA, so the probe is not fully independent of the RNA
  modality; it measures transferable representational structure rather than external ground truth.
- Pre sections are core biopsies, post sections whole surgical specimens; per-patient spread
  reflects this.
- Only 35 (mostly surface) proteins and a hexagonal spot grid are available; single-cell
  resolution is not.
- The 64-hidden / 50-epoch budget is lower than the original 128 / 300; the negative result is
  robust across both the reduced budget and the cross-patient/section designs, but a larger
  model could in principle change the AE numbers (not the PCA baseline).

## 10. Reproducibility

- Seeds: 0,1,2 (primary); 0,1 (secondary). Deterministic methods use seed 0.
- Per-fold checkpoints: `benchmark/folds/*.pkl` (resumable; `--no-resume` to recompute).
- Run provenance: `benchmark/run_record.json`; git commit recorded in the gate report.
- Command: `python benchmark/run_benchmark.py --analyses patient section --seeds 3 --jobs 6`
  (secondary re-run with `--seeds 2`).
