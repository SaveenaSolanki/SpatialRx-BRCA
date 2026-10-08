# 01 — Scientific Claim Ledger

**Audited commit:** `fd3de57` · **Audit date:** 2026-10-08 · **Branch:** `audit/phase0`

Every quantitative or interpretive statement in
`outputs/SpatialRx-BRCA_report.md` and `outputs/SpatialRx-BRCA_next_steps.md` is mapped
below to the script that produces it, its inputs, its output artifact, the method, and the
result of this audit's verification.

**Verification codes**

- `RECOMPUTED ✓` — the number was recomputed from stored artifacts during this audit and matched.
- `MATCHES TABLE` — the report number agrees with its committed output table (traceable).
- `DISCREPANCY` — the report and/or table is internally inconsistent or wrong.
- `NOT RE-EXECUTED` — verified only at the artifact level; the module was not re-run this audit.
- `BLOCKED` — cannot be verified from repository content.

---

## A. Cohort and data claims

| ID | Claim (report) | Script / source | Evidence artifact | Method | Status | Audit finding |
|---|---|---|---|---|---|---|
| C01 | GSE289326; 14 GSM → 12 unique sections + 2 replicates; 6 patients; S1 pre-biopsy, S2 post-pET | `srx.SAMPLE_META`; `data/README.md` | `00_sample_qc.csv` | metadata table | `MATCHES TABLE` | consistent; 14 rows, 45,326 spots sum ✓ |
| C02 | 18,085 genes + 35 antibody features per section | `00_ingest.py` (`feature_types`) | each `*.h5ad`, `protein_var` | feature-type split | `RECOMPUTED ✓` | exactly 18,085 / 35 |
| C03 | 45,326 spots ingested | `00_ingest.py` | `00_sample_qc.csv` | sum | `RECOMPUTED ✓` | sum = 45,326 |
| C04 | Atlas 39,428 spots × 3,000 HVGs | `03_tumor_heterogeneity.py` | `atlas_merged.h5ad` | concat + seurat HVG | `RECOMPUTED ✓` (shape) | shape correct, but `obs["sample"]` corrupt (D1) |
| C05 | Pre sections 294–2,256 spots; post 3,598–4,992 | `03` outputs | `08_embedding_metrics.csv` `n_spots` | per-section QC | `RECOMPUTED ✓` | range correct |
| C06 | Replicate sections excluded from paired tests | `srx.PAIRED` | `PAIRED` = 12 sections | constant list | `RECOMPUTED ✓` | replicates absent from `PAIRED` |

---

## B. Spatial autocorrelation (Module 04)

| ID | Claim | Script | Output | Method | Status | Audit finding |
|---|---|---|---|---|---|---|
| C07 | Our Moran's I correlates with the study's `spatial_enrichment.csv`, r = 0.70–0.94, mean ≈ 0.83 | `04_spatial_stats.py` | `04_morans_I_per_sample.csv` + authors' CSVs | Pearson of per-gene I | `RECOMPUTED ✓` | 0.695–0.938, mean **0.833** |
| C08 | Top SVGs: IGKC 0.496, COL1A1 0.487, COL3A1 0.469, COL1A2 0.441, MUC1 0.380 … | `04` | `04_svg_ranking.csv` | mean Moran's I, 14 samples | `MATCHES TABLE` | exact match |
| C09 | 6-neighbour grid graph, array coords, FDR-BH, permutations for top-500 genes | `04` | code | Squidpy `spatial_autocorr` | `MATCHES TABLE` | documented in code; graph choice not compared to authors' graph |

**Method caveat to record:** the study's `spatial_enrichment.csv` column is named `I` and is
**presumed** to be Moran's I. The authors' graph definition, weighting and normalisation are
not available in the repository, so r≈0.83 establishes *metric agreement*, not an exactly
matched reimplementation. Per-sample protein agreement is far more variable (0.22–0.94).

---

## C. RNA ↔ protein agreement and discordance (Modules 05, 06, 06b)

| ID | Claim | Script | Output | Method | Status | Audit finding |
|---|---|---|---|---|---|---|
| C10 | 27 matched markers; median Pearson r: ACTA2 0.19, VIM 0.14, EPCAM 0.12 … CD4 −0.03 | `06_multimodal.py` | `06_rna_protein_agreement.csv` | per-spot RNA log-CPM vs CLR protein | `DISCREPANCY` | **26** markers analysed (HLA-DRA dropped, D10); reported medians otherwise exact |
| C11 | Discordance (zRNA−zProt) spatially structured, mean Moran 0.219, 99.5% positive | `06` | `06_discordance_morans.csv` | Squidpy Moran on obs | `RECOMPUTED ✓` | 0.219, 0.995 (n=364 rows) |
| C12 | 34–35/35 proteins spatially structured (FDR<0.05) in every section | `06b_protein_spatial.py` | `06b_protein_morans.csv` | CLR + Moran | `DISCREPANCY` | true range **30–35** (P4-S1=30, P5-S1=31); report overstates (D9) |
| C13 | Protein validation vs authors r = 0.69 | `06b` | `06b_protein_morans.csv` + authors' CSV | pooled Pearson | `RECOMPUTED ✓` | pooled **0.694**; per-sample 0.22–0.94 |
| C14 | RNA vs protein NNLS deconvolution agrees best for epithelial/tumour | `05_ecosystem.py` | `05_rna_prot_deconv_agreement.csv` | Pearson per lineage/sample | `MATCHES TABLE` | heuristic NNLS; not an independent deconvolution |
| C15 | Isotype-normalisation robustness check available | `05`/`06` docstrings, `srx.py` | — | — | `DISCREPANCY` | `clr_by_factor` defined but **never called**; as written it is a mathematical no-op (D8) |

---

## D. Endocrine-therapy remodelling (Module 07)

| ID | Claim | Script | Output | Method | Status | Audit finding |
|---|---|---|---|---|---|---|
| C16 | Compartment composition pre→post: Tumour 0.091→0.161, Immune 0.180→0.093, Stroma 0.148→0.235, Mixed 0.580→0.511 | `07_treatment.py` | `07_paired_metrics.csv` | mean spot fraction | `MATCHES TABLE` | exact |
| C17 | No metric survives BH at n=6 | `07` | `07_paired_tests.csv` | sign-flip permutation + BH | `MATCHES TABLE` | all `p_perm_BH` = 1.0 |
| C18 | B-cell activity trend p = 0.094 (uncorrected) | `07` | `07_paired_tests.csv` | exact permutation | `DISCREPANCY` | stored `p_perm` = **0.1875** (doubled, D3); 0.094 is the Wilcoxon value / un-doubled exact p |
| C19 | Tumour–immune contact −0.40 (ES −1.88, n=3) | `07` | `07_paired_tests.csv` | contact fraction | `MATCHES TABLE` | Δ −0.398, ES −1.884, n=3 |
| C20 | Program Moran shifts: fibroblast activation +0.160, EMT +0.114, tumour deconv +0.117, T-cell −0.049 | `07` | `07_program_morans_paired_tests.csv` | paired Δ Moran's I | `MATCHES TABLE` | exact |
| C21 | "22 metrics" tested | `07` | `07_paired_tests.csv` | count | `DISCREPANCY` | **24** rows (D11) |
| C22 | Patient, not spot, is the inference unit | `07` code | code | per-patient pairing / permutation | `RECOMPUTED ✓` | design is correct |

**Method caveat:** the sign-flip test is appropriate for n=6, but `exact_perm_paired` applies a
factor of 2 to an already two-sided |stat|-based count; every stored `p_perm` is therefore
inflated by ≈2× relative to the correct exact value. This does not change the BH conclusion
(all 1.0) but it makes the reported "uncorrected trend" numbers unreliable.

---

## E. Multimodal AI benchmark (Module 08)

| ID | Claim | Script | Output | Method | Status | Audit finding |
|---|---|---|---|---|---|---|
| C23 | Silhouette 0.034 / 0.032 / 0.064 (RNA / +protein / +graph) | `08_ai.py` | `08_embedding_metrics.csv` | within-sample | `RECOMPUTED ✓` | means match |
| C24 | Spatial coherence (Moran of dims) 0.281 / 0.277 / 0.656 | `08` | `08_embedding_metrics.csv` | closed-form Moran | `RECOMPUTED ✓` | means match; **by construction** |
| C25 | Niche separation AUROC 0.765 / 0.766 / 0.752 — no gain from protein/spatial | `08` | `08_embedding_metrics.csv` | logistic regression on embedding | `DISCREPANCY / INVALID` | classifier is **fit and scored on the same spots**, and labels come from the same RNA input → leaky, circular metric (D4) |
| C26 | "Honest negative" for more-modalities-better-niches | report §2/§8 | — | interpretation | `PARTIAL` | the *direction* (no gain) may still hold, but the number cannot support a generalization claim until the metric is fixed with a proper split |
| C27 | Transductive-only, no patient-level holdout, no graph-edge-randomisation control | — | — | — | `NOT_IMPLEMENTED` | required by Phase 4; only 3 self-supervised configs exist |

---

## F. Next-steps analyses (Modules 09–12)

| ID | Claim | Script | Output | Method | Status | Audit finding |
|---|---|---|---|---|---|---|
| C28 | CITEgeist Step-2 QP validates NNLS: epithelial r=0.76 | ad-hoc (no script) | `data/processed/09_citegeist/*.csv` | GPU QP (cuOPT) | `NOT RE-EXECUTED` | no generating script committed → not reproducible (D13); `cg_step2.log` shows an AnnData write error |
| C29 | Discordance is not edge-driven (ρ<0 for 7/8 markers) | `11_discordance_cellstate.py` | `11_discordance_cellstate_summary.csv` | Spearman vs edge exposure | `MATCHES TABLE` | implemented |
| C30 | Discordance structured by niche (epithelial/myeloid-rich, fibroblast-poor) | `11` | `11_hotspot_state_deltas.csv` | hotspot composition deltas | `MATCHES TABLE` | implemented |
| C31 | LR: MIF–CD74 and MDK–LRP1 top; CXCL12–CXCR4 ↓; CSF1–CSF1R/TGFB1–TGFBR1 ↑ (n=4) | `10_lr_analysis.py` | `10_lr_scores/summary/paired_tests.csv` | curated list, product score | `MATCHES TABLE` | no traceable LR resource (curated by hand); score is a heuristic, not a communication test |
| C32 | Tumour–immune contact decrease survives size-matched subsampling | `12_subsampling.py` | `12_subsampling_summary.csv` | random + contiguous BFS subsampling | `MATCHES TABLE` | implemented; correctly conservative |
| C33 | Powered cohort scan: no exact-match public dataset | report / next_steps | narrative | literature/GEO scan | `NOT RE-EXECUTED` | no machine-readable search manifest retained |
| C34 | Power estimates (13–26 pairs etc.) | next_steps | narrative | post-hoc power | `NOT RE-EXECUTED` | no script; arithmetic not reproducible from repo |

---

## G. Hypothesis vs evidence (interpretation audit)

| Statement type | Examples | Verdict |
|---|---|---|
| Measured association | Moran validation, RNA–protein r, discordance Moran, compartment fractions | Evidence-backed |
| Exploratory / underpowered | all pET paired trends; contact change (n=3) | Correctly labelled "hypothesis-generating" |
| Negative result | multimodal niche gain | Honest, but metric must be fixed before publication |
| Mechanistic inference | "post-transcriptional regulation", LR signalling | Correctly hedged in report; **not proven** |
| Causal treatment effect | "therapy remodels" | Observational only; report says so. Headline sentence ("endocrine therapy remodels…") should be softened to "is associated with" |
| Clinical utility / response prediction | none claimed | Good — no overclaim found |

---

## H. Ledger summary

| Verification result | Count |
|---|---|
| Recomputed and matched | 9 |
| Matches committed table | 10 |
| Discrepancy / invalid method | 6 (C10, C12, C15, C18, C21, C25) |
| Not re-executed / not reproducible from repo | 4 (C28, C33, C34 + C13 per-sample) |
| Not implemented (required by brief) | 1 (C27) |

**Bottom line:** the report is largely faithful to its own output tables, but six claims have
provenance or validity problems, and two of them (C18 permutation p, C25 niche AUROC) touch
headline conclusions. All are fixable without changing the underlying data.

---

## I. Remediation updates — branch `fix/repro-correctness-v1`

Baseline for comparison: `47052bd` (frozen Phase-0 audit).

### Gate 1 — reproducibility (no claim changes)

- Canonical environment pinned (`environment.yml`, `requirements.lock`); h5ad writes now work.
- `spatialrx validate-data` passes 31/31 checks (strict, sha256 verified); `reproduce --smoke` passes.
- Found and fixed a latent metadata bug: `patient_of` used dash-field index 3 and returned the
  section token (e.g. `S2_D`) instead of the patient. Unused by the old scripts; now covered by a
  regression test. No reported claim depended on it.

### Gate 2 — statistical correctness (claims changed)

| Claim | Before (committed table) | After (corrected) | Reason |
|---|---|---|---|
| C18 B-cell paired trend p | `p_perm = 0.1875` | **0.09375** | `exact_perm_paired` double-counted both tails (`2·s/N`); correct exact two-sided p is `s/N`. The report prose already quoted 0.094, so the table now matches the text. |
| C17/C21 all `07_paired_tests` p | each `2·s/N` | each `s/N` | same fix; all deltas unchanged; no metric becomes BH-significant (`p_perm_BH` min 0.900). |
| C20 program Moran p | each `2·s/N` | each `s/N` | same fix; BH now uses a proper monotone step-up (`p_BH` min 0.660); no program change is significant. |
| C31 CSF1–CSF1R Δ+0.33 | `p = 0.0` (table) | **0.125** | The committed `0.0` is impossible for an exact two-sided sign-flip test (the observed sign vector is always in the null set); the committed script actually yields 0.25, and the corrected value is 0.125. The stored table was **not reproducible** from the committed code (provenance defect). |
| C31 TGFB1–TGFBR1 Δ+0.28 | `p = 0.0` | **0.125** | same. |
| C31 LR paired tests | no correction | added `p_BH` column | multiple-testing correction was absent; all BH values ≥ 0.81. |

**Unchanged by remediation:** every effect-size / delta in modules 07 and 10, all `07_paired_metrics`
values, all program Moran's I values (to floating precision), and the overall conclusion that no
paired pre/post change survives multiple-testing correction at n = 6. The `10_lr` deltas are
unchanged; only their p-values were corrected.

Pending (post-Gate 4 reconciliation): C10 marker count (26 vs 27), C12 protein-FDR range (30–35),
C21 metric count (24 not 22), C25 niche-AUROC validity, C15 isotype robustness removal.
