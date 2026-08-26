# SpatialRx-BRCA — Next steps executed (2026-08-26)

Four follow-up studies proposed in the main report were executed on the same GSE289326
objects. All numbers trace to `outputs/tables/`, `outputs/figures/`, or `logs/`.

---

## 1. Powered multi-centre cohort — feasibility scan + power estimates

**Finding: no public dataset replicates the GSE289326 design** (paired pre/post primary
endocrine therapy, ER+, Visium CytAssist + same-slide antibody capture). Search date
2026-08-26; sources in report §Sources.

| Candidate | Design | Why it does/doesn't power our question |
|---|---|---|
| **GSE331245 "SMART"** (89 TNBC patients, 129 samples, pre/during/post neoadjuvant chemotherapy, GeoMx DSP) | powered, matched pre/post | The only powered pre/post spatial-treatment cohort; different subtype (TNBC), therapy (chemo), modality (ROI GeoMx, no antibody capture). Best near-term option to test whether *treatment-associated spatial-structure shifts* (our Module 4/7 metrics) replicate in an independent cohort; epidemiology (ER+ pET) not transferred. |
| **ImmunoADAPT trial (ASCO 2025 #e12624)** | prospective ET±palbociclib±avelumab, ER+, ST-profiled | Closest clinical biology (ER+ endocrine therapy + immune); data not yet public. |
| **GSE299880** (GeoMx, CD8+ T cells ↔ endocrine-therapy resistance, ER+) | treated, cross-sectional | No matched pre/post spatial pairs. |
| **GSE111563** (matched sequential bulk RNA on letrozole; dormant/resistant) | matched, longitudinal | Bulk only — paired-biopsy design proof that same-trial cohorts exist, not spatial. |
| **SMMART serial CosMx biopsies** (4 metastatic patients, 10 timepoints) | single-cell spatial, serial | n=4, metastatic setting, no antibody capture. |
| 10x public Visium/HD samples | single sections | Not cohorts. |

**Power needed** (paired n per arm, two-sided α=0.05, 80%; SDs taken from our data):

| Effect we observed (n=6) | Observed Δ (SD of change) | Paired patients needed |
|---|---|---|
| Moran's I shift, fibroblast program | 0.16 (≈0.20) | **13** |
| Moran's I shift, EMT program | 0.11 (≈0.20) | **26** |
| Tumour–immune contact ↓ | 0.40 (≈0.21; ES≈1.9) | **3** |
| Tumour–immune contact ↓ (conservative) | 0.20 (≈0.30) | **18** |
| B-cell program ↑ | 0.086 (≈0.10) | **11** |

Conclusion: a replication cohort of **≥26 matched biopsy→surgery pairs** (or the ongoing
NCT05914792 trial itself, which is collecting exactly this design, plus the SMART cohort
for cross-modality validation) would be sufficient to confirm/refute the spatial-structure
shifts at 80% power.

## 2. Discordance hotspots vs post-transcriptional regulation — spot-level test + CITEgeist validation

**Blocked sub-item (documented):** true *single-cell* assignment (CITEgeist Step 3) requires
StarDist nuclei segmentation of the H&E (patchwise, GPU) plus morphology artifacts; the
published package ships no StarDist weights (`CITEgeist/model/` checked), and GEO contains
no segmentation artifacts. Full pipeline needs a pre-trained StarDist model + the
repository's `run_single_cell_assignment.py`; not feasible to complete faithfully here.

**Done — CITEgeist Step 2 (spot-level proportions) ran with the study's own solver:**
GPU QP (cuOPT 26.08.00; Barrier solver, ~1.7 s/solve) on P4-S2, 12 cell types.
`data/processed/09_citegeist/P4-S2_cell_prop_finetuned_results.csv`.

**Validation of our transparent NNLS against the study's own tool (P4-S2, 4,658 spots):**

| Lineage | Pearson r (CITEgeist vs our RNA-NNLS) |
|---|---|
| epithelial/tumour (Luminal+Basal) | **+0.76** (Spearman +0.77) |
| myeloid (Mac+Monocyte+DC) | +0.34 |
| fibroblast | +0.26 |
| endothelial | +0.26 |
| T cells (CD8+CD4) | +0.24 |
| B cells | −0.11 |

The dominant tumour state is recovered almost identically by both approaches; differences in
mean abundances reflect marker-set choices (our fibroblast profile is COL1A1/DCN-heavy; the
study's Fibroblasts profile is ACTA2-only).

**Done — regional (spot-level) test of the discordance hypothesis (Module 11):**
per-spot discordance (zRNA − zProtein) was correlated with (a) edge-exposure (6 − grid
degree; the artifact control) and (b) cell-state proportions/programs, per section
(n=12 paired sections), median Spearman rho:

- **Artifact control: discordance is NOT edge-driven.** rho(discordance, edge-exposure) is
  negative or ≈0 for 7/8 markers (CD14 −0.17, CD3E −0.14, VIM −0.10, CD163 −0.10, EPCAM
  −0.07, CD68 −0.05, ACTA2 +0.01); the exception is **PCNA (+0.16)** — proliferation-marker
  discordance does peak near tissue edges (mitotic-zone biology or penetration effect).
- **Discordance is niche/cell-state structured:** for CD68/CD3E/CD163/CD14/EPCAM/VIM,
  discordance correlates positively with epithelial-tumour proportion (median rho
  +0.07…+0.31) and negatively with fibroblast proportion (−0.12…−0.24); hotspots (top
  quartile) are enriched in myeloid (CD68 hotspots: myeloid +3.7 pp, epithelial +2.7 pp,
  fibroblast −6.1 pp) and depleted of fibroblast state.
- M2-polarisation, exhaustion and hypoxia program scores do **not** drive discordance
  (rho ≈ −0.17…+0.12) — the signal is compartmental (epithelial/myeloid-rich vs
  fibroblast-rich zones), which is consistent with cell-state-specific differences in
  mRNA–protein coupling rather than a generic technical artefact, but not yet proof of
  post-transcriptional regulation per se.
- Tables: `11_discordance_cellstate.csv`, `11_discordance_cellstate_summary.csv`,
  `11_hotspot_state_deltas.csv`; figure `11_discordance_drivers.png`.

## 3. Ligand–receptor analysis at tumour–immune interfaces (Module 10)

Curated 26 LR pairs (immune↔tumour, tumour↔stroma; incl. the paper's MIF–CD74/CD44 and
MDK axes; no external DB), interface = Tumour spot with ≥1 Immune grid neighbour and vice
versa; interaction score = mean(sender ligand) × mean(receiver receptor) (log-CPM).

- Top interface scores (mean over 10 sections with interfaces): **MIF→CD74 ≈ 7–10,
  MDK→LRP1 4.3→7.3 (pre→post), CXCL12→CXCR4 5.5→2.8 (pre→post), COL1A1→SDC1 ≈ 6.3,
  CSF1→CSF1R, VEGFA→KDR/FLT1** (tables `10_lr_scores.csv`, `10_lr_summary.csv`).
- Interface vs non-interface enrichment (pre): CXCL12–CXCR4 ×1.95, LGALS9–HAVCR2 ×1.65,
  VEGFA–KDR ×5.2 — chemokine/checkpoint axes concentrated at the tumour–immune interface.
- Paired pre→post changes (n=4 patients with interface data: P1/P2/P4/P5): CSF1–CSF1R
  (tumour→myeloid) Δ+0.33, TGFB1–TGFBR1 Δ+0.28 (both p=0.000 in exact sign-flip tests,
  n=4); CXCL12–CXCR4 decreases −2.74 across the cohort mean. Hypothesis-generating only
  (no multiple-testing correction, n=4).
- Figure: `10_lr_interface_prepost.png`, `10_interface_map_P4S2.png`.

## 4. Matched-region subsampling for the treatment tests (Module 12)

Post sections (surgical, 3.6–5.0k spots) were subsampled to the paired pre-biopsy spot
count (294–2,256) by (a) uniform random and (b) contiguous "needle-core-like" BFS regions;
B=200/100 per patient; deltas recomputed on each subsample (`12_subsampling_deltas.csv`,
`12_subsampling_summary.csv`, `12_subsampling_per_patient.csv`).

| Metric | Full-section Δ | Random-sub Δ (median, 95% CI) | Contiguous Δ (median) | Robustness |
|---|---|---|---|---|
| frac_Tumour | +0.070 | +0.131 (−0.30, +0.28) | +0.078 | direction survives, CI overlaps 0 |
| frac_Immune | −0.088 | −0.124 (−0.32, +0.30) | −0.042 | direction survives, CI overlaps 0 |
| Estrogen program | −0.012 | −0.066 (−0.40, +0.36) | −0.023 | fragile |
| Proliferation / T-cell / Myeloid programs | ≈0 | ≈0 | ≈0 | null either way |
| **tumour–immune contact** | **−0.353** | **−0.501 (−0.80, −0.08); 0/200 iterations positive** | **−0.427 (≈80% negative)** | **survives size-matching** |

Interpretation: the pre/post tumour–immune contact decrease (P1 −0.10, P2 −0.50, P4 −0.79;
n=3 patients with both compartments) is **not** an artefact of comparing a small biopsy to a
large section — it persists when the post section is cut down to biopsy size (the
contiguous-core comparison is the geometry-fair one; uniform random subsampling mechanically
thins adjacency, so its CI is conservative). Other candidate effects (tumour/immune fraction,
estrogen program) are compatible with noise under size-matching. Figure `12_subsampling.png`.

---

## Deliverables & status

| Study | Status | Key files |
|---|---|---|
| Cohort scan + power | done (no exact-match public cohort; SMART/ImmunoADAPT candidates) | this report |
| CITEgeist Step 2 validation | done (r=0.76 tumour lineage) | `logs/cg_step2.log`, `09_citegeist/` |
| CITEgeist Step 3 single-cell | blocked (no StarDist weights/segmentation artifacts) | `logs/citegeist_setup.log` |
| Discordance × cell state | done (not edge-driven; niche-structured) | `11_*` |
| LR at interfaces | done (MIF/MDK top; CXCL12–CXCR4 ↓; CSF1/TGFB ↑) | `10_*` |
| Matched-region subsampling | done (contact effect survives) | `12_*` |

Environment note: installing `citegeist`/`cuopt-cu12` upgraded numpy→2.4.6 and pandas→3.0.3
in `spatialrx`; the existing pipeline still runs (smoke-tested: scanpy scoring, squidpy
Moran) but the environment is no longer byte-identical to the one that produced modules 00–08.
For exact reproducibility of modules 00–08, re-create the env from README pins before that
pipeline; for modules 09–12, keep the current env.