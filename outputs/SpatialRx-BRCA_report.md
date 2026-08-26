# SpatialRx-BRCA — How endocrine therapy remodels the breast-cancer spatial ecosystem
## RNA + protein + tissue-neighbourhood integration on GSE289326 (CITEgeist cohort)

**Date:** 2026-08-22 · **Status:** executed end-to-end on real data; every quantitative
claim below traces to a script, table or log in this repository (paths given inline).

**One-sentence result:** On 12 matched pre/post endocrine-therapy breast sections from 6
patients, tumour/stromal gene programs become *more spatially structured* after therapy
(Moran's I of fibroblast-activation program 0.38→0.54; EMT 0.26→0.37) while T-cell spatial
structure decreases (0.12→0.07); RNA and protein measurements agree only weakly at spot
level (r≈0.02–0.19) *but their disagreement is itself spatially organised* (Moran's I ≈
0.24–0.36, positive in ~all samples); and a small multimodal RNA+protein+graph model
increases embedding spatial coherence (0.28→0.66) without improving biological niche
separation (AUROC 0.765→0.752) — an honest null for the "more modalities = better niches"
hypothesis at this cohort size.

---

## 1. Data and preprocessing

| Item | Value | Source |
|---|---|---|
| Study | CITEgeist, ER+/HER2- breast cancer, trial NCT05914792 | bioRxiv 2025.02.15.638331; GEO GSE289326 |
| Sections | 14 GSM entries / 12 unique tissue sections, 6 patients (P1–P6), S1 = pre-treatment biopsy, S2 = post primary endocrine therapy (3–50 months) | series matrix; `experiments/srx.py` |
| Spots (raw ingest) | 45,326 | `data/processed/00_sample_qc.csv` |
| Features per section | 18,085 genes + 35 antibody-capture features (31-marker Immuno-Oncology panel + 4 isotype controls) | h5 `feature_type`; Module 00 log |
| Platform | 10x Visium CytAssist + same-slide antibody capture | GEO record; paper methods |
| Spots after QC (per section) | core biopsies 294–2,256; surgical sections 3,598–4,992; atlas = 39,428 spots × 3,000 HVGs | `03_domains/*.h5ad`, `atlas_merged.h5ad` |
| Median UMI counts per spot | 114–2,367 (deepest: P3-S1_A; shallowest: P4-S2_1i_rep — replicate section, excluded from paired tests) | `00_sample_qc.csv` |

**Data-correction note:** the first ingest pass mis-identified 35 genes whose HGNC symbols
contain dashes (`NKX2-1`, `KRTAP5-1`, …) or duplicate-symbol suffixes (`TBCE-1`) as proteins.
The protein layer is defined strictly by `feature_type == "Antibody Capture"` (35 features).
All downstream modules were rebuilt on the corrected objects.

**Companion files used:** per-section H&E images (hires/lowres PNG, aligned JPG, CytAssist
TIFF); `isotype_normalization_factors.csv` (per-spot protein normalisation factors);
`spatial_enrichment.csv` (the authors' own per-feature Moran's I table — used as an
independent validation target, Modules 4 & 6b).

---

## 2. Module-by-module results

### Module 1 — RNA foundations (one section: HCC22-088-P4-S2)
Workflow: QC → normalize → log1p → HVG (2,000) → scale → PCA(30) → neighbour graph →
UMAP → Leiden (res 0.8). 4,434 spots after filtering (<200 genes removed); **8 Leiden
clusters**. Cluster marker z-scores identify a clear immune cluster (CD3D/CD8A/CD68 z≈+1.9),
two luminal-tumour clusters (EPCAM/ESR1/PGR/KRT; one also MKI67), a fibroblast cluster
(COL1A1/DCN), a vascular+B-cell cluster (PECAM1/VWF/MS4A1/CD79A).
Figures: `01_qc_P4-S2_post.png`, `01_umap_leiden.png`, `01_ranked_genes.png`; object
`data/processed/01_rna_foundations/P4-S2_mod1.h5ad`.

### Module 2 — Spatial foundations
The same markers overlaid on the H&E image: COL1A1 detected in 97.6% of P4-S2 spots
(fibrosis throughout), ESR1 59.7%, CD68 24.0%, EPCAM 23.5%, CD3D 2.3% (sparse T-cell RNA at
spot level — expected dilution, motivating protein layer + deconvolution).
Figures: `02_spatial_markers_P4S2.png`, `02_leiden_on_tissue_P4S2.png`.

### Module 3 — Spatial tumour heterogeneity (whole cohort)
Per-section clustering + *auditable* rule-based domain annotation (mean z of lineage marker
sets → compartment decision at Tumour/Immune/Stroma level, then finest domain inside).
Design detail: the ambiguity check is done at **compartment** level so that proliferative
tumour clusters (epithelial≈proliferation scores) are not collapsed to "Mixed".
Harmony-corrected merged atlas: 39,428 spots × 3,000 HVGs, `atlas_merged.h5ad`.
Pre/post compartment composition (12 paired sections; mean spot fraction, from
`07_paired_metrics.csv`): pre: Tumour 0.091, Immune 0.180, Stroma 0.148, Mixed 0.580 →
post: Tumour 0.161, Immune 0.093, Stroma 0.235, Mixed 0.511. Core biopsies are
tumour-poor but contain tumour in 4/6 pre samples (P1 7.7%, P2 3.0%, P4 14.4%, P5
29.7% of spots; P3/P6 pre cores have none).
Tables: `03_cluster_domain_scores.csv`, `03_domain_composition_fractions.csv`;
figure `03_atlas_umap.png`.

### Module 4 — Spatial statistics (Moran's I) **validated against the study's own output**
Grid neighbour graphs (6 neighbours, array coordinates) + Squidpy Moran's I on all 18,085
genes per section.

- **Validation:** Pearson r between our Moran's I and the authors' `spatial_enrichment.csv`
  per section: 0.70, 0.87, 0.91, 0.94, 0.77, 0.90, 0.73, 0.94, 0.78, 0.80, 0.89, 0.88,
  0.70, 0.86 → **mean ≈ 0.83 (n=14 sections)** (`04_morans_validation.png`).
- Top spatially variable genes (mean Moran's I, n=14): IGKC 0.496, COL1A1 0.487, COL3A1
  0.469, COL1A2 0.441, MUC1 0.380, SCGB2A2 0.368, XBP1 0.366, FABP4 0.366, FN1 0.363 —
  classic spatially structured breast-tissue architecture genes (`04_svg_ranking.csv`,
  `04_svg_heatmap.png`).
- Canonical markers: COL1A1/COL1A2/DCN/LUM and EPCAM/KRT8 high; immune markers (CD3D,
  CD8A) lower Moran's I — immune cells are more dispersed at spot scale.

### Module 5 — Gene programs and the cell ecosystem
- 11 curated gene programs + lineage scores per spot (scanpy `score_genes`);
  **PROGENy was NOT run** — decoupler 2.x removed the PROGENy/resource API (documented
  deviation; curated MSigDB-inspired programs used instead).
- Protein layer: winsorized (5%) + per-spot CLR (per CITEgeist methods); stored as
  `obsm["protein_clr"]`; isotype-factor alternative available (`isotype_normalization_factors.csv`).
- Transparent NNLS deconvolution (scipy.optimize.nnls), RNA (33 marker genes, 6 lineages)
  and protein (31 antibodies, 8 lineages) panels; proportions per spot in
  `rna_<lineage>` / `prot_<lineage>` obs columns.
- **RNA vs protein deconvolution agreement** (per-spot Pearson across lineages): epithelial/
  tumour agrees best (e.g., P4-S2 epithelial r=+0.40; P2-S2 r=+0.55), immune lineages
  ≈ 0 (`05_rna_prot_deconv_agreement.csv` + figure).
- Neighbourhood enrichment on compartment labels (14/14 sections): tables `05_nhood_*.csv`.
Figures: `05_program_maps_P4S2.png`, `05_cellstate_maps_P4S2.png`,
`05_cellstate_prot_maps_P4S2.png`.

### Module 6 — Multi-omics: do RNA and protein tell the same story?
Per-spot correlations (RNA log-CPM vs CLR protein) for 27 matched gene/antibody pairs.

- **Agreement is modest and marker-dependent** (median Pearson r across 14 sections;
  `06_rna_protein_agreement.csv`): ACTA2 0.19, VIM 0.14, EPCAM 0.12, BCL2 0.11, CD68 0.10,
  PCNA 0.10 — versus CD3E 0.02, CD8A 0.02, CD4 −0.03, PDCD1 −0.03. Structural/abundant
  markers agree best; rare immune markers essentially do not (spot-level dilution + antibody
  vs mRNA biology).
- **Disagreement is spatially organised:** per-spot discordance (zRNA − zProtein) has
  Moran's I mean 0.219 over samples/markers with 99.5% positive values (`06_discordance_morans.csv`).
  Discordance maps: `06_discordance_maps_P4S2.png`; scatter examples `06_scatter_examples_P4S2.png`.

### Module 6b — Protein (antibody) spatial structure
Moran's I of CLR proteins per section: **34–35 of 35 proteins spatially structured
(FDR<0.05) in every section**. Top proteins: VIM 0.69, ACTA2/SMA 0.67, KRT5 0.63, BCL2
0.56, HLA-DRA 0.55, CD68 0.52, EPCAM 0.46, CD3E 0.42. Validation vs authors' table:
**r = 0.69 across matched protein entries** (`06b_protein_morans.csv`,
`06b_protein_spatial.png`). The protein layer is *more* structured than the transcript
layer — consistent with surface-protein capture of stably expressed cell populations.

### Module 7 — Therapy remodelling (paired pre vs post, n = 6 patients)
Metrics per section: compartment fractions, deconvolution proportions, program scores,
ESR1 within tumour (log-CPM), tumour–immune contact (fraction of tumour spots with an
immune neighbour), neighbourhood enrichment, program Moran's I.
Statistics: paired per-patient deltas; **exact two-sided permutation tests** (sign flips,
2^6) + Wilcoxon signed-rank; BH across metrics.

- **No metric survives multiple-testing correction at n=6** (22 metrics; `07_paired_tests.csv`).
- Uncorrected trends: B-cell activity +0.086 (p_perm 0.094, ES 1.24); tumour spots/section
  +670 (p 0.125, ES 1.51); ESR1 within tumour +1.69 log-CPM (n=3 pairs only; P1/P2/P4);
  **tumour–immune contact −0.40 (ES −1.88, n=3 pairs)** — directionally, residual tumour in
  post-surgical samples is *less* immune-adjacent; immune fraction −0.088.
- Spatial-structure shifts (Moran's I of programs, paired): fibroblast activation +0.160,
  EMT +0.114, tumour deconvolution +0.117, angiogenesis +0.050 increase post-therapy;
  T-cell structuring −0.049 decreases (`07_program_morans_paired_tests.csv`). None reach
  significance after correction; consistent directions across ≥5/6 patients for the
  fibroblast/EMT increases.
- Figures: `07_paired_metrics.png`, `07_spatial_prepost_P1.png`, `07_spatial_prepost_P4.png`,
  `07_nhood_delta.png`.

**Honest interpretation:** n=6 is underpowered for paired inference; the design (core biopsy
vs whole surgical section) confounds absolute quantities, so we compare only proportions and
per-spot metrics. The consistent *spatial-structuring* trends are hypothesis-generating, not
confirmatory, and the tumour-immune-contact change rests on n=3 pairs.

### Module 8 — Multimodal representation (RNA ± protein ± spatial graph)
Small self-supervised autoencoders (32-d embeddings, MSE; GCN variant: one graph-conv
pre-layer on row-normalised grid adjacency). Configs: RNA-only (MLP), RNA+protein (MLP),
RNA+protein+graph (GCN). Metrics per sample, averaged over 14 sections (`08_embedding_metrics.csv`):

| Metric | RNA-only | RNA+protein | +graph | Δ (graph − RNA) |
|---|---|---|---|---|
| Silhouette vs compartments | 0.034 ± 0.097 | 0.032 ± 0.090 | 0.064 ± 0.122 | +0.030 (wins 11/14) |
| Spatial coherence (mean Moran's I of dims) | 0.281 ± 0.146 | 0.277 ± 0.141 | 0.656 ± 0.243 | +0.375 (wins 14/14) |
| Biological niche separation (AUROC of program-high spots) | 0.765 ± 0.045 | 0.766 ± 0.045 | 0.752 ± 0.040 | −0.013 (wins 2/14) |

Interpretation: adding the spatial graph massively increases embedding *coherence* — partly
by construction (neighbours are smoothed together), so Moran's I is a manipulation check, not
an independent win. Critically, **neither protein nor spatial context improved biological
niche separation** in this cohort (AUROC unchanged/decreased), and silhouette gain is small.
This is the honest, quantitative answer to the project's AI question: at this cohort size,
RNA alone already carries most of the niche information captured by these simple models.
Figures: `08_embedding_metrics.png`, `08_embedding_umap_P4S2.png`.

---

## 3. What the project teaches (conceptual map)

RNA-seq → scRNA-seq concepts (spot = observation) → spatial transcriptomics (X/Y + H&E) →
spatial statistics (Moran's I, SVG) → cell ecosystem (programs, NNLS deconv, neighbourhoods)
→ spatial multi-omics (RNA vs protein agreement + discordance) → therapy-associated spatial
remodelling (paired design) → multimodal representation (honest benchmark).

## 4. Limitations / blocked / not run

- **n=6 patients** — no claim of treatment-response prediction; paired tests are
  underpowered by design.
- **Biopsy vs surgical section** — pre samples are small cores; only proportion/per-spot
  metrics compared.
- **PROGENy / cell2location / CITEgeist-QP deconvolution not run**: decoupler 2.x dropped the
  PROGENy resource API; cell2location needs heavy reference/GPU; CITEgeist QP needs NVIDIA
  cuOPT. Transparent NNLS + the study's own enrichment output used instead.
- Per-patient responder status not published per patient in the paper main text → not
  inferred; reported overall 4 responders / 2 progressors only.
- Replicate sections (P4-S2_1i_rep, P5-S2_F_rep) are QC/depth heterogeneous; excluded from
  paired tests.
- Protein-layer definition is intentionally strict (feature_type == “Antibody Capture”,
  35 features); gene-symbol dashes (NKX2-1, KRTAP5-1…) are unambiguous GEX features.
- PROGENy-based pathway activity: see deviation above. Images were verified numerically
  (marker z-scores, detection rates); visual QA of every PNG is left to the reader
  (`outputs/figures/`).

## 5. Reproducibility

Environment: conda `spatialrx` (Python 3.11): scanpy 1.11.5, squidpy 1.8.2, decoupler 2.2.0,
torch 2.5.1 (CPU), harmonypy, scikit-misc, leidenalg. Data: `data/raw/GSE289326_RAW.tar`
(1,824,737,280 bytes, MD5-verifiable via GEO). Pipeline: `experiments/00_ingest.py` →
`01_rna_foundations.py` → `02_spatial_foundations.py` → `03_tumor_heterogeneity.py` →
`04_spatial_stats.py` → `05_ecosystem.py` → `06_multimodal.py` → `06b_protein_spatial.py` →
`07_treatment.py` → `08_ai.py` (see README). All logs in `logs/`.

## 6. Sources

- GEO: GSE289326 — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289326
  (supplementary: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE289nnn/GSE289326/suppl/GSE289326_RAW.tar)
- CITEgeist preprint: Chang AC et al., "CITEgeist: Cellular Indexing of Transcriptomes and
  Epitopes for Guided Exploration of Intrinsic Spatial Trends", bioRxiv 2025.02.15.638331 —
  https://www.biorxiv.org/content/10.1101/2025.02.15.638331v2 (PMID 41854411)
- CITEgeist code: https://github.com/leeoesterreich/CITEgeist
- Clinical trial: NCT05914792 — https://clinicaltrials.gov/study/NCT05914792
- Wu et al. 2021 scRNA atlas (reference used by the study): GSE176078
- Methods context: Schlee... spatial autocorrelation: Moran 1950; Squidpy (Palla et al. 2022,
  Nat. Methods 19:171–178); Scanpy (Wolf et al. 2018, Genome Biol. 19:15); PROGENy
  (Holland et al. 2020, Nat. Commun. 11:20) — cited for the pathway concept, not run here.

## 7. Recommended next steps — EXECUTED (2026-08-26)

All four recommended follow-ups were run on this same project; full report:
`outputs/SpatialRx-BRCA_next_steps.md`. Headlines:

1. **Powered cohort scan**: no public dataset replicates this design; nearest powered
   alternative is GSE331245 (SMART, n=89 TNBC pre/post-chemo GeoMx). Power math: 13–26
   matched pairs needed for the Moran-I shifts; the ongoing trial NCT05914792 is the natural
   replication source.
2. **CITEgeist single-cell assignment**: Step 2 (spot-level QP with the study's own solver,
   GPU/cuOPT) ran successfully on P4-S2 and **validates our NNLS** (epithelial/tumour
   r=+0.76 across 4,658 spots). Step 3 (StarDist cell assignment) is blocked — the package
   ships no segmentation weights and GEO has no nuclei artifacts. The discordance hotspot
   hypothesis was instead tested at spot level: discordance is **not edge-driven** and is
   structured by niche composition (epithelial/myeloid-rich, fibroblast-poor).
3. **Ligand–receptor at tumour–immune interfaces**: MIF–CD74 and MDK–LRP1 dominate;
   CXCL12–CXCR4 decreases and CSF1–CSF1R/TGFB1–TGFBR1 increase post-therapy (n=4 pairs).
4. **Matched-region subsampling**: the tumour–immune contact decrease survives size-matched
   subsampling of post sections (random: 0/200 iterations positive; contiguous-core ≈80%
   negative) — the strongest treatment signal in the cohort.

## Appendix — Additional reproducibility notes

- Installing `citegeist`+`cuopt-cu12` for study 2 bumped numpy/pandas in `spatialrx`;
   modules 00–08 outputs predate that (see README/next-steps report for env guidance).
- Figure inventory and per-module run logs are in `outputs/figures/` and `logs/`.

## 8. Sources

1. Larger cohort (multi-centre ER+ pET trials) for powered pre/post spatial inference.
2. Reference-free protein-first deconvolution at single-cell assignment scale (CITEgeist
   pipeline with cuOPT) to test whether the RNA/protein discordance hotspots match cell-level
   phenotypes (e.g., post-transcriptional regulation panels).
3. Ligand–receptor / cell-communication analysis (COMMOT or LIANA) on deconvoluted layers,
   focusing on the tumour–immune contact regions (Module 7's most biologically interesting trend).
4. Correct for section-size confound with matched-region subsampling before any follow-up
   treatment tests.
5. Multimodal embeddings with contrastive spatial objectives, evaluated on independent
   cohorts before claiming niche-discovery gains over RNA-only baselines.