# SpatialRx-BRCA — Working plan & task ledger

## Objective
Map how primary endocrine therapy (pET) remodels the breast-cancer spatial ecosystem using
GSE289326 (CITEgeist study): 12 matched pre/post samples from 6 ER+ breast cancer patients,
profiled on 10x Visium CytAssist with same-slide antibody capture (RNA + protein per spot).

Central questions:
1. How do tumour / immune / stromal states change *location* within tissue after therapy?
2. Do RNA and protein measurements tell the same biological story (and where do they disagree)?
3. Does a multimodal (RNA+protein+spatial-graph) representation recover biologically coherent
   spatial tumour states better than RNA alone?  (NOT patient-response prediction — n=6 too small.)

## Data
- GEO: GSE289326 (public 2026-03-19), supplementary `GSE289326_RAW.tar` (~1.7 GB).
- 14 GSM entries → 12 unique tissue sections:
  - P1-S1 (pre), P1-S2 (post); P2-S1, P2-S2; P3-S1_A, P3-S2; P4-S1, P4-S2 (+ P4-S2_1i_rep replicate);
    P5-S1, P5-S2 (+ P5-S2_F_rep replicate); P6-S1, P6-S2_D.
  - S1 = pre-treatment core biopsy; S2 = post-primary-endocrine-therapy surgery (3–50 months pET).
  - Trial NCT05914792; 4 responders, 2 progressors (per imaging + ctDNA).
  - P4-S2 harbours ESR1 D538G (paper case study: MDK signalling).
- Per-sample SpaceRanger-style files: filtered_feature_bc_matrix.h5, matrix.mtx.gz, barcodes,
  features.tsv.gz (RNA genes + antibody features ending "-1"), tissue_positions, scalefactors,
  H&E images (aligned_tissue_image.jpg, cytassist_image.tiff, hires/lowres png),
  isotype_normalization_factors.csv (protein normalisation), spatial_enrichment.csv (CITEgeist??).

## Software
- Conda env `spatialrx` (Python 3.11): scanpy, squidpy, decoupler, leidenalg, python-igraph,
  matplotlib, seaborn, statsmodels, scikit-learn, torch (CPU suffices for module AI).
- Anything GPU (cuOPT/CITEgeist QP) NOT required; we approximate deconvolution with transparent
  NNLS/signature scoring and compare against the study-provided enrichment output.

## Module pipeline (experiments/ scripts, all runnable end-to-end)
1. `00_ingest.py` — extract tar, build per-sample AnnData (RNA + protein layers, spatial,
   images), save canonical merged objects + sample metadata table.
2. `01_rna_foundations.py` — RNA-only workflow on one sample: QC → norm → HVG → PCA → UMAP →
   Leiden (learn "a spot is an observation").
3. `02_spatial_foundations.py` — overlay expression onto tissue coordinates; canonical marker maps.
4. `03_tumor_heterogeneity.py` — cluster all samples; marker DEGs (rank_genes_groups); annotate
   molecular tissue domains.
5. `04_spatial_stats.py` — spatial neighbour graph, Moran's I (squidpy), SVG ranking.
6. `05_ecosystem.py` — gene-program scores (curated programs + decoupler PROGENy), RNA-based and
   protein-based cell-state scoring/NNLS, cellular neighbourhood analysis.
7. `06_multimodal.py` — per-spot RNA vs protein agreement (matched markers), discordance maps,
   spatial structure of discordance.
8. `07_treatment.py` — paired pre/post: domain composition shift, program shift, tumour–immune
   contact, neighbourhood change; paired permutation tests (n=6).
9. `08_ai.py` — RNA-only vs RNA+protein vs RNA+protein+spatial graph embeddings; metrics:
   silhouette vs domains, Moran's I, marker enrichment, treatment-state separation.

## Deliverables
- outputs/figures/*.png
- outputs/SpatialRx-BRCA_report.md (canonical report: summary, methods, results, figures, caveats)
- README.md, CHANGELOG.md, PLAN.md (this file)

## Verification rules
- Every quantitative claim must trace to a script output or figure (logs/ or outputs/).
- Pre/post stats are paired (6 patients); report effect sizes + permutation p-values, not
  overclaimed ML accuracy.
- No fabricated numbers: run everything real.

## Task ledger
- [x] 2026-08-22: env assessment (scanpy broken: numba vs numpy2.5); GEO + paper recon (CITEgeist
      preprint biorxiv 2025.02.15.638331; GSE289326 structure confirmed; antibody suffix “-1”
      caveat found: use feature_type only, not name suffix).
- [x] download GSE289326_RAW.tar (1.7 GB, byte-exact)
- [x] conda env spatialrx ready (py3.11, scanpy 1.11.5, squidpy 1.8.2, torch 2.5.1)
- [x] 00 ingest (corrected 35-antibody panel; companion CSVs; per-sample h5ads + QC table)
- [x] 01–08 + 06b runs + verification (Moran validation vs authors: r≈0.83 RNA, 0.69 protein)
- [x] annotation fix: compartment-level ambiguity (proliferative-tumour clusters no longer collapse to Mixed)
- [x] report + figures + final QA; all artifacts verified on disk

## Known residuals / future
- Paired pET changes: n=6 underpowered; trends only (B-cell ↑, tumour-immune contact ↓, stromal structuring ↑). Contact decrease survived size-matched subsampling (module 12).
- PROGENy/cell2location not run (API/GPU constraints). CITEgeist Step-2 QP now RUN on P4-S2
  (validates our NNLS, r=0.76 epithelial); Step-3 single-cell assignment blocked (no StarDist
  weights/segmentation artifacts) — documented in next-steps report.
- Per-patient responder status not published; do not infer.
- Env drift: citegeist/cuopt install bumped numpy/pandas in spatialrx (modules 00-08 predate).