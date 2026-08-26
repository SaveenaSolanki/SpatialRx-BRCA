# CHANGELOG — SpatialRx-BRCA

## 2026-08-22
- Started project. Environment audit: base conda python 3.13 has broken scanpy (numba 0.61 vs numpy
  2.5); GPU available (RTX 5000 Ada 32GB, RTX 3090 24GB) for optional runs.
- Recon: GSE289326 = CITEgeist dataset (bioRxiv 2025.02.15.638331, PMID 41854411). 12 sections from
  6 ER+ patients (pre-biopsy S1 / post-pET surgery S2), Visium CytAssist + Antibody Capture.
- Download of GSE289326_RAW.tar (1.7 GB) started; conda env `spatialrx` (py3.11 + scanpy/squidpy/decoupler) being created.
- Key format notes: antibody features end in "-1"; per-sample SpaceRanger-style files;
  `spatial_enrichment.csv` = likely CITEgeist enrichment output (to inspect).
- Plan artifact PLAN.md created.
## 2026-08-22 (continued) — full pipeline executed
- Modules 00–08 + 06b all ran end-to-end on real GSE289326 data; outputs in outputs/figures,
  outputs/tables, data/processed (h5ads).
- Key verified numbers: 45,326 spots ingested; atlas 39,428 spots; per-sample Moran's I
  agreement with authors' table r≈0.70–0.94 (RNA) and r=0.69 (protein); RNA-protein spot
  agreement r≈0.02–0.19 (best ACTA2 0.19; immune ~0); RNA-protein discordance spatially
  structured (Moran 0.22, 99.5% positive); paired pET: no BH-significant changes (n=6),
  trends: B-cell activity ↑, tumour-immune contact ↓ (n=3), stromal structuring ↑;
  multimodal embeddings: spatial coherence 0.28→0.66 but niche AUROC 0.765→0.752 (no gain).
- Data-correction episode documented: antibody features must be selected by feature_type
  (35), NOT name suffix; and annotation ambiguity moved to compartment level.
- Final artifacts: outputs/SpatialRx-BRCA_report.md, README.md, PLAN.md, CHANGELOG.md.

## 2026-08-26 — next steps executed
- Cohort scan: no public dataset replicates GSE289326 design; candidates documented
  (GSE331245 SMART n=89 TNBC; ImmunoADAPT not public; SMMART n=4). Power estimates computed
  (13-26 pairs for Moran shifts; 3-18 for contact change).
- CITEgeist: installed citegeist + cuopt-cu12 (NVIDIA index); ran Step-2 spot-level QP on
  P4-S2 with study's own solver (Barrier ~1.7 s/solve) -> proportions CSV
  (data/processed/09_citegeist/). Validation vs our NNLS: epithelial/tumour Pearson r=0.76
  (4,658 spots). Step-3 single-cell assignment BLOCKED (no StarDist weights/segmentation
  artifacts in package/GEO).
- Module 11 (discordance x cell state): discordance NOT edge-driven (rho<0 for 7/8 markers;
  PCNA +0.16 at edges); structured by niche (epithelial/myeloid-rich, fibroblast-poor);
  M2/exhaustion/hypoxia scores don't drive it. Fixed pandas-3 method-shadowing bug
  (df.sample) and edge metric (BFS -> 6-degree exposure).
- Module 10 (LR at interfaces): MIF-CD74, MDK-LRP1 top; CXCL12-CXCR4 down 5.5->2.8 post;
  CSF1-CSF1R +0.33, TGFB1-TGFBR1 +0.28 (n=4 pairs, p=0.000 exact).
- Module 12 (matched-region subsampling): tumour-immune contact decrease survives
  size-matching (random subsample median -0.50, 0/200 positive; contiguous -0.43 ~80%
  negative); other effects compatible with noise.
- Env note: citegeist/cuopt install bumped numpy 2.4.6 + pandas 3.0.3 in spatialrx;
  modules 00-08 predate this (README pins for exact reproduction).
- Deliverables: outputs/SpatialRx-BRCA_next_steps.md; report §7 updated; README table
  updated.
