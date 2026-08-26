# SpatialRx-BRCA

Mapping how endocrine therapy remodels the breast-cancer spatial ecosystem using
spatial transcriptomics + spatial proteomics (GSE289326, the CITEgeist cohort).

## The question

1. How does primary endocrine therapy change *where* tumour, immune and stromal states
   occur inside breast cancer tissue?
2. Do RNA and protein measurements tell the same biological story — and where/how much
   do they disagree, spatially?
3. Does a multimodal (RNA + protein + spatial-graph) representation recover biologically
   coherent spatial tumour states better than RNA alone?
   (Deliberately **not** patient-response prediction: n = 6.)

## Data

- **GSE289326** (public 2026-03-19), supplementary tar `GSE289326_RAW.tar` (~1.7 GB)
- 14 GSM entries / 12 unique tissue sections from 6 ER+/HER2- breast cancer patients in
  trial NCT05914792 (primary endocrine therapy, women ≥ 70 y)
- S1 = pre-treatment core biopsy; S2 = post-primary-endocrine-therapy surgical specimen
  (3–50 months of pET); P4/P5 have replicate sections
- Platform: 10x Visium CytAssist with same-slide antibody capture
  (18,085 genes + 35 protein features: 31-marker Immuno-Oncology panel + 4 isotype controls)
- Companion files per section: SpaceRanger outputs (h5/mtx/tsv), H&E images,
  `isotype_normalization_factors.csv` (protein normalisation), and the authors'
  `spatial_enrichment.csv` (their per-feature Moran's I table — used as a validation target)

## Environment

```bash
conda env create -n spatialrx python=3.11   # or: conda create -y -n spatialrx python=3.11 pip
conda activate spatialrx
pip install numpy\<2.2 scipy pandas scikit-learn statsmodels matplotlib seaborn
pip install scanpy squidpy leidenalg python-igraph decoupler harmonypy scikit-misc
pip install torch --index-url https://download.pytorch.org/whl/cpu   # CPU is enough
```

## Pipeline (experiments/)

| # | Script | What it does | Outputs |
|---|--------|--------------|---------|
| 00 | `00_ingest.py` | extract tar, build per-sample AnnData (RNA counts, protein counts, spatial coords, images) | `data/processed/00_ingest/*.h5ad`, `00_sample_qc.csv` |
| 01 | `01_rna_foundations.py` | RNA-only workflow on one section: QC → norm → log → HVG → PCA → UMAP → Leiden | figures `01_*` |
| 02 | `02_spatial_foundations.py` | overlay expression onto the H&E image | `02_*` figures |
| 03 | `03_tumor_heterogeneity.py` | per-section clustering + auditable domain annotation; harmony-corrected merged atlas | `03_domains/*.h5ad`, `atlas_merged.h5ad`, tables |
| 04 | `04_spatial_stats.py` | grid neighbour graphs, Moran's I for all genes, SVG ranking, **validation vs authors' Moran's I** | `04_*` tables/figures |
| 05 | `05_ecosystem.py` | gene-program scores; protein CLR normalisation; transparent NNLS deconvolution (RNA + protein); neighbourhood enrichment | `05_ecosystem/*.h5ad`, tables, figures |
| 06 | `06_multimodal.py` | per-spot RNA vs protein agreement for 27 matched markers; spatial structure of discordance | `06_*` tables/figures |
| 07 | `07_treatment.py` | paired pre/post analyses (compartments, programs, ESR1, tumour–immune contact, neighbourhoods, Moran's I shifts); exact permutation tests | `07_*` tables/figures |
| 08 | `08_ai.py` | RNA-only vs RNA+protein vs RNA+protein+graph embeddings (small MLP/GCN autoencoders); silhouette / spatial-coherence / niche-separation metrics | `08_*` tables/figures |
| 06b | `06b_protein_spatial.py` | Moran's I of the protein layer; validation vs authors' antibody enrichment | `06b_*` |
| 10 | `10_lr_analysis.py` | ligand–receptor scores at tumour–immune interfaces, pre vs post | `10_*` |
| 11 | `11_discordance_cellstate.py` | discordance hotspots vs cell states + edge-artifact control | `11_*` |
| 12 | `12_subsampling.py` | size-matched subsampling of post sections (random + contiguous core) | `12_*` |

Reproduce end-to-end:

```bash
conda run -n spatialrx python experiments/00_ingest.py
conda run -n spatialrx python experiments/01_rna_foundations.py
conda run -n spatialrx python experiments/02_spatial_foundations.py
conda run -n spatialrx python experiments/03_tumor_heterogeneity.py
conda run -n spatialrx python experiments/04_spatial_stats.py
conda run -n spatialrx python experiments/05_ecosystem.py
conda run -n spatialrx python experiments/06_multimodal.py
conda run -n spatialrx python experiments/07_treatment.py
conda run -n spatialrx python experiments/08_ai.py
```

## Outputs

- `outputs/figures/` — all PNG figures referenced in the report
- `outputs/tables/` — all result tables (CSV)
- `outputs/SpatialRx-BRCA_report.md` — canonical report
- `data/processed/` — intermediate AnnData objects
- `logs/` — run logs

## Key findings (see report for full detail & caveats)

1. **Spatial autocorrelation reproduces the study's own analysis**: our Moran's I
   correlates r = 0.70–0.94 (per sample) with the authors' `spatial_enrichment.csv`;
   top spatially structured genes are classic breast-tissue architecture genes
   (IGKC, COL1A1, COL3A1, MUC1, SCGB2A2...).
2. **RNA and protein agree only weakly at spot level** (median Pearson r ≈ 0.02–0.19;
   structural markers like ACTA2/VIM/EPCAM agree best, immune markers ~0).
3. **Disagreement itself is spatially organised**: per-spot RNA–protein discordance has
   Moran's I ≈ 0.24–0.36 and is positive in essentially all samples/markers.
4. **Post-therapy, tumour/stroma programmes become more spatially structured**
   (mean Moran's I of fibroblast activation 0.38 → 0.54, EMT 0.26 → 0.37) while
   T-cell spatial structure decreases (0.12 → 0.07); none of the paired changes
   survive multiple-testing correction at n = 6 (largest trend: B-cell activity ↑,
   p = 0.094 uncorrected).
5. **Multimodal embeddings**: adding protein + spatial graph strongly increases
   embedding spatial coherence (mean Moran's I of dims 0.28 → 0.66; partly by
   construction) and marginally silhouette, but does **not** improve biological
   niche separation (AUROC 0.765 → 0.752) — an honest negative for the
   "more modalities = better niches" hypothesis at this scale.

## Documented deviations & limitations

- decoupler 2.x removed the PROGENy/resource API → curated MSigDB-inspired programs
  used instead of PROGENy (script documents this).
- No cell2location/CITEgeist QP deconvolution (reference-based/heavy GPU dependencies);
  transparent NNLS with RNA + protein panels used instead; the study's own enrichment
  output was used only for validation.
- n = 6 patients; pre samples are core biopsies (328–2,256 spots), post samples are
  surgical sections (4,351–4,992 spots) → only proportions / per-spot-normalised
  metrics are compared pre/post.
- Per-patient responder status is not published per patient in the paper main text;
  not inferred.