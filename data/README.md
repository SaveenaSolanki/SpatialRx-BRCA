# GSE289326 data — fetch instructions

The raw + processed data are **not committed to this repository** (GitHub hard limit
100 MB/file; the tar is 1.7 GB). Everything needed to reproduce the pipeline:

## 1. Download (one time, ~1.7 GB)

```bash
mkdir -p data/raw
cd data/raw
wget -c "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE289nnn/GSE289326/suppl/GSE289326_RAW.tar"
```

**Checksum (SHA-256):**
`0292ba6990606e52fb8153abfd35ad662c1c266dbe7641a6c51f206b95e75843  GSE289326_RAW.tar`

Verify:
```bash
sha256sum -c <<< "0292ba6990606e52fb8153abfd35ad662c1c266dbe7641a6c51f206b95e75843  GSE289326_RAW.tar"
```

## 2. Ingest

```bash
conda run -n spatialrx python experiments/00_ingest.py
```

This extracts the tar into `data/raw/extracted/` and builds per-sample AnnData objects in
`data/processed/00_ingest/` (RNA counts, protein counts, spatial coordinates, images).
Subsequent modules consume those objects; see README for the full pipeline.

## Contents of the tar (per section, 14 GSM entries)

- `*_filtered_feature_bc_matrix.h5`, `*_matrix.mtx.gz`, `*_barcodes.tsv.gz`,
  `*_features.tsv.gz` — 18,085 genes + 35 antibody-capture features per section
- `*_tissue_positions.csv.gz`, `*_scalefactors_json.json.gz` — spatial layout
- `*_tissue_hires_image.png.gz`, `*_tissue_lowres_image.png.gz`, `*_aligned_tissue_image.jpg.gz`,
  `*_cytassist_image.tiff.gz` — histology
- `*_isotype_normalization_factors.csv.gz` — per-spot protein normalisation factors
- `*_spatial_enrichment.csv.gz` — authors' per-feature Moran's I table (validation target)

## Provenance

- GEO accession: GSE289326 — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289326
- Study: CITEgeist (Chang et al.), bioRxiv 2025.02.15.638331, PMID 41854411
- 14 GSM entries / 12 unique tissue sections from 6 ER+ breast-cancer patients
  (NCT05914792): S1 = pre-treatment biopsy, S2 = post primary endocrine therapy.