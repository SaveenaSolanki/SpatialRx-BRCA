"""
Module 00 — Ingest GSE289326: extract tar, build per-sample AnnData objects.
Creates data/processed/00_ingest/<specimen>.h5ad with:
  - X            : raw RNA counts (Gene Expression)
  - lay['protein']: raw protein (antibody) counts, same spots
  - obsm['spatial']: tissue coordinates (pxl_row/col in fullres image)
  - uns['spatial'] : squidpy-compatible visium schema (hires/lowres images + scalefactors)
  - obs extras   : n_genes, n_counts, n_proteins, protein_umis, treatment, patient, sample
Also saves a sample QC table and dumps format notes about companion files
(spatial_enrichment.csv, isotype_normalization_factors.csv).

Usage: conda run -n spatialrx python experiments/00_ingest.py
"""
from __future__ import annotations

import gzip
import io
import json
import os
import sys
import tarfile

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import PIL.Image as Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from srx import RAW_DIR, EXTRACT_DIR, PROC_DIR, SAMPLE_META

TAR = os.path.join(RAW_DIR, "GSE289326_RAW.tar")


def extract_tar():
    if os.path.exists(EXTRACT_DIR) and os.listdir(EXTRACT_DIR):
        print(f"[ingest] using existing extraction: {EXTRACT_DIR}")
        return
    os.makedirs(EXTRACT_DIR, exist_ok=True)
    print("[ingest] extracting tar (may take a few minutes)...")
    with tarfile.open(TAR, "r") as tf:
        members = tf.getmembers()
        print(f"[ingest] {len(members)} members in tar")
        tf.extractall(EXTRACT_DIR, filter="data")  # py>=3.12 safe filter
    print("[ingest] extraction done")


def find_files(sample_key: str):
    """Return dict of paths for files belonging to one specimen."""
    base = {"sample": None, "h5": None, "positions": None, "scalefactors": None,
            "hires": None, "lowres": None, "tissue": None, "isotype": None,
            "enrichment": None, "features": None, "barcodes": None, "mtx": None}
    for root, _, files in os.walk(EXTRACT_DIR):
        for f in files:
            if not f.startswith(sample_key):
                continue
            if f.endswith("filtered_feature_bc_matrix.h5"):
                base["h5"] = os.path.join(root, f)
            elif f.endswith("tissue_positions.csv.gz"):
                base["positions"] = os.path.join(root, f)
            elif f.endswith("scalefactors_json.json.gz"):
                base["scalefactors"] = os.path.join(root, f)
            elif f.endswith("tissue_hires_image.png.gz"):
                base["hires"] = os.path.join(root, f)
            elif f.endswith("tissue_lowres_image.png.gz"):
                base["lowres"] = os.path.join(root, f)
            elif f.endswith("aligned_tissue_image.jpg.gz"):
                base["tissue"] = os.path.join(root, f)
            elif f.endswith("isotype_normalization_factors.csv.gz"):
                base["isotype"] = os.path.join(root, f)
            elif f.endswith("spatial_enrichment.csv.gz"):
                base["enrichment"] = os.path.join(root, f)
            elif f.endswith("features.tsv.gz"):
                base["features"] = os.path.join(root, f)
            elif f.endswith("barcodes.tsv.gz"):
                base["barcodes"] = os.path.join(root, f)
            elif f.endswith("matrix.mtx.gz"):
                base["mtx"] = os.path.join(root, f)
    return base


def read_gz_df(path, **kwargs):
    with gzip.open(path, "rt") as fh:
        return pd.read_csv(fh, **kwargs)


def load_image(p):
    """Load a .png.gz / .jpg.gz image into an RGB numpy array."""
    with gzip.open(p, "rb") as fh:
        img = Image.open(io.BytesIO(fh.read())).convert("RGB")
    return np.asarray(img)


def build_sample(gsm: str, specimen: str) -> pd.Series | None:
    files = find_files(gsm)
    if not files["h5"]:
        print(f"  !! {gsm}: no h5 found — files: {[k for k,v in files.items() if v]}")
        return None

    adata = sc.read_10x_h5(files["h5"], gex_only=False)
    adata.var["feature_name"] = adata.var_names
    adata.var_names_make_unique()
    ft = adata.var["feature_types"].astype(str)
    # Proteins = ONLY Antibody Capture rows. NOT the "-1" suffix trick: many human
    # gene symbols legitimately contain dashes (NKX2-1, KRTAP5-1, KRT14-1...) or
    # get "-1" disambiguation suffixes for duplicate symbols — those are GEX rows.
    is_ab = ft.str.fullmatch("Antibody Capture|Protein", case=False)
    n_ab = int(is_ab.sum())
    n_gene = int((~is_ab).sum())
    print(f"  {gsm} {specimen}: {adata.shape[0]} spots x {adata.shape[1]} features "
          f"({n_gene} genes, {n_ab} antibodies)")

    # ---- split RNA vs protein ----
    rna = adata[:, ~is_ab].copy()
    prot = adata[:, is_ab].copy()
    if sp.issparse(adata.X):
        prot_X = adata.X[:, is_ab.values].tocsr()
    else:
        prot_X = sp.csr_matrix(adata.X[:, is_ab.values])
    # protein matrix stored in obsm (dims differ from genes); AnnData layers must
    # match (obs, var), so obsm + uns['protein_var'] is the layout
    rna.obsm["protein_counts"] = prot_X
    rna.uns["protein_var"] = pd.DataFrame(
        {"index": prot.var_names.astype(str)},
        index=prot.var_names.astype(str),
    )
    del adata

    # ---- spatial positions ----
    pos_path = files["positions"]
    with gzip.open(pos_path, "rt") as fh:
        first = fh.readline().strip()
    if first.startswith("barcode"):
        pos = pd.read_csv(pos_path, compression="gzip")
    else:  # headerless Visium v1
        pos = pd.read_csv(pos_path, compression="gzip", header=None,
                          names=["barcode", "in_tissue", "array_row", "array_col",
                                 "pxl_row_in_fullres", "pxl_col_in_fullres"])
    pos = pos.set_index("barcode")
    common = rna.obs_names.intersection(pos.index)
    pos = pos.loc[common]
    pos = pos.reindex(rna.obs_names)  # match h5 order
    rna.obsm["spatial"] = pos[["pxl_row_in_fullres", "pxl_col_in_fullres"]].to_numpy(dtype=float)
    rna.obs["array_row"] = pos["array_row"].to_numpy(dtype=float)
    rna.obs["array_col"] = pos["array_col"].to_numpy(dtype=float)
    rna.obs["in_tissue"] = pos["in_tissue"].to_numpy(dtype=int)

    # ---- images + scalefactors (squidpy schema) ----
    scf = json.load(gzip.open(files["scalefactors"], "rt"))
    lib = f"{gsm}_{specimen}"
    uns_spatial = {}
    entry = {"scalefactors": scf, "metadata": {"chemistry_description": scf.get("chemistry_description", "Visium CytAssist")}}
    if files["hires"]:
        entry["images"] = {"hires": load_image(files["hires"])}
    if files["lowres"]:
        entry["images"].setdefault("hires", None)
        entry["images"]["lowres"] = load_image(files["lowres"])
    uns_spatial[lib] = entry
    rna.uns["spatial"] = uns_spatial
    rna.uns["library_id"] = lib
    rna.uns["scalefactors"] = scf
    if files["tissue"]:
        rna.uns["hne_jpg"] = files["tissue"]  # path kept for reference

    # ---- companion files (saved as CSV next to the h5ad; h5ad-uns can't
    # reliably serialize their mixed dtypes) ----
    comp_dir = os.path.join(PROC_DIR, "00_ingest")
    if files["isotype"]:
        iso = read_gz_df(files["isotype"])
        iso.to_csv(os.path.join(comp_dir, f"{specimen}__isotype_factors.csv"), index=False)
    if files["enrichment"]:
        try:
            enr = read_gz_df(files["enrichment"])
            enr.to_csv(os.path.join(comp_dir, f"{specimen}__spatial_enrichment.csv"), index=False)
        except Exception as e:
            print(f"  !! enrichment parse failed: {e}")

    # ---- obs annotations ----
    meta = SAMPLE_META[SAMPLE_META.gsm == gsm].iloc[0]
    rna.obs["sample"] = specimen
    rna.obs["patient"] = meta["patient"]
    rna.obs["treatment"] = meta["treatment"]
    rna.obs["gsm"] = gsm
    rna.obs["n_genes"] = np.asarray((rna.X > 0).sum(axis=1)).ravel()
    rna.obs["n_counts"] = np.asarray(rna.X.sum(axis=1)).ravel()
    rna.obs["protein_umis"] = np.asarray(prot_X.sum(axis=1)).ravel()
    rna.obs["n_proteins"] = np.asarray((prot_X > 0).sum(axis=1)).ravel()

    out = os.path.join(PROC_DIR, "00_ingest", f"{specimen}.h5ad")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    rna.write(out)
    return pd.Series({
        "gsm": gsm, "specimen": specimen, "patient": meta["patient"],
        "treatment": meta["treatment"], "spots": rna.n_obs,
        "genes": rna.n_vars, "proteins": n_ab,
        "median_counts": float(np.median(rna.obs["n_counts"])),
        "median_genes": float(np.median(rna.obs["n_genes"])),
        "median_prot_umis": float(np.median(rna.obs["protein_umis"])),
        "file": out,
    })


def main():
    os.makedirs(PROC_DIR, exist_ok=True)
    extract_tar()

    # inspect tar layout once
    if not os.path.exists(os.path.join(EXTRACT_DIR, "_layout.txt")):
        with tarfile.open(TAR, "r") as tf:
            names = [m.name for m in tf.getmembers()][:12]
            with open(os.path.join(EXTRACT_DIR, "_layout.txt"), "w") as fh:
                fh.write("\n".join(names))
        print("[ingest] tar layout head:", names)

    rows = []
    for _, m in SAMPLE_META.iterrows():
        r = build_sample(m.gsm, m.specimen)
        if r is not None:
            rows.append(r)
    qc = pd.DataFrame(rows)
    qc.to_csv(os.path.join(PROC_DIR, "00_sample_qc.csv"), index=False)
    print("\n[ingest] sample QC table:")
    print(qc.to_string(index=False))


if __name__ == "__main__":
    main()