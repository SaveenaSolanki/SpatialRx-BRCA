"""``spatialrx validate-data`` — early cohort/data integrity checks.

Checks (PASS / FAIL / SKIP):
  1. raw tar present (and optionally checksum-verified)
  2. per-section extracted SpaceRanger files present
  3. processed module-00 h5ads load with expected dimensions and metadata
  4. sample table internally consistent (patient/treatment/replicate structure)
  5. optional: module-00 QC totals

Non-strict mode treats *missing data* as SKIP (so the command is usable on a
machine without the cohort).  Strict mode turns missing data into FAIL.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

from . import data as datamod
from .config import PATHS, Paths

REQUIRED_SUFFIXES = [
    "filtered_feature_bc_matrix.h5",
    "matrix.mtx.gz",
    "barcodes.tsv.gz",
    "features.tsv.gz",
    "tissue_positions.csv.gz",
    "scalefactors_json.json.gz",
    "spatial_enrichment.csv.gz",
    "isotype_normalization_factors.csv.gz",
]


def _check(name: str, status: str, detail: str = "") -> dict:
    return {"name": name, "status": status, "detail": detail}


def sha256_file(path: Path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def check_raw_tar(paths: Paths, do_hash: bool = False, strict: bool = False) -> dict:
    tar = paths.raw / datamod.RAW_TAR
    if not tar.exists():
        return _check("raw_tar", "FAIL" if strict else "SKIP", f"missing {tar}")
    size = tar.stat().st_size
    if size != datamod.RAW_TAR_BYTES:
        return _check(
            "raw_tar", "FAIL",
            f"size {size} != expected {datamod.RAW_TAR_BYTES}",
        )
    if do_hash:
        digest = sha256_file(tar)
        if digest != datamod.RAW_TAR_SHA256:
            return _check("raw_tar", "FAIL", f"sha256 {digest} != expected")
        return _check("raw_tar", "PASS", f"sha256 verified ({digest[:12]}…)")
    return _check("raw_tar", "PASS", f"{size} bytes (set --hash to verify sha256)")


def check_extracted(paths: Paths, strict: bool = False) -> list[dict]:
    checks = []
    if not paths.extracted.exists() or not any(paths.extracted.iterdir()):
        status = "FAIL" if strict else "SKIP"
        return [_check("extracted", status, f"missing/empty {paths.extracted}")]
    for _, r in datamod.sample_meta().iterrows():
        gsm, spec = r["gsm"], r["specimen"]
        missing = []
        for suffix in REQUIRED_SUFFIXES:
            if not (paths.extracted / f"{gsm}_{spec}_{suffix}").exists():
                # some sections may nest files; fall back to a glob
                if not list(paths.extracted.glob(f"*{spec}*{suffix}")):
                    missing.append(suffix)
        checks.append(
            _check(
                f"extracted:{spec}",
                "PASS" if not missing else ("FAIL" if strict else "SKIP"),
                "all files present" if not missing else f"missing {missing}",
            )
        )
    return checks


def check_processed(paths: Paths, strict: bool = False) -> list[dict]:
    from .io import read_h5ad

    checks = []
    ingest = paths.ingest_dir()
    for _, r in datamod.sample_meta().iterrows():
        spec, patient, trt = r["specimen"], r["patient"], r["treatment"]
        p = ingest / f"{spec}.h5ad"
        if not p.exists():
            checks.append(
                _check(f"processed:{spec}", "FAIL" if strict else "SKIP", f"missing {p}")
            )
            continue
        try:
            adata = read_h5ad(p)
            # In ingested objects the 35 antibody features live in
            # obsm['protein_counts'] / uns['protein_var']; var holds genes only.
            n_prot = (
                int(adata.obsm["protein_counts"].shape[1])
                if "protein_counts" in adata.obsm
                else 0
            )
            problems = []
            if adata.n_vars != datamod.N_GENES_EXPECTED:
                problems.append(f"n_vars={adata.n_vars}")
            if n_prot != datamod.N_PROTEINS_EXPECTED:
                problems.append(f"n_proteins={n_prot}")
            if "protein_counts" not in adata.obsm or adata.obsm["protein_counts"].shape[0] != adata.n_obs:
                problems.append("protein_counts shape")
            if "spatial" not in adata.obsm or adata.obsm["spatial"].shape != (adata.n_obs, 2):
                problems.append("spatial coords")
            for col, expect in (("sample", spec), ("patient", patient), ("treatment", trt)):
                if col not in adata.obs or not (adata.obs[col] == expect).all():
                    problems.append(f"obs[{col}]")
            checks.append(
                _check(
                    f"processed:{spec}",
                    "PASS" if not problems else "FAIL",
                    f"{adata.n_obs} spots x {adata.n_vars} genes, {n_prot} proteins"
                    if not problems
                    else "; ".join(problems),
                )
            )
        except Exception as e:  # noqa: BLE001
            checks.append(_check(f"processed:{spec}", "FAIL", f"{type(e).__name__}: {e}"))
    return checks


def check_sample_table() -> dict:
    problems = datamod.validate_sample_table(datamod.sample_meta())
    return _check(
        "sample_table",
        "PASS" if not problems else "FAIL",
        "consistent" if not problems else "; ".join(problems),
    )


def check_qc(paths: Paths, strict: bool = False) -> dict:
    qc = paths.processed / "00_sample_qc.csv"
    if not qc.exists():
        return _check("sample_qc", "FAIL" if strict else "SKIP", f"missing {qc}")
    try:
        totals = datamod.cohort_spot_totals(qc)
    except Exception as e:  # noqa: BLE001
        return _check("sample_qc", "FAIL", f"{type(e).__name__}: {e}")
    return _check(
        "sample_qc",
        "PASS",
        f"{totals['n_sections']} sections, {totals['n_patients']} patients, "
        f"{totals['total_spots']} spots, {totals['total_genes']} genes, "
        f"{totals['total_proteins']} proteins",
    )


def run_validate_data(
    paths: Paths | None = None,
    strict: bool = False,
    do_hash: bool = False,
) -> dict:
    paths = paths or PATHS
    checks: list[dict] = [check_sample_table(), check_raw_tar(paths, do_hash, strict)]
    checks += check_extracted(paths, strict)
    checks += check_processed(paths, strict)
    checks.append(check_qc(paths, strict))

    n_pass = sum(c["status"] == "PASS" for c in checks)
    n_fail = sum(c["status"] == "FAIL" for c in checks)
    n_skip = sum(c["status"] == "SKIP" for c in checks)
    return {
        "ok": n_fail == 0,
        "strict": strict,
        "n_checks": len(checks),
        "n_pass": n_pass,
        "n_fail": n_fail,
        "n_skip": n_skip,
        "checks": checks,
    }


def print_report(report: dict) -> None:
    for c in report["checks"]:
        mark = {"PASS": "✓", "FAIL": "✗", "SKIP": "–"}[c["status"]]
        line = f"  [{mark}] {c['name']:34s} {c['detail']}"
        print(line)
    print(
        f"\n  {report['n_pass']} passed, {report['n_fail']} failed, "
        f"{report['n_skip']} skipped (strict={report['strict']})"
    )
