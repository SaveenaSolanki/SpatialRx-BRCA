# 07 — Stage B report revision + issue closures

**Branch:** `fix/repro-correctness-v1` · **Date:** 2026-08-26
**Scope:** code fixes for the tracked defects, plus the report revision that Stage A deferred
(claim C25), plus the engineering gaps (CI, scripted CITEgeist run).

---

## 1. What was fixed in code

| Issue | Change | Evidence |
|---|---|---|
| **ISSUE-001** (HLA-DRA parser) | New tested helper `spatialrx/proteins.py` (`protein_to_symbol`, `antibody_index`) replacing the `name.split("-")[0]` heuristic at **7 call sites** (`05`, `06`, `06b`, `11`). Handles panel suffixes (`CD68-1`→`CD68`), isoforms (`PTPRC-2`→`PTPRC`, first isoform kept), dashed symbols (`HLA-DRA` preserved) and isotype exclusion. Ambiguous gene symbols such as `NKX2-1` are resolved by passing the gene-symbol universe. | `tests/test_proteins.py` (11 tests) |
| **ISSUE-002** (isotype no-op) | Deleted the unused `clr_by_factor` function (a per-spot factor cancels under CLR). Claim was already withdrawn in Stage A. | `experiments/06_multimodal.py` |
| **ISSUE-003** (NaN z-rows) | `cluster_zscore_table` now guards the per-spot division, **preserves the previous NaN semantics** (spots excluded from cluster means) but reports the affected spot count instead of a silent warning. | `experiments/03_tumor_heterogeneity.py` |
| **ISSUE-005** (CITEgeist unscripted) | New `experiments/09_citegeist_step2.py`: stages the SpaceRanger layout from GEO files, runs the QP solver with the curated panel profile, writes proportions + a CITEgeist-vs-NNLS comparison + a JSON run record, and exits 2 with a clear message if cuOPT/GPU is unavailable. | `logs/rerun_09_citegeist.log`, `data/processed/09_citegeist/09_run_record.json` |
| Stale constant | `srx.PROTEIN_LINEAGE` now describes the **actual** panel (31 antibody features = 30 distinct targets incl. two PTPRC isoforms) instead of a non-existent 44-marker list; isotypes moved to `PROTEIN_ISOTYPES`. | `experiments/srx.py` |
| R1 residual (env drift) | `experiments/srx.py` now opts in to `anndata.settings.allow_write_nullable_strings`, so pandas-3 environments can write h5ad from **every** pipeline script (previously only via `spatialrx.io`). | chain re-run below |
| R2 residual (no CI) | `.github/workflows/tests.yml`: Python 3.11, pinned core stack, CPU torch, `pytest` + coverage, package/CLI smoke checks. | workflow file |

## 2. Re-run of the affected chain

`05 → 06 → 06b → 07 → 10 → 11 → 12` was re-executed with the fixes (logs `logs/rerun_*.log`).

Post-fix effects, measured against the pre-fix tables (`/tmp/pre_hla/`):

| Module | Effect of the fix |
|---|---|
| 05 | **`prot_myeloid` mean +0.031** (HLA-DRA now enters the myeloid protein signature); other protein lineages shifted ≤0.01; RNA outputs unchanged (max shift 0.0). |
| 06 | 26 analysable pairs (unchanged in count — see §3); paired correlations unchanged to <1e-6. |
| 06b / 07 / 10 / 11 / 12 | Numerically identical (max shift 6.7e-07; all p-values byte-identical). |

## 3. The HLA-DRA finding (two independent causes)

Fixing the parser did **not** add an HLA-DRA pair to Module 06, because a second, independent
limitation exists: **this dataset's GEX panel contains no MHC-region genes at all**
(0 `HLA*` features among the 18,085 genes, checked in three sections).
So `HLA-DRA` has a protein feature but no RNA counterpart, and 26 is the maximum number of
analysable pairs. The parser fix still matters: it (a) prevents silent truncation to `HLA`
in any future dataset that *does* have the RNA probe, and (b) repaired module 05's protein
signature, where HLA-DRA *is* used (visible as the `prot_myeloid` shift).

The report's Module 6 text now states both reasons explicitly.

## 4. Report revision (Stage A claim C25)

`outputs/SpatialRx-BRCA_report.md`:

- **One-sentence summary** no longer cites the transductive `0.765 → 0.752` comparison; it
  cites the Gate-4 held-out result (RNA PCA 0.7035 best).
- **Module 8** section rewritten: design (leave-one-patient-out / leave-one-section-out),
  the full held-out table, the three statistical statements (all learned methods below PCA,
  p=0.031; protein −0.0018, p=0.25; graph −0.043, p=0.031), the mechanical-coherence caveat
  (randomised-graph Moran ≈ 0.001) and the label-shuffle controls. The legacy script is
  labelled historical evidence and explicitly "must not be cited as a performance comparison".
- **New §7b** "Why the remaining gaps are *data*-limited, not code-limited" — per-gap
  explanation of what measurement would close it.
- **New §9** revision history mapping each gate to its commit and report.

## 5. Issue register after this stage

| Issue | Status |
|---|---|
| ISSUE-001 HLA-DRA parser | **RESOLVED** (parser fixed + tests; residual data limitation documented) |
| ISSUE-002 isotype no-op | **RESOLVED** (function deleted) |
| ISSUE-003 NaN z-rows | **RESOLVED** (guarded, counted, results unchanged) |
| ISSUE-004 leaky Module-08 metric | RESOLVED (Gate 4, unchanged) |
| ISSUE-005 CITEgeist unscripted | **RESOLVED** (script + run record; reproduces r=0.761) |
| ISSUE-006 no external validation | **OPEN — data-limited** (no public cohort with this design; see report §7b) |

## 6. Verification performed

```bash
python -m pytest tests/ -q                      # 74 passed (11 new protein tests)
python experiments/05_ecosystem.py ... 12_subsampling.py   # full affected chain, exit 0
python experiments/09_citegeist_step2.py --specimen HCC22-088-P4-S2
  # -> epithelial r = 0.761 (reproduces the earlier ad-hoc run)
```

Known upstream quirk (not ours): CITEgeist's `validate_cell_proportions` emits a malformed
logging format string (`'  - %s: %s%'` with two arguments), which surfaces as a non-fatal
logging error. Results are unaffected; recorded here for anyone re-running the script.
