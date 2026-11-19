# SpatialRx-BRCA — Tracked scientific-correctness issues

Status vocabulary: `OPEN` · `MITIGATED` · `RESOLVED`.

Issues are never silently closed: a claim must either be corrected in the report
(`docs/05_CLAIM_RECONCILIATION_STAGE_A.md`) **and** fixed in code, or remain listed here.

---

## ISSUE-001 — `HLA-DRA` silently dropped by the protein-feature parser (**RESOLVED** 2026-08-26)

- **Discovered:** Phase-0 audit (defect D10); confirmed in Stage A.
- **Location:** `experiments/06_multimodal.py` (`sample_agreement`), also affects any code
  using `n.split("-")[0]` to map antibody features to gene symbols.
- **Description:** the matched-marker map is built as
  `pidx = {n.split("-")[0]: i for i, n in enumerate(pnames)}`. The antibody panel contains the
  feature `HLA-DRA` **without** a `-1` suffix; `"HLA-DRA".split("-")[0]` yields `"HLA"`, so the
  lookup for `PAIRS["HLA-DRA"]` (`"HLA-DRA"`) fails and the pair is skipped.
- **Impact:** Module 06 analyses **26 of 27** configured markers; HLA-DRA RNA–protein
  concordance is absent from `06_rna_protein_agreement.csv`. The report previously stated 27.
- **Evidence:** `outputs/tables/06_rna_protein_agreement.csv` has 26 distinct `gene` values;
  `data/processed/05_ecosystem/*_ecosystem.h5ad` `uns["protein_var"]` contains `HLA-DRA`.
- **Status:** `RESOLVED` (2026-08-26) — parser replaced by `spatialrx/proteins.py`
  (`protein_to_symbol` / `antibody_index`) at all 7 call sites; 11 regression tests in
  `tests/test_proteins.py`; Modules 05/06/06b/07/10/11/12 re-run.
- **Residual (data, not code):** `HLA-DRA` is still absent from the Module-06 pair table because
  this dataset's GEX panel contains **no MHC-region genes at all** (0 `HLA*` features), so no RNA
  counterpart exists. 26 is therefore the maximum number of analysable pairs. The fix does change
  Module 05: `prot_myeloid` mean +0.031 (HLA-DRA now enters the myeloid protein signature).
- **Evidence:** `docs/07_STAGE_B_REPORT_REVISION.md` §2–§3; `logs/rerun_*.log`.

---

## ISSUE-002 — isotype-normalisation "robustness check" is dead code / no-op (**RESOLVED** 2026-08-26)

- **Location:** `experiments/06_multimodal.py` `clr_by_factor` (never called);
  `experiments/05_ecosystem.py` docstring.
- **Description:** the branch multiplies every protein in a spot by a single per-spot factor
  before CLR; that factor cancels under CLR, so it cannot change the result. It was never called.
- **Status:** `RESOLVED` (2026-08-26) — the unsupported claim was removed in Stage A (C15) and
  the unused function is now deleted from `experiments/06_multimodal.py`.

---

## ISSUE-003 — `cluster_zscore_table` produces NaN rows for marker-free spots (**RESOLVED** 2026-08-26)

- **Location:** `experiments/03_tumor_heterogeneity.py`.
- **Description:** `X = log1p(X / X.sum(axis=1) * 1e4)` divides by a per-spot marker sum that can
  be zero for spots with no lineage-marker counts; those rows become NaN and are silently dropped
  by `groupby().mean()`. A `RuntimeWarning` is emitted during atlas reconstruction.
- **Status:** `RESOLVED` (2026-08-26) — the division is guarded, the previous NaN semantics are
  preserved (spots excluded from cluster means), and the affected spot count is reported
  explicitly instead of a silent warning. Reported atlas compartments are unchanged.

---

## ISSUE-004 — original Module 08 metric was leaky (RESOLVED in Gate 4)

- **Location:** `experiments/08_ai.py` `niche_separation` (fit and evaluated on the same spots).
- **Status:** `RESOLVED` in Gate 4 — replaced by a leakage-free, patient-held-out benchmark under
  `benchmark/`; the original table is preserved as historical evidence.
- See `benchmark/statistical_analysis.md`.

---

## ISSUE-005 — CITEgeist Step-2 result has no generating script (**RESOLVED** 2026-08-26)

- **Location:** was `data/processed/09_citegeist/` (CSV + log only, ad-hoc run).
- **Status:** `RESOLVED` — `experiments/09_citegeist_step2.py` stages the SpaceRanger layout from
  the GEO files, runs the QP solver, writes the proportions CSV, a CITEgeist-vs-NNLS comparison
  (`outputs/tables/09_citegeist_vs_nnls.csv`) and `09_run_record.json`; exits 2 with guidance when
  cuOPT/GPU is unavailable. Re-run reproduces epithelial r = 0.761.
- **Note:** CITEgeist Step 3 (cell-level assignment) remains blocked — no StarDist weights ship
  with the package and GEO has no nuclei artifacts; see report §7b.

---

## ISSUE-006 — no external validation (**OPEN — data-limited**)

- **Status:** `OPEN` and **not closable with this dataset**. No public cohort shares the design
  (matched pre/post primary endocrine therapy, ER+, Visium + same-slide antibody capture);
  the nearest powered resource is GSE331245 (SMART: 89 TNBC patients, pre/post chemotherapy,
  GeoMx) — different subtype, therapy and assay.
- **What would close it:** a second matched pre/post spatial cohort (≥13–26 pairs; power table in
  `outputs/SpatialRx-BRCA_next_steps.md` §1).
- **See:** main report §7b.
