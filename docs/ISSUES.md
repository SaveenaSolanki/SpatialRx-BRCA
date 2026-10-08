# SpatialRx-BRCA — Tracked scientific-correctness issues

Status vocabulary: `OPEN` · `MITIGATED` · `RESOLVED`.

Issues are never silently closed: a claim must either be corrected in the report
(`docs/05_CLAIM_RECONCILIATION_STAGE_A.md`) **and** fixed in code, or remain listed here.

---

## ISSUE-001 — `HLA-DRA` silently dropped by the protein-feature parser (OPEN)

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
- **Status:** `OPEN` (claim corrected in Stage A; the parser itself is **not** fixed — that is a
  code change requiring re-running Module 06, which is outside the Gate 1–4 documentation scope).
- **Planned fix:** replace the `split("-")[0]` heuristic with an explicit protein→gene alias map
  (strip a trailing `-<int>` only), add a regression test, then re-run Module 06 and re-audit C10.
- **Do not treat as resolved.**

---

## ISSUE-002 — isotype-normalisation "robustness check" is dead code / no-op (MITIGATED)

- **Location:** `experiments/06_multimodal.py` `clr_by_factor` (never called);
  `experiments/05_ecosystem.py` docstring.
- **Description:** the branch multiplies every protein in a spot by a single per-spot factor
  before CLR; that factor cancels under CLR, so it cannot change the result. It was never called.
- **Status:** `MITIGATED` — the unsupported claim was removed in Stage A (C15). The function
  remains for provenance.
- **Planned:** delete the function in the post-Gate-4 code cleanup.

---

## ISSUE-003 — `cluster_zscore_table` produces NaN rows for marker-free spots (OPEN, low)

- **Location:** `experiments/03_tumor_heterogeneity.py`.
- **Description:** `X = log1p(X / X.sum(axis=1) * 1e4)` divides by a per-spot marker sum that can
  be zero for spots with no lineage-marker counts; those rows become NaN and are silently dropped
  by `groupby().mean()`. A `RuntimeWarning` is emitted during atlas reconstruction.
- **Status:** `OPEN` (does not change the reported atlas compartments; flagged for transparency).
- **Planned:** guard the division and count affected spots explicitly.

---

## ISSUE-004 — original Module 08 metric was leaky (RESOLVED in Gate 4)

- **Location:** `experiments/08_ai.py` `niche_separation` (fit and evaluated on the same spots).
- **Status:** `RESOLVED` in Gate 4 — replaced by a leakage-free, patient-held-out benchmark under
  `benchmark/`; the original table is preserved as historical evidence.
- See `benchmark/statistical_analysis.md`.

---

## ISSUE-005 — CITEgeist Step-2 result has no generating script (OPEN)

- **Location:** `data/processed/09_citegeist/` (CSV + log only, ad-hoc run).
- **Status:** `OPEN`.

---

## ISSUE-006 — no external validation (OPEN)

- **Status:** `OPEN` — documented as unavailable; compatibility table not yet produced.
