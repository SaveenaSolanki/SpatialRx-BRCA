# 05 — Claim-to-Evidence Reconciliation (Stage A)

**Branch:** `fix/repro-correctness-v1`
**Stage A baseline:** `7c6d4ca` (Gates 1–3 complete)
**Scope:** documentation-only correction of six claims. No result tables were recomputed, and
no computational code was modified in this stage.

**Rule applied:** every corrected statement is traced to an original data file, analysis output,
or statistical computation already in the repository. No unsupported claim was replaced with an
assumption; where an analysis was never performed (C15), the claim was removed.

---

## 1. Reconciliation table

| Claim | Original claim (source) | Corrected claim | Evidence source | Verification status |
|---|---|---|---|---|
| **C10** | "27 matched gene/antibody pairs" (report §2 Module 6; README pipeline table; `06_multimodal.py` docstring) | **26** analysable matched pairs (27 configured; `HLA-DRA` unmatched) | `outputs/tables/06_rna_protein_agreement.csv` → 26 distinct `gene` values; `06_multimodal.py` matches `n.split("-")[0]`, and the protein feature is `HLA-DRA` (no `-1`), so the lookup key becomes `HLA` | `CORRECTED_AND_VERIFIED` |
| **C12** | "34–35 of 35 proteins spatially structured (FDR<0.05) in every section" (report §2 Module 6b) | **30–35 of 35** per section (median 35/35; 9/14 sections at 35/35; lowest P4-S1 = 30/35) | `outputs/tables/06b_protein_morans.csv`, `pval_norm_fdr_bh < 0.05` grouped by sample | `CORRECTED_AND_VERIFIED` |
| **C18** | "B-cell activity +0.086 (p_perm 0.094, ES 1.24)" (report §2 Module 7); README "p = 0.094 uncorrected" | `p_perm = 0.09375` (= 6/64) after the Gate-2 fix; the prose value 0.094 is now **backed by the table** (previously the table stored 0.1875) | `outputs/tables/07_paired_tests.csv`, row `score_B_cell_activity` (Δ 0.085606, p_perm 0.09375, BH 0.900, ES 1.2388) | `CONFIRMED_CORRECT` (post-Gate-2) |
| **C21** | "22 metrics" (report §2 Module 7) | **24** metrics | `outputs/tables/07_paired_tests.csv` → 24 rows | `CORRECTED_AND_VERIFIED` |
| **C31** | "both p=0.000 in exact sign-flip tests, n=4" (next_steps §3; CHANGELOG) | corrected exact two-sided `p_perm = 0.125` (BH-adjusted 0.81 across all 26 LR pairs) | `outputs/tables/10_lr_paired_tests.csv`, rows `CSF1–CSF1R` and `TGFB1–TGFBR1` (Δ +0.3326 / +0.2840, p_perm 0.125, p_BH 0.8125, n=4) | `CORRECTED_AND_VERIFIED` |
| **C15** | "isotype-factor alternative available" (report §2 Module 5); "robustness check … isotype-normalisation factors" (`06_multimodal.py` docstring); "isotype-factor alternative" (`05_ecosystem.py` docstring) | **Only per-spot CLR was applied.** The isotype-factor branch is not implemented (dead code; a per-spot factor cancels under CLR). No isotype-normalised robustness analysis is claimed | `06_multimodal.py`: `clr_by_factor` is never called; `05_ecosystem.py`: only `clr_protein` (CLR) is computed; no isotype-based result exists in `outputs/` | `UNSUPPORTED_CLAIM_REMOVED` |

---

## 2. Traceability detail

- **C10** — the count is a property of the committed table, not an inference:
  `06_rna_protein_agreement.csv` contains exactly 26 unique `gene` values. The configured
  list has 27 entries; the missing pair is `HLA-DRA`. The matching logic builds
  `pidx = {n.split("-")[0]: i}` from protein names and looks up `PAIRS[gene]`; for the protein
  feature `"HLA-DRA"` this yields the key `"HLA"`, so the lookup `"HLA-DRA"` fails. This is the
  known defect **D10** recorded in `docs/02_IMPLEMENTATION_STATUS.md`. The code is unchanged.
- **C12** — per-section significant counts from `06b_protein_morans.csv`:
  35, 35, 35, 35, 34, 35, 30, 35, 33, 31, 35, 35, 34, 35 (14 sections). Range 30–35;
  9 sections at 35/35. This is the known defect **D9**.
- **C18 / C21** — corrected by the Gate-2 fix (`docs/04_REMEDIATION_GATE_REPORT.md` §4.1).
  All effect sizes are byte-identical to the pre-remediation table; only `p_perm` (and its BH
  step-up) changed.
- **C31** — corrected by the same Gate-2 fix. The committed table's `p=0.0` is impossible for
  an exact two-sided sign-flip test (the observed sign assignment is always in the null set);
  the committed script yields 0.25; the corrected value is 0.125. This provenance defect is
  recorded in `docs/04_REMEDIATION_GATE_REPORT.md` §4.3.
- **C15** — no isotype-normalised result exists anywhere in `outputs/`. The `clr_by_factor`
  function is defined but never invoked, and as written it multiplies each spot by a single
  factor before CLR, which cancels. The claim was therefore removed rather than "verified".

---

## 3. Consistency checks executed

| Check | Result |
|---|---|
| C10: `06_rna_protein_agreement.csv` distinct markers | 26 |
| C12: min / max / sections-at-35 in `06b_protein_morans.csv` | 30 / 35 / 9 of 14 |
| C21: rows in `07_paired_tests.csv` | 24 |
| C18: `score_B_cell_activity` p_perm | 0.09375 |
| C31: CSF1–CSF1R and TGFB1–TGFBR1 p_perm / p_BH | 0.125 / 0.8125 |
| C31: minimum `p_BH` across all 26 LR pairs | 0.8125 |
| Corrected strings absent from report (old values) | `27 matched`, `34–35 of 35`, `22 metrics`, `p=0.000`, `isotype-factor alternative available` all removed |
| Full test suite (`pytest tests/`) | **51 passed** (verifies no computational code changed) |

Evidence log: `logs/remediation/stageA_verification.txt`.

---

## 4. What was preserved

- `outputs/tables/*` — **no result table was recomputed or edited in Stage A.**
- All computational code — the only `.py` edits are **docstrings** (verified in the commit
  diff); no executable statement changed.
- Every original claim is still recorded in `docs/00_REPOSITORY_AUDIT.md`,
  `docs/01_SCIENTIFIC_CLAIM_LEDGER.md`, `docs/02_IMPLEMENTATION_STATUS.md`,
  `docs/03_REPRODUCIBILITY_RISKS.md` and `docs/04_REMEDIATION_GATE_REPORT.md`.

## 5. Files changed in Stage A (documentation only)

| File | Change |
|---|---|
| `outputs/SpatialRx-BRCA_report.md` | C10, C12, C18, C21, C15 |
| `outputs/SpatialRx-BRCA_next_steps.md` | C31 |
| `CHANGELOG.md` | C31 inline correction + dated reconciliation entry |
| `README.md` | C10 marker count |
| `experiments/06_multimodal.py` | docstring C10 + C15 (no logic change) |
| `experiments/05_ecosystem.py` | docstring C15 (no logic change) |
| `docs/05_CLAIM_RECONCILIATION_STAGE_A.md` | this table |

---

## 6. Residual claims (explicitly deferred)

These were **not** changed in Stage A because they require Gate 4 / later work:

- **C25** (Module 08 niche AUROC) — validity depends on the Gate-4 benchmark rebuild; the report
  still quotes `0.765 → 0.752` as a transductive, non-held-out result. It must be re-stated in
  Stage B.
- **CITEgeist step not scripted** (`data/processed/09_citegeist/`) — no generating script.
- **External validation** — not performed.
