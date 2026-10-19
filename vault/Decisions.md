# Decisions Log

> HISTORICAL: September 2026 Track B/C simulation decisions, superseded as project goals by `../START_HERE.md`. Therapy effects are assumed model parameters, not measured patient drug responses. Historical proposals below are not authorization to rerun spent or frozen analyses.

## 2026-09-23 — Adaptive therapy high-rho stratification

**Post hoc simulation finding (corrected 2026-10-04, script 85, `output/rho_threshold_check.json`):** 10 of 61 patients have rho > 0.02/day, and **all 10** progress earlier under adaptive. None of the 51 patients with rho <= 0.02 does (highest rho among them: 0.0191; lowest rho among the 10: 0.0277). The earlier "23 of 61 ... 10/23" matched a cutoff of about rho 0.0088, not 0.02.
The standalone 10/23 statement is superseded by the script-85 stratification above.

**Hypothesis:** The 80% controller setpoint is too tight for fast-growing
tumors — the tumor rebounds from 80% to 100% in days, forcing near-continuous
dosing, so adaptive degenerates to MTD.

**Test:** Rerun with THRESHOLD_OFF=0.50. If high-rho patients now show real
holidays and preserved sensitive clone, artifact. If not, the simulation does not support this threshold-artifact hypothesis; it cannot establish real biology.

**Historical result report:** The author reported that lowering the threshold did not restore holidays and retained 0.80. The exact 0.50 comparison is unverified: no threshold-specific output/key is cited. At 0.80 the saved simulation reports mean drug use 32.36% of MTD and 48.87 holidays (`output/adaptive_cohort_summary.json`, keys `mean_drug_percent_of_mtd` and `mean_num_holidays`). This simulation uses assumed kill rates and does not establish a biological mechanism.

**Decision criteria (original):**
- Artifact → rerun primary with 0.50, update committed result
- Superseded biological interpretation: report only the post hoc simulated stratification (10/10 above 0.02 progress earlier; 0/51 at or below 0.02). This does not establish that adaptive therapy works clinically.

## 2026-09-22 — E_MAX_RATIO = 1000

**Problem:** At 200x, MTD shrinks 0-7% (too weak). At 2000x, saturates
(too strong).

**Historical tuning report (unverified):** Sweep 150x, 500x, 1000x, 2000x on 3 patients. The following shrinkage ranges and the preceding 200x description have no retained run artifact/key cited; they are author-reported model-tuning observations, not measured responses.
- 150x → 0-7% shrinkage
- 500x → 25-58%
- 1000x → 40-62% (model-tuning choice; no clinical comparator source)
- 2000x → 57-62% (saturated)

**Decision:** E_MAX_RATIO = 1000. Same value scales to all patients
via per-patient rho.

## 2026-09-22 — Per-patient E_max from rho

**Historical problem report (unverified spread; parameter-table calculation missing):** Global E_max produced a reported spread in effective drug:growth ratio
across patients — same bug class as the rho scalar.

**Fix:** E_max_patient = patient_rho × E_MAX_RATIO. Every patient gets the
same nominal drug:growth ratio by construction. This is assumed model scaling (`src/44_adaptive_therapy.py`, `E_MAX_RATIO`), not measured patient drug response.
