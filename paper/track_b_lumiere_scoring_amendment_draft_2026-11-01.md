# Track B LUMIERE scoring amendment (DRAFT r3 for review; r1 2026-11-01, r3 2026-11-02)

**Status: DRAFT. Not adopted. Not frozen.** It does not authorize the MRI-QC preflight, voxel reads or scoring. The preflight stays **blocked** until the user explicitly approves its exact read scope: visit-2 and visit-3 atlas CT1 voxels and visit-2 atlas brain-mask voxels, 84 files (`paper/track_b_readiness_audit_2026-11-01.md` §5). LUMIERE results are exploratory, not fresh confirmation. Status documents (`START_HERE.md`, roadmap) are not changed until the user adopts this draft.

**Separate protocol.** This amendment governs LUMIERE only. The MU/Duke 90-day protocol (`paper/track_b_pde_external_protocol_draft_2026-10-28.md`; 76–104-day window in `paper/track_b_fixed_growth_baseline_spec_2026-10-31.md` §1.1) stays unchanged and does not apply to LUMIERE. Results under the two protocols are never pooled.

Revision r2 changes: verbatim pass rule (§3a); corrected fold rationale and a concrete fold table (§3b); time-information mismatch (§6); synthetic PDE stress timing over history gaps (§5a); outlier and minimum-count recommendations (§5b, §7).

## 1. Eligible patients by week-label gap (metadata only)

Source: `output/track_b_lumiere_inventory_20261031/header_check.json`, field `gap_weeks` = week(v3) − week(v2) from folder names. No voxel read.

| Week-label gap v2→v3 | 11 | 12 | 13 | 14 | 15 | Total |
|---|---|---|---|---|---|---|
| Patients | 3 | 11 | 8 | 3 | 3 | **28** |

History gap v1→v2 (same file, week labels): range 1 to 155; 25 of 28 are 6–23. Outside 6–23: Patient-086 at 155 weeks, Patient-070 at 1 week, Patient-078 at 4 weeks.

## 2. Interval rule

- **L1.** The LUMIERE interval is the **published week-label difference**, Δw. Window: 11 ≤ Δw ≤ 15, first usable post-op triplet only.
- **L2.** Where a model needs days, it uses **nominal days = 7 × Δw**. This is an **approximation**, not an exact interval. Week labels are binned; the exact interval is unknown.
- **L3.** No conversion to the MU/Duke 90-day or 76–104-day window. "LUMIERE 11–15 weeks" is its own estimand.

## 3. What carries over from the existing Track B design

| Item | Applies to LUMIERE? | Note |
|---|---|---|
| Primary score: patient-macro plain Dice of the full-grid future core, `dice_with_empty_rule` | Yes | Core = DeepBraTumIA labels 1\|2. Agreement with an automated tool, not truth. |
| Required mask baselines | **Proposed (r3):** copy-last, time-aware fixed growth, PDE-only, FiLM time-aware U-Net (`paper/track_b_lumiere_time_aware_amendment_draft_2026-11-02.md`) | Replaces the gap-blind fixed growth and U-Net. Same patients, origins, targets, support, 2 mm grid. See §6. |
| Count model (Poisson, offset log interval) | Yes, count metrics only | Offset log(7 × Δw), approximate (L2). |
| Volume-trend, geometric train-rate, ablations | Yes, secondary | Rates × nominal days (L2). |
| Pass rule | Yes, verbatim (§3a) | The appendix interval clause "one nominal day-90 window with eligible 76–104-day intervals" is **replaced** for LUMIERE by L1–L3. |
| No model or scorer receives future tumor size | Yes | |
| Outer folds | 5 patient-grouped folds (§3b) | |
| Amendment R (U-Net reuse on full fingerprint) | Yes | Adopted 2026-11-01. |
| Support | Origin-only rule (`paper/track_b_no_contact_dataset_plan_2026-10-31.md` §4) | Released origin brain mask + history labels; not validated. |

### 3a. Pass rule, copied verbatim

Source: `paper/research_plan_complete_project_appendix.md`.

§1.3 (comparison and "better" rule):
> For lower-is-better scores, define improvement as baseline score minus hybrid score. For Track C, define improvement as hybrid score minus comparator score. Calculate paired 95% cluster-bootstrap confidence intervals by resampling whole donors, patients, or virtual parameter scenarios, preserving all repeated cells, visits, noise levels, and replicates within each unit. Use a max-T bootstrap adjustment across all required comparator contrasts within a track. Do not use cells, voxels, or repeated noisy runs as independent bootstrap units.
>
> A track supports a within-study advantage only if all eligibility and validation gates pass and every adjusted interval for the full hybrid versus a required ablation or baseline lies strictly on the favorable side of zero.

§3.3 amendment 2026-10-31 (extract):
> Use one nominal day-90 window with eligible 76–104-day intervals and report each actual interval; the estimand is conditional on that interval. No model or scorer receives future tumor size. [...] If the U-Net is unavailable for any patient, the comparison is incomplete, not passed. The primary outcome is patient-macro Dice for the forecast tumor core. Report patient-bootstrap intervals for paired differences.

§3.5:
> If the PDE does not beat every required mask baseline (§3.3 amendment) on plain patient-macro Dice under the frozen margin rule, report no demonstrated added value.

**Reading.** The appendix has **no numeric margin**. The "frozen margin rule" is the §1.3 rule: every max-T-adjusted paired 95% patient-bootstrap interval (Dice is higher-is-better, so improvement = hybrid − baseline) must lie strictly above zero. For LUMIERE the interval clause of §3.3 is replaced by L1–L3; the rest applies unchanged. **Open for review:** confirm this reading.

**Proposed LUMIERE pass-rule wording (r3, not adopted):** the PDE candidate passes only if, on patient-macro plain Dice, every max-T-adjusted paired 95% patient-bootstrap interval for the candidate versus each of the four required comparators — copy-last, time-aware fixed growth, PDE-only and the FiLM time-aware U-Net — lies strictly above zero. The claim is limited to performance against these four specified comparators (`paper/track_b_lumiere_time_aware_amendment_draft_2026-11-02.md` §4).

### 3b. Fold assignment (deterministic, patient-grouped, balanced by Δw)

Correction to r1: r1 said 3-patient gap groups "cannot be split across 5 folds". That was wrong. Three patients can go to three different folds.

**Proposed rule.** Sort the 28 eligible patients by (Δw, patient number). Assign them in turn to folds 0, 1, 2, 3, 4, 0, 1, … . One patient = one fold (all visits of a patient stay together). No randomness, no outcome input. Source of Δw and IDs: `header_check.json`.

| Fold | n | Δw counts (11/12/13/14/15) | History gaps v1→v2 (weeks) | Patients |
|---|---|---|---|---|
| 0 | 6 | 1/2/2/0/1 | 13, 6, 6, 16, 12, 23 | 003, 042, 068, 013, 074, 023 |
| 1 | 6 | 1/2/2/0/1 | 18, 16, 18, 14, **155**, 22 | 025, 043, 075, 032, 086, 055 |
| 2 | 6 | 1/2/1/1/1 | 12, 12, 13, 14, 6, 15 | 029, 054, 077, 040, 067, 066 |
| 3 | 5 | 0/3/1/1/0 | 14, 13, 13, 11, **1** | 006, 058, 091, 051, 070 |
| 4 | 5 | 0/2/2/1/0 | 14, 15, 14, 18, **4** | 034, 062, 007, 064, 078 |

Each Δw group spreads over at least three folds. Training sizes 22, 22, 22, 23, 23. The three history-gap outliers fall in three different folds. If QC excludes patients, the remaining patients **keep their fold**; no re-folding (§7).

## 4. Primary analysis: pooled or by gap

- **Primary: pooled** over Δw = 11–15 (all QC-passing patients).
- **By gap: descriptive only.** Report per-Δw patient-macro Dice for each arm, with n. Groups with n ≤ 3 (Δw 11, 14, 15) are labelled "descriptive, n = 3". No test, no pass or fail by gap.
- Weak point: pooling assumes the task is similar across 11–15 weeks. 28 patients cannot test that.

## 5. How time-dependent arms use the interval

| Arm | Uses time? | Rule |
|---|---|---|
| Copy-last | No | — |
| Time-aware fixed growth (proposed) | Yes | r = k × 7 × Δw(v2,v3); one k per outer fold from training patients; 20 mm cap. See §6. |
| PDE-only and PDE candidate | Yes | Calibrate on v1→v2 with 7 × Δw(v1,v2); forecast v2→v3 with 7 × Δw(v2,v3). Both approximate. |
| Count model | Yes | Offset log(7 × Δw). |
| Geometric train-rate / volume-trend | Yes | Training rate × 7 × Δw. |
| FiLM time-aware U-Net (proposed) | Yes | z = log(Δw/13) for forecast and history intervals, via FiLM. See §6. |

### 5a. PDE runtime stress over history gaps (synthetic only)

Record: `output/track_b_timing_pilot_20261031/lumiere_history_gap_stress.json` (created 2026-11-02T07:02 UTC; exit code 0). Synthetic volumes on the 91 × 109 × 91 grid, one synthetic patient (same masks for every case), only the history gap changed (7 × weeks nominal days); forecast 91 nominal days; same PDE settings as `lumiere_timing_record.json`; same PC. All 10 runs status "ok". Total about 50,650 s (14.1 h).

**Measured (synthetic, one case per gap):**

| History gap (weeks) | Calibration-only (s) | Full forecast (s) |
|---|---|---|
| 1 | 690.8 | 642.1 |
| 4 | 3,076.7 | 2,872.0 |
| 6 | 974.3 | 1,042.4 |
| 23 | 2,521.3 | 3,183.6 |
| 155 | 19,549.0 | 16,098.3 |

Earlier single case (`lumiere_timing_record.json`, 30-day history): 1,544.3 s calibration-only, 1,250.4 s full.

Run time is **not monotone** in the gap (4 weeks is slower than 6). The 155-week case is about 6× the 23-week case.

**Projected (not measured).** Per analysis = 28 full forecasts + 112 calibration-only fits (each patient calibrates in the 4 folds where it trains). Measured times are used for the six patients outside 11–23 weeks (1, 4, 6, 6, 6, 155 weeks: Patient-070, -078, -042, -068, -067, -086). The other 22 patients (11–23 weeks) have no direct measurement, so two bounds:

| Scenario for the 22 mid-gap patients | PDE h per analysis | PDE h, four analyses | Total with Amendment R (U-Net 2 × 223.3 h) | Total, no reuse (4 × 223.3 h) |
|---|---|---|---|---|
| Low: earlier 30-day case | 80.9 | 323.4 | 770.0 | 1,216.6 |
| High: 23-week case | 116.6 | 466.2 | 912.8 | 1,359.4 |
| Old estimate (all patients = 30-day case) | 57.8 | 231.2 | 677.6 | 1,124.1 |

- Patient-086 (155 weeks) alone: about 26.2 PDE h per analysis.
- The 112 calibration-only calls per analysis already count one inner cross-fitted forecast for each patient in the four folds where that patient trains. These are scenario projections, not bounds; patient geometry and overhead remain uncertain.
- U-Net hours are unchanged (223.3 h per analysis, measured step × projected step count). The §6 option B time-aware U-Net would add about 447 h more.

Weak points: one synthetic case per gap; real tumor shapes may change solver cost; mid-gap patients are bracketed, not measured; the cause of the non-monotone timing is not found. Resource needs only, not a deadline.

### 5b. History-gap outliers

**Recommendation:** keep Patient-086 (155 weeks), Patient-070 (1 week) and Patient-078 (4 weeks) in the **primary** cohort. Flag them in every per-patient table. Predeclare one **sensitivity analysis** now: the same comparison without patients whose v1→v2 Δw is outside 6–23 weeks (n = 25 before QC). The bounds 6–23 are set now from metadata only. Never exclude any patient based on outcome, Dice or forecast quality. If the primary and sensitivity results disagree in direction, report both and claim neither.

## 6. Time-information mismatch (revised r3, 2026-11-02)

**The problem.** The PDE arms see each patient's Δw. The frozen gap-blind U-Net and fixed growth do not, so a win over them could come from the time input.

**Proposed resolution (not adopted): replace, not add.** The required comparators become copy-last, time-aware fixed growth, PDE-only and the FiLM time-aware U-Net. The exact recipes, synthetic checks, pass-rule wording, run-time scenarios and gate list are in `paper/track_b_lumiere_time_aware_amendment_draft_2026-11-02.md`. Prototype code: commit `ff4b9f0`, not wired into any runner.

**Superseded (r2):** option B ("add time-aware baselines in addition", six required comparators, a 5-channel broadcast U-Net, and "+447 h on top of 677.6 h"). It is withdrawn. A broadcast constant channel is removed by InstanceNorm, and the "rules out time input alone" wording overclaimed.

**Run time (wall-clock hours on this PC; projections, not bounds):** 555.4–698.2 h or 787.0–929.8 h for the four analyses, depending on the U-Net step basis. The 112 calibration-only calls per analysis include inner cross-fitting. Details in `paper/track_b_lumiere_time_aware_amendment_draft_2026-11-02.md` §5.

**Limit:** 22–23 training pairs per fold is little data for the U-Net to learn a time effect.

## 7. Registration failures and minimum count

- Run the frozen QC (`paper/track_b_no_contact_dataset_plan_2026-10-31.md` §5b) once per patient, **before any target label is read**.
- A failed patient is excluded from **all arms**. Reason reported. No retry.
- Freeze the passed patient list (with hash) before scoring. Report n passed, n failed, and the §1 and §3b tables for passed patients.
- **Recommendation: adopt the minimum of 25 QC-passing patients.** If fewer than 25 pass: **stop**, report the QC result, do not score, do not re-fold, do not relax QC. Below 25, some folds would have 4 or fewer test patients and the bootstrap intervals would be very wide. Weak point: 25 is a judgment, not a power calculation.

## 8. Exact QC calculations (frozen in `src/fk_longitudinal_prep.py:registration_qc`)

Parameters p = (a_x, a_y, a_z, t_x, t_y, t_z): Euler angles in degrees, shifts in mm. R = Rz(a_z)·Ry(a_y)·Rx(a_x), about the voxel-grid centre.

- **Angle:** `max(|a_x|, |a_y|, |a_z|)` degrees. Pass if ≤ 10.
- **Shift:** `sqrt(t_x² + t_y² + t_z²)` mm. Pass if ≤ 10.
- **NCC:** Pearson correlation of fixed and warped intensities inside the origin brain mask. Pass if ≥ 0.80. Non-finite NCC fails.
- Pass only if all three pass.

Limits: the angle is the largest single Euler angle, not the total rotation angle (total can be up to about √3 times larger). The shift is about the grid centre, not the brain centre. Thresholds came from synthetic phantoms. Any change needs a dated amendment before the preflight.

## 9. Gates (unchanged by this draft)

| Gate | Status |
|---|---|
| License (derived outputs) | OPEN |
| Label map (DeepBraTumIA version, label 0, no cavity) | OPEN |
| School written evidence | OPEN |
| Registration (real patients) | BLOCKED until preflight approval |
| Support (real patients) | OPEN |
| Protocol freeze and hashes | OPEN |
| Scoring | BLOCKED |

## 10. Decisions for the user

| # | Choice | Recommendation |
|---|---|---|
| 1 | Pass rule reading (§3a): no numeric margin; max-T adjusted intervals strictly above zero | Confirm |
| 2 | Fold rule (§3b): sorted by Δw, round-robin | Adopt |
| 3 | Time mismatch (§6): replace the gap-blind U-Net and fixed growth with time-aware versions (`paper/track_b_lumiere_time_aware_amendment_draft_2026-11-02.md`), or keep gap-blind baselines and limit the claim | Replace |
| 4 | Outliers (§5b): keep, flag, predeclared 6–23-week sensitivity | Adopt |
| 5 | Minimum 25 QC passes, stop rule, no re-folding (§7) | Adopt |
| 6 | Preflight read scope (84 files) | Your explicit approval needed; still blocked |

No voxels were read, no one was contacted, nothing was downloaded, and no patient was scored for this draft. The §5a timing uses synthetic volumes only.
