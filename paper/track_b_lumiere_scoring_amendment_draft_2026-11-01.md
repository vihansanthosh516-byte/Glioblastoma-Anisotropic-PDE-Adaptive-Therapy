# Track B LUMIERE scoring amendment (DRAFT for review, 2026-11-01)

**Status: DRAFT. Not adopted. Not frozen.** It does not authorize the MRI-QC preflight, voxel reads or scoring. The preflight stays blocked until the user approves its exact read scope (`paper/track_b_readiness_audit_2026-11-01.md` §5). LUMIERE results are exploratory, not fresh confirmation.

**Separate protocol.** This amendment governs LUMIERE only. The MU/Duke 90-day protocol (`paper/track_b_pde_external_protocol_draft_2026-10-28.md`; 76–104-day window in `paper/track_b_fixed_growth_baseline_spec_2026-10-31.md` §1.1) stays unchanged and does not apply to LUMIERE. Results under the two protocols are never pooled.

Sources: `paper/track_b_no_contact_dataset_plan_2026-10-31.md` (§5a horizon, §5b frozen QC), `paper/track_b_fixed_growth_baseline_spec_2026-10-31.md` (comparators, pass rule, recipes), `src/fk_longitudinal_prep.py`.

## 1. Eligible patients by week-label gap (metadata only)

Source: `output/track_b_lumiere_inventory_20261031/header_check.json`, field `gap_weeks` = week(v3) − week(v2) from folder names. No voxel read.

| Week-label gap v2→v3 | 11 | 12 | 13 | 14 | 15 | Total |
|---|---|---|---|---|---|---|
| Patients | 3 | 11 | 8 | 3 | 3 | **28** |

History gap v1→v2 (same file, week labels): range 1 to 155; 25 of 28 are 6–23. Outside 6–23: one patient at 155 weeks, one at 1 week, one at 4 weeks.

## 2. Interval rule (resolves weeks vs the 90-day protocol)

- **Rule L1.** The LUMIERE interval is the **published week-label difference**, Δw. The window is 11 ≤ Δw ≤ 15, first usable post-op triplet only.
- **Rule L2.** Where a model needs time in days, it uses **nominal days = 7 × Δw**. This is an **approximation**. Week labels are binned, so the exact interval is unknown. Results never state exact days.
- **Rule L3.** No conversion to the MU/Duke 90-day or 76–104-day window. "LUMIERE 11–15 weeks" is its own estimand.

## 3. What carries over from the existing Track B design

| Item | Applies to LUMIERE? | Note |
|---|---|---|
| Primary score: patient-macro plain Dice of the full-grid future core, `dice_with_empty_rule` | Yes | Core = DeepBraTumIA labels 1\|2. Agreement with an automated tool, not truth. |
| Required mask baselines: copy-last, fixed growth, PDE-only, full-grid U-Net | Yes, all four | Same patients, origins, targets, support, 2 mm grid. |
| Count model (Poisson, offset log interval) | Yes, count metrics only | Offset uses log(7 × Δw), approximate (L2). |
| Volume-trend, geometric train-rate, ablations | Yes, secondary | Rates × nominal days (L2). |
| Pass rule: beat **every** required mask baseline on patient-macro Dice under the frozen margin and interval rule | Yes | **OPEN:** copy the exact margin and interval wording verbatim from the frozen appendix at freeze time. Not restated here to avoid drift. |
| No model or scorer receives future tumor size | Yes | |
| Outer folds | 5 patient-grouped folds; training sizes 22, 22, 22, 23, 23 | **OPEN:** fold assignment rule (seeded random; stratify by Δw or not). Recommend seeded random, no stratification (3-patient groups cannot be split across 5 folds). |
| Amendment R (U-Net reuse on full fingerprint) | Yes | Adopted 2026-11-01. |
| Support | Origin-only rule (plan §4) | Released origin brain mask + history labels; not validated. |

## 4. Primary analysis: pooled or by gap

- **Primary: pooled** over Δw = 11–15, all 28 patients (minus QC exclusions).
- **By gap: descriptive only.** Report per-Δw patient-macro Dice for each arm, with n. Groups of 3 (Δw 11, 14, 15) are labelled "descriptive, n = 3". No test, no pass or fail by gap.
- Weak point: pooling assumes the forecast task is similar across 11–15 weeks. 28 patients cannot test that.

## 5. How time-dependent arms use the interval

| Arm | Uses time? | Rule |
|---|---|---|
| Copy-last | No | — |
| Fixed growth | No | One radius per outer fold, chosen on training patients across all Δw. Limit: one radius ignores the 11–15-week spread. |
| PDE-only and PDE candidate | Yes | Calibrate on v1→v2 with 7 × Δw(v1,v2) nominal days; forecast v2→v3 with 7 × Δw(v2,v3). Both approximate. |
| Count model | Yes | Offset log(7 × Δw). |
| Geometric train-rate / volume-trend | Yes | Training rate × 7 × Δw. |
| U-Net | No (see §6) | — |

**OPEN:** the v1→v2 history gap of 155 weeks (one patient) and 1 and 4 weeks (one patient each). Recommend: keep them (no post-hoc exclusion) and report their per-patient scores. Alternative: a predeclared history-gap limit, set now, before any read.

## 6. U-Net across variable gaps

- The frozen U-Net has **three input channels**: visit-1 core, visit-2 core, origin domain. It has **no time input**.
- **Recommendation: keep the recipe unchanged.** Δw spans 11–15 (about ±15% of the middle). 28 patients give 22–23 training pairs per fold; a time channel with 5 levels would be weakly learned. Changing the recipe would also break the frozen comparator.
- **Limitation (state with every result):** the U-Net forecasts the same mask for any Δw. The PDE arms get the interval. So the U-Net is gap-blind, and that can lower its Dice. This can favour the candidate. A win over the U-Net must carry this caveat.
- Option for your review (not recommended): add a constant Δw channel. This would be a recipe change and needs its own dated amendment before freeze.
- The history v1→v2 gap also differs per patient, and the U-Net does not see it either.

## 7. Registration failures before scoring

- Run the frozen QC (plan §5b) once per patient, **before any target label is read**.
- A failed patient is excluded from **all arms**. The reason is reported. No retry.
- Freeze the passed patient list (with hash) before scoring. Report n passed, n failed, and the Δw table of §1 for passed patients only.
- If fewer than 5 patients pass per outer fold on average (fewer than 25 total), **OPEN:** recommend stop and report, not re-fold. Needs your decision.

## 8. Exact QC calculations (frozen in `src/fk_longitudinal_prep.py:registration_qc`)

Parameters p = (a_x, a_y, a_z, t_x, t_y, t_z): Euler angles in degrees, shifts in mm. Rotation R = Rz(a_z)·Ry(a_y)·Rx(a_x), about the voxel-grid centre.

- **Angle:** `max(|a_x|, |a_y|, |a_z|)` in degrees. Pass if ≤ 10.
- **Shift:** `sqrt(t_x² + t_y² + t_z²)` in mm. Pass if ≤ 10.
- **NCC:** Pearson correlation of fixed and warped intensities inside the origin brain mask. Pass if ≥ 0.80. Non-finite NCC fails.
- Pass only if all three pass.

Limits: the angle is the largest single Euler angle, not the total rotation angle; the total rotation can be up to about √3 times larger. The shift is about the grid centre, not the brain centre, so a pure rotation about another point shows as shift. Thresholds came from synthetic phantoms. Any change needs a dated amendment before the preflight.

## 9. Open choices for your review

1. Exact margin and interval wording of the pass rule (copy verbatim at freeze).
2. Fold assignment rule (recommend seeded random, no stratification).
3. History-gap outliers: keep and report (recommended) or a predeclared limit.
4. U-Net: unchanged recipe (recommended) or a Δw channel by amendment.
5. Minimum passed count after QC (recommend stop below 25).
6. DeepBraTumIA version and label 0: still open from documents.
7. Preflight: blocked until you approve the exact 84-file read scope.

No voxels were read, no one was contacted, nothing was downloaded, and no patient was scored for this draft.
