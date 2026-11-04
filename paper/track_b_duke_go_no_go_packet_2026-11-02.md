# Track B Duke go/no-go packet (2026-11-02)

> **Historical status note (2026-11-03):** The user now reports that Duke granted download access and started a local download. This supersedes only the access-status line below. The exact accepted terms are not saved, the archive is incomplete, and all eligibility, approval, protocol, QC, and scoring gates remain open. Metadata-only inventory is authorized after download completion; no MRI, mask, or landmark-coordinate contents are authorized by that scope. See [the current readiness note](track_b_duke_readiness_2026-11-03.md).

**Historical verdict as of 2026-11-02: NO-GO** for restricted Duke access and for any patient work then.
At that time, Duke access had not been confirmed in writing, school approval was not saved, nothing was adopted, and no dataset was downloaded or opened. The current access update appears above; remaining gates still block patient work.

## 1. Source and approval status

Source: `paper/track_b_duke_source_ledger_2026-11-02.md`.

| Item | Status |
|---|---|
| Duke written access | Not granted. Pending. No written reply on file. Duke declined the first request; a resubmission on 2026-10-29 is user-reported only. |
| Dataset (Zenodo record 20548528) | 61 patients, 503 sessions, 3 to 18 sessions each, 136.3 GB. Files are restricted. |
| License | CC BY-NC-ND 4.0, or a custom Duke license. The CC license allows producing/reproducing adaptations for noncommercial use but does not authorize sharing adaptations. Whether forecasts, masks or models are adaptations, and which license applies, are unknown; the custom text has not been seen. |
| School approval | User says the school said OK on 2026-10-31. No written copy. No date, approver or conditions. Not passed. |
| Prior use of Duke data here | No evidence of any. Patient overlap with MU, UCSD-PTGBM, CFB-GBM, PREDICT-GBM or LUMIERE is unknown. |
| Registration | Atlas alignment versus session-to-session registration is unresolved. |

Weakest point: every "open" item rests on a missing file in the repo. The user may hold emails that are not saved here.

## 2. Eligible triplet count

**Unknown (blocked).**
Why: the metadata files are restricted and not on this PC. 61 patients and 503 sessions are totals, not triplets. A triplet needs three usable sessions (two history, one target). Surgery, treatment, date and label checks are not done. The gate code refuses a hand-typed count (`src/fk_precision_gate.py`). No validated metadata report exists.

## 3. Proposed endpoint and contrast (PROPOSED, not adopted)

Both need a dated amendment signed before any target mask is opened (`paper/track_b_next_phase_plan_2026-11-02.md`).

- **Proposed endpoint:** agreement (Dice) with the source's automated core labels, where core = labels 1 or 3. The cavity (label 4) is kept separate. This is not true tumor burden. The target is the next observed visit, conditional on its interval.
- **Proposed primary contrast:** the PDE candidate versus a two-history matched-volume non-PDE control (`paper/track_b_two_history_control_spec_2026-11-02.md`). The deterministic pre-boundary contrast is proposed as primary or co-primary. Reason: the boundary noise model was calibrated on the PDE arm. **Update after the synthetic screen (`output/track_b_two_history_screen_20261102/screen.json`):** the control's default was not stronger than last-mask signed-distance in growth or shrinkage, and the best grid cell was selected on that same screen and is not adopted. So this contrast stays PROPOSED and is not ready to be primary. See the spec's "Interpretation correction" section.
- Weakest point: the strength of the control is unverified. A development sensitivity grid is required first. A null or negative result is plausible and would be reported unchanged.
- No learned residual exists. No positive result is promised.

## 4. Runtime projection

Source: `output/track_b_feasibility_pilot_20261102/pilot.json` and `paper/track_b_compute_audit_2026-11-02.md`.
These are **projections, not bounds**. They are serial hours on this PC, from synthetic cases. They are not a Duke estimate. The Duke grid, gaps and eligible N are unknown.

| N | PDE calls | FiLM steps | Upper-leaning scenario (pilot, scaled; representative gap) | Upper-leaning, extreme 23-week gap | Lower scenario (recorded times) |
|---|---|---|---|---|---|
| 20 | 380 | 129,600 | 308 h | 546 h | 306 h |
| 30 | 570 | 194,400 | 463 h | 819 h | 459 h |
| 40 | 760 | 259,200 | 617 h | 1,092 h | 612 h |
| 61 | 1,160 | 396,000 | 942 h | 1,667 h | 935 h |

Weak points: the PDE time was timed with reduced settings and scaled 20x. Timings varied about 2x between sessions. A 155-week gap was not timed in the new pilot; the old record shows it costs 5 to 8 times more per call. Whether this budget is acceptable is a human decision.

## 5. Open gates

1. Duke written access (open).
2. Duke license and derived-output terms (open; custom text unseen).
3. School written approval saved in the repo (open).
4. Eligible triplet count from validated metadata (blocked by 1 to 3).
5. Patient overlap with other cohorts (unchecked).
6. Compute budget approved by the user (open).
7. Cohort-specific amendment: endpoint, support, interval rule, folds, control, comparator family, failure rules. Signed and hashed before target masks (not written).
8. Control sensitivity grid on development cases (not run).
9. Sealed run: pre-run seal with git state, package versions, eligibility manifest (`src/fk_track_b_run_guard.py`; not yet used on real data). The current git commit must equal the sealed commit exactly. "Clean" means the sealed code and protocol files have no uncommitted change; the rest of the repository is not checked. This is a reproducibility check, not protection against a determined editor.

## 6. Stop rules

- Stop if access is not granted in writing.
- Stop before sharing any forecast, mask or model until the applicable license and whether the output is an adaptation are resolved. Keep restricted access blocked until Duke grants access in writing.
- Stop if the school approval is not saved in writing.
- Stop if the eligible triplet count is too small. The minimum is not set; it must be written in the amendment first.
- Stop if the seal check fails. Do not edit the seal to pass.
- Stop if the compute estimate is above the approved budget. Narrow scope only by a dated amendment before outcomes, and not by weakening only the strongest baseline.
- Report null or negative results unchanged.

## 7. Decisions needed from the user

1. Save the written Duke reply and the school approval in the repo. Or tell me they do not exist.
2. Which Duke license applies, and do derived forecasts count under it?
3. A compute budget in hours.
4. Is the endpoint (labels 1|3, cavity separate, next observed visit conditional on interval) acceptable to draft?
5. Two-history control: keep it as a descriptive diagnostic only (recommended), or select its settings on development-only training folds under a dated amendment before it can be primary or co-primary.
6. Minimum eligible patients for the Duke amendment.
7. Approve running the control sensitivity grid on development cases (synthetic or already-approved data only).
