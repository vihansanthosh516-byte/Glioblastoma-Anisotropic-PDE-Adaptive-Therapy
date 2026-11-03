# Track B LUMIERE go/no-go memo (decision audit, 2026-11-02)

Scope: decision audit only, from `880afcb`. No patient voxels read. No preflight. No contact. No draft adopted. No timing or scoring job. Draft A-2026-11-02-T, the scoring draft and the pre-run freeze stay **unadopted**.

**Recommendation: NO-GO for the full four-analysis run.** Option B below is the only run I would defend, and only after you accept the stated limits in §4.

## 1. Source gates (none relaxed)

### 1a. Reuse terms for results and model weights — OPEN
- Evidence: Figshare items are CC0 (`paper/track_b_no_contact_dataset_plan_2026-10-31.md` §2). The readme last line says "non-commercial use" (same file). Suter et al., *Sci Data* 9:768 (2022), PMC9755255, rechecked 2026-11-02: "This dataset may only be used for non-commercial purposes." The article (not the data) is CC BY 4.0. It states no terms for derived data, results or trained models.
- New evidence 2026-11-02: the NITRC DeepBraTumIA page (https://www.nitrc.org/projects/deepbratumia) now loads. It gives the **tool** license as "Attribution Non-Commercial No Derivatives". No source says whether this license reaches the tool's output masks in LUMIERE. This is not inferred.
- Unknown: may summary results be published; may weights be kept or shared; does the tool's "No Derivatives" term apply to masks derived from its output.
- Endpoint effect: **none on the score.** Pending clarification, this project would share summary numbers only and withhold weights, masks and images. That is a conservative project policy, not a stated license requirement for those outputs.
- Close by: a written statement from the LUMIERE authors or the DeepBraTumIA authors. That needs contact, which is not authorized. Without it, the gate stays OPEN and the conservative rule stays.

### 1b. DeepBraTumIA version — OPEN
- Evidence: not stated on the ARTORG page (checked 2026-11-01), the NITRC page (checked 2026-11-02; the only version string, "v2.1.76-2", is the NITRC site footer), or the *Sci Data* paper (checked 2026-11-02).
- Unknown: which model version made the masks.
- Endpoint effect: **yes, indirect.** The target is the tool's mask. Without a version, the target cannot be reproduced or compared with other DeepBraTumIA outputs. All 84 visits in the release come from one tool run as published, so the comparison between arms inside LUMIERE stays internally fair.
- Close by: a version statement in a published record, or the author contact (not authorized). If neither, report as a stated limit. It cannot become PASS from documents found so far.

### 1c. Label 0 — OPEN
- Evidence: published CSV map 1 Necrosis, 2 Contrast-enhancing, 3 Edema (`paper/track_b_no_contact_dataset_plan_2026-10-31.md` §3). NITRC page and *Sci Data* paper (2026-11-02) name the three compartments and give no numeric labels and no background meaning. No cavity label.
- Unknown: is label 0 a reviewed "no tumor" value, or only "the tool gave no label". Where cavity falls.
- Endpoint effect: **yes, on interpretation.** The reference positive class is labels 1|2 inside support S; label 0 is treated as reference negative. Dice uses true positives, false positives and false negatives, not true negatives. A predicted-positive voxel at label 0 is a false positive against the released annotation, even if the unobserved biology differs. A cavity hidden in label 0 may penalize arms that grow into it. The score remains well-defined as agreement with the automated mask, not true tumor.
- Close by: a published definition (none found) or author contact (not authorized). Otherwise a stated limit: "agreement with automated DeepBraTumIA labels 1|2, label 0 treated as negative, cavity not separated".

## 2. Scientific payoff of the four-analysis run

Facts that cap the payoff:
- **Prior use.** LUMIERE was already used in this project (`src/91_lumiere_volume_forecast.py`, `src/130_h3v_frozen_lumiere.py`). It is not a fresh, unseen test set. Results are exploratory, not confirmation.
- **Automated masks.** The target is DeepBraTumIA output, not expert truth (§1).
- **At most 28 patients** (`output/track_b_lumiere_inventory_20261031/header_check.json`). The primary minimum is 25; QC loss of 4 stops the run.
- **Gates.** Three source gates stay OPEN (§1); label 0 directly affects the endpoint.
- **Cost.** 556.7–699.5 h (S1) or 788.3–931.1 h (S2) wall-clock hours on this PC, a projection, not a bound. The PDE part counts inner cross-fitting but relies on sparse synthetic timings; the central-gap refit is counted at full cohort size (`paper/track_b_lumiere_time_aware_amendment_draft_2026-11-02.md`; timing records in `output/track_b_timing_pilot_20261031/`). These scenarios correspond to about 23–39 days of serial PC time; actual cost is uncertain.

A positive result (every max-T adjusted lower bound > 0) **could support**: "On LUMIERE (exploratory; ≤28 patients; automated DeepBraTumIA masks; label 0 treated as negative), the PDE candidate outperformed the four specified comparators on patient-macro Dice, and this held across the three prior scenarios and the central-gap refit."

It **could not support**: a causal PDE effect; clinical use; agreement with true tumor; generalization to other cohorts or tools; novelty; a confirmatory test (data was used before); sharing of masks or weights.

A negative or incomplete result is likely with n ≤ 28 and four contrasts under max-T. It would show only that this design did not detect a gain.

## 3. Full design vs smaller designs

| Design | Runs | Cost | Claim lost |
|---|---|---|---|
| A. Full (draft §5g) | 3 prior scenarios + central-gap refit | 556.7–931.1+ h | none (within §2 limits) |
| B. Primary only | primary prior run, all four comparators, max-T | not separately projected; must be projected from the timing records before any decision | prior-support robustness; robustness to gap length |
| C. Primary, fewer comparators | drop FiLM U-Net or PDE-only | lower, not projected | "outperformed the four comparators"; it beats only weaker baselines |
| D. No run | none | 0 h | any LUMIERE result; keep the synthetic protocol implementation as a methods result |

C is not defensible: the comparators were set to block weak-baseline claims. B is defensible only if its secondary runs are dropped **before** any score is seen, as a written amendment.

## 4. Decision

**Recommendation: D now; B later only if you want a LUMIERE number.** The full run costs weeks of this PC for an exploratory, reused-data, automated-mask result with an endpoint gate (label 0) that cannot close without contact.

Decisions still needed from you:
1. Pick D (stop the LUMIERE run) or B (primary only).
2. If B: approve a written amendment that drops the three secondary analyses before any patient read, and a runtime projection for B alone.
3. If B: accept §1b and §1c as **stated limits** (they stay OPEN, not PASS), or authorize author contact to try to close them.
4. If B: approve the 84-file preflight (still blocked) and then the pre-run freeze.

Nothing above is adopted by this memo.
