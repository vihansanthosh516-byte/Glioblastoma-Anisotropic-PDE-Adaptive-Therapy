# Track B no-contact dataset plan (2026-10-31, draft, NOT adopted)

**User rule (2026-10-31):** contact no one (no authors, Duke, TCIA helpdesk or school). Use only public data and published documents. MU may be dropped.

No patient data was downloaded, read or scored for this plan. The no-scoring rule, the full-grid U-Net requirement and unadopted amendment A-2026-10-31-R are unchanged.

## 1. Search result: no fresh no-request cohort passes

Earlier screen (`paper/track_b_external_cohort_screen_2026-10-28.md`): UTSW-Glioma, UCSF-ALPTDG, SAILOR, MOTUM, UPENN-GBM, BrainMD, CPTAC-GBM, UCSD-BMETS fail at least one of these: three masked visits, GBM, no data-use agreement.

New checks today:

| Candidate | Finding (source) | Result |
|---|---|---|
| OpenNeuro ds007045 (Sci Data 2026, PMC13498554) | 337 IDH-wildtype GBM; **one preoperative MRI per patient**. Labels 0–3 defined (0 = background). Atlas registration checked by a neuroradiologist. | Fails: not longitudinal. |
| Duke (Zenodo 20548528) | Files restricted; needs a request. | Fails no-contact rule. |
| CFB-GBM v2.0 (MELBA 2026:038) | Longitudinal RANO dataset; license and access not confirmed today; already a separate project side study. | Not fresh; not checked. |
| UCSD-PTGBM held-out subset | Planned release in early 2026; already a project side study. | Not fresh; not checked. |
| LUMIERE (Sci Data 2022, PMC9755255; figshare 10.6084/m9.figshare.c.5904905.v1) | 91 GBM patients, 638 study dates (about 7 per patient). Automated masks only (DeepBraTumIA, HD-GLIO-AUTO). Skull-stripped with HD-BET. Released images are **unregistered**; no visit-to-visit registration. Dates are in weeks from pre-op. Non-commercial use; no data-use agreement mentioned. | Public, no request. Already used in this project, so exploratory only. |

**Conclusion:** no untouched public GBM cohort with three or more masked visits and no access request was found. Any patient result from now on is **exploratory, not confirmation**.

## 2. What the no-contact rule means for the gates

The open gates needed written answers from data owners. With no contact, those answers cannot come. So the gates can only close if the protocol is amended to accept **published documentation plus our own history-only checks**. That is weaker evidence. It must be labelled "assumed from documentation" in every result.

## 3. Proposed plan (needs the user's approval before it applies)

**Step 1. Choose one primary dataset: LUMIERE (recommended).**
- Why: public, no request, many visits per patient (more triplets than MU's 18 in the development set), and its documentation is explicit that images are unregistered. That means we must do registration ourselves, and we can document it fully.
- Weakest point: automated masks, not expert masks; week-level dates; already used in the project.

**Step 2. Define documentation-based gate rules (dated amendment, before any target mask is opened).**

| Gate | No-contact rule | Label in results |
|---|---|---|
| Label meaning | Use the label values in the segmentation tool's published documentation (DeepBraTumIA / HD-GLIO-AUTO). Read label values from file headers and documentation only, never from target voxels. | "Label meaning from tool documentation." |
| Valid support | Use a brain mask made by HD-BET (the tool the source names) on the **history visits only**, at the forecast origin. Never use the target visit's labels. | "Support assumed: origin HD-BET mask." |
| Registration | Register each target-visit MRI to the origin-visit MRI with a fixed rigid method, using MRI only (no masks). Score registration with an MRI-only metric. Freeze the pass threshold before any target mask is opened. Exclude failed cases and report the count. | "Registration done by us; MRI-only QC." |
| Dates | Use week counts. Widen the eligible window to whole weeks (for example 11–15 weeks for a nominal 13-week forecast). Freeze it before scoring. | "Dates at week resolution." |
| School approval | Not applicable while the project is personal research (`START_HERE.md` scope, 2026-10-20). If the project goes to a school or fair, this gate reopens. | — |
| Access/licence | LUMIERE non-commercial terms. No sharing of data or derived images. | — |

**Step 3. Code work with no patient data (can start now).**
- A LUMIERE loader that reads only metadata and history-visit MRI until the gates are frozen.
- An MRI-only registration step and its QC metric, tested on synthetic volumes.
- An HD-BET-style support rule, tested on synthetic volumes.
- A reuse wiring for the U-Net behind the fingerprint check (stays off until amendment R is adopted).

**Step 4. Freeze, then score once.** Freeze the cohort list, folds, gates, window and code hash. Then open target masks and score all required baselines (copy-last, fixed growth, PDE-only, full-grid U-Net). Report the result unchanged, including a negative result.

**Compute (projection, synthetic timing only):** on the current frozen settings, about 399.5 CPU hours with U-Net reuse for four analyses on 18 patients (`paper/track_b_amendment_decision_2026-10-31.md`). LUMIERE has more patients, so the cost grows roughly with the number of training pairs. Recompute after the cohort is frozen.

## 4. Alternative

Keep MU (already coded) and apply the same documentation-based rules to it. Faster to run, but only 18 development patients, and MU was already scored in this project.

## 5. Not decided here

- Adoption of this plan (needs the user's approval).
- Adoption of amendment A-2026-10-31-R.
- The PDE runtime decision.
