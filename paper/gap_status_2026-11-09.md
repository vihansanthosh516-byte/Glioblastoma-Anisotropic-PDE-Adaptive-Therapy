# Open-gap status record (2026-11-09)

Checked against the working tree at commit `efc4ea6`. This record changes no frozen rule, no result file, and no configuration value. No patient data, image, mask, or policy outcome was read or generated for it. No code changed, so no test run was needed for this record.

Status words: **fixed** (closed by a documentation change made with this record), **already recorded** (the files already said it; nothing to change), **blocked** (needs a user decision, a user-supplied document, or work that must not start yet), **irreducible** (a limit of the data or design that no honest fix removes).

## Whole project

| Gap | Status | Evidence and action |
|---|---|---|
| Tracks A and B have no matched patients; the link to Track C is assumed synthetic signals | Irreducible | `tracks/README.md` shared rules; `START_HERE.md` 2026-11-09 prior-art bullet. Every summary must keep the words "simulation" and "assumed". |
| The framework cannot show measured patient drug resistance or real treatment benefit | Irreducible | Same sources. No claim of either is made. |
| No single list says which plan is active | Fixed | Plan-status table added to `tracks/README.md`. One file, `paper/trackABC_integrated_research_plan_revised_2026-11-08.md`, is linked from no hub and is marked "status undecided" until the user decides. |
| Evaluation rules must be fixed before scoring | Already recorded | Track C gate in `tracks/track-c/README.md`; A5 frozen record. No rule was changed here. |

## Track A

| Gap | Status | Evidence and action |
|---|---|---|
| A5 failed its frozen 5-of-6 rule (won 4 of 6 donors) | Irreducible; failure preserved | `output/trackA_yu2020/result.json`; `START_HERE.md` Yu 2020 bullet. Candidate mean AURC 0.1479 vs logistic 0.1764, wins in GS1, GS2, GS5, GS6, losses in GS3, GS13. |
| Candidate did not beat its no-gate version; Brier worse | Irreducible | Same bullet: no-gate AURC 0.1473 vs 0.1479; Brier 0.2335 vs 0.1944; GS13 logistic 0.0109 vs 0.0726. |
| 119-of-121-gene score is a proxy, not a drug-resistance test | Irreducible | `paper/trackA_tmz_response_119of121_amendment_2026-11-08.md`. Paired differences split evenly across the six donors. |
| Labels mark sampled regions, not verified malignant cells | Irreducible | A5 record: mixed-cell, author-supplied MRI-marked point labels. |
| No fresh independent cohort passes all gates | Blocked | `START_HERE.md` Track A status bullet: strongest near-miss GSE286419 has four pairs, below six. Needs a new dataset; scouting is paused by the user. |
| Calibration, within-library, scANVI, attention analyses deferred | Blocked | Same bullet: "not run; deferred pending a future dated amendment". Needs a dated amendment from the user before any run. |

## Track B (Duke)

| Gap | Status | Evidence and action |
|---|---|---|
| Confirm files hold labels, visits, patient matches | Already recorded (stale gap) | Metadata and header-only inventory, 2026-11-08: record lists 61 patients, 503 sessions. A private run on 2026-11-04 read target masks. |
| Eligible patient count and visit-pair list | Already recorded, but spent | A frozen cohort exists in private storage. It was scored on 2026-11-04, so it is exploratory and cannot be fresh confirmation. |
| Already-scored cohort must be labelled exploratory | Already recorded | `tracks/track-b/README.md` and the module header of `src/track_b_transition_estimand.py` both say so. |
| 2-mm prototype versus the scored 10-mm analysis | Already recorded; limitation | The 2-mm prototype is not a replication. Running it on the same cohort would not be independent confirmation and is not planned. |
| Licence interpretation; access and school approval on file | Blocked on the user | See "Items the user must provide". No approval text was written or inferred. |
| Registration, brain support, label meaning, QC frozen before scoring | Blocked | No Duke protocol is frozen. A freeze needs a cohort that has not been scored; none exists in this dataset for the primary endpoint. |
| Duke-specific end-to-end pipeline | Blocked | Prototypes are synthetic-only and not wired to registration. Building further is pointless for confirmation while the only cohort is spent. |
| Two-history control is not a proven strong comparator | Fixed (label added) | `output/track_b_two_history_screen_20261102/screen.json`, synthetic, 27 cases per regime, pre-boundary Dice at true volume: default control 0.8854 vs simpler last-mask control 0.8924 (growth; paired difference -0.0071, SE 0.0065) and 0.7607 vs 0.7965 (shrinkage; -0.0358, SE 0.0244); identical in plateau. Some grid settings scored slightly higher, but the grid was inspected after the screen, so no setting is carried forward. Note added to `tracks/track-b/README.md`. |
| Runtime for a representative Duke case | Blocked | Timing a case means reading patient images. Not run under the current instruction. |

## Track C

| Gap | Status | Evidence and action |
|---|---|---|
| Assumed biology and signal quality | Irreducible | Simulation by design. |
| Progression rules are simulation definitions | Irreducible | `observed_two_scan` is a synthetic proxy, not formal RANO. |
| No validated simultaneous interval method | Blocked | Max-T failed all 24 cells; Student-t with Bonferroni passed 14/24, and 35/44 in the world-count extension (`START_HERE.md`, 2026-11-08 bullets). |
| GE-betting is only a prototype | Blocked | `src/trackc_ge_betting.py`: 10 synthetic checks, unwired, no independent review. It needs a written protocol, an outside review, and an outcome-blind coverage check under the exact world-generation and random-stream rules. |
| 14-day precision target is infeasible at current seed capacity | Irreducible at current capacity; target kept | `paper/trackc_population_mean_feasibility_audit_2026-11-08.md`: at least 207 worlds per family are needed to clear the all-zero-data lower bound; the seed block allows 7 per family. The 14-day target is unchanged. **Confirmatory Track C evaluation is stopped** until the user either authorizes more worlds (with fresh reserved seeds and a runtime cap) or explicitly changes the target. Timing extrapolation for 100 worlds per family is about 79 hours from a two-patient-per-cohort pilot. |
| Seed schedule not materialized or reserved | Blocked | Depends on the world count, which depends on the decision above. Reserving seeds now would fix a count that cannot meet the target. |
| Full run driver, penalty-sensitivity, bootstrap orchestration missing | Blocked | Listed in `not_yet_locked_before_evaluation` in `configs/trackABC_integrated_v2.json`. Not built: with evaluation stopped, a driver would have nothing valid to run. |
| Settings file is not loaded and enforced by a full runner | Blocked | Config status reads "not an execution seal". Commit `efc4ea6` added sampler/config parity checks in `tests/test_trackc_integrated_sampler.py`; those were not duplicated. Full enforcement belongs in the driver. |

## Housekeeping

| Gap | Status | Evidence and action |
|---|---|---|
| Many unrelated modified and untracked files | Already recorded | 130 status entries on 2026-11-09; 14 modified tracked files. None were touched, staged, or reverted for this record. |
| Zenodo token pasted into chat | Blocked on the user | No token string was found in tracked files or in the searched untracked text files. Whether it is still active cannot be checked from here. |

## User report received (2026-11-09)

The user states: the school approved the project, the user obtained approval from the dataset administrator, and "the file said I can share it with anyone". This is recorded as **user-reported**. No approval or licence document is saved in the project, and none was written or reconstructed here. Two points stay unverified until the documents are on file: (a) which file contains the sharing statement and what exactly it covers (raw data, derived numbers, figures, or model weights); (b) how it fits the Zenodo record, which lists CC BY-NC-ND 4.0 or a Duke custom licence, where CC BY-NC-ND does not permit sharing adapted material. Until the file is saved and read, Duke-derived outputs stay out of Git and out of shared reports. Items 1 to 3 below change from "not provided" to "reported, document not on file".

**Second user report (2026-11-09):** the user states that the file says the user may change the version (make adapted material). This is also user-reported; the file is not saved in the project and was not read here. Permission to make adaptations is consistent with CC BY-NC-ND 4.0 and does not by itself grant permission to share adapted material. Local project use was already user-authorized, so this report unblocks no new analysis. The sharing question for derived numbers, figures, and model weights stays open until the file is on file and read.

## Items the user must provide

1. Which Duke access option was accepted (CC BY-NC-ND 4.0 or the Duke custom licence), with the written terms or a screenshot of the grant.
2. The school approval record for use of the Duke dataset.
3. A decision on whether Duke-derived aggregate numbers, figures, or model weights may be shared. Until then they stay out of Git and out of shared reports.
4. Revocation of the pasted Zenodo token and creation of a new one, kept out of chat and out of the repository.
5. A decision on `paper/trackABC_integrated_research_plan_revised_2026-11-08.md`: active, proposed, or superseded.
6. For Track C: authorize more worlds with a runtime cap, or explicitly change the 14-day target. Without one of these, confirmatory evaluation stays stopped.
