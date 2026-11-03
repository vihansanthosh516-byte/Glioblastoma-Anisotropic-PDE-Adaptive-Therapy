# Track B Duke rights and approval ledger (2026-11-02)

Term note: this ledger says **sessions** (a scan date with its images and mask). Other project files say scans or visits; they mean the same unit here, but that is not checked.

**Status:** written-records check only. Nothing downloaded. Nobody contacted. No patient data opened. Public pages were read on 2026-11-02 (WebFetch, read-only).

## 1. Duke access and license

| Item | Status | Source |
|---|---|---|
| Written access grant | **Not on file. Status: pending, not granted.** | No approval file in `paper/`, `data/manifests/`, `vault/` (grep 2026-11-02). `paper/track_b_source_gate_evidence_2026-10-31.md:24` |
| History | Duke declined the first request (missing use and affiliation details). User reports a resubmission on 2026-10-29. No written reply saved. | `paper/track_b_pde_external_protocol_draft_2026-10-28.md:54`; `paper/track_b_source_gate_evidence_2026-10-31.md:44` |
| Public page does not prove access | Right. The page is public. The files are restricted. | https://zenodo.org/records/20548528 |

What the Zenodo page says (read 2026-11-02):
- Title: A Longitudinal Glioblastoma MRI Dataset with Anatomical Landmark Pairs for DIR Validation. Version v1, published 2026-06-12, modified 2026-09-23.
- License: CC BY-NC-ND 4.0. The page also offers a custom Duke license (Office for Translation & Commercialization) for other uses.
- Access: restricted. Request form and login needed.
- Size: 136.3 GB. (An older local note says 119.3 GB, `START_HERE.md:213`. The page changed. Use 136.3 GB.)
- Counts: 61 patients, 3 to 18 sessions each, 503 sessions.
- Six image types per session: T1, T2, FLAIR, ADC, DWI, T1-c. 1 mm, skull-stripped, MNI-aligned NIfTI.
- Label map: 1 necrotic core, 2 non-enhancing T2/FLAIR, 3 enhancing tumor, 4 resection cavity.
- Metadata: per-patient landmark CSVs; patient CSVs with scan dates, demographics and treatment timelines (dates relative to first scan); a mask per session.
- File names and sizes: not listed. The page does not say which files are public. All files appear under "Restricted". So **public vs restricted per file: unknown**; treat all as restricted.
- Derived-output terms: the official CC BY-NC-ND 4.0 legal code permits producing and reproducing Adapted Material for noncommercial purposes, but does not grant permission to share Adapted Material. It does **not** categorically forbid creating derivatives. Whether forecasts, masks, trained models or summary results are Adapted Material, and how the separately offered Duke custom license treats them, is **unknown**. The custom license text is **not seen**. Which license applies to this project is **unknown**. Source: https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode.en (Sections 1 and 2(a)(1)).

What the GitHub README says (https://github.com/deshanyang/GBM-DIR-QA, read 2026-11-02):
- No license stated. No registration step stated. Calls the Zenodo repository "publicly available" with a blank link. This conflicts with the Zenodo "Restricted" label. Zenodo is the authority.
- Scan dates are relative to the first scan, auto-extracted, "may be imperfect or offset by a few days", not individually validated.
- 61 GBM patients from Duke; says pre-op and post-recurrence scans with landmarks. A contact is listed there (not recorded here, not used).
- Atlas alignment versus session-to-session registration: **unresolved**. The page says MNI-aligned (atlas alignment); that does not show the sessions are registered to each other. An earlier local record says sessions are not registered to each other: `paper/track_b_external_cohort_screen_2026-10-28.md:41`. The 2026-11-02 README fetch did not repeat it.

## 2. School approval

**Status: user-reported, not passed.**
- User reported on 2026-10-31 that the school said this use is OK. No written copy saved. No date, approver or conditions on file. `paper/track_b_source_gate_evidence_2026-10-31.md:25,43,57`; `paper/track_b_amendment_decision_2026-10-31.md:115`.
- It is not known whether any approval names Duke, the license, or derived outputs.
- A search of the repo for approval/agreement files found none (find, depth 3, 2026-11-02).

## 3. Eligible triplet count

**Unknown (blocked).**
- No Duke folder in `data/external/` (listing: no duke or GBM-DIR entry). `data/manifests/` has no Duke file (listing: consort_mu, forecast_pairs_mu, patient_manifest, scan_manifest, split_mu, trackA_donors only).
- The 61 patients and 503 sessions are dataset totals, not a triplet count. Each patient has 3 to 18 sessions, but surgery, treatment, date and label checks are not done.
- The precision gate code refuses a hand-entered count and needs a validated metadata report: `src/fk_precision_gate.py:113-144`. No such report exists.

## 4. Prior use of Duke data in this project

**No evidence of Duke patient data being used.**
- No Duke data files found locally (see 3).
- Records say no Duke file was read: `START_HERE.md:205`; `paper/track_b_source_gate_evidence_2026-10-31.md:5`.
- Duke appears only in documents, protocol drafts and the gate code (`src/fk_precision_gate.py:25`), tested on synthetic inputs.
- Patient overlap of Duke with MU, UCSD-PTGBM, CFB-GBM, PREDICT-GBM or LUMIERE is **unchecked** and cannot be checked without authorized metadata: `paper/track_b_pde_external_protocol_draft_2026-10-28.md:22`. Independence is not established.

## Gate result

| Gate | Status |
|---|---|
| Duke written access | Open (pending) |
| Duke license and derived-output terms | Open (CC BY-NC-ND 4.0 or custom; custom text unseen) |
| School written approval | User-reported, not passed |
| Eligible triplet count | Unknown (blocked) |
| Prior-use overlap | No local use found; patient overlap unknown |

Decision: **no-go for any restricted-file access** until written Duke access and written school approval are saved here. Weakest point: every "open" item rests on absence of a file in the repo, and the user may hold emails not saved to the repo.
