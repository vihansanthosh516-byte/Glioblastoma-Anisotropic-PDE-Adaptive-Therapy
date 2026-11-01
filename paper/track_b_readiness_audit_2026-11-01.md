# Track B readiness audit, no voxel reads (2026-11-01)

**Scope.** Documents, code and header records only. No data was downloaded. No one was contacted. No patient MRI or mask voxel was read. No patient was scored. Synthetic tests are not counted as real-data gate passes.

Sources checked: `START_HERE.md`, `paper/trackABC_execution_plan_2026-10-31.md`, `paper/track_b_no_contact_dataset_plan_2026-10-31.md` (the "LUMIERE plan"), `paper/track_b_source_gate_evidence_2026-10-31.md`, `paper/track_b_amendment_decision_2026-10-31.md`, `paper/track_b_fixed_growth_baseline_spec_2026-10-31.md`, `paper/track_b_pde_external_protocol_draft_2026-10-28.md`, `src/fk_longitudinal_prep.py`, `src/fk_lumiere_headers.py`, `output/track_b_lumiere_inventory_20261031/header_check.json`.

## 1. Gate table

| Gate | Status | Exact evidence | What would close it |
|---|---|---|---|
| License and derived outputs | **OPEN** (use allowed; derived-output terms unknown) | LUMIERE plan §2: Figshare items CC0; `LUMIERE-readme.pdf` and Suter et al. *Sci Data* 9:768 (2022) say non-commercial use only. No source states rules for derived masks, images or weights. | Published text. Until then: non-commercial use, summary numbers only, no shared derived masks or images. |
| Label map | **OPEN** | LUMIERE plan §3: map from published `LUMIERE-pyradiomics-deepbratumia-features.csv` — 1 Necrosis, 2 Contrast-enhancing, 3 Edema; core = 1\|2. DeepBraTumIA version not stated. Label 0 not defined. No cavity label. Masks are automated. | Version and label-0 meaning from published documents. If none exist, keep as stated limits; the gate cannot become PASS from documents found so far. |
| Support rule | **OPEN** | LUMIERE plan §4; `src/fk_longitudinal_prep.py:origin_support` (no target input; synthetic test only). Support = released origin-visit atlas `brain_mask` + history label voxels. Mask is HD-BET output, not validated. | Freeze in the protocol. Real-data check (mask present, non-empty, on grid) needs MRI-access approval. |
| Registration method and QC | **Method: OPEN (draft, synthetic only). Real-patient QC: BLOCKED** | LUMIERE plan §5; `register_rigid_mri`, `registration_qc`; thresholds NCC ≥ 0.80, rotation ≤ 10°, shift ≤ 10 mm, set on synthetic phantoms only. | User authorizes MRI voxel reads (§3 proposal), then run QC once with frozen settings. |
| School approval | **OPEN** | Source-gate record contact log: user report 2026-10-31 that the school said OK. No written copy on file. | Save the written approval in the project with date, approver and conditions. |
| Protocol freeze and hashes | **OPEN** | No freeze file exists. No code, cohort or config hashes recorded. | Dated freeze amendment with SHA-256 of code, cohort list, thresholds and record files, before any target mask is opened. |
| Header inventory | **PASS for headers only** | `header_check.json`: 28 eligible, 63 excluded, 84 visits, 504 headers read, one grid 182×218×182 at 1 mm, 0 problems. | Voxel content stays **unchecked**. That part needs MRI access. |
| U-Net reuse (amendment R) | Adopted; wired; synthetic tests only | Amendment memo §5. | Not a data gate. |

## 2. What can close without patient data, and what cannot

**From published documents only:**
- Derived-output license terms.
- DeepBraTumIA version; label-0 meaning.
- Writing and hashing the protocol freeze (rules, code, cohort list).

**Needs your written school evidence:**
- School approval gate.

**Needs your OK to read patient MRI voxels:**
- Real-patient registration QC.
- Support-mask content check (non-empty, inside grid).
- Voxel-level check that files open and match their headers.

**Needs target-mask access (not asked for now):**
- Any scoring. This stays blocked until the freeze is done and you approve.

## 3. Stale or conflicting claims found

| File | Line | Problem |
|---|---|---|
| `paper/track_b_source_gate_evidence_2026-10-31.md` | 3 | Says amendment R "is still **not adopted**". Stale: adopted 2026-11-01. |
| `paper/track_b_source_gate_evidence_2026-10-31.md` | 23 | School gate: "No local record found." Conflicts with its own line 41 (user-reported OK). Status OPEN is right; wording is stale. |
| `paper/track_b_source_gate_evidence_2026-10-31.md` | whole | Gate table covers the MU/Duke path only. It has no LUMIERE rows. The LUMIERE gates live in the LUMIERE plan §8. A reader may think MU gates apply to LUMIERE. |
| `paper/track_b_fixed_growth_baseline_spec_2026-10-31.md` | 3 | "not frozen until the support gate passes" — does not say which cohort's support gate. Consistent but ambiguous. |
| `paper/track_b_pde_external_protocol_draft_2026-10-28.md` | 3 | MU/Duke candidate protocol, 90-day horizon. LUMIERE uses an 11–15-week window. No line says which protocol governs LUMIERE. The freeze must state it. |

No conflict found in `START_HERE.md`, the roadmap, or the LUMIERE plan after commit `4e5073e`. None of them marks a synthetic check as a real-data pass. I did not edit the files above. `START_HERE.md` has your uncommitted Track A change; I did not touch it.

## 4. Proposal: read-only MRI-QC preflight (NOT RUN; needs your OK)

**Purpose.** Run the frozen registration QC once on real MRI. Decide which of the 28 patients pass, before any target mask is opened.

**Order.**
1. Freeze the QC rule first: thresholds, `downsample=2`, Powell settings, code hash. No change after the run.
2. Run the preflight once.
3. Freeze the cohort (passed patient list) in the protocol.

**Files it would read** (inside `data/external/lumiere/Imaging-v202211.zip`, 28 eligible patients, visits 2 and 3 of the first triplet only):
- Visit 2 (origin): `DeepBraTumIA-segmentation/atlas/skull_strip/ct1_skull_strip.nii.gz` and `.../skull_strip/brain_mask.nii.gz`.
- Visit 3 (target visit): `.../skull_strip/ct1_skull_strip.nii.gz` only.
- Total: 84 files.

**Files it would NOT read:** any `seg_mask` (no label voxels, history or target), T1, T2, FLAIR, native-space files.

**Voxel data inspected:**
- CT1 intensities at visits 2 and 3, downsampled to 2 mm, for the rigid NCC fit.
- Origin brain mask as the fit region.
- Sanity checks: shape equals header, all values finite, brain mask non-empty.
- Limit: support for QC is the brain mask only. History label voxels are left out because reading them needs `seg_mask`. They are added later for scoring support only.

**What it would save** (`output/track_b_lumiere_mri_qc_preflight_<date>/qc.json`, refuses to overwrite):
- Per patient: NCC before and after, rotation (degrees), shift (mm), PASS or FAIL, failure reason, sanity-check flags.
- Totals: passed, failed, excluded.
- SHA-256 of the code, the zip and each file read; run environment.
- No voxel arrays, no images, no figures, no transforms as images. Nothing to share outside the project.

**Weak points to keep beside the result:**
- Thresholds come from synthetic phantoms, not LUMIERE.
- Tumor change between visits lowers NCC, so real growth can cause exclusion. Report the excluded count with every result.
- CT1 at visit 3 shows the target tumor. The preflight reads it but saves only QC numbers. It must not be used to tune anything.

**Your decision needed before it can run:** written OK to read the 84 CT1 and brain-mask files listed above.
