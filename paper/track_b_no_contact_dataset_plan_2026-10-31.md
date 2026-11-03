# Track B LUMIERE exploratory plan (2026-10-31; ADOPTED 2026-11-01 for exploratory work only)

**User decisions (2026-11-01):** (1) adopt this plan for exploratory Track B work only; (2) adopt amendment A-2026-10-31-R, U-Net reuse only, under the full-fingerprint match rule, with the central-gap refit separate; (3) keep the full-grid U-Net required; the crop is not adopted. These decisions do **not** authorize patient scoring, target-mask access or reading MRI/mask voxels. Ask the user before reading any patient MRI or mask voxel data. LUMIERE results stay exploratory, not fresh confirmation or proof of novelty.

**Working constraint for this step (2026-10-31):** no outreach to authors, Duke, TCIA or the school during this step. This is a temporary working constraint for this step, not a standing project rule. (Correction: an earlier version of this file called it a standing user rule. That was wrong.)

**Status.** The user chose LUMIERE for **exploratory** Track B development, subject to a metadata-only cohort check. LUMIERE is not fresh confirmation: this project used it before (`src/91_lumiere_volume_forecast.py`, `src/130_h3v_frozen_lumiere.py`). **Exposure list (added 2026-11-01 audit, `paper/trackAB_audit_and_gap_closure.md` B-G1):** LUMIERE was part of the spent H-2 test set (58 LUMIERE patients, pre-op inputs; `vault/PAPER_FINDINGS_LEDGER.md` 2ax); script 134 reused H-2 LUMIERE logits post hoc as development data; scripts 104 and 106 compared LUMIERE segmentations (ledger 2ab, 2aj); script 91 forecast LUMIERE volumes (ledger 2l). PREDICT-GBM also includes 61 LUMIERE patients (`paper/track_b_prior_art_followup_2026-10-31.md`). Any LUMIERE result here is an exploratory re-analysis of an exposed public cohort, not external or independent validation. LUMIERE does not raise novelty by itself.

No target-label voxel, mask voxel or image voxel was read for this plan. Nothing was downloaded: the LUMIERE files in `data/external/lumiere/` were already on disk (dated 2026-10-04/05). No patient was scored. The no-scoring rule and the full-grid U-Net requirement are unchanged. Amendment A-2026-10-31-R was adopted 2026-11-01 (reuse only).

## 1. Why LUMIERE

The fresh-cohort search found no untouched public GBM cohort with three or more masked visits and no access request (`paper/track_b_external_cohort_screen_2026-10-28.md`; new 2026-10-31 check: OpenNeuro ds007045 has one preoperative MRI per patient, PMC13498554). LUMIERE is public and many-visit. Any result is exploratory.

## 2. License and reuse conditions (verified; conflict recorded, not resolved)

| Source | Stated terms |
|---|---|
| Figshare API, collection 10.6084/m9.figshare.c.5904905 (v1, 2022-12-13) | Collection-level license field: none. |
| Figshare items 21249516 (MRI and segmentations), 21195556 (expert rating), 21195559 (demographics), 21266241 (readme), 21195430, 21195562, 21187030, 21187033, 21251643 | Each item: **CC0**. |
| `LUMIERE-readme.pdf` (local copy of item 21266241), last line | "This dataset is provided for non-commercial use." |
| Data descriptor, Suter et al., *Sci Data* 9, 768 (2022), PMC9755255 | Article text CC BY 4.0. Usage note: "This dataset may only be used for non-commercial purposes." |

**Reading:** the item metadata says CC0, and both the readme and the paper add a non-commercial condition. This project is non-commercial, so it meets the stricter reading. **No source states rules for derived outputs** (trained weights, derived masks, figures). This plan does not infer any. Until a source states them, do not publish or share derived masks or images; report summary numbers only.

## 3. Segmentation tool and tumor-core label mapping

- **Tool for the target:** DeepBraTumIA (https://www.nitrc.org/projects/deepbratumia), atlas-space output `DeepBraTumIA-segmentation/atlas/segmentation/seg_mask.nii.gz` (readme). HD-GLIO-AUTO is not used for the target.
- **Version:** **not stated** in the paper or readme. The only version given is PyRadiomics 3.0.1. Record this as a gap; do not guess.
- **Official-documentation check (2026-11-01):** ARTORG MIA software page (https://www.artorg.unibe.ch/research/mia/software/index_eng.html) names the sub-compartments "necrotic, contrast-enhancing, and edema" but states **no version and no label values** (no label 0, background or cavity). The NITRC project page (https://www.nitrc.org/projects/deepbratumia) refused the connection twice on 2026-11-01, so it is **not checked**. Version and label-0 meaning stay **OPEN**. Nothing is inferred.
- **NITRC recheck (2026-11-02):** the page loaded. It states no tool version (the "v2.1.76-2" footer is the NITRC site) and no label values. It gives the tool license as "Attribution Non-Commercial No Derivatives"; whether that reaches LUMIERE masks is not stated. Gates stay **OPEN** (`paper/track_b_lumiere_go_no_go_2026-11-02.md`).
- **Label mapping** (from the published `LUMIERE-pyradiomics-deepbratumia-features.csv`, columns `Label`/`Label name`; no mask voxel read): 1 = Necrosis, 2 = Contrast-enhancing, 3 = Edema. HD-GLIO-AUTO: 1 = Non-enhancing, 2 = Contrast-enhancing.
- **Tumor core (draft):** DeepBraTumIA labels 1 or 2. Edema (3) is not core. Label 0 is not defined in the sources; by the tool output, it is "no tumor label", not a reviewed negative.
- **Resection cavity:** DeepBraTumIA has **no cavity label** in the published mapping. Cavity is therefore not separated from background. Record as a limitation.

**Automated-mask limitations (keep prominent in every result):**
- The target is an automated DeepBraTumIA mask. There is no expert voxel mask. A forecast scored against it measures agreement with the tool, not with true tumor.
- The tool was built for pre-operative scans; post-operative performance is not reported for LUMIERE.
- No cavity label; cavity and blood products may be mislabelled as necrosis or enhancement.
- The paper notes that segmentation outliers were deliberately kept.

## 4. Valid support: origin-only rule (draft)

The released images were already skull-stripped with HD-BET, and the release includes the brain mask used (`DeepBraTumIA-segmentation/atlas/skull_strip/brain_mask.nii.gz`, readme). This plan does **not** rerun HD-BET and does **not** claim the mask is validated.

**Rule:** support S = the released atlas-space brain mask of the **origin visit (visit 2)**, plus any history-visit label voxel outside it (so a history core is never clipped). No target-visit file enters S. Code: `src/fk_longitudinal_prep.py:origin_support` (its signature has no target input; tested).

**Label in results:** "Support = released origin-visit skull-strip mask; not independently validated."

Target core outside S counts as a miss and is reported per patient.

## 5. MRI-only registration and QC (draft, frozen before any target mask is opened)

- DeepBraTumIA already maps each visit to the atlas (readme; per-visit `.tfm`). Visits are not registered to each other.
- **Step:** rigid refinement of the target-visit atlas-space CT1 onto the origin-visit atlas-space CT1, by NCC inside the origin support. Inputs are MRI only; no label map is an input. Code: `src/fk_longitudinal_prep.py:register_rigid_mri`.
- **QC (frozen draft constants):** pass only if NCC after ≥ 0.80, maximum rotation ≤ 10°, and shift ≤ 10 mm.
- **Failure handling:** a failed patient is excluded from scoring for all arms. Every exclusion and its reason is reported (`exclusion_report`). A failure is never retried with other settings.
- **Weakest point:** the threshold 0.80 is set on synthetic phantoms only. It is not tuned on LUMIERE. A tumor change between visits lowers NCC, so true growth can cause exclusion. Report the excluded count beside every result.

## 5a. LUMIERE horizon rule (2026-11-01)

- The gap is the **difference in published week labels**: week(v3) − week(v2), from the `week-NNN` folder names.
- Primary window: 11 ≤ week(v3) − week(v2) ≤ 15 (week labels). First usable post-op triplet only.
- LUMIERE dates are released only as week bins. The exact day gap is **unknown**. This rule does **not** claim exact days.
- This rule is **not** the same as the MU/Duke 90-day nominal horizon or any 76–104-day window (`paper/track_b_pde_external_protocol_draft_2026-10-28.md`). Do not convert between them. Results from the two rules are not pooled.
- **Proposed for the protocol freeze (not yet wired or frozen):** the PDE and fixed-growth arms use 7 × (week difference) days as the model time step. This is a **modelling convention**, not a measured gap. All arms would use the same value.

## 5b. Registration QC: frozen settings (frozen 2026-11-01, before any MRI voxel read)

Code: `src/fk_longitudinal_prep.py`, SHA-256 `6588f511d9e9fcc88e26591aea19db38c5673acc6d618188c0b24791ea45d4bf` (working-tree file at commit time). Change only by a dated amendment, before the preflight runs.

| Setting | Frozen value |
|---|---|
| Fixed image | Visit-2 (origin) atlas `skull_strip/ct1_skull_strip.nii.gz` |
| Moving image | Visit-3 atlas `skull_strip/ct1_skull_strip.nii.gz` |
| Fit region | Visit-2 atlas `skull_strip/brain_mask.nii.gz` > 0. History label voxels are not added for QC (no `seg_mask` is read). |
| Grid | Released 1 mm atlas grid, 182 × 218 × 182 (from headers). No resampling before QC. |
| Intensity | Raw values, no normalisation (NCC is invariant to linear intensity change). |
| Transform | Rigid, 6 parameters: 3 Euler angles (deg, x-y-z order, R = Rz·Ry·Rx) about the grid centre, 3 shifts (mm). Start at zero. |
| Interpolation | Linear (`scipy.ndimage.affine_transform`, order 1, constant 0 outside). |
| Optimiser | Powell. Stage 1 on every 2nd voxel (`downsample=2`): xtol 1e-2, ftol 1e-4, maxiter 2000. Stage 2 full grid from stage 1: xtol 1e-3, ftol 1e-5, maxiter 2000. |
| Metric | NCC inside the fit region. |
| Pass rule | NCC after ≥ 0.80 **and** max |angle| ≤ 10° **and** shift norm ≤ 10 mm. |
| Failure | Excluded from all arms. Reason reported. No retry, no other settings. |
| Sanity checks | Shape equals header; all values finite; brain mask non-empty. A failed check is an exclusion with that reason. |

Gate status: settings **frozen**; real-patient QC **BLOCKED** until the user approves the exact preflight (`paper/track_b_readiness_audit_2026-11-01.md` §5).

## 6. Metadata-only inventory

Record: `output/track_b_lumiere_inventory_20261031/inventory.json`. Source: the central directory (file names) of `Imaging-v202211.zip`; no image or mask bytes read.

| Item | Count |
|---|---|
| Patients in the zip | 91 |
| Study dates | 638 |
| Usable studies (CT1, T1, T2, FLAIR + DeepBraTumIA atlas seg + atlas brain mask) | 599 (matches the paper's 599 four-sequence study dates) |
| Patients with ≥ 3 usable post-op visits | 62 |
| All consecutive usable post-op triplets | 325 (not independent; not used) |
| **Primary eligible: first triplet, v2→v3 gap 11–15 weeks** | **28 patients** |

Window: 11–15 in published week labels (see §5a). The exact day gap is unknown because dates are week-binned; an earlier version of this line gave a day range, which is withdrawn. The window is applied to each patient's **first** triplet only, never to a later one that happens to fit. Patients excluded by the window: 34.

**Weakest point:** 28 patients is small. Fold sizes would be about 6, 6, 6, 5, 5. Inventory rules are draft; the file-name check does not prove a file opens or has the right content.

## 6a. Header check (2026-11-01)

`src/fk_lumiere_headers.py` opens each NIfTI file inside the zip and reads only its 348-byte header. Record: `output/track_b_lumiere_inventory_20261031/header_check.json`. Files per visit: DeepBraTumIA atlas `seg_mask`, atlas `brain_mask`, and the four atlas skull-stripped images.
- 28 eligible patients confirmed (63 excluded), 84 visits, 504 headers read, 0 problems.
- Every file is on one grid: 182 × 218 × 182 at 1.0 mm (MNI atlas). At 2 mm: 91 × 109 × 91.
- Tests: `tests/test_fk_lumiere_headers.py` (a spy confirms each read is exactly 348 bytes; synthetic zip).
- Weak point: a header check does not prove the voxel content is valid.

## 7a. Re-timed runtime at the LUMIERE grid (measured synthetic, 2026-11-01)

Record: `output/track_b_timing_pilot_20261031/lumiere_timing_record.json`. Synthetic volumes on the 2 mm protocol grid 91 × 109 × 91 (padded to 112 × 112 × 96 for the U-Net). Source headers are 182 × 218 × 182 at 1 mm, same PC, torch 2.14.0+cpu, 8 threads.

| Measured (synthetic) | Value |
|---|---|
| Full-grid U-Net step | 7.973 s mean (6.948–10.492), 40 steps |
| U-Net prediction, 3 seeds, 1 patient | 5.27 s |
| PDE full forecast, 1 patient | 1,250.4 s |
| PDE calibration-only, 1 patient | 1,544.3 s |

| Projected (not measured) | Wall-clock hours on this 8-thread PC |
|---|---|
| U-Net, one analysis (112 pairs × 300 × 3) | 223.3 (max step 293.8) |
| PDE, one analysis (28 full + 112 calibration-only), **lower bound** | 57.8 |
| Four analyses, no reuse | 1,124.1 |
| **Four analyses, amendment R reuse (adopted)** | **677.6** |

Unit note (2026-11-01 audit, B-G4): these are wall-clock hours with 8 Torch threads, not single-CPU hours. 677.6 h is about 28 days of continuous running; 1,124.1 h is about 47 days. The max-step figure gives about 819 h with reuse. Analyses proposed in the unadopted scoring draft (history-gap sensitivity, time-aware baselines) are not included.

Weak points: the PDE full forecast was faster than calibration-only here, the reverse of the first pilot; PDE time depends on the synthetic case, so one case is a weak estimate. The U-Net step was slower than in the earlier pilots (5.34–6.94 s). The central-gap U-Net is counted at full size. These are resource needs, not a deadline. **Section 7 below is superseded by this section.**

## 7. Runtime projection (HISTORICAL; superseded by §7a — do not cite)

Assumptions: 28 patients, 5 outer folds, so 112 training pairs per analysis (each patient trains in 4 folds). U-Net 112 × 300 epochs × 3 seeds = 100,800 steps at the measured 5.340 s full-grid step (`crop_timing_record.json`). PDE 28 full + 112 calibration-only forecasts at 3,128.6 s and 1,808.7 s (`timing_record.json`); a lower bound.

| Scenario | U-Net h | PDE h (lower bound) | Total CPU h, four analyses |
|---|---|---|---|
| Frozen, no reuse | 598.1 | 322.4 | 920.5 |
| Frozen, with reuse (amendment R; draft status at the time) | 299.0 | 322.4 | 621.5 |

**Weakest point:** the steps were timed on a synthetic 120×120×78 grid. The LUMIERE atlas grid at 2 mm may differ in size, so step time may differ. Re-time on a synthetic volume of the true grid shape before freezing. These are resource needs, not a deadline.

## 8. Gates

| Gate | Status |
|---|---|
| License | Recorded; CC0 items plus non-commercial condition; derived-output rules **not stated**. |
| Label meaning | Mapping from published CSV; tool version **not stated**; label 0 and cavity undefined. Official-doc check 2026-11-01: still open (§3). |
| Valid support | Draft origin-only rule; source mask not validated. |
| Horizon | §5a: 11–15 week-label difference; not a day window. |
| Registration | Settings frozen (§5b); synthetic tests only. Real-patient QC **BLOCKED** until the user approves the preflight. **Not passed.** |
| U-Net reuse (amendment R) | Adopted; wired with full-fingerprint match; central-gap separate. Synthetic tests only. |
| Inventory | Done, metadata and headers only (28 eligible; one 182×218×182 grid). |
| School-use approval | **User-reported approval (2026-10-31): the user says the school said this use is OK.** Not yet a pass under the written-evidence rule (`paper/track_b_source_gate_evidence_2026-10-31.md`): save the written approval (email or form, with date, approver and any conditions) in the project, then mark it passed. |
| Protocol freeze and hashes | **Not done. Open.** |
| Patient voxel reads and scoring | **Blocked.** Ask the user first. |

**No scoring** until the gates above are frozen by a dated amendment, and the code and cohort are hashed.

## 9. Code done in this step (synthetic tests only)

`src/fk_longitudinal_prep.py`: `Visit`, `primary_triplets`, `origin_support`, `register_rigid_mri`, `apply_rigid`, `registration_qc`, `exclusion_report`. Tests: `tests/test_fk_longitudinal_prep.py` (5 passed, scratch copy, Python 3.14). `src/fk_lumiere_headers.py`: header-only reader (no voxel arrays). U-Net reuse wired in `src/fk_learned_baselines.py` and `src/fk_mu_patient_runner.py`. Not built: a voxel reader (needs the user's OK first). Protocol freeze and hashes: not done.
