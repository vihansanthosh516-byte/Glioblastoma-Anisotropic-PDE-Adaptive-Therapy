# Track B amendment draft: visit-1 alignment for LUMIERE (2026-11-10)

**Status: DRAFT, NOT ADOPTED.** Must be decided before any coverage count, history-input build, or target read. Written without any patient coverage count; no setting here was chosen from a cohort size.

## Problem

The frozen plan registers only visit 3 to visit 2 (`paper/track_b_no_contact_dataset_plan_2026-10-31.md` section 5b). Visit 1 enters the history inputs (FiLM channel 1, PDE start state, transition history) and the coverage common support, but it is aligned to visit 2 only through the publisher's per-visit native-to-atlas registration. That alignment has no QC in this project.

## Options

| | 1. MRI-only rigid registration of visit 1 to visit 2 | 2. Accept publisher atlas alignment for visit 1 |
|---|---|---|
| Method | Same as visit 3: `register_rigid_mri`, fixed = visit-2 atlas CT1, moving = visit-1 atlas CT1, fit region visit-2 atlas brain mask, frozen settings and QC (NCC >= 0.80, max Euler angle <= 10 deg, shift <= 10 mm) | Identity in atlas space; no new step |
| Reads | +28 visit-1 atlas CT1 files (MRI only; no labels) | none |
| History inputs | Visit-1 label map resampled with the visit-1 rigid transform (nearest neighbour) before the 2 mm step | Unchanged |
| Coverage chain (visit 1) | includes the rigid factor, as for visit 3 | no rigid factor |
| QC failure | patient excluded from all arms before any target read; reason reported; no retry | not applicable |
| Limits | Rigid only; thresholds came from synthetic phantoms; more exclusions possible, which can push the cohort below the frozen minimum of 25 (then stop) | Unverified; treats visit 1 differently from visit 3; residual misalignment biases history inputs and coverage silently |

## Recommendation (not adopted)

Option 1. It applies one rule to every non-origin visit and makes visit-1 alignment checkable. Its cost (possible extra exclusions) is measured and reported, not hidden. The decision must be made before any count, and must not be revisited after the count is known.

## Text to adopt (if Option 1 is chosen)

> Visit 1 is rigidly registered to visit 2 with `register_rigid_mri` using the frozen section 5b settings and QC, with visit-1 atlas CT1 as the moving image. Patients failing visit-1 or visit-3 QC are excluded from all arms before any target read. The visit-1 history label map and the visit-1 coverage chain use the visit-1 rigid transform. The read scope adds the 28 visit-1 atlas CT1 files. No tumor label is read for registration.
