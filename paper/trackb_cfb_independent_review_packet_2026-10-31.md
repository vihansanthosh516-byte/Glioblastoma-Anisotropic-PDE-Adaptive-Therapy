# Track B CFB-GBM independent MRI review packet

**Status:** preparation only. No independent review has occurred. No registration pass is recorded. This packet does not authorize opening GTV masks or scoring forecasts.

## Purpose

Record an independent, target-blind assessment of whether the selected CFB-GBM MRI visits are sufficiently aligned for the planned two-history-to-next-visit PDE forecast. Review gross anatomical alignment only. This is a technical dataset-quality gate, not a clinical read or a claim that tumor boundaries are correctly registered.

## Prepared material

- Image aid: `output/cfb_mri_paginated_review_full_20261031/`
- Manifest: `output/cfb_mri_paginated_review_full_20261031/manifest.json`
- Blank pair-level form: `paper/trackb_cfb_registration_review_template_2026-10-31.csv`
- Candidate set: 38 pseudonymous patient IDs, each with t0, t1 and t2 T1-Gd images; review t0-to-t1 and t0-to-t2 separately (76 rows).

The image aid contains 760 PNG pages and 5,668 sampled slices. It displays T1-Gd only, at 320 pixels per visit, sampled at approximately 5-mm spacing. It does not contain full-volume overlays or landmarks. Use it for navigation and gross triage only. **It cannot establish a formal registration pass.** Formal review needs the full-resolution 3D images in a viewer that preserves physical coordinates and supports side-by-side or checkerboard inspection.

## Blinding and access limits

The reviewer will not see GTV arrays or contours, forecast outputs, model scores, or patient-level outcomes. The reviewer will use only the three T1-Gd MRI volumes and non-tumor anatomical landmarks. Do not open GTV or brain-mask voxel arrays as part of this review. The existing aid-generation manifest records that neither GTV nor brain-mask arrays were opened.

The reviewer should be independent of model development and scoring, and should have relevant MRI or image-registration experience. Record a reviewer code and qualifications without placing personal contact details in the project output.

## Freeze before review

The current protocol draft does not yet specify a validated numeric tolerance. Therefore, do not mark any row “pass” or “fail” until the reviewer and project lead have frozen:

1. Landmark definitions and the minimum number and spatial distribution needed per pair.
2. Landmark visibility rules and physical-coordinate measurement method.
3. Maximum tolerable landmark displacement in millimeters, with justification tied to the image resolution and intended PDE endpoint.
4. Rules for postoperative anatomy change, missing coverage, distortion, motion and artifacts.
5. Pass, fail and indeterminate rules; reviewer disagreement and adjudication; and treatment of unreviewed cases.

Pilot the frozen rubric on a small set chosen without target masks or model outcomes. Revise it only during the pilot, then freeze it before reviewing the remaining pairs. Preserve every pilot decision and protocol version. If no defensible tolerance or qualified reviewer is available, leave the registration gate unresolved.

## Review procedure

For each row in the CSV, open the pair's full-resolution t0 and comparison visit (t1 or t2) T1-Gd volumes. Confirm orientation and physical coordinates. Assess predeclared non-tumor landmarks across the brain, including central and peripheral structures where visible. Record landmark locations and displacement in millimeters under the frozen rubric. Note systematic translation or rotation, local distortion, missing coverage and artifacts. Assign only the frozen rubric's outcome. Use “indeterminate” when visibility, coverage or evidence is insufficient; do not convert missing evidence to a pass.

Review t0-to-t1 and t0-to-t2 independently. The 5-mm page aid may help locate regions but must not replace full-volume review. Retain the reviewer's full case-level form, exclusions, disagreements and reasons.

## Separate target-support gate

Even a gross-alignment pass does not establish that GTV labels have valid negative support throughout the scoring grid, preserve the intended target definition, or have acceptable provenance. Those questions require separate source evidence and a separately authorized mask-access decision. Do not use MRI appearance, a brain mask, or an empty label region as proof of valid negative tumor labels.

## Current decision

No patient mask will be scored yet. After the independent review and separate target-support/provenance audit, update the Track B protocol and decide whether any CFB cases meet the locked eligibility rules. If no cases pass, report the gate failure and keep CFB out of the PDE patient-scoring analysis.
