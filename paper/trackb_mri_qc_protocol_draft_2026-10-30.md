# Track B CFB MRI-only registration review: draft

Status: proposed, not frozen and not a registration pass. This targets gross anatomical alignment only and does not authorize GTV access or scoring.

## Inputs and cohort

Use the 38 unique patients in the existing candidate crosswalk whose week-derived t1-to-t2 intervals fall in the 76–104-day filter. Their table-derived values are 77, 84, 91 and 98 days; exact acquisition-day eligibility is unverified. Hash the exact crosswalk, each t0/t1/t2 T1-Gd image and each accompanying brain mask before review. Record missing, corrupt or unusable files. Do not read a GTV array, use a tumor contour to select slices or landmarks, or calculate performance on target masks. The existing header check finds all 114 image/mask pairs and all 38 MRI triplets share header geometry; this alone does not demonstrate alignment.

## Proposed review

1. Freeze the software version, normalization, interpolation, physical-coordinate handling and image presentation before new reviews. Inspect t0-to-t1 and t0-to-t2 for every candidate, using only T1-Gd images and brain masks. Show whole-volume views in axial, coronal and sagittal planes, with checkerboards or overlays and fixed slice spacing. Clearly label visit pairs and preserve all case-level review records.
2. Calculate brain-mask overlap and MRI image-similarity measures only as descriptive triage. The earlier 107-case distribution and its IQR-based flags can prioritize attention; they are not clinically validated acceptance thresholds. Mask source ambiguity and treatment-related anatomy change limit their interpretation.
3. Have independent reviewers assess visible non-tumor anatomical landmarks and systematic translation, rotation, distortion, missing coverage and image artifacts. Specify landmark definitions, visibility requirements, physical-distance tolerances and review instructions in advance. Proposed tolerances must be justified and piloted independently; a numerical cutoff is not yet validated. Record pass, fail or indeterminate and reason codes separately for each visit pair, without access to tumor labels or prediction outcomes. Predefine adjudication for disagreement and unavailable reviewers.
4. Retain the full review flow and count failures and indeterminate cases. A gross-alignment pass does not certify tumor-boundary registration. Unreviewed cases do not pass by default. If the intended review staff or a defensible rubric is unavailable, report the gate unresolved.

## Separate support decision

Determine annotation coverage, valid negative labels, mask provenance and the relationship of binary GTV to the original Track B target independently of MRI registration. A brain mask alone cannot establish those properties. Only after the reviewed registration and support gates and any necessary protocol amendment are complete can a separate decision authorize target access. The separate CFB annotation-change endpoint is not the Track B PDE growth outcome.
