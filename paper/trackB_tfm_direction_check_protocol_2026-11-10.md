# Track B: frozen protocol for the LUMIERE `.tfm` direction check (2026-11-10, frozen before any mask read)

**Purpose.** Decide which way the publisher's `DeepBraTumIA-segmentation/native/transformation/ct1.tfm` maps points, using non-tumor brain masks only. No eligible-patient count, CT1 intensity, or tumor label is used. This protocol was committed before any mask voxel was read.

## Inputs (84 scans = 28 planned triplets x 3 visits; all present, checked by file name)

- Native mask: `DeepBraTumIA-segmentation/native/skull_strip/ct1_brain_mask.nii.gz` (voxels > 0), with its own header affine.
- Atlas mask: `DeepBraTumIA-segmentation/atlas/skull_strip/brain_mask.nii.gz` (voxels > 0), on the 182 x 218 x 182 atlas grid.
- `ct1.tfm` and header affines as parsed by `src/fk_lumiere_transform_audit.py` (ITK LPS converted to RAS with diag(-1,-1,1)).
- Never read: CT1/T1/T2/FLAIR voxels, any `seg_mask`, any HD-GLIO segmentation.

## Method (identical for both candidates)

Let `T` be the stored transform in RAS and `A`, `N` the atlas and native affines. Two candidate atlas-to-native point maps:

- Candidate F ("stored T maps atlas to native"): `M_F = T`.
- Candidate I ("stored T maps native to atlas"): `M_I = inv(T)`.

For each candidate, the native mask is resampled onto the atlas grid with one output-to-input voxel map `inv(N) @ M @ A`, using `scipy.ndimage.affine_transform`, order 0 (nearest neighbour), constant 0 outside. The result is compared with the atlas mask by Dice.

## Decision rule (frozen)

Direction is **resolved** in favour of candidate X only if all three hold:

1. Median Dice of X across the 84 scans is at least 0.80.
2. X has the higher Dice in at least 80 of 84 scans.
3. The median of (Dice_X - Dice_other) is at least 0.10.

Otherwise the direction stays **unresolved** and the coverage count stays blocked. No threshold, input or method is changed after reading. Per-scan Dice values for both candidates are reported unchanged.

## Limits

- Brain masks come from HD-BET on different spaces; they will not match perfectly even with the right transform.
- One direction applies to all 84 scans (one tool, one release); per-scan direction is not chosen.
- The check validates the transform direction for the CT1 brain mask only, not for the other sequences.
