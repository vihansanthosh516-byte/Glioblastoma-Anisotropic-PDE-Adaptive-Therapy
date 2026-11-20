# Track B: image-feature inputs (VIG-F) — written before any run (2026-11-19)

Question, target, score unchanged. Development only; no test is spent. Same six cohorts, leave-one-cohort-out, bootstrap 20,000, seed 20261119.

## Why image features, and what the earlier work says
Level-reversion variants gave no clean gain (`paper/trackB_vig2_noise_result_2026-11-19.md`). The next honest source of new information is the scan itself. Features are chosen **from the segmentation labels only, and fixed here, before any model sees them**. They do not come from the CFB/WAW result.

## Features (all at the INPUT scan; volumes in mm3 so voxel size does not matter)
- `log_edema` = log(volume of label 2 (FLAIR abnormality) + 1)
- `nec_frac` = label 1 volume / core volume
- `enh_frac` = label 3 volume / core volume
- `log_cavity` = log(volume of label 4 (resection cavity) + 1); missing for LUMIERE (no cavity label there)
Sources: MU `mu_mask_counts.csv`; UCSD and UCSF released masks; WAW released masks; CFB our frozen tool's masks; LUMIERE volumes per label from its segmentation table. Matching to pair rows is by (patient id, input core volume); rows that do not match exactly one scan keep NaN features and the count of such rows is printed.

## Candidates (all absolute loss; fitted on the five training cohorts of each fold)
- **GBR-std (reference):** gradient boosting, absolute loss, on [y_in, dt, r_hist, has_prior], target = y_out - y_in (as the comparator in the locked protocol).
- **GBR-img (reference, an existing method with image features):** same, plus the four image features.
- **VF1 (VIG + boosted residual):** `VIG(dt, history) + GB(residual)`; GB: absolute loss, max_depth 2, 100 iterations, learning rate 0.05, min_samples_leaf 40; inputs = the four image features + y_in.
- **VF2 (VIG + linear residual):** `VIG(dt, history) + median regression (L1) on the four standardized image features + y_in`; missing values set to the training median.
VIG is `fit_vig` unchanged, refitted in each fold. No other settings are tried.

## Selection rule (fixed now)
A candidate (VF1 or VF2) is selected only if ALL hold:
1. cohort-equal mean difference to persistence over six cohorts is lower than VIG's by at least 0.010, and lower than GBR-img's;
2. its mean over the four earlier cohorts is not worse than VIG's;
3. no cohort is worse than persistence by more than +0.03;
4. it beats persistence in at least 5 of 6 cohorts.
If none passes: nothing is frozen; result "no image-feature gain". Nothing is added or retuned after the run.

## What a pass would mean
Only a candidate to freeze before one future untouched cohort. Not an improvement claim.

## Weak points stated first
- Six cohorts; features come from different mask makers (expert-edited, FeTS-style, our tool). Mask style is itself a cohort signature, so a "gain" could be a cohort-label artifact.
- LUMIERE has no cavity label.
- CFB and WAW were already seen as a test; these folds are informed by that for VIG itself, not for the choice of features.
