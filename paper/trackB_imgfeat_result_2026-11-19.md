# Track B: image-feature inputs, development result (2026-11-19)

Pre-registration (committed first, `c202b6e`): `paper/trackB_imgfeat_prereg_2026-11-19.md`. Scripts `src/track_b_imgfeat_build.py`, `src/track_b_imgfeat_dev.py`; numbers from `output/track_b_imgfeat_dev_20261119/imgfeat_result.json`. All 1,271 pairs matched a scan (no missing features except LUMIERE cavity, by design). Development only. Six cohorts, leave one cohort out.

## Result: **no candidate passes. Nothing is frozen.**

Difference to persistence (patient-level MAE, negative = better):

| Method | MU | UCSD | UCSF | LUMIERE | CFB | WAW | Mean (6) | Mean (4 earlier) | Better in |
|---|---|---|---|---|---|---|---|---|---|
| VIG (reference) | -0.034 | -0.021 | -0.024 | -0.022 | +0.010 | -0.043 | -0.0224 | -0.0253 | 5 of 6 |
| GBR-std (existing method) | -0.002 | +0.042 | +0.243 | -0.035 | +0.068 | -0.044 | +0.0452 | +0.0620 | 3 of 6 |
| GBR-img (existing method + image features) | -0.027 | +0.110 | +0.214 | -0.036 | +0.022 | -0.006 | +0.0464 | +0.0655 | 3 of 6 |
| VF1 (VIG + boosted image residual) | -0.038 | -0.023 | +0.005 | -0.042 | +0.023 | -0.075 | -0.0252 | -0.0247 | 4 of 6 |
| VF2 (VIG + linear image residual) | -0.046 | -0.028 | +0.071 | -0.055 | +0.018 | -0.110 | -0.0249 | -0.0145 | 4 of 6 |

Why they fail the rule: the six-cohort mean of VF1 and VF2 is only 0.003 lower than VIG's (needed 0.010); neither beats persistence in 5 of 6 cohorts; VF2 is worse than persistence on UCSF by +0.071 (limit +0.03).

## What this tells us
- The four mask-derived features (edema extent, necrosis share, enhancing share, cavity size) add **no measurable forecasting information** on top of VIG in these cohorts. The gain over VIG is 0.003.
- Plain gradient boosting is bad here: worse than persistence on average, with or without image features (+0.045 / +0.046), and very bad on UCSF (+0.24 / +0.21). So "learn it from the image numbers" is not a way out with 1,271 pairs.
- VIG stays the most stable method (5 of 6 cohorts, worst +0.010), with a small gain (-0.022).

## Weakest points
- Features are simple volume summaries. Texture, location, or raw-image deep features are not tested. A negative here does not rule them out; with six cohorts and mask styles that differ by cohort, they would be hard to test cleanly.
- Mask style is a cohort signature (expert-edited, FeTS-style, our tool).
- LUMIERE has no cavity label.
- WAW (23 patients) can swing the mean.

## What this means for the project
Three development families in a row (level-reversion, noise-aware level-reversion, image features) gave no clean gain over VIG. The data we hold is too small and too mixed to show a better method. **The bottleneck is data, not the model:** a clean, large, uniformly labelled cohort is needed. See `paper/trackB_next_data_plan_2026-11-19.md`.
