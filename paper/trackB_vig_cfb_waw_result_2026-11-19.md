# Track B: VIG hybrid on CFB-GBM + WAW — result of the one locked run (2026-11-19)

Protocol: `paper/trackB_vig_cfb_waw_protocol_2026-11-15.md`. Every number below is from `output/track_b_vig_cfb_waw_20261115/result.json` (SHA-256 of the sealed pair table 2f4ef74a…; seed 20261115; 50,000 patient-level resamples). Run once. Nothing was changed after it.

## Plain result
**Primary rule failed.** VIG did not beat all eight comparators. It beat four and did not beat four.

| Comparator | Patient MAE | VIG minus comparator (mean) | 97.5% CI | Passes (upper limit < 0)? |
|---|---|---|---|---|
| persistence | 1.1538 | -0.0029 | [-0.051, 0.034] | **No** |
| population rate | 1.1542 | -0.0034 | [-0.052, 0.035] | **No** |
| empirical-Bayes random slope | 1.2247 | -0.0738 | [-0.213, 0.037] | **No** |
| Gompertz, absolute loss | 1.1190 | +0.0319 | [-0.015, 0.078] | **No** (VIG is worse) |
| two-point exponential | 3.4735 | -2.32 | [-6.93, -0.080] | Yes |
| Gompertz | 1.2285 | -0.0777 | [-0.148, -0.007] | Yes |
| gradient boosting | 1.4651 | -0.3143 | [-0.446, -0.184] | Yes |
| gradient boosting, absolute loss | 1.2285 | -0.0777 | [-0.154, -0.0003] | Yes |

VIG patient MAE 1.1508 (181 patients, 295 pairs). VIG without the history term: 1.1484. Best single method on this test: Gompertz with absolute loss, 1.1190.

## The four report questions
1. **Positive?** No. VIG is not better than persistence (difference -0.003, CI includes 0, sign-flip p = 0.92).
2. **How much?** About 0.3% lower error than persistence. That is no measurable gain. The exploratory effect (0.027) did not reproduce.
3. **More unique?** The method is still a new combination (`paper/trackB_vig_prior_art_2026-11-14.md`). On this test it is not better than simple baselines, so it adds no usable forecasting gain here.
4. **Fix?** Do not tune on these cohorts: both are spent. See "what to try next".

## Secondary results (reported, no claim; `secondary_no_claim` in `result.json`)
- CFB alone (158 patients, tool-made targets): persistence 0.883, VIG 0.891 (VIG minus persistence +0.009, CI [-0.016, 0.036]).
- WAW alone (23 patients, released masks): persistence 3.014, VIG 2.933 (CI [-0.43, 0.10]). Error is more than three times larger than on CFB.
- VIG minus VIG without history: +0.002, CI [-0.045, 0.037]. The history term adds nothing.

## Weakest points
- Power was about 0.5 (protocol section 8). The fail is "not confirmed, underpowered", not "disproved". But the point estimate is also near zero, so a larger sample is unlikely to turn it into a clear win.
- CFB targets come from our own tool (G1 mean error 0.25 log units, `gate_g1.json`). This noise lifts every method.
- WAW has only 23 patients and a MAE near 3.0. Its pairs are early post-operative scans, where resection changes the core. The result is dominated by this unstable subset.
- Two mask sources in one pooled test.
- This is the second outside run of VIG; the first (LUMIERE) also did not pass.

## What this means for Track B
Two locked outside tests (LUMIERE, CFB+WAW) did not confirm VIG. The exploratory gain on spent cohorts (0.027 vs persistence) did not hold on new data. The honest status: **VIG is not a confirmed better forecaster than persistence.** Both test cohorts are now spent for this method.

Options for the user (not decided here): (a) a new question with new data, since the scan-interval signal was not visible at the intervals of these cohorts; (b) report the pair of negative results as the finding (a new combination method that fails two locked tests, plus the audit trail). Closing Track B is not the default (user rule).
