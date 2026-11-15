# Track B problem review (2026-11-14) — exploratory, post hoc

Scope: next-scan tumor-core volume forecast. Question, target and score are unchanged.
All new numbers come from read-only scripts in `output/track_b_review_diag_20261114/` (`diag1`..`diag4`, each with a `_stdout.txt`). They read data and saved predictions. They fit no new method and write no result file.
Label: every number here is post hoc and descriptive. UCSD and UCSF are spent test sets; nothing here can be used as a new test.

## 0. Corrections first

- My earlier reply said "11 pairs fall to zero, 15% of error" was not found. That was wrong. It reproduces: 11 MU pairs have outcome volume 0 and hold 15.5% of the persistence error (`diag1_stdout.txt`).
- The prompt says MU masks are automatic and UCSF masks are expert. The TCIA page for MU-Glioma-Post says the automatic masks were "manually validated and refined by neuroradiologists" (https://www.cancerimagingarchive.net/collection/mu-glioma-post/). `paper/trackB_mu_tcks_protocol_2026-11-14.md` line 47 calls them automatic. So "automatic vs expert" is not an established difference. How UCSD masks were made was not checked here.
- The prompt says later surgery dates are missing. UCSF has columns for the 2nd and 3rd surgery and an on-treatment flag. UCSD has last-surgery, last-radiation and Avastin (bevacizumab) dates. MU has a "Multiple surgeries" flag (51 yes) but no later surgery date (`diag2_stdout.txt`, `diag3_stdout.txt`).

## 1. The real cause (three findings)

1. **v8 forecasts a level, not a change.** Its regression has a slope of 0.54 on input size and an intercept of 8.835 learned on MU (`diag2_stdout.txt`). On the 137 UCSF patients its forecast is +0.81 log units above the input scan on average (`diag1_stdout.txt`). The true mean change is +0.05 (7.922 to 7.969, `diag2_stdout.txt`). v8 pulls every tumor toward the MU average size. This alone can explain why v8 is worse than persistence on UCSF. Weakest point: I did not re-score v8 with the bias removed, so the size of this effect on the MAE is not measured.
2. **The signal the MU methods learned is regression to the mean, and it does not transfer.** Mean log change by input-size third (small / middle / large): MU +0.587 / -0.108 / -0.411 (`diag3_stdout.txt`); UCSF -0.180 / +0.171 / +0.152; UCSD -0.141 / -0.106 / -0.524 (`diag4_stdout.txt`). On MU, small tumors grow and large tumors shrink. On UCSF the pattern is absent. Weakest point: the cause of the MU pattern is unknown. It can be mask noise, or the MU scan schedule (MU includes scans soon after surgery; UCSF scans are a median of 496 days after surgery, `diag2_stdout.txt`).
3. **Small lesions make most of the error in every cohort.** Pairs with outcome volume 0: MU 11 pairs = 15.5% of error; UCSD 2 pairs = 26.6%; UCSF 5 pairs = 23.1% (`diag1_stdout.txt`). Pairs with input < 500 voxels: MU 41 pairs hold 21.9% of error (mean error 2.34 vs 1.17 for the rest, `diag1_stdout.txt`); UCSF mean error 1.66 vs 0.63 (all 186 eligible, `diag2_stdout.txt`). The log(voxels + 1) target makes a small lesion that disappears cost 4 to 7 log units.

## 2. Problem list

Severity: C = critical, H = high, M = medium, L = low.

| # | Sev | Problem | Evidence | Fix | Weakest point of the fix |
|---|---|---|---|---|---|
| 1 | C | Learned methods predict the level `y_out`. Persistence is not a special case. | Finding 1 | Model the change `y_out - y_in`. Add a weight that can go to 0 (= persistence). Limit the size of the correction. | A method that can fall back to persistence can tie it, not beat it. |
| 2 | C | The MU gain is regression to the mean. It is absent on UCSF. | Finding 2 | Set the shrink strength per cohort from input scans only (two masks of the same scan). No shrink if measured noise is low. | If the MU pattern is biology, a noise-based weight is wrong. |
| 3 | H | "Automatic vs expert masks" is not established. | Section 0 | Correct the project text. Measure mask noise per cohort directly (Section 3). | The public record gives no per-scan record of what the radiologists changed. |
| 4 | H | Outside success rule compared v8 with population rate only. Persistence was the best baseline on UCSD (1.342 vs 1.441) and UCSF (0.832 vs 0.851). The question says beat ALL. v8 minus persistence on UCSD is -0.123, not -0.222 (`diag2_stdout.txt`). | result.json files; `diag2_stdout.txt` | Rule: the upper CI limit must be below 0 against each of the six baselines. All six must pass. | Harder to pass. |
| 5 | H | Zero and small volumes dominate the score. | Finding 3 | Keep target and score. Add a size-dependent noise term and a "lesion disappears" part driven by input size. Report error by size third as a secondary table. | v8 collapse AUC was only 0.61 (`output/track_b_mu_v8_20261114/result.json`). |
| 6 | M | "Collapse = surgery" is not supported. MU: 5 of 22 collapse pairs have a recorded therapy or progression event between the scans; collapse rate is lower with an event (0.031 vs 0.099). UCSF: 0 of 6 zero-outcome pairs have a surgery between the scans under either sign reading; their inputs are small (log size 3.5 to 6.5). | `diag2_stdout.txt`, `diag3_stdout.txt` | Treat collapse as small-lesion drop-out, not as surgery. | MU has no later surgery dates. UCSF sign convention is unclear. |
| 7 | M | Treatment data exist but were replaced by assumed windows. UCSF flag says 77 of 186 patients were on chemo at scan 1; the 180-day window labelled 16 of those 77 as adjuvant. Bevacizumab, TTF, brachytherapy and surgery fields were not used in any cohort. | `diag3_stdout.txt` | One treatment state at the input scan, from fields all cohorts have: on chemo (yes/no), days since last radiotherapy, bevacizumab before the scan, days since last surgery. | UCSD has no chemo end date. UCSF has no bevacizumab field. |
| 8 | M | Future information in the phase feature. It uses the midpoint of the two scan dates and therapy END dates. | `src/track_b_mu_tcks.py` lines 61-65 | Define the state at `day_in` with events on or before `day_in` only. | Loses some real signal. The scan interval `dt` stays (it is the forecast horizon). |
| 9 | M | UCSF table has bad date cells. 22 of 298 rows (14 of 186 eligible) hold a day value above 15,000 (largest 328,372). The data dictionary and the column names disagree on sign. No range check in code. | `diag3_stdout.txt`, `diag2_stdout.txt` | Fixed valid range; out-of-range becomes missing; assert before the run. | The true dates stay unknown. |
| 10 | M | Baselines are fitted with squared loss and scored with absolute error (Gompertz least squares, default boosting). Candidates use median loss. | `src/track_b_mu_tcks.py` lines 92, 121 | Fit each trainable baseline with absolute loss. | Small effect: the best baseline is persistence, which has no fit. |
| 11 | M | Unpenalised median regression on collinear inputs (`y_in`, `x_hat`, `log_whole`, `log_edema`). Coefficients +1.08 and -0.82 offset each other. | `diag2_stdout.txt` | Fewer inputs. Use ratios, not levels. Penalty chosen inside training folds. | Less flexible. |
| 12 | M | Bootstrap: 5,000 resamples for a 99.4% to 99.6% interval. About 16 resamples define each tail. Percentile method only. | `src/track_b_mu_v3.py` lines 46-49 | 50,000 resamples plus a paired sign-flip test. | None of note. |
| 13 | M | Power. SD of the per-patient paired difference: MU 0.86, UCSD 1.24, UCSF 0.73. | `diag2_stdout.txt` | Section 5. | Needs more patients than exist. |
| 14 | L | Volume is in voxels; spacing is never checked. All sampled masks are 1 mm isotropic (MU and UCSF 240x240x155; UCSD 256x256x256), so no error was found. | `diag1_stdout.txt` (MU every 25th mask, UCSD every 10th, UCSF 1 patient) | Assert spacing or convert to mm3. | Only a sample was checked. |
| 15 | L | Label 4 (resection cavity) is in all three cohorts but is never an input. Cavity can be confused with necrosis (label 1), which is in the target. | `diag1_stdout.txt` | Add cavity volume at the input scan as an input. | Not tested. |
| 16 | L | Folds are patient-level (no patient in two folds). No leak found. Baseline predictions are joined by row order with a check on patient ID only. | `diag1_stdout.txt`; `src/track_b_mu_v3.py` line 55 | Join on patient and timepoint keys. | None. |
| 17 | L | Selection. A pair exists only if a next scan exists. Input must be above 0 but the outcome can be 0. 29 MU pairs with empty input are excluded. | `diag1_stdout.txt` | Report it. Not fixable inside the locked score. | No fix. |
| 18 | L | Scan intervals differ: MU median 77 days, UCSF 59 days. Baseline errors are not comparable across cohorts. | `diag1_stdout.txt`, `diag2_stdout.txt` | Report error by interval bin. | Smaller groups. |
| 19 | L | Code risks. Column lookup by first substring match. Input path in a Downloads folder. No input file hash in results. | `src/track_b_v8_ucsf_alptdg.py` lines 25, 45 | Exact column names. SHA-256 of each input in `result.json`. | None. |
| 20 | H | Reuse and selection: 17 locked tries; all cohorts spent. | `paper/trackB_mu_method_tries_log_2026-11-14.md` | A new untouched cohort only. | None is open to the user now (Section 5). |

Not changed: the locked scripts. They are records of spent tests. Fixes go into a new script under a new protocol.

## 3. Mask-source shift

- The stated shift is not proven (Section 0). The measured differences are: input size (MU mean 8.32, UCSF 7.92, UCSD 9.10 log voxels), prior-scan availability (MU 60%, UCSF 0%), and time since surgery (`diag2_stdout.txt`).
- Useful check: run ONE fixed public tool on the INPUT scans of all cohorts. Use its mask as a second reading of the same scan. The disagreement between the dataset mask and the tool mask gives a noise value per cohort and per size, with no outcome data.
- Tool: FeTS / FL-PoST, version pinned. Reason: it is public, made for post-treatment glioma, and the MU paper says it made the MU masks before refinement (https://www.nature.com/articles/s41597-025-06011-7). Assumption to verify before locking: it runs offline with fixed weights and gives labels 1 to 4. HD-GLIO is not suitable if it gives no necrosis label (assumption from memory; ledger 2aj also records that its LUMIERE maps were unusable).
- Keep the dataset mask as the target. Do not replace the target with the tool mask on cohorts that have a dataset mask.
- "Untouched" status: re-segmentation does not restore it. The outcomes of these patients were already seen. MU, UCSD and UCSF stay development data.
- Weakest point: the tool was possibly trained on scans from these sites. It will then agree too well there and the noise value will be too low. Ledger 2ai shows this effect: AI mask vs radiologist-corrected mask has median Dice 1.000.

## 4. One hybrid method: anchored two-reader shrinkage with Gompertz drift

Families joined: measurement-error (empirical-Bayes) model + mechanistic Gompertz growth + a disappear (hurdle) part.

Inputs for one pair: `m1` = log(core + 1) from the dataset mask of the input scan; `m2` = the same from the fixed tool on the same scan; `dt`; treatment state `g` at the input scan (problem 7 definition); cohort `c`.

1. Noise, from inputs only. For cohort `c` and size third `b`: `s2[c,b] = (1.4826 * MAD(m1 - m2))^2 / 2`. Assumption: the two readings have equal, independent error.
2. Cohort prior, from inputs only: `mu[c] = median(m1)`, `tau2[c] = max(Var(m1) - mean(s2), 1e-6)`.
3. Latent size: `k = s2 / (s2 + tau2)`; `x = m1 - k * (m1 - mu[c])`.
4. Gompertz drift in log volume: `x_next = K - (K - x) * exp(-a[g] * dt)`, `a[g] >= 0`. `K` and `a[g]` are fitted on development data with absolute loss.
5. Disappear part: `p0 = logistic(b0 + b1*x + b2*log(dt) + b3*sqrt(s2))`. If `p0 > 0.5` the model forecast is 0, else `x_next`.
6. Anchor: `forecast = m1 + L * (model forecast - m1)`, with `L` in {0, 0.25, 0.5, 0.75, 1}. `L` is chosen by leave-one-cohort-out to minimise the worst cohort's (method MAE minus persistence MAE). If no `L > 0` gives a worst case below 0, then `L = 0`, the method is declared failed in development, and no outside test is spent.

Why it can beat persistence and the population rate: both ignore measurement error in the input scan. Where measured noise is high, shrinking the input toward the cohort centre removes part of that error. The cohort centre comes from the TEST cohort's own input scans, so the +0.81 offset of v8 cannot occur.

Direct judgement: I rate an outside pass as unlikely. On UCSF there is no regression to the mean to use (finding 2), so the method becomes persistence there and ties. The largest gain seen on MU in 14 tries is 0.079. Novelty is not verified: each part is standard, the noise-aware latent size was already tried with a history prior (try 2, MAE 1.555), and no prior-art search was done for this combination.

## 5. Locked plan

- Development: MU (334 pairs, 135 patients), UCSD (46 pairs, 33 patients), UCSF (186 eligible patients). All exploratory. Leave-one-cohort-out.
- Gate 0, before any modelling: the two-reader check of Section 3. Written prediction: measured noise is higher on MU than on UCSF. If it is not, noise is not the cause and the method stops.
- Outside test (untouched): Burdenko-GBM-Progression, TCIA. 180 patients, 645 studies, 1 to 8 follow-up scans each. Images are under the NIH Controlled Data Access Policy (https://www.cancerimagingarchive.net/collection/burdenko-gbm-progression/). The project record says the request needs a faculty-level sponsor through dbGaP, and that the user cannot ask anyone for data (`paper/data_request_drafts_2026-11-12.md` line 76; `paper/data_request_plan_2026-11-12.md` line 53). It has no tumor masks, so target masks must come from the fixed tool. The number of eligible pairs is unknown.
- Other candidates checked: UCSD-PTGBM held-out part, 40 patients and 60 timepoints (https://pmc.ncbi.nlm.nih.gov/articles/PMC12886825/), so at most 20 pairs, below the project minimum of 25. LUMIERE is spent. SAILOR and BrainMD need a data-use agreement (`paper/track_b_external_cohort_screen_2026-10-28.md`).
- Comparator, fixed in advance: all six baselines, each fitted on development data with absolute loss, then frozen.
- Success rule: one run. For each of the six baselines, the upper limit of the 95% CI of (method minus baseline), per-patient paired, must be below 0. All six must pass. 50,000 bootstrap resamples; a paired sign-flip test is reported with it.
- Minimum sample size (two-sided 5%, 80% power, `n = 7.85 * SD^2 / effect^2`):

| Assumed effect | Assumed SD | Patients needed |
|---|---|---|
| 0.079 (MU, post hoc, so optimistic) | 0.86 (MU) | about 930 |
| 0.123 (UCSD, v8 vs persistence) | 1.13 (UCSD) | about 670 |
| 0.222 (UCSD, v8 vs population rate; likely inflated, n = 33) | 1.24 (UCSD) | about 245 |
| 0.222 | 0.86 | about 120 |

  With 180 patients and SD 0.86, the smallest effect with 80% power is about 0.18. These are my calculations from the SDs in `diag2_stdout.txt`. They assume the effect has the same sign outside, and UCSF shows the opposite sign.
- Multiple testing: one method, one cohort, one run. "All six must pass" needs no split of the 5% level. A second outside run of any kind halves the level. No claim is made from development data.

## 6. Plain statement

No fair positive result is possible with the data that exists and is open to the user. All three usable cohorts are spent. The one untouched cohort needs a sponsor. The honest effect size (about 0.08) needs about 930 patients; no public cohort has that.

Data that makes it possible: about 250 (optimistic) to 930 (honest) untouched patients; two or more dated post-treatment scans each; T1, T1 with contrast, T2 and FLAIR, so one tool can segment all scans; treatment dates that include later surgery and bevacizumab; and repeat scans a few days apart to measure mask noise directly (QIN-GBM-Treatment-Response and RIDER Neuro MRI are candidates from memory; access not verified).
