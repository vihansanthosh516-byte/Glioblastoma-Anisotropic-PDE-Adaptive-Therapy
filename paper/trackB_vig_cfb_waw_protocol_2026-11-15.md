# Track B: VIG hybrid on CFB-GBM + WAW — protocol, locked before any outcome is opened (2026-11-15)

Question, target and score are unchanged: from a post-treatment scan plus scan dates, forecast y = log(core volume in mm3 + 1) at the next scan; core = necrosis label + enhancing label; score = mean absolute error, averaged per patient first. (1 mm3 = 1 voxel in the earlier cohorts.)

This is the **second** outside run of VIG. The first (LUMIERE) did not pass (`output/track_b_vig_lumiere_20261114/result.json`).

## 1. Status of the data at the lock

- Spent, used for fitting only: MU (335 pairs), UCSD (55), UCSF (186), LUMIERE (400); 976 pairs, 442 patients (`output/track_b_code_audit_20261114/fix_and_loco4_stdout.txt`).
- Test set: CFB-GBM v4 (TCIA, CC BY 4.0, on disk) and WAW (Zenodo 10.5281/zenodo.22205631, CC BY 4.0).
- Seen before the lock, CFB: availability tables, scan-interval quantiles, headers, brain masks, T1Gd intensities at t1 for registration checks. **Zero GTV voxels. No tumor label. No volume.** (`START_HERE.md`, CFB section; `output/track_b_cohort_screen_20261115/screen_stdout.txt`.)
- Seen before the lock, WAW: four metadata tables (scan days, grade, file list). **No mask, no image, no radiomics table.** The mask archive is not on disk.
- Upper limits before the "input core > 0" rule: CFB 158 patients (155 t0-to-t1 pairs, 107 t1-to-t2 pairs), WAW 28 patients (48 pairs) (`screen_stdout.txt`).

## 2. Method (unchanged)

VIG as in `paper/trackB_vig_lumiere_protocol_2026-11-14.md`, section "Method": `y_hat = y_in + g(dt) - h(c)`, code `fit_vig` and `vig_curve` in `src/track_b_vig_lumiere.py`, not edited.

One change, fixed now: the parameters are fitted on all four spent cohorts (976 pairs) and not on three. Reason: LUMIERE is spent and adds pairs with an earlier scan. The pair builders are those of `output/track_b_code_audit_20261114/fix_and_loco4.py` (with the two data fixes).

## 3. Lock chain (order is part of the protocol)

| Step | Action | What it may read | Output, committed before the next step |
|---|---|---|---|
| 1 | Commit this protocol | — | this file |
| 2 | Write `src/track_b_vig_cfb_waw.py`. Fit VIG and all comparators on the 976 pairs. | Spent cohorts only | `frozen_params.json` and its SHA-256 |
| 3 | Train the segmentation model (section 4). Run gate G1. | Spent cohorts only | weights SHA-256, `gate_g1.json` |
| 4 | Run gate G2. | CFB t0 images and t0 human GTV only | `gate_g2.json` |
| 5 | Label all eligible CFB scans. Download WAW `3_2_postop_masks.zip` (2.6 MB; user approval needed). Build the pair table. | CFB images; WAW masks | sealed `test_pairs.csv`; nobody prints or opens it |
| 6 | Inventory mode: prints counts of eligible pairs and patients per cohort. Nothing else. | sealed table | `inventory.json` |
| 7 | One scoring run. | sealed table, frozen parameters | `result.json`, `predictions.csv` |

Stop rules. If G1 or G2 fails: no CFB label is made, WAW masks are not downloaded, no test runs, both cohorts stay unspent. If step 6 shows fewer than 100 eligible patients in total: no scoring run; both cohorts stay unspent for this method only if step 7 never ran. The script refuses to run step 7 without steps 2 to 6 on disk and refuses to overwrite `result.json`.

## 4. Core labels for CFB (the target mask is made by one frozen tool)

Reason: CFB gives one binary GTV label. GTV is the radiotherapy target, not necrosis + enhancing. Using GTV would change the target. Standard four-sequence tools (FeTS / FL-PoST) cannot run: T2 exists in 144 of 527 CFB scans (`screen_stdout.txt`).

- Model: nnU-Net v2, `3d_fullres`, default planning, two input channels (T1Gd, FLAIR). Labels: 1 necrosis / non-enhancing core, 2 FLAIR abnormality, 3 enhancing, 4 resection cavity.
- Training data: every MU and UCSF-ALPTDG scan with a mask (spent cohorts). UCSD is held out for gate G1.
- Training: one model on all training scans (fold `all`), 250 epochs, seed fixed, last checkpoint. No model selection. Training may resume from a checkpoint across GPU sessions. The nnU-Net version is pinned and written to `gate_g1.json`.
- Assumptions to check before training (a failed check is recorded as a deviation before step 3): label numbers are identical in MU, UCSF and UCSD masks (the loaders in `src/track_b_v8_*.py` treat them so); a 250-epoch trainer exists in the pinned version, else a subclass that only sets the epoch count.
- CFB inference: T1Gd and FLAIR of the same visit, multiplied by the released brain mask, at native grid. Volume = label voxels x voxel volume in mm3. Scans without both sequences are not eligible.

**Gate G1 (cross-site, spent data only).** On all UCSD scans with a released mask, d = |log(tool core + 1) - log(released core + 1)|. Pass: median d <= 0.20 AND mean d <= 0.50. (Assumed thresholds, fixed before training. The mean limit is under half of the pooled persistence error 1.148, `loco4_all_baselines_stdout.txt`.)

**Gate G2 (CFB images, input side only).** On t0 scans with a human GTV (191 files) and tool core > 0: the share of tool-core volume that lies inside the GTV dilated by 5 mm. Pass: median share >= 0.70. (Assumed threshold. Assumption: the t0 GTV covers residual enhancing tumor and the cavity.) The t0 GTV is never a target in this protocol.

## 5. Test pairs

- CFB: visits t0, t1, t2. Time = `time_diff_t0 (weeks)` x 7. t0 is after surgery or biopsy and before radiotherapy; it is taken as a post-treatment scan (surgery is a treatment; LUMIERE had 99 such pairs).
- WAW: post-operative exams only. Time = `days_postop`. Core from the released masks (labels 1 + 3). Pre-operative exams are never used. All grades (the development cohorts also mix grades).
- Pairs: consecutive scans of one patient. Eligible: interval 14 to 365 days and input core > 0.
- Last change `c`: from the earlier post-operative scan if its core > 0, else 0.
- No patient is removed after step 6 for any reason.

Details fixed while writing the script (2026-11-15, step 2, before any outcome existed; input side only):
- CFB pairs are between adjacent visits only (t0-to-t1, t1-to-t2). A visit without both sequences breaks the chain.
- WAW exams with no `days_postop` value (12 of 160 in the metadata table) are dropped before pairing.
- Builders were tested with random fake volumes and a fake mask archive only.

## 6. Comparators (all fitted in step 2 on the 976 pairs, then frozen)

Persistence, two-point exponential, population rate, empirical-Bayes random slope, Gompertz, gradient boosting, Gompertz with absolute loss, gradient boosting with absolute loss. Definitions as in `output/track_b_code_audit_20261114/loco4_all_baselines.py`. Gradient boosting uses the inputs present in all cohorts (input size, interval, earlier rate, has-earlier-scan).

## 7. Primary success rule

For EACH of the eight comparators: the upper limit of the two-sided **97.5%** CI of (VIG error minus comparator error) is below 0. Patient-level paired bootstrap, resampled inside each cohort (CFB, WAW), 50,000 resamples, seed 20261115. All eight must pass.

Multiple testing: "all eight must pass" needs no split. This is the second outside run of VIG, so the 5% level is halved (rule set in `paper/trackB_problem_review_2026-11-14.md`, line 92).

Reported with it, no claim: the 95% CI; paired sign-flip test; each cohort alone; CFB t0-to-t1 and t1-to-t2 apart; WAW grade 4 only; VIG without the history term; error by interval bin; tool-noise note from G1. After the primary run only: the tool run on WAW images against the released WAW masks (needs `2_2_postop_scans_coreg.zip`, 3.3 GB, user approval; if not approved it is recorded as not done).

## 8. Power (stated before the run)

`output/track_b_vig_cfb_waw_20261115/power_stdout.txt`, against persistence only, assuming effect 0.027 and SD 0.154: at the 2.5% level, power is 0.56 with 186 patients, 0.49 with 158, 0.31 with 100. 80% power needs 310 patients. The "input core > 0" rule will lower the count. The eight-comparator rule lowers power more.

Reading of the outcome, fixed now: pass = "confirmed on one pooled outside set, second attempt, tool-made targets for CFB". Fail = "not confirmed; the test was underpowered", not "disproved".

## 9. Weakest points

- Power is about 0.5 or less. A fail is the more likely outcome even if the exploratory effect is real.
- CFB targets come from our own two-sequence tool. Tool error adds noise to every method. Noise in the input size also makes "pull back toward the mean" methods look better than persistence. VIG has such a term; so do Gompertz, random slope and gradient boosting. A win over persistence alone would be weak; the rule needs all eight.
- Two mask sources in one pooled test (tool for CFB, released masks for WAW). Training cohorts use released masks.
- Gate thresholds (0.20, 0.50, 0.70, 100 patients) are judgement values, not derived.
- Median intervals are near 84 to 98 days. The interval part of VIG is near 0 there (curve in the LUMIERE protocol). The run mostly tests the drift and history parts, as on LUMIERE.
- CFB t0 is before radiotherapy; radiotherapy lies inside every t0-to-t1 pair. CFB times are whole weeks.
- The effect size used for power was chosen post hoc on spent data (`paper/trackB_code_audit_2026-11-14.md`).
- VIG is a combination of known parts; the scan-interval effect is published (`paper/trackB_vig_prior_art_2026-11-14.md`). A pass would confirm a small forecasting gain, not new biology.
- Training cost is an estimate: one nnU-Net model, 250 epochs, more than one free GPU session (assumption).
