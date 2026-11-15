# Track B new question: 1-year progression from one post-treatment MRI — locked protocol (2026-11-14)

Written before any model in this file is fit or scored. This is a NEW Track B question (user request "find something else"); the volume-forecast question stays closed as negative (`paper/trackB_mu_method_tries_log_2026-11-14.md`).

## Question
From one post-treatment MRI scan plus basic clinical facts, predict whether the tumor progresses within 365 days after that scan.

## Data
- Development: MU-Glioma-Post input scans of the primary pairs (unique patient + timepoint), 5 patient folds (`data/manifests/split_mu.csv`).
  - Event: first-progression day - scan day in (0, 365]. Non-event: progression > 365 days later, or no progression and (death day or last scan day) - scan day > 365. Otherwise excluded. Scans on/after first progression excluded.
- Test (run once): UCSF-ALPTDG scan 1. Event: "Days from 1st scan to next progression" <= 365. Non-event: that value > 365, or no progression and death > 365 days after scan 1. Otherwise excluded (follow-up unknown).
  - Primary: patients not in the "UCSF-PDGM Mapping" sheet / H-3 list. Secondary: all eligible. The progression outcome was never used before in this project.

## Features (same definitions both cohorts; input scan only)
log core volume (labels 1+3), log enhancing volume (label 3), v5 composition (log_edema, nec_frac, enh_frac, log_whole, compact), age, GBM flag, log(days since diagnosis/first surgery + 1), on-treatment at scan (MU: inside RT [start,end] or additional-therapy window; UCSF: "On tx at 1st scan" not "n").

## Baselines (all existing kinds)
B1 log core volume alone; B2 log enhancing volume alone (RANO-style size); B3 clinical logistic (age, GBM, log days since dx, on-tx); B4 clinical + volume logistic; B5 composition logistic (all features, L2, C=1); B6 gradient boosting (HistGradientBoostingClassifier defaults, random_state 0) on all features.

## New hybrid (fixed now)
Mechanistic growth-to-threshold features joined with a classifier:
- From MU training pairs only: phase-specific median log-growth rates (as v1 TCKS: on-treatment vs off) and a size-dependent Gompertz curve.
- For each scan: projected log-volume change over 365 days under the Gompertz + on-treatment rate model, and projected days to 40% volume growth (capped at 730).
- Classifier: L2 logistic (C=1) on standardized [all baseline features + the 2 mechanistic features + their interaction with log core volume].

## Score and success
- Score: AUC. Dev: out-of-fold AUC on MU. Comparator = best baseline by MU dev AUC (chosen before the test).
- Test: all models fit once on all MU dev rows, applied to UCSF.
- Success: lower bound of the 95% patient-bootstrap CI of (hybrid AUC - comparator AUC) on UCSF primary > 0. Bootstrap 5000, seed 20261114.

## Weak points (stated in advance)
- Different label sources (MU automatic masks, UCSF expert masks).
- Progression definitions differ between centres; follow-up for non-progressors is partly unknown (exclusions).
- UCSF patients are mixed grades; MU mostly GBM.
- I already saw UCSF time-2 volumes (not progression); the progression column was not examined beyond counting non-missing values.
