> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Project overview for an independent rating

You know nothing about this project. Please read it and rate it honestly. I want hard feedback, not praise.

## What it is
A single-author science-fair project in computational glioblastoma (GBM, a deadly brain tumour). Working title: "A Multi-Scale Computational Framework for Glioblastoma: From Single-Cell Drug Screening to Patient-Level Tumor Growth Forecasting". Inventory checked 2026-10-19: 159 numbered top-level Python files, 208 top-level source files, 217 recursively. Scope of this historical GBM overview is simulation and retrospective public-data analysis; this is not a repository-wide activity audit. Nothing was tested in a lab or on live patients. Paper status was author-reported at this snapshot; current dated status is in ../../START_HERE.md.

## Three tracks
**Track A: single cell to clinical.** Data: the multiomic-gbm single-cell atlas (140,355 cells, 15,000 used, 21 patients, three zones: tumour Core, Periphery, Healthy) and TCGA-GBM (n = 518).
- A graph neural network (C-GAT) classifies the zone of a cell: 78.7% vs scVI 73.0%, on a random cell-level split with patient overlap. Patient-ID lookup reaches 82.8%; do not interpret 78.7% as a biology result (ledger A1/2o).
- Historical latent-landscape claim is unverified here without an exact method/metric artifact; internal SPIB checks do not prove biological states.
- Age predicts survival in TCGA (HR 1.031 per year, p = 8.88e-16 (`../../output/survival_stats_summary.json`, ledger A3)). Age is already well known.
- A virtual drug screen is in silico only.

**Track B: tumour growth model and adaptive therapy.** Data: MU-Glioma-Post MRI (154 patients fitted, 133 scored), UCSF-PDGM diffusion MRI (62 patients).
- A reaction-diffusion model (anisotropic Fisher-Kolmogorov PDE) is used to forecast tumour shape at the next scan. Dice 0.259 vs 0.233 for "no change" (+0.026, 95% CI 0.008-0.044, Holm p = 0.0015, n = 133). The median gain is only +0.005. It loses on tumours that shrink. Target choice was post hoc after a volume probe; low absolute Dice and different splits/windows prevent direct comparison to later ladder results. Current paired manifest result is in ledger 2y.
- Adaptive dosing vs maximum tolerated dose, 61 simulated patients: 54/61 have less drug-resistant tumour. Adaptive uses 32.36% of the dose at assumed E_MAX ratio 1000 (`../../output/emax_sweep_summary.json`, ledger 2n). Resistance is model-defined; 42 negative-rho cases were excluded (2bi). But the final tumour is larger, and 10/61 progress earlier.
- All 10 patients with growth rate > 0.02/day progress earlier under adaptive. 0 of 51 below. The cutoff was chosen after seeing the data; controller sensitivity later gives only 7/10 at OFF=.50 (ledger 2bn), so the default 10/10 is not robust.
- Simulated anisotropic tumours are more elongated (1.21 vs 1.03). This compares a simulation with a simulation.
- No DTI forecast gain was observed in the older -0.001 Dice comparison (p=.97; ledger B2); exact solver 2ak is a separate later contrast.

**Track C: digital twin and reinforcement learning (RL).** Data: real fitted growth rates, assumed drug kill rates, three synthetic patient sets.
- Oracle-conditioned PPO wins on day-90 volume in 20/20 simulated cohort64 cases (100%; `../../output/rl_kill_conditioned/evaluate.json:sets.cohort64.arms.ppo_cond.win_rate_vs_stupp_pct`). Blind ppo66 wins only 25% there; 20/21 on fitted real_test growth rates still assumes kill rates. But a one-line rule does as well, and the policy never improves time to progression.
- A simple "paced" dosing rule extends time to progression: +48 to +53 days on synthetic sets, +8 to +9 days on 91 real MU patients, +16 to +19 days on 34-38 LUMIERE patients, with kill rates drawn at random. These are assumed-parameter simulations, not measured gains. The quoted ranges and all-positive CI assertion remain unverified without each exact arm/cohort/CI key in `../../output/probe_paced_policies/evaluate.json`; do not pool them.
- "Probe-then-commit" learns kill rates from the first days of data: 100% in the stated noise-free model; about 62% is the cohort64 rule arm at sigma=.10, three-seed mean (45/75/65%). Other noisy cohorts differ; no universal “realistic-noise” performance (ledger C2).

## What I found wrong with my own work (checked this week)
1. Track A's accuracy used a random split by cell, so the same patient is in train and test. A model that only knows the patient ID scores 82.8%, above the C-GAT. With a split by patient, simple models score 60-66%. Completed five-fold means: clean PCA graph 63.8% (SD13.3), cVAE/no-label-edge 74.2%, published-edge 75.8%; only the first avoids both named leaks. Donor/composition confounding remains (ledger 2o, cgat_leak_free.json).
2. Track A does not feed Track B. The score was a default 1.0; IDs lacked a measured linkage (ledger 2ap). Zone CSVs are likely synthetic with unverified exact lineage (ledger 2m); this does not establish measured expression-to-MRI biology.
3. In Track C, the "real patient" set uses real growth rates but the same assumed drug kill rates for everyone.
4. The second-drug explanation is supported only by the distinct modeled arm comparison `../../output/dual_mtd_baseline.json`, ledger 2i; it does not establish clinical rescue.
5. A growth-only forecast loses to "no change" on LUMIERE too: HD-GLIO-AUTO n=40 and DeepBraTumIA n=50 are different tool-specific volume endpoints (`../../output/lumiere_volume_forecast.json`, ledger 2l), not H-2 recurrence-test counts.
6. CGGA check of the TCGA age result (n = 374): the effect is much weaker (HR 1.009, Holm p = 0.082). Sex/MGMT intervals are inconclusive. IDH-wild-type primary subgroup age is an exception (HR1.018, Holm p=.016; ledger 2k).
7. The MGMT test is underpowered. It could only detect HR 1.74 or more.

## Please answer
1. Rate the project 1-10 for: novelty, scientific rigour, validity of the claims, honesty about limits, clarity of the story. Give one sentence for each.
2. Give an overall rating as a science-fair entry (high school level, regional to international).
3. Which finding is most publishable? Which is weakest?
4. Does the title overpromise? Suggest a more honest one if so.
5. What are the three most important things to fix or add in the next two weeks?
6. What would a skeptical judge ask first?
