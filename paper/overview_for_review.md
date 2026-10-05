# Project overview for an independent rating

You know nothing about this project. Please read it and rate it honestly. I want hard feedback, not praise.

## What it is
A single-author science-fair project in computational glioblastoma (GBM, a deadly brain tumour). Working title: "A Multi-Scale Computational Framework for Glioblastoma: From Single-Cell Drug Screening to Patient-Level Tumor Growth Forecasting". About 97 numbered Python scripts. Everything is simulation or retrospective analysis of public data. Nothing was tested in a lab or on live patients. The paper is not written yet.

## Three tracks
**Track A: single cell to clinical.** Data: the multiomic-gbm single-cell atlas (140,355 cells, 15,000 used, 21 patients, three zones: tumour Core, Periphery, Healthy) and TCGA-GBM (n = 518).
- A graph neural network (C-GAT) classifies the zone of a cell: 78.7% vs scVI 73.0%, on a random cell-level split.
- An energy-landscape model found two metastable states.
- Age predicts survival in TCGA (HR 1.031 per year, p = 9e-16). Age is already well known.
- A virtual drug screen is in silico only.

**Track B: tumour growth model and adaptive therapy.** Data: MU-Glioma-Post MRI (154 patients fitted, 133 scored), UCSF-PDGM diffusion MRI (62 patients).
- A reaction-diffusion model (anisotropic Fisher-Kolmogorov PDE) is used to forecast tumour shape at the next scan. Dice 0.259 vs 0.233 for "no change" (+0.026, 95% CI 0.008-0.044, Holm p = 0.0015, n = 133). The median gain is only +0.005. It loses on tumours that shrink.
- Adaptive dosing vs maximum tolerated dose, 61 simulated patients: 54/61 have less drug-resistant tumour. Adaptive uses 32% of the dose. But the final tumour is larger, and 10/61 progress earlier.
- All 10 patients with growth rate > 0.02/day progress earlier under adaptive. 0 of 51 below. The cutoff was chosen after seeing the data.
- Simulated anisotropic tumours are more elongated (1.21 vs 1.03). This compares a simulation with a simulation.
- DTI fibre orientation adds nothing to the forecast.

**Track C: digital twin and reinforcement learning (RL).** Data: real fitted growth rates, assumed drug kill rates, three synthetic patient sets.
- An RL policy that knows the kill rates beats the standard schedule on day-90 tumour volume (100%). But a one-line rule does as well, and the policy never improves time to progression.
- A simple "paced" dosing rule extends time to progression: +48 to +53 days on synthetic sets, +8 to +9 days on 91 real MU patients, +16 to +19 days on 34-38 LUMIERE patients, with kill rates drawn at random. All CIs are above 0.
- "Probe-then-commit" learns kill rates from the first days of data: 100% with no noise, 62% with realistic noise.

## What I found wrong with my own work (checked this week)
1. Track A's accuracy used a random split by cell, so the same patient is in train and test. A model that only knows the patient ID scores 82.8%, above the C-GAT. With a split by patient, simple models score 60-66%. A leak-free C-GAT is still running (first two folds average 68-74% for clean variants; the published setup gives about 82%).
2. Track A does not feed Track B. The link used a score that is 1.0 for every patient (IDs do not match), and the file behind it is probably synthetic.
3. In Track C, the "real patient" set uses real growth rates but the same assumed drug kill rates for everyone.
4. A "dual-drug rescue" result came from the second drug, not from adaptive control.
5. A growth-only forecast loses to "no change" on new LUMIERE patients too (n = 40 and 50).
6. CGGA check of the TCGA age result (n = 374): the effect is much weaker (HR 1.009, Holm p = 0.08). Sex and MGMT are null.
7. The MGMT test is underpowered. It could only detect HR 1.74 or more.

## Please answer
1. Rate the project 1-10 for: novelty, scientific rigour, validity of the claims, honesty about limits, clarity of the story. Give one sentence for each.
2. Give an overall rating as a science-fair entry (high school level, regional to international).
3. Which finding is most publishable? Which is weakest?
4. Does the title overpromise? Suggest a more honest one if so.
5. What are the three most important things to fix or add in the next two weeks?
6. What would a skeptical judge ask first?
