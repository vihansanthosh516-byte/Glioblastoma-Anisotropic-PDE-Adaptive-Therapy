# New project options (research of 2026-10-18)

Goal set by the user: a completely new, strong project in computational biology with substantial mathematics, for the Greater Austin regional fair (GARSEF; February 2027). Old GBM work may be discarded. Nothing here is a promise of a win.

How this was checked: about 17 web searches and 4 page fetches, snippet level unless stated. "Not found" means not found in those searches, not proof that it does not exist.

## Provisional fair and ISEF notes
These snippet-derived rules, scoring weights and award listings lack a retained authoritative page/version check here. Treat them as unverified, especially all 2027 details.

- GARSEF is ISEF-affiliated. Senior = grades 9-12. For 2027 Computational Biology (CBIO) and Biomedical Engineering (ENBM) were added. 1st and 2nd place go to the Texas fair; six Senior Best of Fair go to ISEF. Source: https://austinscience.org/secondary (via search snippet). The 2027 handbook PDF could not be read: category list, dates and forms are UNVERIFIED. Check them: https://austinscience.org/images/pdfs/JSHandbook.pdf
- ISEF judging (science): research question 10, design 15, execution 20, creativity and impact 20, presentation 35 (interview 25, includes "interpretation and limits of results"). Source: https://www.societyforscience.org/isef/grand-award/criteria/
- ISEF 2025 Mathematics winners (titles only, https://www.societyforscience.org/press-release/regeneron-isef-2025-full-awards/): 1st Prisoner's Dilemma strategy dynamics; 2nd quantum-sound factorization, tensor decomposition for PDEs; 3rd coded federated learning rates, Make24 bounds, model compression theory; 4th p-adic generating functions, X-ray image statistics, Shapley value for win ratio, topological data analysis. Interpretation of these selected titles: bounds, theorems and dynamics appear prominent. Titles alone do not establish that data-fitting is absent.
- A brain tumour recurrence project (RecuNet) placed 4th in Biomedical Engineering at ISEF 2025. Not read. Brain tumour recurrence and growth are already a student topic.

## What the search found about the old project (and why it is not enough)
- The recurrence-map idea is not new (PREDICT-GBM, a 2025 Nature-family paper, student projects).
- Growth models fail on real data in published work too (rat C6 study; DSC falling from 0.85 to 0.46 as the interval grows in a diffusion-model preprint, https://arxiv.org/pdf/2509.10824). I found NO paper that compares tumour-growth forecasts with a no-change baseline or reports a forecast skill score.
- "Limits of predictability" paper (https://pmc.ncbi.nlm.nih.gov/articles/PMC12564088/, read): simulation on virtual tumours from one rat, empirical thresholds, no theorem, no competing baseline. So a theorem-level limit is not covered by it.
- I found no source that derives an expected Dice under boundary displacement or noise, or connects it to a persistence forecast.
- Monotone finite-difference schemes for anisotropic diffusion exist in numerical analysis (https://arxiv.org/pdf/1503.08177). Glioma papers found do not report positivity. Application-level gap only.

## Option 1 (recommended): A theory of forecast skill for tumour growth from segmentations
Question: when can ANY growth forecast beat "no change", given that tumour masks are measured with error?
Mathematics (core; items 2–5 are unproved project hypotheses, not established bounds or results):
1. Exact identity: Dice(A,B) = 1 - |A symmetric-difference B| / (|A|+|B|). For positive-volume sets with a smooth regular boundary and a small normal displacement field u, the first-order approximation is: 1 - Dice ~ (mean |u|) x S / (2V) (S surface area, V volume; a sphere gives 3/(2R) x mean |u|).
2. Measurement noise on both scans gives a Dice ceiling for every forecast, however good (folded-normal moments).
3. Skill score SS = (Dice_model - Dice_persistence) / (1 - Dice_persistence), as in weather forecasting. Closed form for the best possible skill SS_max as a function of displacement and noise.
4. Fisher-KPP front speed c = 2 sqrt(D rho) turns displacement into days: a closed-form forecastability horizon, and a per-patient Forecastability Index = expected displacement / noise.
5. Error budget: split forecast error into noise, speed error, direction error; prove when persistence must win.
Computation (checks the theory): simulate Fisher-KPP on brain geometry with the existing monotone solver you already built (tests numerical consequences; it cannot prove an unproved theorem); measure noise floor and displacement on real pairs (MU 334, LUMIERE, UCSF); report what share of real intervals are forecastable at all.
Why it fits: real mathematics (geometry, probability, PDE asymptotics), data present locally; CPU feasibility and GARSEF CBIO/Math eligibility remain to be checked, uses your measured segmentation noise floor.
Novelty status: no prior work found in about 6 searches. UNVERIFIED. Could be wrong.
Risks: the theory may turn out to be known in forecasting or image analysis; real intervals may be mostly forecastable (then the finding is the opposite and still valid); the theorems must be correct and proved cleanly.

## Option 2 (component, not a project alone): rigorous monotone solver for anisotropic Fisher-KPP
Prove the discrete maximum principle for the Selling-based scheme; show the standard central scheme with clamping creates spurious anisotropy effects. Strong math, low novelty (application of known numerics). Proposed for use inside Option 1, after the stated proof and numerical checks.

## Option 3: start in a different area (not recommended now)
Examples that win in Math: evolutionary game dynamics, hypergraph epidemic thresholds. I found no verified open question in these that you could finish by February, and they would not use your data or code. Needs a faculty mentor to pick a real open problem. Not researched deeply.

## Open items (UNVERIFIED)
- GARSEF 2027 categories, dates, forms and whether the project must be in Math or CBIO.
- Whether Option 1's results are already known.
- Whether any UT Austin faculty could mentor (optional, may help the interview).

## Sources used
https://www.societyforscience.org/press-release/regeneron-isef-2025-full-awards/ ; https://www.societyforscience.org/isef/grand-award/criteria/ ; https://austinscience.org/secondary ; https://pmc.ncbi.nlm.nih.gov/articles/PMC12564088/ ; https://arxiv.org/pdf/2509.10824 ; https://arxiv.org/pdf/1503.08177 ; https://arxiv.org/pdf/2206.09516 (open problems survey, not read in full)
