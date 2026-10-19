> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Research prompt: methods and data to make the GBM digital twin work

Paste everything below the line into a new Claude session (deep research mode if available).

---

You are a research assistant for a high-school science fair project on glioblastoma (GBM). Find published methods, public datasets and code that could make three parts of the project work. Report only what you can cite with a working link (paper DOI/arXiv, dataset page, GitHub). Never invent a source, a number or a dataset. If you are not sure, say "not verified".

## Historical context (corrected 2026-10-19; current goal is ../../START_HERE.md)

- Data already used: PREDICT-GBM (TUM 142 dev; LUMIERE 61 + RHUH 40 former test: 101 total, 98 primary-evaluable; now SPENT), MU-Glioma-Post (334 recurrence-training pairs; old 153 source-patient total needs a manifest; separate forecast cohort is 133/328), UCSF-PDGM + UCSF-ALPTDG (59 total former test, 49 enhancing-evaluable; completed and now SPENT), CGGA and TCGA clinical tables.
- Models: monotonicity-checked Fisher-Kolmogorov simulator (2 mm; synthetic numerical checks, notable speed/resolution error, ledger 2q/2ac), isotropic and DTI-atlas anisotropic tensors; 3D U-Net with soft top-k ranking loss; released PREDICT-GBM U-Net (rebuilt; exact reproduction is unverified without a hash/comparison artifact).
- Compute: Author-reported Kaggle/CPU availability; T4 quota/no-paid-cluster assertions lack dated account/environment evidence.
- Rules: prospective declaration and no tuning-until-positive are intended rules, not complete historical access certification. Amendment 12 reused H-2/H-3 as development; both sets are spent and cannot be future locked tests. Negative results stay reported.

## Part 1 (historical primary null, positive secondary): predict where GBM recurs, as a radiotherapy target

Result so far: at equal volume to the 15 mm standard plan, all-recurrence coverage +0.081 (CI +0.065 to +0.097) on 101 unseen patients; enhancing-recurrence coverage +0.009 (n.s.); both are endpoint-specific: n=101 secondary all recurrence and n=98 primary enhancing (ledger 2ax). UCSF completed: secondary +0.0849 [0.0588,0.1150], n=59; primary +0.0044 [-0.0182,0.0369], n=49 (2bd). Neither primary met its rule; both former test sets are spent.
Find:
1. Methods published 2022-2026 that improve ENHANCING-recurrence coverage over a 15 mm (or 20 mm) CTV at equal volume (deep learning, physics, hybrids, diffusion MRI, perfusion, PET, radiomics). Give the reported gain and the test design.
2. Public datasets with pre-op MRI + segmentation + a later recurrence segmentation that are NOT LUMIERE, RHUH, TUM, MU, UCSF (for a third external test). State the access type (open / data use agreement / restricted).
3. Known reasons why enhancing recurrence is hard to predict (distant recurrence share, pseudoprogression, registration error) and how papers handle them.

## Part 2 (failed): predict GBM growth over time (digital twin)

Result so far: on MU-Glioma-Post (2-3 scans per patient, months apart) a Fisher-Kolmogorov forecast did not beat "no change" (Dice +0.013, CI includes 0); no paired PDE advantage over the population rule was demonstrated (ledger 2y PDE-minus-geometric -0.0019 [-0.0173,0.0130], not equivalence). Older +0.026 script-81 and baseline-ladder values use different splits/windows. DTI comparison needs the exact older -0.001,p=.97 context or later monotone 2ak. Fitted kill multiplier remains exploratory (ledger 2at/2av); no patient-specific drug efficacy follows.
Find:
1. Papers where a patient-specific growth model beat a persistence (no-change) or volume-extrapolation baseline on held-out later scans of REAL GBM patients. Give the baseline, metric, n and data.
2. Methods for calibration with only 2-3 time points: Bayesian / hierarchical priors, physics-informed neural networks, neural operators, data assimilation, reduced-order models, mixed-effects growth models.
3. Public longitudinal GBM datasets with 3+ MRI time points per patient, ideally with segmentations and dates.
4. What a fair "digital twin" claim needs (definitions in recent reviews, e.g. NASEM 2023 report on digital twins), so the project uses the term honestly.

## Part 3 (failed): per-patient treatment effect (drug and radiation suppression)

Result so far: radiation and temozolomide effects use assumed population kill rates (linear-quadratic alpha/beta, TMZ kill); the fitted multiplier was exploratory and did not demonstrate a patient forecast advantage (ledger 2at/2av; source119). Yes/no treatment fields use assumed timing rather than detailed observed dose/response.
Find:
1. Public GBM imaging datasets that include per-patient treatment records: radiotherapy dose and dates, TMZ cycles and dates, extent of resection, steroids, bevacizumab.
2. Papers that estimated per-patient radiosensitivity or TMZ response from imaging and validated it on later scans (e.g. proliferation-saturation index, LQ-model fits, ADC change, perfusion).
3. Imaging or molecular biomarkers that predict treatment response and are available in public data (MGMT methylation, IDH, ADC, rCBV, amino-acid PET).
4. Whether "virtual treatment" comparison (simulating different schedules for one patient) has been validated anywhere, or remains simulation only.

## Output format

For each part: a table (method or dataset | citation link | key number | test design | access | fits our compute? yes/no). Then a ranked list of the 3 most promising next steps with a one-sentence reason each. Then a list of claims you could NOT verify. Keep it under 2,500 words.
