> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Research prompt: methods and data to make the GBM digital twin work

Paste everything below the line into a new Claude session (deep research mode if available).

---

You are a research assistant for a high-school science fair project on glioblastoma (GBM). Find published methods, public datasets and code that could make three parts of the project work. Report only what you can cite with a working link (paper DOI/arXiv, dataset page, GitHub). Never invent a source, a number or a dataset. If you are not sure, say "not verified".

## Project context (facts, do not change)

- Data already used: PREDICT-GBM (TUM 142 dev; LUMIERE 61 + RHUH 40 locked test), MU-Glioma-Post (153 patients, 334 scan pairs), UCSF-PDGM + UCSF-ALPTDG (59 grade-4 patients, new locked test), CGGA and TCGA clinical tables.
- Models: exact monotone Fisher-Kolmogorov solver (2 mm), isotropic and DTI-atlas anisotropic tensors; 3D U-Net with soft top-k ranking loss; released PREDICT-GBM U-Net (rebuilt, reproduces released predictions exactly).
- Compute: Kaggle free GPU (T4, 30 h/week), CPU PC. No paid cluster.
- Rules: every test is declared before the data is read; a locked test set is used once; no re-running until positive; negative results are reported.

## Part 1 (works): predict where GBM recurs, as a radiotherapy target

Result so far: at equal volume to the 15 mm standard plan, all-recurrence coverage +0.081 (CI +0.065 to +0.097) on 101 unseen patients; enhancing-recurrence coverage +0.009 (n.s.); a second locked test on UCSF is pending.
Find:
1. Methods published 2022-2026 that improve ENHANCING-recurrence coverage over a 15 mm (or 20 mm) CTV at equal volume (deep learning, physics, hybrids, diffusion MRI, perfusion, PET, radiomics). Give the reported gain and the test design.
2. Public datasets with pre-op MRI + segmentation + a later recurrence segmentation that are NOT LUMIERE, RHUH, TUM, MU, UCSF (for a third external test). State the access type (open / data use agreement / restricted).
3. Known reasons why enhancing recurrence is hard to predict (distant recurrence share, pseudoprogression, registration error) and how papers handle them.

## Part 2 (failed): predict GBM growth over time (digital twin)

Result so far: on MU-Glioma-Post (2-3 scans per patient, months apart) a Fisher-Kolmogorov forecast did not beat "no change" (Dice +0.013, CI includes 0); a population growth-rate rule was as good; DTI direction did not help; one per-patient kill multiplier fitted on the previous pair was null.
Find:
1. Papers where a patient-specific growth model beat a persistence (no-change) or volume-extrapolation baseline on held-out later scans of REAL GBM patients. Give the baseline, metric, n and data.
2. Methods for calibration with only 2-3 time points: Bayesian / hierarchical priors, physics-informed neural networks, neural operators, data assimilation, reduced-order models, mixed-effects growth models.
3. Public longitudinal GBM datasets with 3+ MRI time points per patient, ideally with segmentations and dates.
4. What a fair "digital twin" claim needs (definitions in recent reviews, e.g. NASEM 2023 report on digital twins), so the project uses the term honestly.

## Part 3 (failed): per-patient treatment effect (drug and radiation suppression)

Result so far: radiation and temozolomide effects use assumed population kill rates (linear-quadratic alpha/beta, TMZ kill); the fitted per-patient effect was null; datasets record treatment only as yes/no with assumed timing.
Find:
1. Public GBM imaging datasets that include per-patient treatment records: radiotherapy dose and dates, TMZ cycles and dates, extent of resection, steroids, bevacizumab.
2. Papers that estimated per-patient radiosensitivity or TMZ response from imaging and validated it on later scans (e.g. proliferation-saturation index, LQ-model fits, ADC change, perfusion).
3. Imaging or molecular biomarkers that predict treatment response and are available in public data (MGMT methylation, IDH, ADC, rCBV, amino-acid PET).
4. Whether "virtual treatment" comparison (simulating different schedules for one patient) has been validated anywhere, or remains simulation only.

## Output format

For each part: a table (method or dataset | citation link | key number | test design | access | fits our compute? yes/no). Then a ranked list of the 3 most promising next steps with a one-sentence reason each. Then a list of claims you could NOT verify. Keep it under 2,500 words.
