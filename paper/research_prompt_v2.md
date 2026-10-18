# Research prompt v2 (paste everything below the line into a new Claude session, deep research mode)

---

You are a research assistant for a high-school science fair project on glioblastoma (GBM). Find published methods, public datasets and code that could IMPROVE the results below. Report only what you can cite with a working link (DOI, arXiv, dataset page, GitHub). Never invent a source, number or dataset. If unsure, write "not verified". Say if you read the full paper or only the abstract.

## Facts (do not change)
- Task: from a PRE-OPERATIVE MRI, predict where GBM recurs, then build a radiotherapy plan at the SAME VOLUME as the standard 15 mm margin around the tumor core.
- Model: 3D U-Net, 2 mm grid, loss = BCE + Dice + soft top-k ranking. Inputs: pre-op MRI (T1c, FLAIR), tumor core mask, tissue maps. A Fisher-Kolmogorov PDE model was also tried.
- Training data: 636 samples (TUM 142, MU-Glioma-Post 334 pairs, LUMIERE 61, RHUH 40, UCSF 59).
- Compute: Kaggle free T4 (30 h per week), one CPU PC. No paid cluster.
- Rules: each test is declared before data is read; a locked test set is used once; negative results are reported.

## Results so far (all coverage = fraction of recurrence inside the plan, equal volume, minus the 15 mm plan)
- All recurrence: +0.081 (CI +0.065 to +0.097, n=101) and +0.085 (CI +0.059 to +0.115, n=59). Works, replicated twice.
- ENHANCING recurrence: +0.009 (CI -0.012 to +0.031, n=98) and +0.004 (CI -0.018 to +0.037, n=49). Null, twice.
- Known limit: in 32 of 49 UCSF patients the 15 mm plan already covers at least 0.999 of enhancing tissue, so there is little room to gain.
- Growth forecast (FK PDE vs "no change"): no gain. Dice +0.013, CI includes 0. DTI direction did not help.
- Per-patient treatment effect (radiation, TMZ): fitted effect was null.
- No eligible new locked test set found. CFB-GBM v2.0 (no pre-op MRI, one merged GTV label) and Burdenko-GBM-Progression (post-op baseline, restricted licence, no enhancing label) were rejected.

## Questions
A. Enhancing recurrence (main goal)
1. Why is enhancing recurrence hard to predict from pre-op MRI? Cite studies on distant recurrence share, pseudoprogression, registration error, residual tumor after surgery.
2. Methods 2022-2026 that raise enhancing-recurrence coverage above a 15 mm or 20 mm CTV at equal volume. Give gain, test design, n.
3. Does adding the post-operative MRI (resection cavity, residual tumor) help? Which papers show this, and by how much?
4. Which extra inputs help: diffusion (ADC, DTI), perfusion (rCBV), amino-acid PET, MGMT/IDH, extent of resection? Which are in public data?
5. Training tricks for small data: pretraining, foundation models, multi-task loss, label smoothing for registration error, test-time augmentation, ensembling. Gains with CI.
6. Better loss or metric for equal-volume plan coverage (e.g. direct optimisation of coverage, differentiable top-k, conformal or risk-controlled plans).

B. New locked test set
7. Public datasets with pre-op MRI + segmentation + later recurrence segmentation (with a separate enhancing label), at least 40 patients, NOT LUMIERE, RHUH, TUM, MU, UCSF, CFB-GBM, Burdenko. Check BraTS 2025, UCSD-PTGBM, GLIS-RT, IvyGAP, CPTAC-GBM, TCGA-GBM, ReMIND, others. State access type (open / DUA / restricted) and what labels exist.

C. Growth forecast and digital twin
8. Papers where a patient-specific growth model beat a persistence or volume-extrapolation baseline on held-out later scans of real GBM patients. Give baseline, metric, n.
9. Calibration with only 2-3 time points: hierarchical Bayes, mixed-effects, PINNs, neural operators, data assimilation.
10. What a fair "digital twin" claim needs (NASEM 2023 report, recent oncology reviews). Quote the definition and the validation requirements.

D. Treatment effect
11. Public GBM imaging sets with per-patient radiotherapy dose, TMZ dates, resection extent, steroids.
12. Per-patient radiosensitivity or TMZ response estimated from imaging and validated on later scans.
13. Has any "virtual treatment schedule" comparison been validated on real patients, or is it simulation only?

## Output (max 3000 words)
For each section A to D: a table (method or dataset | link | key number | test design | access | read full paper? | fits Kaggle T4? yes/no). Then the 3 best next steps, ranked, one sentence each. Then a list of claims you could NOT verify. Flag any paper that tests on its own training centre only.
