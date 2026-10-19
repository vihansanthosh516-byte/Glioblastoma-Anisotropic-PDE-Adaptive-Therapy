# Problem-first brief: are the RANO 2.0 volume thresholds above the measurement noise? (DRAFT 2026-10-18)

Status: CURRENT explanatory brief, updated 2026-10-19. Section 4 retains the original proposal as history, superseded by the completed analyses below. Its tests were later declared and run (Amendments 15 to 21; ledger 2br to 2bx); see START_HERE.md for current state. Novelty = not found in about 12 searches (snippet level); UNVERIFIED.

## 1. The problem (from the literature)
- RANO 2.0 volumetric rules: progression = volume increase of 40% or more; partial response = decrease of 65% or more; reference = nadir for progression, baseline for response. The 0.5 cm3 floor is our operational choice, NOT in the RANO 2.0 text (checked 2026-10-18). Sources: https://pmc.ncbi.nlm.nih.gov/articles/PMC11630866/ (guide), https://pmc.ncbi.nlm.nih.gov/articles/PMC12986764/ (RANO 2.0 vs mRANO vs RANO, 137 patients).
- Segmentation variability motivates this work. The earlier ICC range 0.22 to 0.96 and quoted expert-supervision requirement are unverified here; exact source table/section is missing (https://pmc.ncbi.nlm.nih.gov/articles/PMC10513769/, previously attributed to a LUMIERE study with 80 patients, 502 time points and two tools; these source-specific counts remain unverified, not the dataset inventory).
- A 2026 narrative review states that no automated tool had been formally validated against RANO 2.0 (as of Nov 2024; a historical review statement, not evidence of current absence) (https://www.frontiersin.org/journals/neurology/articles/10.3389/fneur.2026.1797487/full).
- NOT found in the approximately 12 snippet-level searches completed (no retained search log): any study that computes a minimal detectable change or within-subject variability for glioblastoma volumes, or that estimates how often segmentation error flips a RANO category near the 40% / 65% lines.

## 2. Established method families proposed for use (not invented here)
The linked application and exact source equations have not been checked in full here. The formula below assumes independent repeat errors with equal variance and a normal-error approximation.

- Bland-Altman limits of agreement and minimal detectable change MDC = 1.96 x sqrt(2) x SEM (used for glioma segmentation comparisons, https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0164891).
- Error propagation for the log volume ratio; reliable-change logic; Bayesian classification with measurement error.

## 3. What the repo already holds (verified in files, exploratory)
`output/seg_noise_floor.json` (script 104, LUMIERE, two automated tools HD-GLIO-AUTO and DeepBraTumIA, contrast-enhancing label): 599 scans, 419 with both positive, 89 patients.
- Median |ln volume ratio| between the two tools: 1.37.
- Noise floor in ln volume by tumour-size tertile (median volume): 1.03 (0.9 cm3), 0.88 (4.0 cm3), 0.35 (15.7 cm3).
- Change disagreement (|ln change in tool A - ln change in tool B|, removes a constant tool bias): median 1.90 over 330 pairs; 1.20 over 181 pairs of RATED follow-ups (58 patients).
- For scale: ln(1.4) = 0.34 (progression line), ln(0.35) = -1.05 (response line).
Also: `output/upenn_rater_noise.json` (script 109) = expert correction of an AI mask on 220 scans (an anchored-correction sensitivity estimate; the lower estimate is zero/degenerate in many reported strata).
Strongest caveat: two tools with different label definitions are not a scan-rescan pair. Disagreement overstates pure measurement noise; expert correction understates it. These are sensitivity brackets, not formal confidence bounds on pure measurement error; definition effects contribute to tool disagreement.

## 4. Original proposal: probabilistic RANO (superseded status)
Current status: A15–A21 analyses have run (ledger 2br–2bx); the calculator and timeline viewer exist (2bv, 2bw), with arithmetic/display checks only. A1 proof and clinical validation remain unfinished. The following records the original proposal.

Question: given measured segmentation noise, how often does the RANO 2.0 category depend on the tool or the noise, and what size and change are needed before a call is reliable?
Mathematics: variance of ln(V1/V0) from per-scan noise that depends on volume; MDC curve versus volume; probability of each RANO class as a function of true change and noise; optimal threshold and "indeterminate band" under a stated loss; hierarchical measurement model with two tools as replicates; identifiability of tool bias versus noise.
Computation: LUMIERE two tools plus expert RANO ratings (data on disk), UPenn expert-correction noise, simulations to check the formulas; a 4D / timeline viewer showing volume with an uncertainty band and RANO class probabilities for a patient.
Declared tests (each written before it runs): flip rate of RANO class between the two tools on rated follow-ups; agreement of each tool with the expert rating; the volume below which a +40% call cannot be made at 95% confidence, for the lower and upper noise bounds.
Risks: tool disagreement may be mostly definition mismatch; real scan-rescan data are controlled-access; LUMIERE ratings cover selected follow-ups; the finding may be "known to radiologists".
Link to the old work: measurement noise may constrain forecasts, but no causal explanation of model performance was established. PDE forecasts did not significantly beat persistence; the population-rate geometric baseline had positive GBM-only gain (ledger 2w).
