# Problem-first brief: are the RANO 2.0 volume thresholds above the measurement noise? (DRAFT 2026-10-18)

Status: research brief only. Nothing here is declared or run. Novelty = not found in about 12 searches (snippet level); UNVERIFIED.

## 1. The problem (from the literature)
- RANO 2.0 volumetric rules: progression = volume increase of 40% or more; partial response = decrease of 65% or more; measurability floor 0.5 cm3. Sources: https://pmc.ncbi.nlm.nih.gov/articles/PMC11630866/ (guide), https://pmc.ncbi.nlm.nih.gov/articles/PMC12986764/ (RANO 2.0 vs mRANO vs RANO, 137 patients).
- Volumes come from segmentation, which is noisy: manual and automated inter-rater ICC from 0.22 to 0.96 across studies; automated tools "need expert supervision" (https://pmc.ncbi.nlm.nih.gov/articles/PMC10513769/, LUMIERE, 80 patients, 502 time points, two tools).
- A 2026 narrative review states that no automated tool had been formally validated against RANO 2.0 (as of Nov 2024) (https://www.frontiersin.org/journals/neurology/articles/10.3389/fneur.2026.1797487/full).
- NOT found in searches: any study that computes a minimal detectable change or within-subject variability for glioblastoma volumes, or that estimates how often segmentation error flips a RANO category near the 40% / 65% lines.

## 2. Solution tools that exist in the literature (not invented here)
- Bland-Altman limits of agreement and minimal detectable change MDC = 1.96 x sqrt(2) x SEM (used for glioma segmentation comparisons, https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0164891).
- Error propagation for the log volume ratio; reliable-change logic; Bayesian classification with measurement error.

## 3. What the repo already holds (verified in files, exploratory)
`output/seg_noise_floor.json` (script 104, LUMIERE, two automated tools HD-GLIO-AUTO and DeepBraTumIA, contrast-enhancing label): 599 scans, 419 with both positive, 89 patients.
- Median |ln volume ratio| between the two tools: 1.37.
- Noise floor in ln volume by tumour-size tertile (median volume): 1.03 (0.9 cm3), 0.88 (4.0 cm3), 0.35 (15.7 cm3).
- Change disagreement (|ln change in tool A - ln change in tool B|, removes a constant tool bias): median 1.90 over 330 pairs; 1.20 over 181 pairs of RATED follow-ups (58 patients).
- For scale: ln(1.4) = 0.34 (progression line), ln(0.35) = -1.05 (response line).
Also: `output/upenn_rater_noise.json` (script 109) = expert correction of an AI mask on 220 scans (a LOWER bound on noise).
Strongest caveat: two tools with different label definitions are not a scan-rescan pair. Disagreement overstates pure measurement noise; expert correction understates it. The honest result is a bracket, not a single number.

## 4. Candidate project: probabilistic RANO
Question: given measured segmentation noise, how often does the RANO 2.0 category depend on the tool or the noise, and what size and change are needed before a call is reliable?
Mathematics: variance of ln(V1/V0) from per-scan noise that depends on volume; MDC curve versus volume; probability of each RANO class as a function of true change and noise; optimal threshold and "indeterminate band" under a stated loss; hierarchical measurement model with two tools as replicates; identifiability of tool bias versus noise.
Computation: LUMIERE two tools plus expert RANO ratings (data on disk), UPenn expert-correction noise, simulations to check the formulas; a 4D / timeline viewer showing volume with an uncertainty band and RANO class probabilities for a patient.
Declared tests (each written before it runs): flip rate of RANO class between the two tools on rated follow-ups; agreement of each tool with the expert rating; the volume below which a +40% call cannot be made at 95% confidence, for the lower and upper noise bounds.
Risks: tool disagreement may be mostly definition mismatch; real scan-rescan data are controlled-access; LUMIERE ratings cover selected follow-ups; the finding may be "known to radiologists".
Link to the old work: the same noise floor limits growth forecasts, which explains why growth models did not beat "no change".
