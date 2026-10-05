# Limits: every way I found to fix each one

Written 2026-10-04. Two source types. **[S]** = a web search this session (link given; I read search summaries, not full papers). **[K]** = my general knowledge, not searched, so check before citing. **[R]** = checked in this repo today.
Cost key: **now** = existing data and code, hours. **download** = needs a public dataset. **compute** = long run. **no** = cannot be done in a science-fair project, so state it as future work.

## New limit found in this research (not in the earlier list)
**L0. C-GAT and all Track A classifiers use a random cell-level split, not a patient-level split.** [R] `src/04_export_for_attention_model.py` makes a stratified 80/20 split with `random_state=42` over 15,000 cells. `src/13_gat_train.py` reuses it. The graph also carries a `same_patient` edge feature (`src/12_gat_build_graph.py`), and test nodes sit in the training graph. Cells from one patient land in train and test, so accuracy (78.7%, and the 94% in the README) may be inflated. Other fields saw large drops under patient-level evaluation: accuracy 98.6% to 95.9% in one case and 96.9% to 78.0% in another [S](https://www.biorxiv.org/content/10.1101/2025.05.12.653389.full.pdf) (those were not GBM studies).
- **Fix 1 (now):** leave-one-patient-out or grouped k-fold on `obs["ID1"]` for LR, RF, NMF and scVI. Cheap.
- **Fix 2 (now):** retrain C-GAT with test patients removed from the graph (inductive). Small model (41,499 parameters).
- **Fix 3 (now):** 5 seeds plus bootstrap CIs for every row of the leaderboard.
- **Fix 4 (compute):** use all 140,355 cells, not 15,000.
- **Fix 5 (now):** test for batch effects: batch-correction methods often leave residual batch signal that classifiers exploit [S](https://pmc.ncbi.nlm.nih.gov/articles/PMC10914288).
- **If accuracy drops:** report the patient-level number as the headline and keep the old one as "cell-level, optimistic".

## Data limits
**L1. Forecast: Dice only 0.26, median gain +0.005.**
1. Add treatment information. Treatment-aware diffusion models (TaDiff) condition on MRI plus treatment type and day [S](https://pubmed.ncbi.nlm.nih.gov/40031286/). Caveat: they encode treatment as global variables, not voxel-wise dose [S](https://arxiv.org/pdf/2606.05113). Cost: compute, plus training data.
2. Add voxel-wise radiotherapy dose. Burdenko-GBM-Progression has 180 patients, pre-RT MRI, follow-ups (1 to 7 time points) and RTDOSE/RTPLAN for a subset [S](https://wiki.cancerimagingarchive.net/x/KgWwC). Cost: download.
3. Hybrid mechanistic-learning: a growth ODE plus a guided diffusion model, with radiotherapy effects in the ODE [S](https://arxiv.org/abs/2509.09610v1). Cost: compute.
4. Spatially varying parameters from Bayesian inference [S](https://arxiv.org/pdf/2209.12089). Cost: compute.
5. Stronger baselines before claiming a gain: linear volume extrapolation, per-patient growth-rate carry-forward [K]. Cost: now.

**L2. Forecast loses on shrinking tumours (80/152 shrank or stayed).**
1. Model treatment kill by dose (L1 fix 2). The simple kill term already made the forecast worse (`NEGATIVES_REVISITED.md`), so use dose maps, not a global term.
2. Separate pseudoprogression and cavity collapse. LUMIERE has expert RANO ratings for 91 GBM patients, 638 study dates [S](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9755255/). Cost: download.
3. Better stratification. Script 79's growth classifier failed (AUC 0.46) using only scan history and timing. Add MGMT, IDH, extent of resection. RHUH-GBM has molecular data, volumetric assessments and survival for patients with total or near-total resection [S](https://arxiv.org/pdf/2305.00005). Cost: download.
4. Report shrinking tumours as out of scope in the abstract. Cost: now.

**L3. UCSF-PDGM: grade 2-3, one scan per patient, no replication.**
1. Replicate the tensor study on a second cohort. I did not find a second public longitudinal GBM set with DTI in these searches [S], so this may need a new data source. Cost: download, unknown.
2. Clean the tensors first. Free-water elimination corrected FA separated later-recurrence areas from recurrence-free edema (AUC 0.9 vs 0.77 uncorrected) [S](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12696142/). Apply to UCSF-PDGM, then rerun the three modes. Cost: now plus compute.
3. State "grade 2-3, shape test only" in the paper. Cost: now.

**L4. Only 21 real patients (`real_test`) in Track C.**
1. Fit growth rates on LUMIERE (91), RHUH-GBM and Burdenko (180) and make a larger real set. Cost: download plus compute.
2. Real kill rates are not observed anywhere in the vault. Fit them from response during therapy in patients with dated treatment and volumes [K]. Cost: now to compute.
3. Report intervals for every real_test number. Cost: now.

**L5. D cannot be identified; 124/154 growth rates sit at the lower bound.**
1. Report posteriors, not points. Bayesian calibration (GPHMC) returns the posterior and parameter correlations; the field names non-identifiability as a main difficulty [S](https://hal-univ-tlse3.archives-ouvertes.fr/INRIA/hal-01324849). Cost: compute.
2. Use informative population priors, e.g. on the D/rho ratio, then update [K].
3. Use spatial information, not volume only: contours at two thresholds (T2/FLAIR vs T1Gd) [K]. The repo already has labels 1, 2, 3.
4. Physics-informed networks can estimate D and rho from a single 3D snapshot [S](https://arxiv.org/pdf/2412.05330). Cost: compute.
5. Use 3 or more scans per patient. Cost: now.
6. Fix the lower-bound problem: 64 of the 124 are shrinking tumours (`output/inverse_est_metrics.json`). Fit a signed model with a response term, not a pure growth bound. Cost: now.

**L6. 15,000-cell subsample, one seed, no CI.** See L0 fixes 3 and 4.

**L7. TCGA: n = 150 with expression; age is a known factor.**
1. External validation on CGGA. A published study used TCGA (n = 518) for training and CGGA (n = 350) for validation [S](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2022.900911/pdf). Cost: download.
2. Add IDH, MGMT and extent of resection. A review found only 5.9% of studies used extent of resection, 17.6% IDH and 29.4% MGMT [S](https://medrxiv.org/node/464632). Restrict to IDH-wildtype GBM. Cost: now to download.
3. Score the inflammation signature as a survival predictor on CGGA. Cost: download.

## Model limits
**L8. Adaptive results rest on an assumed kill scale (E_MAX = rho x 1000).**
1. Calibrate the kill scale to response curves in MU-Glioma-Post patients on therapy [K]. Cost: now.
2. Report outcomes across a sweep of the scale, not one value. The vault has a 3-patient sweep only. Cost: compute.

**L9. Script 48 (3D) and script 60 used D about 10x below the literature.** Script 60 was re-run (scripts 66, 77, 80). **Script 48 was not.**
1. Rerun script 48 at D_white 0.13 and D_gray 0.013 (the Swanson-type values in `Script-78-Horizon-Crossover.md`). 8 patients, 50^3 grid, 3600 steps. Cost: compute, probably hours. This is the cheapest fix of the whole list.

**L10. Script 47 dual-agent uses 2.5-3.3x more drug than MTD.**
1. Compare against a dual-agent MTD, not a single-agent MTD. Cost: now.
2. Put drug exposure in one unit for all arms. Cost: now.
3. Add a dose-matched arm, as in the equal-budget RL work. Cost: compute.

**L11. No clonal competition; resistance is a fixed fraction and partly plastic.**
1. Adaptive therapy needs competition between sensitive and resistant cells. In the Moffitt prostate trial (17 patients) the aim was to keep sensitive cells to suppress resistant ones [S](https://elifesciences.org/articles/76284). The trial used PSA as a trackable biomarker. GBM has no such blood marker [K], so the main barrier is measurement. State this plainly.
2. Add history-dependent phenotype variables, calibrated to spheroid data with temozolomide [S](https://www.biorxiv.org/content/10.1101/2023.11.24.568421.full.pdf). Cost: compute. The repo's script 83 already has a two-state switching model, so this would be a richer version.
3. Add spatial competition (resistant and sensitive cells in the PDE, not a well-mixed Lotka-Volterra model) [K]. Cost: compute.
4. Serial ctDNA or biopsy data to calibrate resistance would be ideal. It is not publicly available for GBM [K]. Cost: no.

**L12. Anisotropy never matched to real tumour shape; tract alignment is lower under aniso.**
1. Check the alignment metric. Script 74 reports tract alignment as 0.4535 (aniso) vs 0.5011 (iso) [R]. I did not read its definition. Confirm what it measures and why the isotropic value sits near 0.5. Cost: now.
2. Test on recurrence location. Standard planning is an isotropic margin, and 52% of patients progressed outside the 2 cm margin in one study; fibre-density-weighted maps improved cell-migration prediction [S](https://www.mate.polimi.it/biblioteca/add/qmox/17-2018.pdf). Needs DTI plus recurrence masks. Cost: download, unknown.
3. Compare to the Painter-Hillen formulation, which has a specific diffusion-tensor model for white-matter tracts [S](https://www.math.ualberta.ca/~thillen/paper/PainterHillen2013.pdf). Cost: now to compute.

## Evaluation limits
**L13. RL wins only on day-90 volume; never wins on time to progression.**
1. Train and score on TTP. Script 68 started this (`dagger_ttp`). Extend it to every policy. Cost: compute.
2. Use survival-aware off-policy evaluation. A 2026 paper (arXiv 2603.22900) tailors OPE for right-censored outcomes [S](https://arxiv.org/pdf/2603.22900). Needs real treatment variation. Cost: download.
3. Offline RL work for chemotherapy dosing stresses safety constraints and heterogeneous test conditions [S](https://arxiv.org/pdf/2508.17212). Add a safety gate and report it. Cost: now.

**L14. Pacing gains only +2.1 d on real parameters.**
1. Find out why. Compare real_test kill rates and growth rates to the synthetic sets [K]. Cost: now.
2. Enlarge the real set (L4). Cost: download.
3. Report the synthetic sets as model checks, not evidence. Cost: now.

**L15. Probe-then-commit needs precise measurement.**
1. Automated volumetry is much more reproducible than manual RANO: ICC 0.97 vs 0.42 (1D) and 0.61 (2D) [S](https://arxiv.org/pdf/2209.01402). Use the measured volume CV as the noise level, not 5-10% [K].
2. Infer kill rates in a Bayesian way with an explicit noise model. Cost: compute.
3. Longer probes (L = 3 and more) average noise. Cost: now (already simulated).
4. A power analysis: noise level at which the rule's win rate falls below the blind policy. Cost: now.

**L16. The 0.02 threshold has 10 patients and was chosen after seeing the data.**
1. Fix the cutoff before testing on a new cohort (LUMIERE, Burdenko) [K]. Cost: download.
2. Fit a continuous model (adaptive earlier vs rho) with a confidence band, not a cutoff. Cost: now.
3. Report the sharp 10/10 vs 0/51 split as descriptive only. Cost: now.

**L17. Forecast target chosen after a volume probe.**
1. Run the unrun secondary: whole tumour without cavity (`--target wt`, script 81). Cost: compute (about 12 h per the vault).
2. Report all targets in one table. Cost: now.
3. Write future targets into a script header before any run (already the repo's habit). Cost: now.

## Negative-result limits
**L18. D_f is not valid.** Already replaced by elongation. To keep D_f, run box counting on larger tumours or finer grids [K]. Cost: compute.

**L19. MGMT does not predict TTP (n = 130, p = 0.60).**
1. Power analysis: what n would detect the 169 vs 181 d gap? [K] Cost: now.
2. More patients (L4 sets with MGMT). Cost: download.
3. Use overall survival, where MGMT is known to matter (Hegi 2005 [VERIFY]), not TTP. Cost: now where OS exists.
4. Adjust for age and extent of resection. Cost: now to download.

**L20. Resistance-driven adaptive gain fails on real data (scripts 76, 82, 83).** This is a finding. Fixes are the L11 model upgrades. If three model classes keep failing, say "does not occur under these models" and stop.

**L21. Growth classifier failed (AUC 0.46).** Add dose, MGMT, IDH, extent of resection (L2 fix 3) or drop it. A failed classifier is a valid result.

## Not tested
**L22. No prospective validation.** Not possible here (no). Closest: temporal validation (train on early-scanned patients, test on later) and a pre-registered analysis plan for the next cohort [K]. State it as future work.

**L23. 3D at literature diffusivity.** See L9.

**L24. Origin of the inflammation-score values (scripts 43-44).** [R] `get_patient_inflammation` reads a per-zone expression table by `patient_id`, and returns 1.0 (no effect) when the patient is absent. Check whether the MU-Glioma patient IDs match any expression rows. If few or none do, the score does nothing and the A-to-B link in the schematic should say so. Cost: now.

**L25. Drug screen is in silico only.**
1. Compare predicted targets to DepMap CRISPR essentiality in 60 patient-derived GBM lines [S](https://era.ed.ac.uk/items/2b20d86b-08e2-4ee8-9fee-951475f046a2) and to the wider DepMap resource [S](https://forum.depmap.org/t/validation-of-crispr-ko/1040). Cost: download.
2. **A caution [K, not searched]:** the top hit MT-CO2 is encoded in mitochondrial DNA. Nuclear CRISPR screens do not knock out mitochondrial genes, so DepMap may not contain it. Check how the "virtual knockout" treats it before claiming a target.
3. Check toxicity to normal tissue with healthy-brain expression. Cost: download.

## Suggested order (cheapest and most important first)
1. **L0** patient-level split for Track A (cheap; the biggest risk to a headline number).
2. **L24** check the inflammation-score link (minutes).
3. **L9** rerun script 48 at literature D.
4. **L10** fair dual-agent baseline.
5. **L17** run the whole-tumour target.
6. **L7** external TCGA-style validation on CGGA, plus added covariates.
7. **L16, L19** power analyses and continuous models.
8. **L1-L4** new cohorts (LUMIERE, RHUH-GBM, Burdenko) for a real external set.
Everything else is either future work or a model upgrade.

## Open reference item
"Weidner 2023": two web searches found no paper by that name. The nearest 2023 hits were "Predictive Digital Twin for Optimizing Patient-Specific Radiotherapy Regimens under Uncertainty in High-Grade Gliomas" [S](https://arxiv.org/pdf/2308.12429) and a Frontiers in AI paper from October 2023 [S](https://www.frontiersin.org/journals/artificial-intelligence/articles/10.3389/frai.2023.1222612/pdf). I do not know the authors. Do not cite "Weidner 2023" until you have the title or DOI.
