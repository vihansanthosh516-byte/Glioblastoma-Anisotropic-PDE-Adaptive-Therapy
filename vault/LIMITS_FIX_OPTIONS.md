> **HISTORICAL BACKGROUND (corrected 2026-10-19):** Older GBM proposals, not the current goal or permission to run analyses. See START_HERE.md for probabilistic RANO. All [S] external numbers/design/availability claims remain unverified without primary-source tables/sections; [K] suggestions are not results. Cost labels and old run-status claims are historical estimates, not tested feasibility. All treatment/RL comparisons below are simulations with assumed kill rates.

# Limits: every way I found to fix each one

Initial draft 2026-10-04; audit/status corrections 2026-10-19. Three source types. **[S]** = a web search this session (link given; I read search summaries, not full papers). **[K]** = my general knowledge, not searched, so check before citing. **[R]** = originally checked in the repo on 2026-10-04; later corrections are dated above.
Cost key: **now** = existing data and code, hours. **download** = needs a public dataset. **compute** = long run. **no** = cannot be done in a science-fair project, so state it as future work.

## New limit found in this research (not in the earlier list)
**L0. C-GAT and all Track A classifiers use a random cell-level split, not a patient-level split.** [R] `src/04_export_for_attention_model.py` makes a stratified 80/20 split with `random_state=42` over 15,000 cells. `src/13_gat_train.py` reuses it. The graph also carries a `same_patient` edge feature (`src/12_gat_build_graph.py`), and test nodes sit in the training graph. Cells from one patient land in train and test, so the old random-cell accuracy 78.7% is optimistic (ledger A1/2o); the former README 94% claim is withdrawn here because its exact run/source is missing. Unverified non-GBM context (exact source table not checked): other fields reportedly saw large drops under patient-level evaluation: accuracy 98.6% to 95.9% in one case and 96.9% to 78.0% in another [S](https://www.biorxiv.org/content/10.1101/2025.05.12.653389.full.pdf) (those were not GBM studies).
- **Fix 1 (now):** leave-one-patient-out or grouped k-fold on `obs["ID1"]` for LR, RF, NMF and scVI. Cheap.
- **Fix 2 (now):** retrain C-GAT with test patients removed from the graph (inductive). The historical 41,499-parameter count lacks an exact model/config count record here and is unverified.
- **Fix 3 (now):** 5 seeds plus bootstrap CIs for every row of the leaderboard.
- **Fix 4 (compute):** use all 140,355 cells, not 15,000 (array inventory: `output/cgat/cvae_latent_full.npy`, shape 140355 by 32; a proposal, not an authorized run).
- **Fix 5 (now):** test for batch effects: batch-correction methods often leave residual batch signal that classifiers exploit [S](https://pmc.ncbi.nlm.nih.gov/articles/PMC10914288).
- **If accuracy drops:** report the patient-level number as the headline and keep the old one as "cell-level, optimistic".

## Data limits
**L1. Historical selected-target forecast summary, superseded by ledger 2w/2y.** Cellular-core OOF mean Dice was 0.2591 vs persistence 0.2332 (`output/forecast_labels_core/results.json`); the later manifest PDE gain was inconclusive and did not beat the population-rate baseline.
1. Add treatment information. Treatment-aware diffusion models (TaDiff) condition on MRI plus treatment type and day [S](https://pubmed.ncbi.nlm.nih.gov/40031286/). Caveat: they encode treatment as global variables, not voxel-wise dose [S](https://arxiv.org/pdf/2606.05113). Cost: compute, plus training data.
2. Add voxel-wise radiotherapy dose. Burdenko-GBM-Progression has 180 patients, pre-RT MRI, follow-ups (1 to 7 time points) and RTDOSE/RTPLAN for a subset [S](https://wiki.cancerimagingarchive.net/x/KgWwC). Cost: download.
3. Hybrid mechanistic-learning: a growth ODE plus a guided diffusion model, with radiotherapy effects in the ODE [S](https://arxiv.org/abs/2509.09610v1). Cost: compute.
4. Spatially varying parameters from Bayesian inference [S](https://arxiv.org/pdf/2209.12089). Cost: compute.
5. Stronger baselines before claiming a gain: linear volume extrapolation, per-patient growth-rate carry-forward [K]. Cost: now.

**L2. Historical forecast/shrinking subgroup.** The old 80/152 summary lacks its exact pair manifest here; later mask subgroup and volume-only negative results are documented separately in ledger 2w/2l.
1. Model treatment kill by dose (L1 fix 2). The historical kill-term worsening claim is unverified here without its exact saved contrast (the old note alone is insufficient), so use dose maps, not a global term.
2. Separate pseudoprogression and cavity collapse. LUMIERE has expert RANO ratings for 91 GBM patients, 638 study dates [S](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9755255/). Cost: download.
3. Better stratification. The reported script-79 growth-classifier AUC 0.46 remains unverified here without its exact output/key; the reported inputs were scan history and timing. Add MGMT, IDH, extent of resection. RHUH-GBM has molecular data, volumetric assessments and survival for patients with total or near-total resection [S](https://arxiv.org/pdf/2305.00005). Cost: download.
4. Report shrinking tumours as out of scope in the abstract. Cost: now.

**L3. Historical UCSF tensor study, no independent replication established.** The grade 2-3 eligibility claim needs its clinical manifest; result CSVs alone do not establish it. Retained 62-case tables and 42/20 split are described in `archive_old_gbm/vault/UCSF-PDGM-Tensor-Study.md`.
1. Replicate the tensor study on a second cohort. I did not find a second public longitudinal GBM set with DTI in these searches [S], so this may need a new data source. Cost: download, unknown.
2. Clean the tensors first. Free-water elimination corrected FA separated later-recurrence areas from recurrence-free edema (AUC 0.9 vs 0.77 uncorrected) [S](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12696142/). Apply to UCSF-PDGM, then rerun the three modes. Cost: now plus compute.
3. State "tensor-assumption shape test only"; verify grade eligibility before using it. Cost: now.

**L4. Only 21 real patients (`real_test`) in Track C.**
1. Fit growth rates on LUMIERE (91), RHUH-GBM and Burdenko (180) and make a larger real set. Cost: download plus compute.
2. No measured patient-specific kill rates were identified in the files reviewed for this note; an exhaustive search inventory is missing. Fit them from response during therapy in patients with dated treatment and volumes [K]. Cost: now to compute.
3. Report intervals for every real_test number. Cost: now.

**L5. In this volume-only fit/model, D is not identifiable; 124/154 rho fits sit at the lower bound.** Source: `output/inverse_est_metrics.json:summary` (`n_rho_at_lower_bound`, `n_patients_with_2plus_scans`, `n_D_at_bound=147`); this is not a universal identifiability theorem.
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
2. Report outcomes across a sweep of the scale, not one value. The original note reported only a three-patient sweep (unverified exact experiment key). Later cohort sensitivity at 500/1000/2000 is summarized in ledger 2j; none validates assumed patient kill rates. Cost: compute.

**L9. Legacy and 10x-higher diffusivity regimes differ.** The old claim of literature calibration is unverified without the source/units. Script 48 alternate-regime results now exist (`output/3d_extension_summary_swanson.json`, ledger 2h); the earlier not-run statement is superseded.
1. Rerun script 48 at D_white 0.13 and D_gray 0.013 (the Swanson-type values in `Script-78-Horizon-Crossover.md`). 8 patients, 50^3 grid, 3600 steps. Cost: compute, probably hours. Historical cost estimate only, not verified timing or current run authorization.

**L10. Script 47 dual-agent exposure exceeds the single-agent MTD reference.** In the selected eight simulations, total model-unit AUC is about 944-1280 vs 385; summing two drugs does not establish clinical dose equivalence. The matched two-drug MTD baseline has mean capped TTP 358.9 vs 305.6 days for dual adaptive (`output/dual_mtd_baseline.json`, ledger 2i).
1. Compare against a dual-agent MTD, not a single-agent MTD. Cost: now.
2. Put drug exposure in one unit for all arms. Cost: now.
3. Add a dose-matched arm, as in the equal-budget RL work. Cost: compute.

**L11. No clonal competition; resistance is a fixed fraction and partly plastic.**
1. Adaptive therapy needs competition between sensitive and resistant cells. In the Moffitt prostate trial (17 patients) the aim was to keep sensitive cells to suppress resistant ones [S](https://elifesciences.org/articles/76284). The trial used PSA as a trackable biomarker. A comparably useful GBM blood marker was not established by the sources reviewed here [K]; do not infer universal absence or a single main barrier.
2. Add history-dependent phenotype variables, calibrated to spheroid data with temozolomide [S](https://www.biorxiv.org/content/10.1101/2023.11.24.568421.full.pdf). Cost: compute. The repo's script 83 already has a two-state switching model, so this would be a richer version.
3. Add spatial competition (resistant and sensitive cells in the PDE, not a well-mixed Lotka-Volterra model) [K]. Cost: compute.
4. Serial ctDNA or biopsy data to calibrate resistance would be ideal. Public availability was not established in this limited search [K]; a dataset inventory is missing. Cost: no.

**L12. Anisotropy never matched to real tumour shape; tract alignment is lower under aniso.**
1. Check the alignment metric. Script 74 reports tract alignment as 0.4535 (aniso) vs 0.5011 (iso) [R], `output/fractal_aniso_vs_iso.json:secondary.tract_alignment_150d`. I did not read its definition. Confirm what it measures and why the isotropic value sits near 0.5. Cost: now.
2. Test on recurrence location. Standard planning is an isotropic margin, and 52% of patients progressed outside the 2 cm margin in one study; fibre-density-weighted maps improved cell-migration prediction [S](https://www.mate.polimi.it/biblioteca/add/qmox/17-2018.pdf). Needs DTI plus recurrence masks. Cost: download, unknown.
3. Compare to the Painter-Hillen formulation, which has a specific diffusion-tensor model for white-matter tracts [S](https://www.math.ualberta.ca/~thillen/paper/PainterHillen2013.pdf). Cost: now to compute.

## Evaluation limits
**L13. Historical tested RL policies improved some day-90 volume scores without a corresponding general TTP gain.** This is limited to the modeled cohorts/policies with assumed kill rates, not every possible RL method; ledger 2i documents the source-linked endpoint-specific comparisons.
1. Train and score on TTP. Script 68 started this (`dagger_ttp`). Extend it to every policy. Cost: compute.
2. Use survival-aware off-policy evaluation. A 2026 paper (arXiv 2603.22900) tailors OPE for right-censored outcomes [S](https://arxiv.org/pdf/2603.22900). Needs real treatment variation. Cost: download.
3. Offline RL work for chemotherapy dosing stresses safety constraints and heterogeneous test conditions [S](https://arxiv.org/pdf/2508.17212). Add a safety gate and report it. Cost: now.

**L14. Pacing gains +2.143 days in the real-growth-input simulation (21 patients, assumed kill rates).** Source: ledger 2i; exact policy output path/key must be recovered from its cited source before a new numeric claim is made here.
1. Find out why. Compare real_test kill rates and growth rates to the synthetic sets [K]. Cost: now.
2. Enlarge the real set (L4). Cost: download.
3. Report the synthetic sets as model checks, not evidence. Cost: now.

**L15. Probe-then-commit needs precise measurement.**
1. Automated volumetry is much more reproducible than manual RANO: ICC 0.97 vs 0.42 (1D) and 0.61 (2D) [S](https://arxiv.org/pdf/2209.01402). Use the measured volume CV as the noise level, not 5-10% [K].
2. Infer kill rates in a Bayesian way with an explicit noise model. Cost: compute.
3. In the assumed-noise simulation, longer probes (L = 3 and more) average noise; this is not empirically measured performance. Cost: now (already simulated).
4. A power analysis: noise level at which the rule's win rate falls below the blind policy. Cost: now.

**L16. The 0.02 threshold has 10 patients and was chosen after seeing the data.**
1. Fix the cutoff before testing on a new cohort (LUMIERE, Burdenko) [K]. Cost: download.
2. Fit a continuous model (adaptive earlier vs rho) with a confidence band, not a cutoff. Cost: now.
3. Report the sharp 10/10 vs 0/51 split as descriptive only. Cost: now.

**L17. Forecast target chosen after a volume probe.**
1. Historical proposal (current execution status/timing unverified without a run record): run the secondary: whole tumour without cavity (`--target wt`, script 81). Cost: compute (about 12 h per the vault).
2. Report all targets in one table. Cost: now.
3. Write future targets into a script header before any run (already the repo's habit). Cost: now.

## Negative-result limits
**L18. D_f is not valid.** Already replaced by elongation. To keep D_f, run box counting on larger tumours or finer grids [K]. Cost: compute.

**L19. MGMT does not predict TTP (n = 130, p = 0.60).**
1. Power analysis: what n would detect the 169 vs 181 d gap? [K] Cost: now.
2. More patients (L4 sets with MGMT). Cost: download.
3. Consider overall survival; the historical Hegi 2005 clinical claim remains unverified here without its exact title/DOI and endpoint. Do not cite this note as clinical evidence. Cost: now where OS exists.
4. Adjust for age and extent of resection. Cost: now to download.

**L20. Resistance-driven adaptive gain was not demonstrated on real growth inputs under the tested assumed-kill models (scripts 76, 82, 83).** This is a finding. Fixes are the L11 model upgrades. If three model classes keep failing, say "does not occur under these models" and stop.

**L21. Historical growth-classifier AUC 0.46 (unverified exact artifact/key; see L2).** Add dose, MGMT, IDH, extent of resection (L2 fix 3) or drop it. A failed classifier is a valid result.

## Not tested
**L22. No prospective validation.** Not possible here (no). Closest: temporal validation (train on early-scanned patients, test on later) and a pre-registered analysis plan for the next cohort [K]. State it as future work.

**L23. 3D alternate-diffusivity sensitivity (literature calibration unverified).** See L9.

**L24. Origin of the inflammation-score values (scripts 43-44).** [R] `get_patient_inflammation` reads a per-zone expression table by `patient_id`, and returns 1.0 (no effect) when the patient is absent. Check whether the MU-Glioma patient IDs match any expression rows. Later ledger A24 reports all 61 inflammation scores equal 1.0 and no effective A-to-B link. This is a recorded result, but its exact score artifact/key is missing here; it does not establish measured cross-track biology. Cost: now.

**L25. Drug screen is in silico only.**
1. Compare predicted targets to DepMap CRISPR essentiality in 60 patient-derived GBM lines [S](https://era.ed.ac.uk/items/2b20d86b-08e2-4ee8-9fee-951475f046a2) and to the wider DepMap resource [S](https://forum.depmap.org/t/validation-of-crispr-ko/1040). Cost: download.
2. **A caution [K, not searched]:** the top hit MT-CO2 is encoded in mitochondrial DNA. Nuclear CRISPR screens do not knock out mitochondrial genes, so DepMap may not contain it. Check how the "virtual knockout" treats it before claiming a target.
3. Check toxicity to normal tissue with healthy-brain expression. Cost: download.

## Original suggested order (2026-10-04; superseded, not current next steps)
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
"Weidner 2023": a matching paper was not located in two web searches; this does not establish that none exists. The nearest 2023 hits were "Predictive Digital Twin for Optimizing Patient-Specific Radiotherapy Regimens under Uncertainty in High-Grade Gliomas" [S](https://arxiv.org/pdf/2308.12429) and a Frontiers in AI paper from October 2023 [S](https://www.frontiersin.org/journals/artificial-intelligence/articles/10.3389/frai.2023.1222612/pdf). I do not know the authors. Do not cite "Weidner 2023" until you have the title or DOI.

## Status of fixes (initial snapshot 2026-10-04; corrections 2026-10-19)
| Limit | Fix | Script / file | Result |
|---|---|---|---|
| L0 Track A split | Patient-level 5-fold CV for LR, RF, scVI-LR; patient-lookup baseline | script 86, `output/patient_level_cv.json` | Keys `results.{random,patient}.{logreg_genes,rf_genes,logreg_scvi,patient_lookup_baseline}.acc.mean`; random split: LR 0.710, RF 0.714, scVI-LR 0.726, lookup 0.8275. Patient-level: LR 0.603, RF 0.620, scVI-LR 0.662, lookup 0.157. About 10 pp of inflation. 21 patients; Healthy has only 3 patients, so class and patient are confounded. RF used 150 trees (script 06 used 500). |
| L0 C-GAT | Leak-free C-GAT, 3 variants, patient folds, final epoch | script 88, `output/cgat_leak_free.json` | Completed retained result: clean PCA graph mean 63.8% (SD 13.3), cVAE without label edges 74.2%, published edges 75.8%; `variants.<variant>.acc.mean/sd`, ledger 2o. cVAE features retain upstream leakage/confounding; not all variants are leakage-free. |
| L9 3D diffusivity | Re-run at D 0.13 / 0.013 | script 48 with `GBM_D_REGIME=swanson` | Alternate-regime sensitivity changes selected high-rho sparing (0.678 to 0.346; 0.109 to 0.051), not general validation; ledger 2h and `output/3d_extension_summary_swanson.json`. |
| L10 dual-agent baseline | MTD plus the same secondary drug | script 89, `output/dual_mtd_baseline.json` | Matched two-drug MTD mean capped TTP 358.9 vs dual-adaptive 305.6 days (ledger 2i; model-unit AUC 1152.5 vs 1015.5). |
| L14 real vs synthetic | Compare parameters | `output/limit_checks.json:L14_real_vs_synthetic` | real_test kill/rho 77.7 vs 2.7-12.5; kill rates are the same for all 21 patients (assumed). |
| L16 threshold | Exact CI, permutation, leave-one-out | `output/limit_checks.json:L16_threshold` | 10/10 vs 0/51; p 5e-5; LOO 98.4%; still post hoc, n=10. |
| L19 MGMT power | Schoenfeld | `output/limit_checks.json:L19_mgmt_power` | Minimum detectable HR 1.74; underpowered. |
| L7 external survival check | CGGA, WHO IV, two batches | script 90, `output/cgga_validation.json` | Age direction replicates, effect weaker (HR 1.009 vs TCGA 1.031; Holm p 0.08). Sex and MGMT null. See ledger 2k. |
| L1/L3 new real patients | LUMIERE volume forecast | script 91, `output/lumiere_volume_forecast.json` | Growth-only model loses to no-change again (ledger 2l). Volumes only, no Dice. |
| L24 inflammation link | Checked IDs | ledger | Score is 1.0 for all 61 patients. A-to-B link = none. Schematic fixed. |

### Original not-done snapshot (superseded; not run authorization)
- **L17** script 81 whole-tumour target: about 12 h compute.
- **External data:** CGGA and LUMIERE tables are downloaded and used (`data/external/`, git-ignored). RHUH-GBM clinical CSV and NIfTI (2.9 GB) are downloaded, not yet analysed (n=40, no MGMT column). Burdenko is restricted, skipped. LUMIERE masks (30 GB zip) not downloaded.
- **L1.5** stronger forecast baselines, **L8** E_MAX sweep, **L15** power analysis, **L5** Bayesian calibration: not started.

Later status: all five C-GAT folds, E_MAX sensitivity, stronger forecast baselines, H-2/RHUH and H-3/UCSF scoring were completed (ledger 2o/2j/2w/2ax/2bd). H-2/H-3 are spent. Burdenko clinical CSV was checked (2az); image access remains separate. Remaining historical status/cost claims above require dated inventories before reuse.

### Extra notes to keep
- **L25 caution:** MT-CO2 is a mitochondrial gene. Coverage by a particular DepMap assay remains unverified without its target list; the mitochondrial-gene caveat alone does not check coverage. Check script 31 before claiming it as a target.
- **L12 alignment metric (script 74):** tract alignment is lower under aniso (0.45 vs 0.50). The cause is NOT checked (see [[Script-74-Fractal-Elongation]]). Do not say aniso "aligns with tracts". Say elongation rises and alignment does not.
