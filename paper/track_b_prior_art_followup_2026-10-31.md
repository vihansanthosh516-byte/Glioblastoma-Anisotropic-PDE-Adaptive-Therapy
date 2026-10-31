# Track B B2.1 prior-art follow-up (2026-10-31)

Follow-up to `paper/track_b_prior_art_2026-10-31.md`. Desk research only. No patient data was scored, no restricted file was accessed, no model was run.

## Evidence levels used

- **FT** = full text read (arXiv HTML), numbers taken from tables.
- **AB** = abstract page only.
- **SN** = search snippet only.
- **API** = GitHub API metadata (`gh api repos/<repo>`), 2026-10-31.

"Not found" means these searches did not find it. It does not mean it does not exist. Search was a general web index plus arXiv HTML. PubMed and Google Scholar were not queried directly. That is the main limit of every check below.

## Corrections to the first report

1. **arXiv 2509.18998 is not a patient study.** It is Schenk, Ayensa Jiménez and Romero, a Bayesian calibration of cell-migration models in **microfluidic devices** (AB, https://arxiv.org/abs/2509.18998). The first report called it "glioblastoma progression". It is in-vitro. Its discrepancy result is not quantified on the abstract page.
2. **PREDICT-GBM cohort is 243, not 255.** 243 = TUM-GBM 142 + LUMIERE 61 + RHUH 40 (FT, https://arxiv.org/html/2509.13360). No "255" appears in the text. The U-Net was trained on an extended set of 327 (243 + 84 from CPTAC-GBM, IVY GAP, TCGA-GBM, UPENN-GBM). Unresolved: the nnU-Net section gives "169 training and 84 test cases per fold" (sum 253).
3. **SoC-DT is closer than the first report said.** It has real patient results and a PDE-plus-neural-corrector baseline (see check 3).
4. **TaDiff Dice is not plain Dice.** It is computed "with an optimized threshold" on slices with tumor ≥ 11 cm² only (FT, https://arxiv.org/html/2309.05406).

## Check 1. PREDICT-GBM: are U-Net and GliODIL on a comparable metric?

- Search: arXiv HTML 2509.13360, main text and Supplement 9.1. Evidence: FT.
- **Metric:** recurrence coverage = |V_target ∩ V_rec| / |V_rec|. The paper says it is "equivalent to the sensitivity". All plans are "iso-volumetric to the standard plan". The standard plan is the tumor core dilated 15 mm inside the brain mask.
- So **U-Net (79.37 ± 2.08) and GliODIL (78.91 ± 2.08) are on the same metric** (pooled enhancing coverage, N = 243, mean ± SE). SOC is 77.34 ± 2.17.
- **This is fixed-volume sensitivity, not Dice.** Plan volume is set equal to the SOC volume, so the metric cannot reward a correct size. It must not be compared with forecast Dice.
- **Supplement 9.1 Dice:** U-Net mean 10.67 ± 0.78 vs SOC 10.34 ± 0.77 (printed as given; likely percent). GliODIL Dice is described only as "slightly higher" median and mean; no number. Dice is low because the plan is a large radiation volume, not a tumor forecast.
- Unresolved: the text does not say outright that the U-Net and GliODIL test sets are identical. The U-Net predictions are out-of-sample on held-out splits.
- Published version: npj Digital Medicine, DOI 10.1038/s41746-026-03194-0 (cited in `paper/track_b_external_cohort_screen_2026-10-28.md`; not re-read here).

## Check 2. Plain Dice against a copy-last-mask baseline, by horizon

Search terms: "glioma longitudinal tumor growth prediction Dice persistence / no-change / copy baseline time interval" (extended web search), plus full-text reads.

| Work | Evidence | Persistence baseline? | Plain Dice? | By horizon? |
|---|---|---|---|---|
| Jung & Shin, arXiv 2609.12435 (OASA) | FT | **Yes** | No: Dice on thresholded smoothed proxy maps (τ = 0.2, 0.5) | No; one aggregate over 15 test triplets |
| TaDiff, arXiv 2309.05406 | FT | No (own ablations only) | No: optimized threshold, slices ≥ 11 cm² | Yes, by treatment-day bin, model only |
| Petersen 2021, arXiv 2106.12917 | FT | Only as a figure line (mean true overlap 0.67) | No: Query Volume Dice picks the sample closest to true future volume | No; binned by true overlap, not time |
| SoC-DT, arXiv 2510.03287 | FT | No | Yes (DSC) | No |
| PREDICT-GBM, arXiv 2509.13360 | FT | No (SOC margin only) | Supplement only, on radiation plans | No |

OASA numbers (15 held-out test triplets, 3 seeds): Dice@0.2 persistence 0.6070, OASA 0.6071 ± 0.0025, dilation 0.6139. Dice@0.5 persistence 0.3981, OASA 0.4269 ± 0.0079. The abstract says CIs vs persistence and dilation include zero for Dice. Weakest point: n = 15, proxy target.

**Finding:** no glioma forecaster found reports plain annotated-mask Dice next to a copy-last-mask baseline, split by horizon. Not found ≠ does not exist. A useful fact: in OASA, simple **dilation of the last mask beat persistence and the learned model** at τ = 0.2. So a fixed dilation baseline belongs in the B2.1 comparison.

## Check 3. Cohort overlap with this project

- **PREDICT-GBM:** all three cohorts overlap prior work here. `START_HERE.md:218` and `paper/track_b_external_cohort_screen_2026-10-28.md:19-21` record that the local PREDICT-GBM work used the 142 TUM patients and that LUMIERE and RHUH were examined. Local outputs exist (`output/predictgbm/final/test_results.json`; scripts `src/111`–`src/123`). **PREDICT-GBM cannot serve as fresh confirmation.**
- **OASA (2609.12435):** uses MU-Glioma-Post from TCIA (FT). MU was examined here (`START_HERE.md:50`). It is a proxy-map forecaster, not our method. It does not overlap our method, but its cohort overlaps our data.
- **SoC-DT (2510.03287):** uses UCSF-ALPTDG (FT). UCSF was the spent H-3 test set (`src/133_h3_ucsf_test_score.py`). Patient count not stated.
- Weakest point: overlap is judged by dataset name. No patient-ID match was done, and none is needed, because whole datasets are already exposed.

## Check 4. Training-only count baseline and a fair mask for Dice

- **A count has no location.** A predicted number of new additions cannot be turned into a Dice mask without a location rule. Any rule that uses the true future volume is banned.
- **Defensible rule (proposal, not tested):** "training-only spatial baseline" = last observed mask ∪ voxels whose training-fold addition frequency, mapped by distance from the last core (for example 1 mm shells, optionally inside a WM/GM mask), passes a threshold. Both the shell frequencies and the threshold are fit on training folds only, per horizon bin. The mask size then comes from training data, not from the test future.
- A simpler version of the same idea: a fixed dilation radius tuned on training folds. OASA shows such a rule can be competitive (check 2).
- **Limit:** this baseline is a spatial rule, not a count model. Its Dice must not be called a count result. The count model is scored on count metrics only (for example Poisson or negative-binomial deviance, MAE of counts, calibration of P(count ≥ 1)). In MS, negative binomial with a log follow-up-time offset is the standard count model (SN, NCT04879628 SAP; PMC2631439; PMC3197685).
- If the shell rule cannot be fit stably on the training folds (too few additions), say so and drop it. Do not invent a location rule after seeing test outcomes.

## Check 5. Count prediction of new glioma lesions

- Search: "glioma 'new lesion' count prediction Poisson OR 'negative binomial' MRI follow-up" (extended). Evidence: SN.
- Hits were MS trials and MS methods papers. No glioma count model was found.
- A glioma review notes that treatment change vs true progression is hard to separate (SN, PMC9367286). So a glioma "new lesion" label is noisier than in MS.
- Not done: direct PubMed and Scholar queries; citation sweeps of LUMIERE, MU-Glioma-Post and UCSF-ALPTDG.

## Check 6. PDE plus learned correction: full texts and new hits

| Work | Evidence | Data | Key fact | Code (API) |
|---|---|---|---|---|
| SoC-DT, Bhattacharya, Singh, Prasanna, arXiv 2510.03287 | FT | Synthetic + UCSF-ALPTDG; 3 timepoints (baseline, first post-treatment, predict second) | DSC, 5-fold patient-wise CV: U-Net 53.09 ± 2.96, Fisher PDE 52.75 ± 3.47, PINN 51.71 ± 4.61, **Hybrid PDE (PDE prior + neural correctors) 54.29 ± 3.14**, SoC-DT 55.45 ± 8.12 | "upon acceptance"; none |
| Liu, Trimboli, Li, arXiv 2607.13877 | AB | 387 synthetic trajectories only | −84.3% masked MSE, +43.5% Dice vs RD; no residual penalty stated | none on page |
| TumorFlow, arXiv 2603.04058 | SN | Pre-op GBM MRI | 75% Dice is vs the biophysical model, not vs follow-up scans | github.com/valentin-biller/lgm (not checked) |
| UDE / neural ODE search | SN | — | arXiv 2511.22292 is Gompertz ODE, not glioma or spatial | — |
| GLIA (Subramanian, Scheufele, Biros) | SN | Single mpMRI inversion | Solver repo github.com/ShashankSubramanian/GLIA; no perturbation-robustness study found | not checked |

**Most important new fact:** SoC-DT already reports a **PDE-prior-plus-neural-corrector baseline on real patients** that edges a U-Net by 1.2 Dice points. The SoC-DT gain over that hybrid (1.16 points) is smaller than its own fold SD (8.12). No significance test on UCSF. No persistence baseline. So "PDE prior + learned correction beats a U-Net" is **not open** as a claim. Any B2.1 gap must lie in the narrower parts: training-only expected additions, stability-shrunk residual, and plain-Dice-by-horizon scoring against persistence.

"TumorGrowthToolkit" was not found under that name (search and `gh api repos/KitwareMedical/TumorGrowthToolkit` → 404). Jaroudi BMCMC and Rockne PI uncertainty were not reached in this pass.

## Check 7. Repos and licenses

From `gh api repos/<repo>` on 2026-10-31 (API):

| Repo | License | Last push |
|---|---|---|
| m1balcerak/GliODIL | MIT | 2026-06-24 |
| BrainLesion/PredictGBM | Apache-2.0 | 2026-10-28 |
| MIC-DKFZ/deep-glioma-growth | MIT | 2021-11-30 |
| jsudg436/longitudinal-proxy-forecasting (OASA) | MIT | 2026-07-02 |

PREDICT-GBM data is on Hugging Face (`LZimmer/PREDICT-GBM`); its data license is not stated in the paper text. The paper itself is CC BY-NC-SA 4.0.

## What this changes for B2.1

- The bar stays fixed: same held-out patients and horizon bins; compare persistence, PDE-only, fold-trained U-Net, and the training-only baseline; report plain Dice by horizon; never use true future size.
- Add a training-tuned **fixed dilation** baseline. OASA shows it can beat persistence and learned models.
- Treat **SoC-DT's Hybrid PDE** as prior art for "PDE + neural correction". Do not claim that combination as new.
- A count model is scored on count metrics only. Its spatial twin (check 4) is a separate baseline.
- No public cohort found here is fresh for this project. PREDICT-GBM, MU and UCSF are exposed.
- Any remaining gap is conditional on searches not yet done (PubMed, Scholar, citation sweeps). No "first" claim.

## Still open

- PubMed and Scholar direct queries for glioma count models.
- Citation sweeps of LUMIERE, MU-Glioma-Post, UCSF-ALPTDG and PREDICT-GBM.
- Full text of 2607.13877 (residual penalty?) and TumorFlow.
- Jaroudi BMCMC, Rockne PI uncertainty, GLIA robustness docs.
- GliODIL Dice value in PREDICT-GBM Supplement 9.1 (only a qualitative statement was found).
