# Integration protocol: Biology-Anchored MRI Decomposition (BAMD) — locked 2026-11-12

Written and committed before any survival value, expression value, or image intensity is read for this study. Data inventory: `paper/integration_data_inventory_2026-11-12.md`. Simulation is not involved; this is a retrospective analysis of public data. No clinical claim.

## Question

Can routine pre-operative MRI carry tumor-cell-state biology into patients who have no RNA, and does the biology-linked part of the MRI predict survival better than the rest of the MRI and the usual baselines?

## Method (BAMD)

1. **Anchor (TCGA, outcome-free).** In TCGA-GBM patients with BraTS 2021 MRI + expert masks + U133A RNA (expected n about 91), compute four Neftel cell-state scores from RNA and learn a ridge map from fixed MRI features to those scores. No survival data is used in this step.
2. **Transfer.** Apply the frozen map to UPENN-GBM (no RNA). This gives MRI-predicted cell states.
3. **Decompose.** Split the MRI features into two parts: the **biology-anchored part** (the predicted states, i.e. the projection onto the map's output) and the **residual part** (the first 3 principal components of the MRI features after removing their linear projection on the predicted states, with the PCA fitted on TCGA and frozen).
4. **Test.** In UPENN, ask which part adds survival information beyond clinical factors and tumor volumes.

What is new, after the prior-art check (`paper/integration_prior_art_2026-11-12.md`): MRI-to-Neftel-state mapping learned from bulk RNA, transferred to an RNA-free external cohort, plus a split of MRI into biology-anchored and residual parts with a pre-registered survival test of each. Closest prior work: Matsumoto et al. 2025 (Neuro-Oncology; MRI predicts a methylation-derived stem-like cell proportion and survival, UPenn data), Sun et al. 2021 (Radiology; MRI to RNA pathways), Braman et al. 2021 (orthogonal fusion in learned embeddings). Not verified: full texts of several related papers were blocked; a human expert has not reviewed the novelty claim.

## Data and inclusion

- **TCGA anchor set:** BraTS 2021 training cases mapped to TCGA-GBM (`data/external/brats2021_map/BraTS2021_MappingToTCIA.xlsx`) that are on disk with masks and present in `data/external/tcga_gbm_u133a/HT_HG-U133A.gz` (first sample per patient, primary tumor `-01` preferred).
- **UPENN test set:** BraTS 2021 cases mapped to UPENN-GBM, on disk with masks, joined to `data/external/upenn_meta/clin.csv` by patient ID (`UPENN-GBM-xxxxx`). Exclude IDH1 mutant. Exclude overall survival missing or not positive. Censoring from `Survival_Censor` (column semantics checked from the UPENN data descriptor before the run and recorded).
- External MRI+RNA sets (CPTAC-GBM, IvyGAP) are not used in the primary analysis.

## Features (fixed list, computed from BraTS images and masks)

Regions: enhancing tumor (label 4), necrotic core (label 1), edema (label 2), whole tumor. Images: T1, T1ce, T2, FLAIR, each z-scored within the brain mask (nonzero voxels).
- Volumes (cm^3, log1p): ET, NCR, ED, WT. Ratios: NCR/(NCR+ET), ET/WT, ED/WT.
- Shape (whole tumor and tumor core): surface-to-volume ratio, sphericity.
- Intensity per region (ET, NCR, ED) x image (4): mean, SD, 10th and 90th percentile = 48.
Total: 4 + 3 + 4 + 48 = 59 features. Missing region (zero volume): intensity features set to the training median, plus the volume itself carries the zero.

## Cell-state scores (TCGA RNA)

Gene lists: Neftel et al. 2019 IDH-wt meta-modules (`data/external/neftel2019/IDHwt.GBM.MetaModules.tsv`, from the bardylab repository, SHA-256 046905b6...). States: MES = MESlike1 + MESlike2, AC, OPC, NPC = NPClike1 + NPClike2 (cell-cycle modules not used). Score = mean z-score (across the TCGA anchor set) of module genes present on the array, minus the mean z-score of all genes. Weakest point: modules come from single cells; bulk scores mix tumor and non-tumor cells.

## Gate 1 (outcome-free; decides whether to continue)

Ridge map, alpha chosen by inner 5-fold CV from {0.1, 1, 10, 100, 1000}, features standardized on training folds. Outer: 10 repeats of 5-fold CV. For each state, Spearman rho between out-of-fold prediction and RNA score (averaged over repeats); permutation p from 1,000 label shuffles of the whole procedure (states shuffled jointly by patient). A state **passes** if mean rho >= 0.25 and Holm-adjusted p < 0.05 over the 4 states.
- If no state passes: stop; report "MRI does not recover cell states at this sample size". No survival analysis of predicted states.
- Passed states only enter the biology-anchored part.

## Primary test (UPENN; survival read once)

Base model B0: Cox with age, sex, gross total resection (yes / no / unknown indicator), log1p ET, NCR, ED volumes.
- **H1 (primary):** B0 + passed predicted states vs B0.
- **H2 (secondary):** B0 + 3 residual PCs vs B0.
- **H3 (secondary):** B0 + both vs B0.
Score: Harrell's C from 10 repeats of 5-fold CV in UPENN (Cox coefficients fitted on training folds only; map and PCA frozen from TCGA). Delta C = C(model) - C(B0), paired within folds. 95% CI by bootstrap over patients (1,000 resamples of the whole CV). Also the likelihood-ratio test on the full UPENN set.
**H1 passes** if the delta-C 95% CI lower bound > 0 **and** LR-test p < 0.05. Everything else is reported, with no re-choice of features, states or models.

## Benchmarks (same UPENN CV and score; all reported)

1. Clinical Cox (age, sex, resection).
2. B0 (clinical + volumes).
3. LASSO Cox on clinical + all 59 MRI features (penalty by inner CV).
4. Ridge Cox on the same.
5. DeepSurv-style MLP (torch; 1 hidden layer of 32, dropout 0.2, Cox partial-likelihood loss, 200 epochs, weight decay 1e-3; fixed, not tuned).
6. Random survival forest and gradient-boosted survival: require the scikit-survival package; if it cannot be installed, these are listed as **not run**.
7. BAMD (B0 + predicted states + residual PCs).
BAMD is called "better" only if its delta-C versus every baseline that ran has a 95% CI above 0. Otherwise the report says which baseline matched or beat it.

## Secondary (TCGA, exploratory, exposed)

TCGA survival was analysed before in this project (clinical covariates, 4 genes). In the TCGA anchor set: Cox C-index (5-fold CV) for RNA states vs MRI-predicted states vs both, labelled exploratory. No pass/fail.

## Order of work

1. Commit this protocol. 2. Feature extraction and RNA scores with synthetic tests. 3. Gate 1 (no survival opened). 4. Commit the frozen map and PCA with hashes. 5. Open UPENN survival once; run primary, secondary and benchmarks. 6. Report everything, including failure.

## Known weaknesses

TCGA anchor set is small (about 91); the MRI-to-state map may be weak. TCGA and UPENN differ in scanners and eras (BraTS harmonized geometry, not intensity; z-scoring per image only partly fixes this). OS is affected by treatment that differs by cohort. Neftel modules are from single cells. One external cohort only.
