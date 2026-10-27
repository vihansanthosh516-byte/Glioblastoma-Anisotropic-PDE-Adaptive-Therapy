# Track B UCSD exploratory amendment — 2026-10-27

**Status:** Freeze this amendment and its runner in Git before reading any further UCSD target-mask voxel arrays. Keep `track_b_ucsd_external_test_protocol.md` unchanged. This amendment defines one exploratory analysis after the prior locked run failed.

## Why this analysis is exploratory

The locked runner created `output/track_b_ucsd_locked_test/analysis_started.json`, read target labels for at least five patient pairs, and stopped during pair 5 when origin T2 normalization failed. Four earlier pairs had been processed. No aggregate result was written, and partial in-memory metrics were not saved. Preserve that marker and all partial-run evidence. The UCSD outcomes are now partially exposed, so this analysis cannot serve as a fresh confirmatory test, regardless of the code freeze below.

## Question and endpoint

Keep the original question and endpoint: from one forecast-origin scan, rank newly annotated NETC+ET voxels at the next available scan. This measures change in source annotations. It does not measure biological tumor growth, recurrence, clinical progression, treatment response, or a fixed-day forecast.

## T2-quality failure policy

Use a complete-case technical-failure policy. Before reading any further UCSD segmentation voxel array, screen each selected origin using its four MRI arrays alone. Require finite, usable robust intensity normalization for every standard modality over the union of nonzero MRI support. Exclude an origin if any modality fails. Do not use current-mask voxels to rescue an origin that fails this MRI-only screen. Do not impute a T2 feature or add a fallback model.

Report all selected aligned pairs, each origin-quality failure, each modality failure, each pair scored, and each metric-evaluable pair. The five origins failing the prior preflight remain technical failures; analyze only origins that pass the frozen screen. This rule depends on forecast-origin MRI inputs, not future annotations. It can still select a nonrepresentative subset if image-quality failure relates to patient or acquisition factors. State this limitation beside all results.

A missing-T2 fallback would keep more cases but change predictor inputs and require a separately developed and validated model on MU. The complete-case rule keeps the existing model fixed and avoids outcome-guided imputation, but reduces the cohort and may select a biased subset. Do not choose between policy options using UCSD outcomes.

## Frozen analysis

- Use the first complete, aligned adjacent pair per patient under the unchanged protocol.
- Use the same origin-only graph-diffusion prior, seven MRI/distance features, fixed boundary-sensitivity gate, MU-Glioma-Post training, and five-fold patient-grouped development comparator selection as the locked runner.
- Keep the next available scan as target, regardless of gap. Do not supply the gap to the model.
- Do not tune or select any model, baseline, threshold, or exclusion using UCSD outcomes.
- Report each evaluable patient’s paired AUC difference, aggregate patient-macro AUCs, a 10,000-draw patient bootstrap interval with seed 20261126, and descriptive strata at ≤90, 91–180, and >180 days. Report positive and negative risk-set voxel counts, zero-addition patients, additions outside the risk set, and technical failures. Voxels are not inferential units.
- Treat the original 0.02 AUC margin only as a descriptive reference. Do not claim that this exploratory analysis meets or fails a confirmatory success criterion.
- Write the marker and result only under the ignored `output/track_b_ucsd_exploratory_20261027/` directory. Never alter the original locked-test marker or output. Do not write subject IDs or voxel predictions to the result.

## Published-method comparison and limits

Related published methods predict longitudinal segmentation, future MRI, or tumor growth, but their inputs and endpoints do not match this one-origin annotation-addition ranking task. TaDiff conditions future MRI and masks on sequential scans and treatment information ([Liu et al.](https://arxiv.org/abs/2309.05406)). A feature/label-fusion study predicts later segmentations using growth-density, texture, and label-fusion components and evaluates nine longitudinal BraTS 2015 patients ([Chen et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8336640/)). Earlier patient-specific reaction-diffusion forecasting used 17 intervals from 10 patients ([Miller et al.](https://pubmed.ncbi.nlm.nih.gov/35603364/)); a prior MRI prediction method used five timepoints from four GBM patients ([Weizman et al.](https://pubmed.ncbi.nlm.nih.gov/23286047/)).

These methods do not provide a directly comparable, ready-to-run baseline for a single origin, an unknown next-visit interval, this post-treatment cohort, and patient-macro AUC of source-label additions. Their reported Dice, volume overlap, or future-MRI metrics are not comparable to this endpoint. The current runner therefore retains its MU-selected comparator and internal ablations. Do not claim superiority over these published methods. A fair reimplementation would require a separate development-only protocol and method-specific data, inputs, and evaluation.

## Weakest points

The UCSD outcomes were partially exposed before this amendment. The MRI-only screen is stricter than normalization over the full support that also includes current labels, so it may exclude cases that could otherwise normalize successfully. Only quality-passing cases will be scored, and the quality screen may select a nonrepresentative subset. The endpoint remains a change in source annotations; variable intervals, treatment, registration, and annotation variation can all affect it. The bootstrap interval may be wide. The published comparisons are contextual, not head-to-head.
