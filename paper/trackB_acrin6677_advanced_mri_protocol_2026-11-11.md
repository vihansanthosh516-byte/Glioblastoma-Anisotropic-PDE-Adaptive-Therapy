# ACRIN 6677 advanced-MRI increment protocol

**Status: frozen for a metadata-only feasibility audit.** Freeze date: 2026-11-11 (America/Chicago). The area forecast has already been evaluated on this cohort, so any advanced-MRI analysis will be exploratory and will not be an independent confirmation. No advanced-MRI measurement magnitudes are authorized by this protocol until the metadata audit meets the frozen feasibility rule below.

## Question

Among patients with three valid central-reader enhancing-area visits, does standardized relative cerebral blood volume (sRCBV) measured at the second area visit improve prediction of the third area beyond an area-only model?

This asks about incremental prediction of a radiographic measurement. It does not test treatment response, drug resistance, causal benefit, or clinical utility.

## Prior-art boundary

ACRIN 6677 publications have already studied DSC-MRI rCBV changes at weeks 2, 8, and 16 in relation to progression-free and overall survival, including pseudoresponse under bevacizumab. The proposed endpoint is next-visit enhancing-area forecast error, not survival or progression classification. This distinction is a narrow difference in endpoint; it does not establish novelty. The published DSC subset was small (37 of 123 had DSC imaging; 21 had interpretable baseline and at least one post-baseline scan), so feasibility is uncertain. Sources: Schmainda et al. (2015), https://pmc.ncbi.nlm.nih.gov/articles/PMC4490871/; Boxerman et al. (2023), https://pmc.ncbi.nlm.nih.gov/articles/PMC9909012/.

The primary predictor will be the SC form's standardized rCBV fields (`sce20`, `sce28`, `sce36`), selected from the data dictionary before inspecting SC rows. The paper describes sRCBV as a standardized intensity scale intended to improve comparability across scanner vendors, models, and field strengths. It is not an absolute blood-volume measurement, and the project will not assign it a physical unit or biological threshold.

## Metadata-only feasibility audit

The audit may read only these fields from the clinical tables:

- S0: patient key `cn`, nominal visit `s0e1`, day `s0e2d`, image-quality code `s0e5`, lesion-presence code `s0e13`, and nonempty/missing status of area field `s0e23`.
- SC: patient key `cn`, nominal visit `sce1`, day `sce2d`, reader key `sce3deid`, image-quality code `sce4`, and nonempty/missing status only for `sce20`, `sce28`, and `sce36`.

The audit will not parse, display, summarize, or model area or rCBV magnitudes. It will not open images, survival fields, progression fields, or treatment-arm fields. The result will contain aggregate counts, reader-count categories, file/schema checks, and a hash of the eligible patient-key set; it will not save patient-level rows or feature values.

## Frozen linkage and eligibility rules

1. Apply the area amendment's existing visit rules without reading area magnitudes: exclude S0 rows with `s0e5 = 1` (uninterpretable); a visit is present when `s0e23` is nonempty or `s0e13 = 1` (no lesion, therefore zero under the area amendment). Group same-patient reader rows by distinct day and order visits by day. Require at least three visits with strictly increasing days and use the first three. The second is the forecast-origin visit; the third is the target. Retain its nominal visit code for matching. If the valid reader rows on that day disagree on the visit code, that day cannot link to SC.
2. Link DSC to the second S0 visit only by exact patient key, exact nominal visit code, and exact integer day. No nearest-date substitution or time tolerance is allowed. A mismatch is missing linkage.
3. A DSC reader row is eligible only when `sce4 = 1` (adequate). It must have at least one nonempty sRCBV field among `sce20`, `sce28`, and `sce36`. Suboptimal or uninterpretable SC rows are not used. No nRCBV or other CBV field will substitute for missing sRCBV.
4. If more than one adequate reader row is available at the exact linked visit, each reader's sRCBV feature will be the arithmetic mean of that reader's nonmissing values across the three sRCBV lesion slots; the patient feature will be the arithmetic mean across eligible readers. The audit will record only whether each slot is present and how many readers qualify, not the values.
5. The feasibility threshold is **at least 25 unique patients** with three valid S0 visits and an exact, adequate, nonempty sRCBV link at visit 2. This reuses the project's pre-existing minimum cohort size; it is a feasibility gate, not a power calculation. If fewer than 25 qualify, stop before reading any area or advanced-MRI measurement magnitudes. Do not relax the gate or change the matching rule after seeing the count.

If the required columns, time units, or form-to-field mapping do not match the dictionaries, stop and record the discrepancy rather than infer a mapping.

## Analysis, only if the feasibility gate passes

- Use one patient and one forecast window per person: the first three valid S0 area visits under the area amendment. Use `y = log(area + 0.1 cm²)` as already frozen there.
- Compare two patient-level models on the same eligible patients: (1) ridge regression with intercept, using only the second-visit log area to predict third-visit log area; (2) the same model plus the second-visit sRCBV feature. The ridge penalty is fixed at `alpha = 1`; there is no hyperparameter search. Standardize predictors using training patients only.
- Generate predictions with leave-one-patient-out cross-validation. No row from a held-out patient may enter preprocessing or fitting. Report the existing copy-last baseline separately as a reference.
- Primary metric: median absolute log error. Report the paired per-patient error difference (area+sRCBV minus area-only), its median and range, and a paired patient bootstrap interval with 10,000 resamples and seed `20261112`. The interval describes the observed cross-validated error differences only; because leave-one-out training sets overlap, it is not a formal population-generalization interval. Do not use a p-value or a pass/fail claim. Report the exact n and all exclusions.
- There is no imputation, no feature substitution, and no threshold tuning. If an sRCBV value is missing, malformed, or non-finite, that patient is not in the paired model set; report the count and reason. The model is not run if fewer than 25 patients remain.
- DCE (`DE`) and spectroscopy (`MS`) are out of scope. They will not be inspected as a way to replace a failed DSC feasibility gate.

## Interpretation limits

Any observed difference will be exploratory because the same cohort's area outcomes have already been examined. It will not be called external validation, a new treatment-response biomarker, drug resistance, evidence that perfusion causes area change, or evidence of clinical benefit. The ACRIN publications already establish that rCBV has been studied as a survival/progression marker in this cohort; a favorable result here would support only a hypothesis for a separate, untouched cohort.

**Weakest point:** the optional DSC sample is small and selected, and bevacizumab can change enhancement through vascular effects. The one-window, low-dimensional model can still be unstable; leave-one-patient-out estimates do not create an independent test set.
