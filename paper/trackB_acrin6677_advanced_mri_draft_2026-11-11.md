# Exploratory advanced-MRI forecasting protocol draft (2026-11-11)

**Draft only. No advanced-MRI rows inspected and no evaluation run.** The area result already informed this follow-up, so reuse of this cohort is exploratory. [Area result](trackB_acrin6677_result_2026-11-11.md)

Question. Does DSC perfusion add information to the last enhancing area when forecasting the next measured area? Use SC as the primary form for this draft. Consider DCE permeability from DE and spectroscopy from MS as separately labeled descriptive extensions, subject to metadata coverage and documented units.

Before fitting, audit the dictionaries for identifiers, units, QC codes, and whether acquisition dates can link to the S0 history visits. Audit overlap using presence and timing only. Freeze the exact predictor column, matching rule, QC exclusions, minimum eligible count, and missing-data rules before reading advanced-MRI values. If a usable linkage or unit is missing, stop and document the gap.

Use the area amendment's window and transform. Use only features acquired on or before the second measured visit. Compare an area-only model and an otherwise identical area-plus-perfusion model on the same eligible patients, with B1 as a reference. Fit preprocessing and model coefficients within patient-level training folds only. Record paired absolute error changes and uncertainty, along with overlap and exclusions. Treat findings as exploratory. A later independent cohort would be needed for a confirmatory test.

Weakest point. The advanced-imaging subset may be small and selected. Perfusion and enhancement both reflect vascular effects under bevacizumab. [Same-trial perfusion study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4490871/)
