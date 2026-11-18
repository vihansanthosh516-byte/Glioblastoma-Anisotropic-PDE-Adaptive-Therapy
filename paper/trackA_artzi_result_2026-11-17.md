# Track A Artzi benchmark result

**Exploratory. MP-RPV failed the frozen pass rule.**

The benchmark completed for all 10 paired patients and all 28 eligible comparators. MP-RPV had mean paired concordance of 0.50, equal to the always-tie reference. Twelve comparators reached 1.00, including paired conditional logistic regression and ridge logistic regression. Every Holm-adjusted p-value for MP-RPV superiority was 1.00. The result file records all comparisons, eligibility decisions, and patient scores. [Result](../output/trackA_artzi_region_benchmark_20261116/result.json) | [Patient scores](../output/trackA_artzi_region_benchmark_20261116/patient_scores.tsv) | [Frozen plan](trackA_artzi_preanalysis_2026-11-16.md)

The worst-patient filter required each selected pair to have a strictly positive signed contrast for every training patient. No pair qualified in any fold. The frozen fallback therefore returned equal core and infiltrative scores for every patient. This explains the candidate failure. The comparator results describe whole-tissue regional expression within the Artzi study. [Fold details](../output/trackA_artzi_region_benchmark_20261116/model_details.json)

Scikit-learn emitted a deprecation warning for the logistic regression penalty argument and flagged the combination of L1 penalty and the default L1 ratio. The benchmark reported no arm failures. Results reflect the locked implementation and software recorded in the manifest. [Manifest](trackA_artzi_benchmark_manifest.json) | [Result](../output/trackA_artzi_region_benchmark_20261116/result.json)

Artzi is now spent. The next candidate can replace the all-patient veto with a prespecified support threshold, develop it on simulation, and evaluate it on an untouched cohort with region labels and an exact patient linkage. Harwood, Gill, and Nam remain leads pending their metadata gates. [Cohort follow-up](trackA_fresh_data_followup_2026-11-13.md)
