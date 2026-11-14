# Track A benchmark readiness

## Verified components

- [Specimen manifest audit](../src/trackA_margin_manifest.py) detects duplicate specimens, repeated origins across studies, development/test overlap, unpaired regions, and declared expression-derived region labels. It checks the donor floor and the best possible Holm result before scoring. Real label provenance and biological origin mapping still require source verification.
- [Shared evaluator](../src/trackA_margin_evaluation.py) requires the same donor IDs for every frozen arm. It computes concordance, per-donor differences, exact sign tests, Holm adjustment, and the frozen pass rule.
- [Baseline adapters](../src/trackA_margin_baselines.py) apply fixed weights with a gene coverage gate and fit a paired centroid on development donors. They do not select the published signatures, transform expression, or establish valid assay coverage.

Synthetic checks verified these components individually and in a combined run from a simulated manifest through candidate and baseline scoring. In the combined fixture, the candidate tied the fixed score and centroid and correctly failed the overall pass rule. This is an implementation check. The separate [frozen simulation result](trackA_margin_simulation_result_2026-11-13.md) remains negative. No real outcomes were opened.

## Pending before the patient run

The [source metadata audit](trackA_source_metadata_2026-11-13.md) records the pending Harwood reply. Verify that the available label set comes from pathology and joins to coded donors and malignant-cell annotations. Confirm that the released assay, state coverage, and independent donor count can support the frozen benchmark.

Freeze the assay universe, independent development studies, pseudobulk adapter, published comparator definitions, scoring transforms, tuning grid, and run hashes before test access. The current adapters cover only part of the [comparator inventory](trackA_regional_prior_art_2026-11-13.md). The weakest-state proposal still lacks evidence of added value over its mean-state ablation. Further method work must use development data or a separately prespecified simulation.
