# Track A development implementation and results

2026-11-23. Implemented the two candidates and ran generated cases to check their mechanics. Both candidates fell short of their strongest comparison in these cases. No real expression or assay outcomes were opened under the [roadmap](trackA_completion_plan.md). The numerical settings below belong to this development run.

## Regional prediction

[Regional model](../src/trackA_regional_model.py) fits paired logistic margins using one canonical specimen contrast per patient. It fits composition, density and depth, malignant expression, their combined nonspatial model and separate spatial baselines. The candidate uses a shared fit across the nominal and alternative reconstructions. It uses the nominal margin only when the independent audit qualifies the bound and both uncertainty intervals certify the same sign. Otherwise it uses the full nonspatial margin. The agreement-only ablation uses the same reconstructions and fitted coefficients.

The development run used 160 generated training pairs, 160 independent calibration pairs and 160 evaluation pairs. The bound used the maximum absolute error over both reconstructions for each calibration pair. The generated reference gave direct access to a clean margin. Real reference coverage still needs to be established. The strongest baseline column records the best observed noncandidate arm in each scenario.

| Generated condition | Reference coverage | Spatial use | Candidate concordance | Nonspatial | Agreement only | Strongest baseline |
|---|---:|---:|---:|---:|---:|---:|
| Clean | 0.9563 | 0.9250 | 0.9125 | 0.6813 | 0.9250 | 0.9313 |
| Divergent error | 0.9750 | 0.0000 | 0.6813 | 0.6813 | 0.6875 | 1.0000 |
| Shared error | 0.9625 | 0.0563 | 0.6438 | 0.6813 | 0.2188 | 0.8000 |
| Unqualified shift | 0.9625 | 0.0000 | 0.6813 | 0.6813 | 0.2188 | 0.8000 |
| Mismatched calibration | 0.0000 | 0.8875 | 0.2437 | 0.6813 | 0.2188 | 0.8000 |

**Failure diagnosis.** With shared error, the gate used spatial evidence on nine pairs and worsened six relative to fallback, with zero improvements. Marginal calibration coverage of 96.25% was insufficient to protect the accepted spatial decisions. A clean calibration set also missed the shared shift entirely. This is a calibration and selection problem under the simulated conditions. The test used a generated clean reference and stylized offsets. An actual follow-up needs a reference that can measure the supported error on real development sections and a prospectively specified assessment of errors among decisions that change after adding spatial evidence.

## Measured TMZ viability

[TMZ model](../src/trackA_tmz_model.py) fits an origin-balanced loss on mean viability and a shared-slope contrast term for independently derived sister cultures. It checks origin and culture disjointness, fits transformations on development rows, and can predict before receiving evaluation viability. Evaluation gives each biological origin equal weight. Baselines include source and matched random signature scores, MGMT, culture-weighted and origin-balanced ridge, the loss with the within-origin term removed, and a shuffle that changes sister links only in the contrast term.

| Generated condition | Candidate origin MAE | Origin mean only | Shuffled links | Origin-balanced ridge | Strongest baseline |
|---|---:|---:|---:|---:|---|
| Signal | 0.5072 | 0.5465 | 0.5018 | 0.5103 | Culture-weighted ridge, 0.5013 |
| Assay shift | 0.4844 | 0.5166 | 0.4804 | 0.4853 | Shuffled links, 0.4804 |
| Null | 0.5399 | 0.5616 | 0.5664 | 0.5421 | Matched signature, 0.4992 |

**Failure diagnosis.** The contrast term improved on the origin-mean-only fit, but the shuffle did better in both signal runs. Separately, expanding the squared residual around each origin mean shows that ordinary origin-balanced ridge equals the joint objective when the contrast weight is the fraction of origins with repeated cultures. The test suite verifies this equality. The candidate at other weights reweights the contrast component of the same loss. A novelty claim or a claim that real sister linkage drives gains requires an advantage over origin-balanced ridge and the shuffled-link control under a design that can distinguish these explanations. The single generated fixture family uses constructed outcomes and cannot estimate a biological effect or power.

## Metadata gate and requests

[Cohort preflight](../src/trackA_cohort_preflight.py) reads a coded CSV manifest containing only provenance. Its common fields cover study, accession, canonical origin and item IDs, split role, assay condition, source evidence, reuse terms, origin lineage and prior exposure. Regional manifests additionally require the modality and independent pathology region. Spatial manifests include the section and registration, spatial-input and independent-audit status. The TMZ manifest records availability of untreated RNA, matched measured viability, assay definition, independent culture derivation and MGMT. Each verification status is `verified`, `missing` or `unknown`. See the module's column constants for the exact header. The gate reports missing fields, independent origin counts, support for calibration or sister contrasts, and whether a shift has enough origins under a distinct condition. A regional sign-test result here is an attainability calculation. Sample size and joint power still belong in the experiment protocol.

The [GBM-Space section inventory](trackA_gbmspace_visium_sections_2026-11-23.csv) was run through the spot-level preflight using explicit unknowns for missing annotation and registration. All 97 sections from 12 identified donors lacked at least one required field, leaving zero verified complete pairs. An additional TMZ screen found a public baseline expression accession, [GSE214966](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214966), linked to a pharmacoscopy study. The paper reports TMZ measurements, but describes baseline RNA for only four patients and the assay measures a separate image-based endpoint. Hold it as a limited development lead under its own endpoint and linkage review. [Study](https://doi.org/10.1038/s41591-024-03224-y)

The metadata-only requests went to the GBM-Space, Petti, Park and Moniz Garcia study teams on November 23. The first Petti address bounced and the request was resent to her current Massachusetts General Hospital address listed by Harvard. The most useful inputs are an independent pathology-to-section/spot/cell crosswalk and the count of origins with at least two independent cultures with matched baseline RNA and TMZ assay. See the [execution ledger](trackA_execution_ledger_2026-11-23.md).

## Reproduction and disposition

Use a scratch copy as required by the repository. With Python 3.14.4, NumPy 2.5.1, SciPy 1.18.1 and scikit-learn 1.9.1, the 13 regression tests passed. The generated benchmark passed its seven implementation checks and reproduced the saved output byte for byte. The metrics above are from the [development output](../output/trackA_development_benchmark/result_with_origin_balanced_baseline.json). Its SHA-256 is `DDE309905E86CCC7F9EBC292F90531030DBE7E7382091612D0BD8BE247A442C1`. The [earlier output](../output/trackA_development_benchmark/result.json) preserves the run before the origin-balanced baseline was added.

```text
python -m unittest discover -s tests -p test_trackA_development.py
python -m src.trackA_development_benchmark
python -m src.trackA_cohort_preflight --help
```

Both claims remain blocked for external confirmation. The next scientific milestone is a qualified development source for independent spatial calibration or sister-culture learning, followed by a separate, sufficiently large untouched evaluation. Source-specific feature extraction, comparison adapters and the signed-off executable external protocols depend on those data. Earlier spent cohorts and negative results remain in the record.
