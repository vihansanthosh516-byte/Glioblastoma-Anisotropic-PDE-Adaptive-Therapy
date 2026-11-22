# Track A research recovery plan

Keep the question fixed: does tumor biology distinguish core from infiltrating tissue, and can a fixed biological score predict measured TMZ response? Keep the regional and drug response experiments separate. Use the patient or verified biological origin as the independent unit.

## Regional experiment

Develop a block pooling score that can retain localized signal. Fit paired logistic regression on training patients, score individual tissue blocks, then pool scores with a bounded smooth maximum. Freeze the clipping limits and pooling temperature using development data only. Compare against pooled paired logistic regression, mean and median block pooling, unbounded smooth maximum, and a composition-only control. Use the same patient splits and paired concordance. Bounded pooling is a candidate, with novelty unresolved. Smooth max and learned pooling have direct precedents in [Pinheiro and Collobert](https://openaccess.thecvf.com/content_cvpr_2015/html/Pinheiro_From_Image-Level_to_2015_CVPR_paper.html) and [Ilse et al.](https://proceedings.mlr.press/v80/ilse18a.html).

Freeze a fresh simulation with diffuse and localized signal, isolated technical artifacts, spatially correlated artifacts, and corruption present only at test time. Include a null with depth and batch imbalance. Compare median pooling directly to determine whether any gain comes from retaining a small patch. The prior median result came from [diffuse signal with matched corruption](trackA_outlier_result.md).

Before opening a new cohort, verify independent block IDs, pathology and patient linkage, assay compatibility, and overlap with development cohorts. Set a minimum meaningful improvement over the strongest locked comparator and size the independent cohort for that paired difference. Freeze the precision criterion, multiplicity handling, eligibility, and missing-data rules. Include independent composition estimates where available. If a cohort only supports whole tissue claims, report that scope explicitly. Artzi and Yu remain spent. [Artzi result](trackA_artzi_result_2026-11-17.md) | [Cohort leads](trackA_fresh_data_followup_2026-11-13.md)

## Measured drug response

Keep TMZ viability or another prespecified measured TMZ endpoint. Require an untreated baseline expression profile, a matched assay outcome, verified biological origin and an untouched cohort. Compare a locked model against MGMT where available, regularized expression, and matched random signatures. Group related cultures by origin and keep an entire assay cohort outside development when feasible. Record assay and dose differences before scoring. The previous proxy and HGCC signature results stay negative. [Failure audit](trackA_failure_audit_2026-11-13.md)

This is a research roadmap. Each experiment needs its own executable plan and pass rule before its test data are accessed.
