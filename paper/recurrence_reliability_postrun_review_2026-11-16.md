# Recurrence reliability: saved-result verification and next-study requirements

Date: 2026-11-16. Read-only verification of the completed locked run, followed by planning. No biological analysis was rerun; no frozen setting was changed.

## Verified on disk

Baseline checkout at review: `838d8de`. The reported result is `output/recurrence_reliability/run_20261116_locked_a1/`.

The saved patient table contains 21 unique patient rows and three opposite determinate directions. Recomputing the binomial interval from the saved event count gives 0.0304889676 to 0.3634239870, matching `result.json`. The saved config SHA-256 matches the current locked config, and all six recorded source-file hashes match the current implementation. This verifies artifact consistency and count arithmetic, not a complete implementation or scientific-validity audit. The fixed result remains inconclusive against the 5% tolerance.

The run is spent. Do not rerun it as fresh confirmation, add patients retrospectively, replace excluded libraries, or select additional scorers because of this outcome.

## Interpretation should remain narrower than “error bars work”

The observed result supports a within-cohort statement: three patients received opposite determinate directions under two fixed operational scoring procedures. Their within-procedure nucleus-bootstrap uncertainty did not describe the difference between procedures. Neither procedure establishes the true biological direction, so this is disagreement rather than a measured misclassification/error rate against ground truth.

The 5.5%/6.3% split-half rates are a diagnostic consistent with the specified cell-sampling model in these specimens. They are not proof of universal calibration. Splits reuse specimens and cells, and the bootstrap conditions on already-computed per-cell scores, reference pooling, control genes, and annotations. It does not integrate every source of biological, technical, or method-selection uncertainty.

The implementation estimates reference bins/control genes/centering from the shared final-cohort nucleus pool. Consequently, patient-level discordance indicators need not be unconditionally independent Bernoulli draws. The Clopper-Pearson arithmetic is reproduced exactly, but the advertised exact population coverage relies on a binomial model that should be assessed and disclosed. This observation does not authorize changing the frozen interval or rerunning the analysis. Any follow-up inference or sensitivity study needs a dated, prespecified analysis with the original result retained.

## Sample-size statements need the zero-event condition

The protocol uses a two-sided 95% interval. Its approximately 72-patient minimum for confirming discordance below 5% assumes **zero observed discordances**, not the current 14.3% point estimate. The one-sided 59-patient calculation is a different interval convention and cannot replace the locked rule.

- With zero events in 59 independent patients, the two-sided exact upper limit is `1 - 0.025**(1/59) = 0.0606089`, above 5%.
- Even if 59 new patients had zero events, simply adding them to the observed 3/21 would yield 3/80 and a two-sided upper limit of about 10.57% under a common binomial model.
- If exactly three events were retained and all additional observations were zero, at least 173 total observations would be needed for that upper limit to fall below 5%. This is arithmetic illustration, not a recommendation to keep recruiting until significance or to pool heterogeneous cohorts.

More patients may establish a discordance rate above tolerance, below tolerance, or remain inconclusive. Do not plan as if a positive reliability claim is inevitable. Wider scorer coverage and an independent cohort answer different questions; neither is required to exhaust every published scorer before a focused methods contribution can be assessed.

## CARE is the preferred next eligibility audit

An independent cohort is valuable for checking transportability of the observed disagreement, regardless of whether it can certify reliability below 5%. CARE's published 59-patient count is not yet an eligible-pair count. Verify patient pairing, IDH status, primary/first-recurrence identity, human clinical origin, overlap with Wang, author malignant-cell labels, library/barcode joins, processed count availability, source terms, and prior project exposure. Do not infer independence solely from accession or institution names.

Keep this step metadata-only. No CARE recurrence scores, state-abundance results, or survival outcomes should be opened to decide methods. If annotations encode published cell-state outcomes, extract only the fields needed for specimen identity and operational malignant/nonmalignant eligibility, and record any incidental exposure.

## Before an external scoring run

Write a separate external-validation protocol and config. Do not edit the completed Wang protocol/config. Resolve:

1. **Estimand:** external same-workflow reproducibility versus deployment of the specific fitted Wang scoring reference. These are different tests.
2. **Reference pool:** freeze whether to transport Wang gene centers/control genes/bins or rebuild the reference according to the same algorithm on CARE. Rebuilding is an algorithm replication, not a frozen fitted-reference validation. Specify feature alignment and missing-gene rules before scores. The original result stores control-gene counts rather than a complete deployable reference artifact, so do not assume that artifact already exists.
3. **Platform/QC/annotation mapping:** define permitted differences prospectively; source-specific author filters and annotation procedures are not automatically interchangeable.
4. **Inference:** assess the shared-reference dependence issue and the independent-unit assumptions. Prespecify intervals, multiplicity for repeat attempts/claims, and external acceptance criteria. Avoid outcome-dependent pooling.
5. **Scope:** same two scorers, program, and primary endpoint unless a distinct exploratory extension is explicitly declared. Do not select methods or thresholds using CARE outcomes.
6. **Run protection:** reserve the eligible cohort, freeze hashes/settings, validate adapters on synthetic inputs in the required scratch copy, then run the authorized external evaluation once.

No independent dataset, including CARE, is guaranteed to improve precision in the desired direction. If CARE fails access, pairing, annotation, or overlap gates, document the failure and keep the original result as a limited within-cohort methods study.
