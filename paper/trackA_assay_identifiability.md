# Track A assay gate novelty and identification audit

The [initial simulation](trackA_assay_gate_result.md) showed conditional protection against a reversing artifact. The [transport stress test](trackA_assay_transport_result.md) showed that the gate fails when the controls miss the artifact. This audit explains why that limit cannot be repaired merely by tuning the threshold.

## Two indistinguishable training mechanisms

Use the scalar model in the [frozen stress plan](trackA_assay_transport_plan.md). Write the spatial contrast as S = B + A + E, with regional biology B, assay contrast A, and measurement noise E. Consider a biological world (B=0.80, A=0 in development and test) and an artifact world (B=0, A=+0.80 in development and A=-0.80 in test). In both worlds the development spatial contrast has the same distribution. Set the other development features equal across worlds and assume the control assay is blind, so its observations also have identical distributions. All constants are copied or derived from the frozen stress design.

Any gate using only that development information has the same selection probability in both worlds. It cannot reliably keep spatial evidence in the first world and reject it in the second. More replicates of a blind control leave this equality intact. This is a mathematical counterexample within the stated model, not a claim that these two mechanisms explain the patient data.

## What a distinctive experiment would add

Cross the same independently pathology-labeled tissue blocks with randomized acquisition conditions. Register adjacent sections when the same material cannot be measured twice. Include same-condition replicates to quantify section-to-section variation. Deliberate positive controls must verify that the chosen contrasts respond to each targeted technical perturbation. Reserve patients and an acquisition condition before model development. A separate control panel for that condition can establish coverage before patient predictions are scored, with no access to held-out region outcomes.

Compare the control-gated spatial predictor with identical predictors with spatial always included or always omitted, plus replicate-based nuisance correction and a robust predictor using the same acquisition labels. Evaluate patient paired concordance and report coverage and failure by acquisition condition. A claim of benefit requires improvement on an independently held-out condition where the spatial feature is informative, and preservation of the expression baseline in a condition where it is misleading. Set the margins and sample size from development evidence before the test. No patient data were used in this audit.

## Prior art and claim boundary

[Risso et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC4404308/) use control genes and replicates to estimate unwanted variation. [Anchor regression](https://arxiv.org/abs/1801.06229) supplies robustness against a specified class of shifts. [Spatial platform benchmarking](https://doi.org/10.1038/s41467-025-64990-y) already uses serial sections across assays. [Xenium contamination and specificity work](https://www.nature.com/articles/s41592-026-03089-8) already evaluates technical reproducibility. These precedents rule out broad novelty claims about controls, robust prediction, or serial-section assays. A narrow candidate contribution is a coverage-qualified decision about whether spatial organization adds predictive information for independent paired GBM regions across acquisition conditions. No exact-match claim or positive patient result is established.
