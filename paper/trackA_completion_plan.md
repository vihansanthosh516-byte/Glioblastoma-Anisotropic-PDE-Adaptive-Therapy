# Track A completion plan

Status: research roadmap. Each confirmatory run needs a separate signed-off executable protocol. No test data are opened under this roadmap.

## 1. Definition of done

| Box | Required evidence | Fails the box |
|---|---|---|
| Positive | Locked candidate achieves a meaningful patient-level gain over the strongest applicable baseline on untouched patients or origins. Prespecified uncertainty and multiplicity checks pass. | Simulation alone, a selected subgroup, a weak comparator, or performance estimated after tuning on the test set. |
| Unique | A claim chart shows exactly what differs from the closest published methods and experimental designs. A same-data ablation demonstrates that the claimed new component contributes. | Renaming known components, unsupported first-ever claims, or gains attributable to a different input set. |
| Substantial | The gain survives a prespecified shift relevant to the claim. Report effect size, uncertainty, coverage, and what inference the result supports. | An artifact shortcut, low coverage, or a benefit that disappears under the intended measurement conditions. |

Track A retains two distinct endpoints: independently labeled core versus infiltrating tissue, and measured temozolomide response. Mark the three boxes for each claim separately. A regional success does not validate drug response. A patient-level result that meets all three boxes supports the bounded claim tested. Clinical benefit needs separate evidence.

## 2. Critical path and parallel lanes

| Order | Workstream | Deliverable required to advance |
|---|---|---|
| A | Regional data | Metadata ledger with accession, license, patient and section identifiers, explicit independent pathology labels, transforms, modality, provenance and overlap flags. Inventory missing cells and complete pairs by cohort. |
| B | Drug response data | Verified untreated expression, MGMT where measured, measured TMZ response, dose and assay details, and culture-to-origin crosswalk. Identify source signature overlap and flag related cultures. |
| C | Method and novelty | Register one primary candidate per endpoint, its failure hypothesis, closest precedents, input requirements, strongest baselines, and the component removed in each ablation. |
| D | Technical development | Demonstrate that candidate checks respond to the artifacts they claim to detect. Use only simulation or declared development data. Compare fixed alternatives on identical splits. |
| E | Lock and feasibility review | Freeze protocol, code, environment, eligible cohorts, model and feature rules, margins, power, multiplicity, missing data, fallback and reporting rules. Commit and push before opening test expression or outcomes. |
| F | Blinded external evaluation | Run the locked pipeline once. Save all patient-level paired scores, changes from baseline, uncertainty, failures, excluded records and reasons. Record deviations and spend the test cohort. |
| G | Paper and reproducibility | Produce evidence table, complete final benchmark, limitations, linked development record, source and output hashes, commands and a reproducibility audit. Mark each box passed, failed or blocked. |

Run A, B and C concurrently. Use technical development to answer specific design questions while waiting for data. Do not treat additional synthetic wins as substitutes for independent validation. The [external design](trackA_external_design.md) explains why the minimum attainable p-value and usable sample size are different questions.

## 3. Regional lane

- **First candidate:** compare paired predictions estimated with nominal segmentation and with a frozen alternative on the same section. Use a validated independent contamination check or an error bound for common errors. Omit spatial input when the required evidence is missing. Agreement between methods alone is insufficient. [Shared-error falsification](trackA_segmentation_result.md)
- **Same-input comparison:** paired expression and composition regression, density and depth controls, the same nonspatial model plus nominal spatial features, alternative segmentation alone, segmentation-free spatial features, and the full candidate. Test shared-error controls and remove each proposed safeguard in turn. Retain one patient paired score across arms.
- **Claim boundary:** require spatial gain beyond composition, density and malignant expression. Require robustness and coverage under an unseen segmentation condition. Claim assay transport only after a genuinely held-out acquisition test. Whole-tissue expression is a valid separately scoped fallback for the same region question, with a separate protocol and source-appropriate baselines.
- **Data order:** resolve the [Petti pathology join](trackA_petti_acquisition_audit.md) for technical development. Audit [GBM-Space](trackA_gbmspace_screen.md), additional CosMx sources and potential overlap. Request [Larroquette metadata](trackA_parallel_eligibility.md) for an expression-only route. Do not substitute fields or replicate sections for patient count.

## 4. Measured TMZ lane

- **First candidate:** learn from within-origin response contrasts only if independent sister cultures exist, with a separate origin-mean term. Freeze the estimator and endpoint before evaluation. Tissue region, treatment receipt, and survival are not substitutions for measured viability.
- **Same-input comparison:** MGMT alone, regularized expression, fixed source signature, matched random signatures, and the same candidate without the within-origin term. Shuffle sister links during development as a falsification control. Keep related cultures together. Reserve independent origins and preferably a separate lab.
- **Data order:** prioritize the [organoid acquisition gate](trackA_petti_acquisition_audit.md). Request authorized sequencing access and coded baseline-to-assay joins. Document the published response-group exposure and keep scoring blinded if this cohort qualifies. Exclude GSE232173 from independent validation of the signature derived from that study. [TMZ screening record](trackA_parallel_eligibility.md)

## 5. Lock checklist

- Provenance manifest plus exact patient and specimen join. Identify which records and outcomes developers have seen.
- Sample size based on a minimum meaningful effect, plausible discordance, and the actual joint pass rule. Freeze paired estimand, confidence intervals, comparator eligibility, multiplicity and handling of ties.
- Predeclare the primary candidate, all eligible comparator adaptations and exclusions, feature data available at prediction, coverage floor, and what triggers the fallback.
- Preserve raw outcomes from the model development team until the code, plan and configuration are hashed. Use a blinded operator when practical.
- Predefine the analysis for adverse assay or cohort shifts, negative controls and ablations. Persist results for every arm and every predeclared check.

## 6. Decision rules after each gate

| Event | Action |
|---|---|
| Missing mapping or insufficient patients | Record a no-go, request the exact missing metadata, and continue parallel cohort work. No expressions or outcomes from the candidate test are opened for tuning. |
| A development test fails | Record whether the failure is implementation, measurement, modeling or identifiability. Fix implementation defects transparently. A new scientific hypothesis gets its own protocol and fresh draws. |
| External test fails | Publish the negative result with limitations. Freeze a distinct candidate or test a new cohort. Do not salvage the verdict by changing the rule or removing difficult patients. |
| Positive result but novelty unclear | Check close published methods, run the locked ablation if specified, and narrow the claim or advance a distinct design. |
| Positive and distinct but fragile | Hold the substantial box open. Use a separate independent shift test with its own lock. |
| All boxes supported | Prepare the final paper and reproducibility record. Lead with the final evidence. Disclose relevant development decisions and retain unsuccessful experiments in the supplement or repository. |

## 7. Reporting and ownership

Maintain a ledger by lane with status, evidence path, next action, dependency, and reason for any stop. Commit and push verified steps. Keep routine updates in the repository and reserve a completion report for supported claims. Flag a concrete external dependency when it prevents meaningful progress. Artzi, Yu, and HGCC remain spent; all prior failures remain recorded.
