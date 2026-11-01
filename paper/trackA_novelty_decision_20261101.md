# Track A novelty decision — 2026-11-01

## Recommendation

Report external confirmation as unavailable. Preserve the completed A5 analysis and all frozen results. Do not re-tune it or reframe a result from exposed data as a new biological discovery. No current option raises Track A's scientific novelty while preserving the six-donor and independence gates.

This decision does not change the current figure. It does not authorize data download, expression review, scoring, or contact.

## Research questions

**Question answered by the completed A5 analysis:** Among six eligible Yu 2020 donors, does the frozen candidate improve donor-held-out Core-versus-Periphery cell classification over ordinary L2 logistic regression under class-balanced AURC?

**Observed answer:** The candidate's mean AURC was lower (0.1479 vs. 0.1764), but it was lower for only four of six donors. The prespecified rule requires at least five. The candidate therefore failed its descriptive success rule. Cell-type-only balanced accuracy was 0.6546, within 0.05 of the best mapped-cell expression model at 0.6789; the frozen reading is comparable prediction from measured cell-type composition. The Brier score was worse for the candidate than logistic regression (0.2335 vs. 0.1944). The eight-donor sensitivity cohort was not run. See `output/trackA_yu2020/result.json` and `paper/trackA_execution_audit_20261031.md`.

A5 includes all eligible tumor/peritumoral point cells; it is not a malignant-cell-only analysis. Its six held-out donors are folds from one previously used cohort. They do not provide an independent external replication.

**Unanswered external question:** Does a frozen, donor-stable Core/Periphery signature in defensibly labeled malignant cells transfer to at least six new, independent IDH-wildtype patients after accounting for normal-brain content, cell-type composition, and established tumor-cell states? No currently screened cohort passes the source, sample-map, access, reuse, and independence gates needed to answer it.

## Comparison of paths

| Path | What it can add | Decision |
|---|---|---|
| Strengthen the approved A5 analysis | Clear donor-level reporting and presentation of its already-frozen logistic, alternative-model, ablation, cell-type, point-size, and within-cell-type controls. The primary result remains negative under its frozen rule. | Use for reporting only. No post-score analysis can turn this into independent confirmation or a new biological claim. |
| Reframe results around data already used | A descriptive observation from an exposed cohort. | Reject as a novelty strategy. The Yu cohort has prior project use; SNUH donors appear in prior project audits; Ivy GAP is already used. Reframing after outcomes were seen would be exploratory and would not satisfy independence. |
| Report external confirmation unavailable | A transparent account of the failed external-data gates and the limits of current evidence. | Recommend. Keep the six-donor and independence requirements intact; reopen only for a specific new accession that plausibly clears every gate. |

## What is and is not new

Core/periphery tumor-cell differences are established. Darmanis et al. reported infiltrating neoplastic cells in peripheral tissue and a shared infiltration signature across four GBM patients ([2017 primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC5810554/)). Wang et al. reported regional malignant-cell abundance and expression differences between MRI-defined HBI and LBI peritumoral regions ([2024 primary paper](https://doi.org/10.1038/s41467-024-55243-5)). Pai et al. compared tumor-Core and tumor-Margin cell populations, derived a shared margin signature, and compared regional tumor-cell states in four grade-4 HGG; the cohort includes mixed molecular drivers ([2025 primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12603948/)). Sonpatki et al. reported six matched Core/Periphery single-cell samples, with major cell populations including neoplastic cells, and spatially mapped GBM niches ([2026 primary paper](https://doi.org/10.1038/s41467-026-69716-2)).

The 2026 atlas is especially close to the regional single-cell question and its six paired donors overlap prior project use recorded for SNUH. It cannot serve as independent confirmation. The narrow A5 contribution is a frozen, donor-held-out benchmark of one candidate classifier against conventional baselines in one existing cohort. Its candidate did not meet the success rule. This is not evidence of a new malignant-cell program, causal mechanism, spatial transfer, or generalizable superiority. Do not claim “first.” Do not raise the internal novelty estimate on current evidence.

## Controls and donor-level validation

**A5 controls already completed:** six-fold leave-one-donor-out evaluation; training-only preprocessing and nested tuning; L2 logistic regression as the primary comparator; three additional model baselines and three candidate ablations; cell-type-only and point-size diagnostics; and within-cell-type checks where donor and class counts met the frozen minimum. Preserve the reported limitations: six donors support descriptive estimates only, point is nested within donor, labels are author-supplied MRI-marked point labels, and no result establishes a cell-intrinsic effect. Do not add post hoc tests.

**Controls required before any future external test:** use the same prespecified patient-level endpoint for every model; compare with ordinary regularized pseudobulk or equivalent baseline; control for normal-brain content and cell-type composition; compare against established cell-state programs and prior regional signatures; include size-matched gene-set controls where applicable; and keep feature selection, preprocessing, and tuning within training donors. Do not treat cells or spatial spots as independent donors. Spatial spots can support a separate spatial description but cannot supply malignant-cell labels.

**Donor validation gate:** require at least six IDH-wildtype patients with both regions mapped to the same donor in public metadata. Keep all samples and cells from each donor in one fold. Aggregate or compare at donor level and give each donor equal weight. Report paired donor effects and uncertainty; with six donors, call results descriptive unless the frozen precision and success criteria support a stronger claim. Require defensible cell-level malignant labels with documented validation, exact region/sample joins, clear no-request access and reuse terms, exact sample-ID cross-checks against local training and prior-use records, and verified cross-study independence.

## Open gates and strongest near-miss

- **No eligible independent cohort passes.** The latest metadata-only screen and exact-ID audit are recorded in `paper/trackA_alternate_dataset_plan_20261101.md`.
- **Strongest label-ready near-miss: GSE286419.** Its paper retains four Core/Margin pairs after pathology/spatial checks, below six. Sample-level IDH-wildtype eligibility, record-specific reuse terms, and biological independence remain unresolved. It fails even before those open gates.
- **GSE318564:** three explicit IDH-wildtype Core/Edge pairs, below six; the join from paired matrices to ready-to-use malignant-cell annotations remains unverified.
- **HRA009021:** the paper reports ten donors, but the public GSA sample-to-region crosswalk is absent and data reuse terms are unresolved. Do not infer regions from suffixes or expression. No contact has been made.
- **SNUH-GliomaAtlas:** six paired donors, but all six appear in the project donor manifest and prior exploratory audit. Fail independence.

Any missing donor map, ambiguous region, unvalidated malignant-cell label, unclear access/reuse term, donor count below six, or unverified overlap is a stop. A new lead must be named before any further search. Broad scouting stays paused.

## Stop rules and effort

Stop current A5 analysis at its frozen result. Do not run the eight-donor sensitivity cohort, add a new gene set, retune a model, search for a favorable subgroup, or change the figure. Do not use spatial spot annotations as cell-level malignant labels. Stop a candidate before download if any required gate above fails.

**Effort:** about 4–8 focused hours (up to one workday) to finalize this evidence-bounded report and keep the roadmap current; no compute or new data needed. If a specific candidate later appears, allow 1–3 days for metadata, terms, donor mapping, and overlap gates. Only after review and approval could a ready-to-use cohort support a frozen analysis, likely 1–3 weeks depending on annotation readiness and compute. HRA009021's roughly 765 GB raw FASTQ route is not a realistic fallback while its map and terms remain unresolved.
