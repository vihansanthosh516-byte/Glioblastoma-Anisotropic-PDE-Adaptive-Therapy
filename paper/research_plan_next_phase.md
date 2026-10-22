# Research plan: next phase

Updated 2026-10-21. This is the current work plan. `START_HERE.md` remains the quick resume guide.

## What this research asks

The current project asks how disagreement between tumor-segmentation tools and assumed measurement noise affect volumetric progression categories and simulated progression-free survival under simplified RANO-inspired rules. It does not measure scan-rescan noise with its own pipeline. Its current evidence is an exploratory real-data tool comparison plus sensitivity simulations anchored to one published repeatability study through an assumed conversion. Novelty and clinical validity are unestablished.

## Where the work stands

The A15-A21 analyses, calculator, and viewer already exist. The technical note is an audited draft, not a finished paper or clinical result. Recent source checks corrected its ranges, citations, reference-scan wording, and limits on execution-history claims. Some primary tables and supplements remain unread. No new analysis is currently authorized by the existing plan.

The user reports that their mentor thinks the current work is good. The mentor's detailed comments and clinical or imaging qualifications are not recorded here.

### Dataset feasibility screen, 2026-10-21

The closest scientific fit found so far is QIN-GBM Treatment Response: its TCIA page describes double-baseline scans typically 2–5 days apart, 54 patients, multiple structural and advanced MRI sequences, and 33.5 GB of DICOM images. A published test-retest analysis used 11 newly diagnosed GBM patients with scans 2–6 days apart. The 11-patient analysis does not establish that all 54 patients have usable scan pairs for this project's segmentation question.

QIN is not cleared for use. Its collection page says users must sign a TCIA Restricted License Agreement. TCIA's newer policy says it stopped hosting controlled-access datasets on July 8, 2025 and moved access and hosting to dbGaP and the NCI Cancer Research Data Commons. The QIN page still shows the older TCIA route, so current request eligibility and hosting route are unresolved. The restricted-license language also limits disclosure to third parties unless permitted; Kaggle use is not confirmed.

The desk-only screen found open longitudinal datasets such as LUMIERE and UCSD-PTGBM, but their follow-up scans do not provide a clean short-interval scan-rescan estimate. No fully open tumor scan-rescan dataset that meets the current measurement objective has been verified in this screen. This is a preliminary search, not an exhaustive review.

Sources checked: [QIN collection page](https://wiki.cancerimagingarchive.net/display/Public/QIN%2BGBM%2BTreatment%2BResponse); [QIN test-retest study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9207284/); [TCIA controlled-data policy](https://www.cancerimagingarchive.net/nih-controlled-data-access-policy/); [TCIA restricted-license text](https://www.cancerimagingarchive.net/wp-content/uploads/TCIA-Restricted-License-20220519.pdf); [LUMIERE dataset paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9755255/); [UCSD-PTGBM dataset paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12886825/).

## What happens next

### Phase 1 — Close the current note review

1. Keep the note's distinction between observed tool disagreement, assumed-noise simulations, and published repeatability converted by assumption.
2. Ask a clinician or imaging physicist to review the RANO framing, 500 mm³ floor, post-treatment scan context, and interpretation of the repeatability study. Record comments as review feedback; do not turn them into new results.
3. Resolve only source gaps that can be checked lawfully and directly. If a supplement or table remains inaccessible, keep the claim marked unresolved instead of delaying the project indefinitely.

**Gate:** The note can move to a final internal technical report when every material claim is traceable or explicitly labeled unresolved, and the clinical framing has received expert review or is clearly marked as not reviewed. It must not claim measured noise, clinical validity, novelty, or a theorem.

### Phase 2 — Test whether a stronger data contribution is feasible

Run a desk-only search for an alternative to QIN and RIDER. Do not download images or run segmentation during this phase. A candidate must have documented lawful access for this project, repeat tumor MRI suitable for estimating measurement variation, enough sequence information for a defensible segmentation method, and a usable license or written access route. Check data counts and series on primary collection pages; do not rely on search snippets alone.

**Gate:** Proceed only if one dataset meets all access, repeat-scan, and measurement-method requirements. RIDER and QIN currently require restricted access. QIN remains the strongest technical fit, but its current access route and permitted compute environment must be confirmed in writing. No qualifying open alternative has been verified yet.

If no dataset passes, stop the new-code path. Keep the existing work as a limited exploratory and simulation-based technical note, or choose a separately justified research question with its own evidence and plan. Do not invent a dataset or force a new analysis just to keep coding.

### Phase 3 — If a dataset passes the gate, preregister before analysis

Write a new amendment before code runs. Define one primary estimand, inclusion and exclusion rules, segmentation inputs, noise model, uncertainty method, and a failure criterion. Make an explicit resource and license plan for desktop or Kaggle use. Keep access-controlled images off any service unless its use is allowed by the dataset terms. Have a mentor review the amendment, then commit the amendment and code before the planned run, following `START_HERE.md` rules.

### Phase 4 — Analyze, test failure modes, and write

Run the declared analysis on the eligible data. Report uncertainty and all planned failures. Compare the result with the existing simulation only where the estimands align. Do not describe tool disagreement as pure noise. Do not describe a sensitivity result as a measured effect. Write a methods-and-results report, then get an independent technical review before public release or submission anywhere.

## Working timeline

Use gates, not a fixed ten-month prerequisite block. Phase 1 can take one to three weeks, depending on reviewer availability. Phase 2 can take two to four weeks. If a dataset passes, allow several months for amendment, segmentation validation, analysis, and independent review. The full path may take six to nine months at roughly 20 hours per week, but it may end earlier if no suitable dataset exists. The current note can still be completed without opening a new code project.

## Codex and Claude handoff

Use Claude for one bounded, read-only audit or source-check task at a time. Codex checks its claims against the repository, primary sources, and project rules, then makes any supported document corrections. Keep the human researcher responsible for understanding every claim, reviewing the analysis, and discussing scientific decisions with a mentor. No agent may claim work, source access, or analysis that it did not perform.

## Current next action

Complete Phase 1 review and begin Phase 2 as a desk-only dataset feasibility search. Do not start new code until a candidate passes the gate and a new amendment is committed.
