# Claude handoff: strengthen the GBM project for a research paper

Workspace: `C:\Users\vihan\20206 science fair`
Prepared: November 16, 2026
Baseline commit: `dbb53b4`; verify the current checkout because another assistant may have made newer commits.

## 1. Your assignment

Continue the existing science-fair research and develop a credible research-paper extension. The selected next direction is a **patient-paired recurrence measurement-reliability study**, conditional on its data and novelty gates. Complete the necessary research, source verification, reproducible implementation, protocol, and permitted analyses when the user asks you to execute this handoff. Do not merely repeat the candidate list or ask the user to choose routine implementation details.

The immediate objective is a decisive answer: can this extension support a distinct, adequately defined study with usable data? If yes, produce and follow a frozen protocol. If no, document why and recommend the strongest remaining option. A transparent no-go is a valid outcome. A positive biological result, novel method, or publication is not guaranteed.

This document supplies the entire workflow, not an already-frozen statistical protocol. Some scientific choices cannot honestly be finalized until the source, annotation, batch, and closest-prior-art checks finish. Resolve those choices using the decision rules below before computing recurrence program outcomes. The user's current request to write this handoff is not itself a request to start experiments in this chat.

## 2. Start here and preserve the workspace

Run `git log -5 --oneline` and `git status --short` first. Read `START_HERE.md`, then applicable `AGENTS.md`, and `CLAUDE.md` for environment/build information. Newer dated records take precedence over historical summaries.

Follow project workflow: inspect status before staging; stage exact relevant paths; commit after verified changes; push immediately after each commit. Never use blind `git add -A`, overwrite unrelated work, reset another assistant's changes, or restore old generated outputs. Several `output/patient_split/` files and `paper/trackc_prior_art_full_2026-11-04.md` were already modified before this handoff. They are not yours to discard or include automatically.

Use `graphify query "<focused question>"` before grepping project content. The graph exists but may lag newer documents; inspect the cited current files afterward. Apply installed skills where relevant. Use normal, precise prose in saved documents; keep chat updates concise. No email, author contact, access-request submission, or upload to a third-party service without the applicable user authorization.

Keep raw data and workbooks local, outside commits. Some report/notes directories are gitignored; force-add only exact intended text artifacts when necessary. Do not change ignore rules globally to include all generated material. Keep a dated handoff current so a different assistant can resume.

## 3. Current scientific state

- **Track A:** molecular/spatial research. Frozen A5 failed its donor-level success rule. The frozen expression-to-TMZ response bridge did not show an association. New recurrence work does not erase those failures or validate the existing regional classifier.
- **Track B:** longitudinal MRI forecasting. MU, UCSD, UCSF, and LUMIERE are spent for the existing forecasting question. CFB/WAW are reserved under the locked test chain. Do not open their protected labels/outcomes, tune on them, or reuse them for this extension.
- **Track C:** treatment-policy simulations. Drug-effect assumptions and simulator validity remain limitations. Current run, sampling, inference, and 14-day precision gates remain unchanged. Do not launch policy evaluation through this handoff.
- **Integration/EARF:** implemented and locked but not evaluated in the current record. Real molecular/MRI bridge tests failed. EARF partly shares reliability anchors between its controller and generating worlds. No validated patient-level molecular–MRI–treatment join exists.

This extension adds longitudinal tissue evidence and a reliability question. It does not create a validated digital twin, estimate patient-specific drug benefit, infer treatment causality, or identify state-transition rates from destructive tissue snapshots.

## 4. Read these files; do not redo completed work

1. `reports/GBM paper strengthening options.md` — broad research and later corrections.
2. `paper/recurrence_reliability_plan_2026-11-16.md` — selected direction, prior art, precision limits.
3. `paper/recurrence_metadata_gate_2026-11-16.md` — official GEO metadata audit.
4. `paper/recurrence_specimen_crosswalk_2026-11-16.md` — supplementary Table 1 crosswalk.
5. `paper/paper_strengthening_handoff_2026-11-16.md` — chronological handoff.
6. `paper/spatial_protein_metadata_gate_2026-11-16.md` — alternative spatial direction and its limitations.

Machine-readable artifacts are in `output/paper_strengthening_metadata/`:

- `recurrence_series_metadata.txt`, `recurrence_sample_manifest.json`, `recurrence_pair_manifest.json`.
- `recurrence_table1_crosswalk.json` — selected sample/assay fields only; excludes survival.
- `recurrence_file_listing.html`, `recurrence_filelist.txt`.
- `recurrence_processed_coverage.json` — file triples for the conservative candidate subset.
- `care_series_metadata.json` — external-cohort leads only.
- `spatial_manifest.json` — alternate spatial resource file listing.

The local source workbook is `43018_2022_475_MOESM2_ESM.xlsx`, SHA-256 `58e42239d4e105cfbbda6441bd640afbd64e0c55f87e07b1e0aaab0d1eff1c42`. It and `PMC9767870_supplementary_media.zip` remain uncommitted. Only selected identity, stage, pairing, diagnosis, IDH, and assay-presence columns from Table 1 were extracted. Do not inspect survival or unrelated gene-result sheets merely because the workbook is present.

## 5. Verified facts and what they do not establish

Discovery/reanalysis source: Wang et al., Nature Cancer 2022, DOI `10.1038/s43018-022-00475-x`, GEO `GSE174554`.

- Paper describes 86 specimens from 49 patients overall and a 36-patient paired subgroup.
- Official GEO series summary advertises 40 primary and 40 matched recurrent specimens. This statement really is in official metadata; an interim note attributing it only to secondary sources was superseded.
- Current GEO lists 113 GSM records: 111 human and two mouse. Only 81 titles explicitly contain `_snRNA`.
- Supplementary Table 1 contains 111 specimen rows; 86 are marked snRNA-seq. Of these, 83 are IDH-wildtype and three IDH-mutant.
- Table 1 has 33 simple groups with one primary and one recurrent specimen, both IDH-wildtype. Of those, **27 pairs / 54 specimens** have unique exact SF-ID matches to GEO with agreeing stage labels. This is the current conservative candidate set, not a final biological-QC cohort.
- Eight Table 1 RNA specimens lack an exact GEO title-ID match. SF9798/SF9494 have conflicting stage labels. A separate group has two recurrence specimens without an explicit first-recurrence field. Do not guess corrections or choose the more favorable recurrence.
- Pair numbers are largely renumbered between Table 1 and GEO. Join by specimen identity and compare group membership; a different numeric pair label is not itself biological disagreement.
- All 54 candidates have a complete listed matrix/features/barcodes triple. There are **61 libraries**, because seven specimens have multiple libraries: SF10099, SF12243, SF4449, SF4810, SF7307, SF10433, SF4209.
- Listed candidate files total **713,265,211 compressed bytes**. Enclosing `GSE174554_RAW.tar` is **1,655,029,760 bytes**. These are file-listing facts, not verified file integrity or a proven minimal download route.
- Raw sequencing reads are cited at EGA `EGAS00001004909`. GEO's `RAW.tar` name does not prove public raw-read access; its listing contains processed matrix files.

No candidate expression matrices, gene lists, barcode contents, recurrence scores, or survival values have been analyzed in this work. Published literature findings were read; that is distinct from opening patient-level study outcomes.

Official metadata endpoints can avoid the GEO browser CAPTCHA:

```text
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174554&targ=self&form=text&view=brief
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174554&targ=gsm&form=text&view=brief
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174554/suppl/filelist.txt
```

Earlier Firecrawl full-text retrieval stalled. If a retrieval route fails, use a bounded-time primary-source fallback and record its limit; do not spend repeated long waits on the same failed call.

## 6. Phase A: close the novelty gate first

Create `paper/recurrence_reliability_prior_art.md` with a compact closest-method comparison. Use primary papers, accessible methods/supplements, and author code. Record the exact question, cohorts, biological program, scoring, patient-level uncertainty, annotation checks, and claimed novelty of each close study. Distinguish full-text verification from abstract-only evidence.

At minimum examine:

1. Wang et al. 2022: https://doi.org/10.1038/s43018-022-00475-x
2. Spitzer et al. 2025: https://doi.org/10.1038/s41588-025-02168-4
3. Nomura et al. 2025: https://doi.org/10.1038/s41588-025-02167-5
4. Alexopoulos and Karampelas, online October 25, 2026: https://doi.org/10.1177/15578100261499487

The last paper already analyzes Wang and CARE patient-level recurrence trajectories, including uncertainty and malignant-annotation checks. Its abstract was accessible; full methods were restricted during the screen. Therefore patient pairing, a larger Wang subset, comparing recurrence programs, using CARE, or rediscovering mesenchymal change is not sufficient novelty.

Search specifically for technical reliability of patient-level recurrence conclusions under sampling, batch handling, scoring, malignant-cell annotation, and state definitions. Do not rely on a search that merely fails to find an identical title. State the contribution narrowly and identify its closest competitor.

**Decision:** proceed as a research-methods/measurement study only if a defensible unresolved question remains. If the same question and evaluation already exist, classify the work as replication/benchmarking or stop that novelty claim. Do not keep adding architectures or signatures until something appears unique.

## 7. Phase B: complete metadata, annotation, and batch eligibility

Create an explicit candidate-patient table derived reproducibly from the existing JSON artifacts, not a handwritten list. Preserve source specimen IDs, GSM, stage, diagnosis, IDH, pair membership, library IDs, and every exclusion reason.

Use the 27 clean cross-source pairs as the conservative starting subset. Resolving additional specimens is allowed using source metadata only. Any expansion must have documented provenance before outcome scoring. Do not force a target count of 36 or 40.

For every specimen, establish:

- Which processed matrix and cell/barcode metadata belong to it.
- Whether author-provided malignant/nonmalignant labels exist and can be joined unambiguously by library and barcode.
- Whether malignancy labels are supported by the authors' CNV/genomic procedure. Such labels are operational annotations, not perfect independent ground truth.
- Whether recurrence order and treatment context support the chosen primary-to-first-recurrence question. If only generic recurrence can be established, change the proposed endpoint transparently before scoring, or exclude the uncertain specimen.
- Whether multiple libraries are independent captures, resequencing, technical duplicates, or different material. Never concatenate, deduplicate, or sum blindly. Preserve library provenance in cell IDs.
- File-level terms, attribution, access route, sizes, checksums where supplied, and intended-use compatibility. Do not bypass controlled access or claim that public visibility establishes permission to redistribute.

Do not infer malignant labels using the same MES/state score that will be tested. If author annotations cannot be joined, evaluate an independent CNV-based operational approach and its limitations before proceeding; do not quietly switch the population to all cells.

Write `paper/recurrence_reliability_eligibility.md` and `output/recurrence_reliability/metadata_manifest.json`. Include gate states, source URLs, hashes, excluded candidates, batch interpretation, and remaining uncertainty. No gene-state outcome analysis yet.

## 8. Phase C: separate technical QC from biological outcome inspection

If expression data must be inspected to establish usable cell support, first freeze a **QC-only plan**. Specify matrix/count checks, barcode joins, gene-ID mapping, library relationships, fixed technical cell filters, minimum malignant-cell support, and exclusions. Derive thresholds from methods/technical justification or synthetic assessment; do not choose them to obtain a desired donor count.

When the user has instructed execution and public-data download is within scope, obtain only the smallest authorized useful files after access/terms checks. Avoid a full raw-read download. Keep large data local and hashed. If another approval rule genuinely applies, complete the concrete download manifest first and explain that exact rule rather than inventing a new approval step.

During QC, do not compute, plot, rank, or test primary-versus-recurrence program differences. Do not inspect survival, response labels, candidate gene associations, or external-cohort outcomes. Freeze the final eligible cohort after QC and record any incidental outcome exposure.

Required QC checks include matrix dimensions, count validity, unique cell identifiers, annotation join rates, gene-ID collisions, specimen/library completeness, duplicate captures, and malignant-cell support. Different source files may use different feature identifiers or contain processed rather than raw counts; verify rather than assume.

## 9. Phase D: write and freeze the actual study protocol

Create `paper/recurrence_reliability_protocol.md` and a machine-readable config. Resolve every analysis choice before recurrence program scoring; no unresolved placeholders in the locked protocol.

### Scientific question

How often does a patient's inferred primary-to-recurrence program change depend on prespecified technical processing choices, when the intended biological construct is held fixed?

This is a reliability study, not a recurrence classifier or a treatment-causality model.

### Biological construct and scoring

Choose one published malignant-cell program through the source/novelty review. Freeze the exact genes, gene-ID mapping, module handling, normalization, control-gene procedure, missing-gene rule, and software versions. Do not select the program by examining recurrence associations.

Prespecify a simple reference score and one justified alternative procedure for the same program. Explain whether they estimate the same mathematical quantity or different operational summaries; do not compare their raw numerical magnitudes as if scales were identical. Retrieve the original method definitions before implementation. Freeze random seeds for any control-gene or cell-sampling steps.

Use different published biological state definitions only in a clearly separate secondary construct-sensitivity analysis. Different constructs can legitimately disagree; that is not automatically technical error.

### Patient-level endpoint

For patient `i` and procedure `k`, define `delta[i,k] = score(recurrence[i], k) - score(primary[i], k)`.

Choose exactly one primary reliability endpoint before scoring. A candidate is the fraction of eligible patients whose directional conclusion differs between the two fixed procedures. Explicitly define increases, decreases, indeterminate results, denominator, missingness, and ties. Calibrate any indeterminate region using independent technical evidence or a declared model; never tune it on recurrence results.

Report continuous patient-level changes alongside categories. Report indeterminate coverage so a procedure cannot appear reliable merely by withholding every conclusion. Distinguish opposite-direction results from one procedure being indeterminate. A small change around zero must not be promoted to strong biological reversal.

### Controls and secondary analyses

- Malignant-cell analysis is primary; whole tissue is a composition-sensitive secondary comparison.
- Compare equal-patient and pooled-cell summaries descriptively, acknowledging their different estimands.
- Prespecify matched-depth cell sampling and library-aware technical checks only after batch relationships are known.
- Use a finite procedure list, frozen before scoring. Report every procedure and failure.
- Do not add new gene panels, predictors, or clinical endpoints to rescue a failed primary result.

### Inference and precision

The independent unit is the patient. Cells, libraries, slides, and resampling draws never increase the patient count. Nested cell resampling describes within-specimen sampling uncertainty; patient-level resampling/inference addresses between-patient uncertainty. State the assumptions of the selected interval or test and address multiple primary comparisons if any.

With zero failures among 27 independent patients, the exact one-sided 95% binomial upper bound is `1 - 0.05**(1/27) = 0.105019...`. Thus 27 patients cannot establish failure below 10% at that confidence level, even in that best observed case. The cohort may shrink after QC. If dependence or shared preprocessing invalidates a simple binomial model, use an appropriate alternative rather than quoting this bound as universally valid.

Choose a scientifically justified reliability tolerance before outcomes, then calculate required precision/sample size. Do not choose the tolerance to fit available `n`. If meaningful confirmation is infeasible, predeclare descriptive/exploratory analysis and report intervals without a manufactured pass/fail win.

Freeze the sample list, analysis plan, software/config/gene hashes, seed schedule, multiplicity, uncertainty method, runtime limit, and failure rules in a commit before biological scoring. Do not describe a Git commit alone as independent preregistration.

## 10. Phase E: implement and verify before scoring

Prefer a small transparent pipeline over a new neural architecture. Use existing repository helpers where suitable after checking what they actually implement. Avoid mutating historical runners or outputs.

Suggested isolated layout:

```text
src/recurrence_reliability/
    metadata.py
    io.py
    qc.py
    scoring.py
    inference.py
    run.py
configs/recurrence_reliability.json
tests/test_recurrence_reliability_*.py
output/recurrence_reliability/<run_id>/
```

Adapt names to repository conventions. Implement separate metadata/QC/score commands so running QC cannot accidentally generate outcomes. Default to development/synthetic inputs; require an explicit locked config for biological scoring. Refuse missing annotations, unresolved duplicate IDs, hash drift, or mismatched specimen/barcode joins instead of silently filling them.

Verify meaningful failure cases with synthetic data: correct specimen/timepoint joins; duplicate-barcode/library handling; count/gene-ID checks; correct recurrence-minus-primary sign; unchanged results under cell-row reordering; no patient counted twice; ties/missing/indeterminate denominator behavior; and inference on hand-computable examples. Check scorer behavior against its published definition. Do not write tests that merely restate implementation.

Store run provenance, environment versions, input/config/code hashes, seeds, exclusions, warnings, and all declared analyses. Do not regenerate historical tests or consume external data to debug.

## 11. Phase F: execute once and interpret honestly

After authorization and all applicable gates pass, run the locked primary analysis. Read output artifacts back. Confirm fields, patient counts, unique IDs, score direction, finite values, interval bounds, denominators, exclusions, and hashes are sane. Compare outputs with the frozen protocol, not expectations of a positive result.

Keep any failed or inconclusive result. An execution bug requires a dated correction, regenerated affected outputs, and disclosure of outcome exposure; it does not restore untouched status. Do not relabel an exploratory result as confirmatory or repeat variants until one wins.

Required outputs: patient-level change table; primary endpoint and interval; indeterminate coverage; technical/batch diagnostics; complete exclusion flow; secondary results marked as such; limitations; reproducible commands; and an updated claim/evidence ledger. Figures should show patients and uncertainty, not thousands of cells as independent evidence.

## 12. Phase G: external validation is a separate gate

Potential lead: CARE, associated with Spitzer et al. and the recent trajectory paper. Official metadata currently lists 121 GSMs at `GSE274546` (10X) and one GSM at `GSE274548` (Smart-seq2). Published study count is 59 patients; that is not a verified eligible-pair count in this project. Two accessions do not imply two independent cohorts.

Before opening external expression/outcomes, audit paired identifiers, assays, malignancy annotations, file terms, source origins, overlap with Wang and all prior project data, and prior exposure. Classify external data as fresh, previously used, or uncertain. Reserve one eligible independent source before tuning. Published aggregate findings may be known; disclose them and distinguish that knowledge from project outcome exposure.

Apply the frozen gene set, scoring, QC, primary endpoint, and inference without outcome-driven adjustment. Prespecify allowed platform mapping. If calibration is required, separate calibration patients from the final test. No external source means no external-replication claim; an honest within-cohort reliability benchmark can still be the output.

## 13. Publication package and relation to A/B/C

Organize one paper around one claim. Candidate framing: reliability of patient-level GBM recurrence-program measurements under technical variation. Introduction motivates the measurement problem; methods define the patient-level target; results report cohort flow, reliability, uncertainty, and any independent replication; discussion separates technical limitations, construct disagreement, and biological interpretation.

Do not force the MRI forecast and simulated controller into the same validated pipeline. They may motivate the broader research program or appear as separate studies. No inferred molecular-to-MRI patient join, clinical treatment recommendation, causal resistance claim, or confirmed state-transition dynamics is justified here.

Prepare a draft manuscript only after the actual evidence is known, with retained negative results and properly limited novelty. A reproducibility package should contain scripts, frozen configs, gene-list provenance, aggregate/allowed outputs, environment instructions, source citations, and an analysis history. Respect source data-use and redistribution terms.

## 14. Alternatives if the selected extension fails

Do not launch every alternative automatically. Give a ranked recommendation based on the failed gate.

- **Independent-mechanism control falsification:** test EARF using genuinely held-out tumor/observation mechanisms, beyond existing parameter sensitivities. Strong feasible simulation option, but separate protocol and current Track C run/precision gates still apply. No run authorized by this handoff alone.
- **Spatial RNA–protein validation:** GSE237183/Zenodo current release 12624860 is open CC BY 4.0. RNA and CODEX use near-adjacent sections. Five donor-code candidates overlap, but no complete section-to-image crosswalk is verified. Keep exploratory/conditional; the old regional benchmark's six-donor rule is not a universal rule for a new endpoint.
- **Empirical MRI repeatability:** valuable for Track B forecast interpretation. RIDER has 19 repeat-scan patients in the official description but currently points to controlled image access. Do not call it unrestricted or divert CFB/WAW outcomes.
- **Experimental perturbation:** a longer-term biological anchor, conditional on independent origins and measured assays; expression perturbation is not clinical efficacy or a drug kill rate.

## 15. Stop rules and completion checklist

Stop the proposed primary study before scoring if no distinct question remains, annotations cannot be linked without circularity, independent patients cannot be established, library duplication cannot be resolved, access/terms fail, or the planned estimand cannot be supported. Where feasible, finish a clearly labeled smaller descriptive question through a dated plan rather than silently changing the endpoint.

The handoff execution is complete when you deliver either:

1. A gated, verified study with frozen protocol, reproducible implementation, declared results, limitations, and a manuscript/reproducibility package; or
2. A documented feasibility/novelty no-go with exact evidence, completed safe independent work, and the strongest concrete next recommendation.

At every stopping point, leave the current status, commands, hashes, findings, unresolved gates, and exact next action in a dated handoff. Commit and push verified relevant changes. Report what you actually found, not what the proposal hoped to find.
