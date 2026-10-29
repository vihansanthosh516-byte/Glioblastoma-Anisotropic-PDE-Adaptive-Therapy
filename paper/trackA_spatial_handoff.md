# Track A handoff: spatial core-vs-edge test (option 2)

Written 2026-10-28 by Claude for a new chat. **Track A only.** Do not touch Track B (Codex owns it) or Track C (another chat runs it).

## Start here
1. Run `git log --oneline -5` and `git status --short`. Codex may have changed files.
2. Read this file, then the 9 team notes in `research_notes/trackA_spatial/`. That folder is Git-ignored, so the notes are local only.
3. Read the Track A lines of `START_HERE.md`: lines 13-31, 91 and 113-132.

## User rules for this work
- Report every result with 4 answers: (1) Is it positive? (2) How positive or negative, with numbers and the source file? (3) Does it make the project more unique, and how was that checked? (4) What fixes it?
- If a result is negative, first start a subagent that looks for legitimate ways to improve it. Put its ideas in the report.
- **No data download without the user's explicit OK.** For each file, give the file name, source URL and size.
- CPU jobs run on this PC, not on Kaggle.
- Do not overclaim. Label results exploratory or post hoc where that applies. Every number cites a file.
- Write in simple ASD-STE100 English.
- Commit format: "Script N: ...". Run `git status` first. Never `git add -A`. Push after each commit. Stage only your own lines in `START_HERE.md`, because Codex often has uncommitted edits there.

## Goal
Test this question: inside ONE tissue section, for the SAME cell type, does gene expression differ between tumor core and tumor edge? The region labels must come from histology or a pathologist, not from expression. Run it as a sealed test with the donor as the unit.

Why this test: the old script 157 result (genes AUC 0.904 vs cell type only 0.713) is confounded, because each region was its own sequencing library (`START_HERE.md` line 91). The old classifier result (0.638 vs a 0.646 baseline, line 21) is negative. **Track A has no positive result today.**

## What the 9-agent science team found (2026-10-28)
Notes are in `research_notes/trackA_spatial/`. The files are: sci-literature, sci-hypothesis, sci-experiment, sci-data, sci-stats, sci-coder, sci-debugger, sci-visualizer, sci-ethics.

**Facts the verifier agent confirmed on disk and on the web:**
- Ravi et al. 2022, *Cancer Cell* 40:639. Dryad doi:10.5061/dryad.h70rxwdmj, license CC0. `10XVisium_2.zip` is 7.87 GB. The current README is `README.md` (9.55 KB).
- Greenwald et al. 2024, *Cell*. Zenodo record 12624860, license CC BY 4.0.
  - 12 `GBM_*.tar.gz` files (380.8 MB) plus `GBM_zh1019inf.zip`.
  - `visium_dataset_description.csv` (10.7 kB) and `HE_scans.zip` (131.2 MB).
- Exact sign-flip test, two-sided p floor = 2/2^n. So n = 5 cannot reach p < 0.05.
- Smallest effect (d_z) with 80% power: n=6 → 1.8; n=8 → 1.2; n=10 → 1.05; n=15 → 0.8. Source: the stats agent's simulation, saved in the old chat's scratchpad, not in the repo.
- **Error found:** `START_HERE.md` line 91 said "13-80 periphery tumor cells". `output/trackA_composition_audit/audit_results.json` gives 14, 13, 14, 56, 9, 26 for the six donors (range 9-56), and 629 for SNU46. Line 91 is now corrected.

**Skeptic verdict: NO-GO for scoring and for downloads now.** The flaws, from worst:
1. **No dataset has verified pathologist core/edge labels inside the same section.**
   - Ravi's TC/TI sample codes suggest that core and edge are separate samples. That brings back the region = library confound.
   - Greenwald's layers and niches look expression-derived. If so, the test is circular.
2. **The 10-donor gate is probably out of reach.** The draft needs 3 development donors plus 10 sealed donors.
3. **The notes give 3 different primary endpoints.**
   - **Fix:** use one primary. It is the donor-level, signed, exact sign-flip test on ONE pre-set published signature, restricted to pure spots or cells. Drop the LOO-AUC endpoint.
4. **One-sided vs two-sided p is mixed up.** `sign_flip_p` in `src/trackA_composition_audit.py` (about line 98) is one-sided. Pick one convention.
5. **The six-gene myeloid set (SPP1, APOE, CD163, C1QA, HLA-DRA, CD74) came from the same atlas.** It sits next to stress and tumor genes, so the set is circular and possibly ambient RNA. Use a published signature that is frozen before any data are seen.
- **Other conflicts:** the composition thresholds (0.2 vs 0.8 purity) and the sign convention (edge − core vs core − edge) disagree across the notes.
- **Skeptic's odds** (a judgment, not computed): about 10-15% for a positive result, and about 5% for a result that is both positive and unique.

## Best lead found at the end (not yet checked)
**Pai et al., bioRxiv 10.1101/2024.11.05.621879,** "Spatial multi-omics defines a shared glioblastoma infiltrative signature at the resection margin".
- It studies contrast-enhancing core vs distal FLAIR margin from supratotal resections.
- It has snRNA-seq and snATAC-seq from the same Core and Margin dissociates in 6 patients, plus spatial transcriptomics. It reports EGFR as a top margin marker.
- This matches our question closely. It is also prior art, so it lowers uniqueness.
- Its data availability and license are **not checked**.

## Next steps (in order; none needs a download)
1. **Zero-download reading.** Read the Methods of Ravi 2022, Greenwald 2024 and Pai 2024. Answer 3 questions:
   - (a) Do core and edge share one section?
   - (b) Are the region labels drawn by a pathologist or from histology?
   - (c) How many donors are there?
2. **If (a) means separate slides, or (b) means expression-derived:** stop Track A spatial. Tell the user plainly. Suggest pausing Track A.
3. **If a dataset survives:** ask the user to approve only the small metadata files first:
   - the Ravi `README.md` (9.55 KB);
   - the Greenwald `visium_dataset_description.csv` (10.7 kB).
4. **Then fix the plan.** Use one primary endpoint, one p convention, a frozen published signature, one purity threshold, and one sign convention. Seal it in code with a SHA-256 hash before any scoring (the pattern is in `src/170_trackc_sealed_missed_scans.py`). Script number 172 is free.
5. **Before citing any count:** periphery cell-type counts sum to only 2,775-2,940 of the 3,000 sampled cells per donor (writer agent); explain the gap. Team summary: `research_notes/trackA_spatial/SUMMARY.md`.
6. **Before coding:** read `research_notes/trackA_spatial/sci-coder.md` (file plan) and `sci-debugger.md` (12 failure-mode tests). Write synthetic tests first.

## Open gates (from START_HERE)
- The Cell Browser reuse terms are unknown.
- The donor mapping is not source-verified (`START_HERE.md` lines 116-131).
- CC0 and CC BY on Ravi and Greenwald do not cover the Cell Browser export.
- Never join Track A data to the Track B MRI data.
