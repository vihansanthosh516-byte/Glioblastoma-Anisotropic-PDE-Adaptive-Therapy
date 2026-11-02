# Track A and Track B audit and gap closure (2026-11-01)

**Scope.** Track A and Track B only. Track C was not inspected. The two tracks use different cohorts; no patient-level link between them exists or is implied.

**Method.** Lead agent (Claude Opus 5.5) plus five read-only subagents (Sonnet 5.5): Track A repository/protocol, Track A novelty/evidence, Track B repository/protocol, Track B methods/code, and a red team. The lead checked each finding used below against files on disk. Tags: **[O]** observed in a file or command output, **[R]** reported by a project document or web page, **[I]** inferred, **[P]** proposed.

**Limits on this audit.** No expression matrix, MRI, mask, or image was downloaded or opened. The local Yu 2020 Table S1 workbook was opened for A-G16. An initial structure listing displayed the first six rows from both sheets, including non-target clinical and point/sample metadata; no expression values were read. The subsequent IDH crosswalk scanned the `PatientID` column to locate GS2, GS3, GS5, and GS6, then read `IDH1 R132 mutation` only for those four rows. No patient was scored. No one was contacted. No figure was changed. One metadata-only web call was made (Zenodo record API, A-G10). One synthetic test file ran in a scratch copy (B-G6). NCBI GEO and BioStudies pages were not readable by the subagent, so GEO/ArrayExpress counts stay [R].

**Running work from another session (not touched).** A synthetic-only PDE timing job (`gap_stress.py`, another Claude session's scratchpad, started 2026-11-01 10:58) is running the last history gap (155 weeks). It will write `output/track_b_timing_pilot_20261031/lumiere_history_gap_stress.json`. The uncommitted r2 of `paper/track_b_lumiere_scoring_amendment_draft_2026-11-01.md` holds `TIMING_PLACEHOLDER` in §5a for that result. This audit did not edit, stage or commit that draft. **Interim, not for citation:** the job's log shows calibration-only PDE times of 691–3,077 s per synthetic case across gaps of 1–23 weeks, against 1,544 s in the pilot used for the 677.6 h projection. PDE time is therefore not a stable lower bound.

---

## 1. Current status, established independently

### Track A
- **[O]** A5 (Yu 2020, GSE117891; six donors; nested leave-one-donor-out) finished. Candidate mean class-balanced AURC 0.1479 vs L2 logistic 0.1764; lower in 4/6 donors; the frozen 5/6 rule failed (`output/trackA_yu2020/result.json`, `frozen_descriptive_reading`). Brier was worse (0.2335 vs 0.1944). The no-gate ablation matched the candidate (0.1473). Cell-type-only balanced accuracy 0.6546 vs 0.6789 for the best expression model: comparable. 2,519 primary cells scored.
- **[O]** Three sealed bulk tests failed their adjusted rules: SenMayo p 0.0035 → adjusted 0.37; OXPHOS 0.0024 → 0.67; TAM 0.0039 → 0.44 (`output/trackA_senmayo|oxphos|tam_5ala/result.json`).
- **[R]** No public cohort passes the six-donor, label, map, terms and independence gates. Closest: GSE286419 with four Core/Margin pairs (`paper/trackA_alternate_dataset_plan_20261101.md`).
- **[O]** Seal hashes in `paper/trackA_*/SHA256SUMS` and the A5 `source_hashes` match the current working-tree files.
- **Plain reading:** Track A has no positive or independent result. It is a well-documented negative and methods record.

### Track B
- **[O]** The LUMIERE exploratory plan and amendment R (U-Net reuse) are adopted (`paper/track_b_no_contact_dataset_plan_2026-10-31.md`, line 3). The scoring amendment is a draft, not adopted.
- **[O]** Header-only check: 28 eligible patients, 84 visits, one grid, 0 problems (`output/track_b_lumiere_inventory_20261031/header_check.json`). The draft's fold table and gap counts reproduce from this file.
- **[O]** Timing arithmetic is correct: 223.25 h U-Net and 57.77 h PDE per analysis; 1,124.1 h and 677.6 h for four analyses (`output/track_b_timing_pilot_20261031/lumiere_timing_record.json`). These are 8-thread wall-clock hours.
- **[O]** No LUMIERE loader, runner, freeze file or gate code exists. The MU runner is hard-wired to 18 patients, day 90 and a 76–104-day window. Every LUMIERE real-data gate is open or blocked.
- **[O]** LUMIERE is an exposed cohort: spent H-2 test set (58 LUMIERE patients, ledger 2ax), scripts 91, 104, 106, 130, 134, and PREDICT-GBM.
- **Plain reading:** Track B has no new patient result. The only patient evidence is the older MU run, which was inconclusive against persistence and geometric train-rate (`output/pde_manifest_monotone/results.json`).

---

## 2. Gap register

Class: **FIX** = fix now (done or doable), **DEC** = needs user decision, **EXT** = needs outside evidence, **NO** = cannot be resolved under current constraints. Effort: S ≤ 2 h, M ≤ 1 day, L > 1 day. Status: **Done** (commit given), **Open**, **Blocked**.

### Track A

| ID | Gap and evidence | Impact | Fix | Class | Depends on | Acceptance | Status | Effort |
|---|---|---|---|---|---|---|---|---|
| A-G1 | A5 headline lacked its weakest points: Brier worse, no-gate ablation equal, GS13 logistic AURC 0.011 vs 0.073 [O result.json] | Result looked closer to success than it was | Added to `START_HERE.md` and ledger 2by | FIX | — | Numbers match result.json | Done (this commit set) | S |
| A-G2 | "4,139 cells passed QC" read as 4,139 scored; 2,519 were scored [O] | Count misread | Stated both in `START_HERE.md` | FIX | — | Sum of per-donor n = 2,519 | Done | S |
| A-G3 | `code_commit_before_run` in result.json was corrected by hand after generation; runner records write-time HEAD [O `paper/trackA_execution_audit_20261031.md` line 17] | Breaks the "never hand-edit outputs" rule; disclosed | Keep the disclosure; do not rerun (spent). Runner `src/189…` is seal-locked, so its HEAD logic is not changed. Future runners must record the commit at run start | FIX (doc) / NO (rerun) | — | Disclosure present in audit and ledger 2by | Done (disclosure) | S |
| A-G4 | Recorded hashes are of CRLF working-tree bytes. Committed LF blobs hash differently: `src/189…` blob f34fa7d5… vs recorded e20a74f9…; CRLF conversion reproduces e20a74f9… [O] | A clean clone on Linux or with `autocrlf=false` fails hash checks | Record here. Verify on Windows with `autocrlf=true`, or convert LF→CRLF before hashing. Add `.gitattributes` only by user decision, since it changes checkouts | FIX (doc) / DEC (.gitattributes) | — | This note exists | Done (doc) | S |
| A-G5 | `paper/trackA_null_result.md` title and "both shifts vanish" overclaim; says "two tests" [O] | A reviewer would flag it | Correction box added; sealed protocols untouched | FIX | — | Box present | Done `698c17b` | S |
| A-G6 | Sealed gene lists were chosen after the script-157 look; no single register of all Track A hypotheses and screens [O START_HERE 113–118] | Forking-path attack | Register in §3 below; called "hypothesis-informed" in null_result note | FIX | — | §3 lists every test and screen | Done | S |
| A-G7 | TAM test scored 5 of 6 TAM genes and 6 of 8 myeloid genes; not disclosed [O `n_genes`] | Frozen list ≠ scored list | Disclosed in `paper/trackA_composition_finding.md`, START_HERE, ledger. Which genes dropped is not recorded; finding out needs the matrix row IDs | FIX / DEC (gene identity) | Data-read approval | Note present | Done (disclosure) | S |
| A-G8 | Execution audit said the Ivy matrix "remains unopened", but tests 1–2 scored it [O] | Contradiction | Clarified | FIX | — | Line reworded | Done `698c17b` | S |
| A-G9 | Alternate plan Aim cites an "already frozen malignant-cell signature" that does not exist [O] | Misstates what exists | Corrected | FIX | — | Correction present | Done `698c17b` | S |
| A-G10 | Zenodo 20709378 licence: one doc "blank", another "CC BY 4.0" [O] | Conflict | API checked 2026-11-01: `cc-by-4.0`; doc updated | FIX | — | Doc cites API | Done `698c17b` | S |
| A-G11 | Roadmap Track A next action and dataset-audit "no result.json" were stale [O] | Fresh session could redo done work | Updated | FIX | — | Lines updated | Done `698c17b` | S |
| A-G12 | Script-157 shuffle control described as "about 0.5"; cell-type-only shuffle is 0.427 (donor range 0.248–0.758) [O `audit_results.json`] | Wrong control value | Corrected in START_HERE with exact values | FIX | — | Values match JSON | Done | S |
| A-G13 | No Track A ledger entries for A5 and the sealed tests [O] | Numbers lacked a ledger citation | Ledger 2by added | FIX | — | Entry cites result files | Done | S |
| A-G14 | No unit tests for sealed runners 173, 174, 185, 189; runner 189 does not use the tested `src/donor_cv.py` [O] | Success-rule and fold code untested | Added synthetic-only tests for shared `173` sign-flip/scale helpers used by 174/185, the 185 loader with synthetic CSVs, and 189 metric/tie helpers. No real inputs used; runner files stay unchanged | FIX | — | Focused tests pass in scratch copy; full runner coverage remains open | Partial | M |
| A-G15 | Shared appendix is pinned in the A5 seal, and Track B edits already broke it once (amendment A4) [O] | Track B edits can break Track A verification | Rule: do not edit `paper/research_plan_complete_project_appendix.md` without re-running the A5 hash check, or freeze a Track A copy | DEC | — | Rule recorded in START_HERE | Open | S |
| A-G16 | Yu 2020 main text identifies GS1 and GS13 as IDH-wildtype ([Oxford Academic](https://academic.oup.com/nsr/article/7/8/1306/5849029); [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8289159/)). Local Table S1 (`data/external/yu2020_supp/Table S1. ClinicalInformation.xlsx`, sheet `Clinical_Info`, column F `IDH1 R132 mutation`) records `-` for CGGA_GS2 row 3, GS3 row 4, GS5 row 6, and GS6 row 7; workbook SHA256 `0676607ff4c57348fc043ab4bace1c3985617d0d90507742344ce19642cb772a` matches the A5 `result.json` input hash | The table records no reported IDH1 R132 mutation for those four donors; comprehensive IDH1/2 status is not established by this field | Report the exact marker and do not use these rows alone to claim full IDH1/2-wildtype status | FIX | — | Per-donor cells and the source limitation are recorded | Done (source-limited) | S |
| A-G17 | The Track A gap register lists calibration curve, within-library, scANVI, and attention analyses beyond A5; `result.json` has no outputs for them [O] | Plan-to-result gap | Keep all four marked **not run; deferred**; do not run or remove them now. Any change requires a future dated amendment | DEC | User direction (resolved) | Audit and execution report state the deferred status and amendment rule | Done (deferred) | S |
| A-G18 | SVM is a constant-score control (AURC 0.25 in 5/6 donors; BA 0.5 in all 6); point-size is constant-score in 4/6 donors [O result.json] | Inflates the model count | Label controls and report point-size fold variation | FIX | — | Execution audit gives donor-specific control status | Done (this commit) | S |
| A-G19 | Script-157 result came from the pre-identity-gate code; SNU17/SNU33 merge unverified [O] | Composition numbers rest on unverified donor IDs | Label as such (done in START_HERE); rerun needs donor evidence | FIX (label) / EXT | Source donor map | Label present | Done (label) | S |
| A-G20 | No independent cohort with ≥ 6 paired IDH-wt donors passes gates [R] | No confirmation possible | Wait for a named accession; do not relax the six-donor gate | NO | — | — | Blocked | — |
| A-G21 | HRA009021 (10 donors) lacks a public sample-to-region map; only the submitter has it [R] | Only 10-donor lead is blocked | Needs contact, which the no-contact rule forbids | DEC | User | — | Blocked | — |
| A-G22 | Cell Browser export reuse terms unknown [R] | Sharing risk | Keep derived data local; cite only | EXT | — | — | Open | — |
| A-G23 | GSE286419 counts use different denominators: GEO lists 6 snRNA-seq donors/12 records; paper Methods reports 8 surgical patients (6 snRNA-seq, 2 spatial); paper Results retain 4 matched pairs after HGG3/HGG5 exclusion [official GEO and published paper] | Near-miss count looked inconsistent | Record all three denominators and exact exclusion in the Track A plan | FIX | — | Plan states 8 total, 6 GEO donors, 4 retained pairs | Done (this commit) | S |

### Track B

| ID | Gap and evidence | Impact | Fix | Class | Depends on | Acceptance | Status | Effort |
|---|---|---|---|---|---|---|---|---|
| B-G1 | Plan listed only scripts 91 and 130 as prior LUMIERE use. LUMIERE was also the spent H-2 test set (ledger 2ax), used by 104/106/134, and in PREDICT-GBM [O] | Independence overstated by omission | Exposure list added to plan line 7 and START_HERE | FIX | — | List present | Done `08a6bd3` | S |
| B-G2 | Appendix §1.2 still says LUMIERE is "exposure disclosure only; not a planned analysis"; no LUMIERE terms row; no frozen exposure manifest [O] | Plan contradicts appendix | Appendix amendment. Appendix is seal-locked by Track A (A-G15), so edit only with a Track A re-verify | DEC | A-G15 | Dated amendment | Open | S |
| B-G3 | No document is named as the LUMIERE scoring protocol; MU protocol draft names a different primary contrast (hybrid vs matched-volume no-PDE) than appendix §3.3 (four required mask baselines) [O] | Two "primaries" | Adopt the scoring draft (or a revision) that names the governing text | DEC | B-G12 | Freeze file names one protocol | Open | S |
| B-G4 | "CPU hours" are 8-thread wall-clock hours; 677.6 h ≈ 28 days [O arithmetic] | Feasibility hidden | Relabelled in plan, START_HERE, roadmap | FIX | — | Units stated | Done `08a6bd3` + this commit | S |
| B-G5 | 677.6 h omits analyses the draft adds: history-gap sensitivity (~+250 h, [I]) and option-B time-aware baselines (~+447 h). Max-step reuse total ≈ 819 h [I] | Real cost ≈ 1,300 h if all adopted | Shown in §6 below; draft owner should add a cost table | FIX (doc) | Running timing job | Cost table in draft | Open (draft owned by other session) | S |
| B-G6 | MU result provenance did not hash `src/fk_fixed_growth.py` or `src/fk_learned_baselines.py` [O] | Baseline code could change after freeze and still pass hash checks | Added both; test asserts them; 14/14 synthetic tests passed in scratch copy (Python 3.14.4) | FIX | — | Test passes | Done `4170221` | S |
| B-G7 | No LUMIERE loader/runner/gate code; MU runner fixed at 18 patients, day 90, 76–104 days, central gap 78–102 days [O] | LUMIERE gates enforced nowhere | Build a `LumiereCohort` loader that refuses to load masks without QC list hash, freeze file and ≥ 25 passes | FIX (after B-G3) | B-G3, B-G12 | Synthetic refusal tests | Open | L |
| B-G8 | Scorer raises if target lies outside support; spec says count as a miss [O `src/fk_mu_scoring.py:298`] | Safe for MU (full-grid support) but would abort a LUMIERE run (origin-only support) | Implement miss-counting in the LUMIERE path only; MU behaviour stays | FIX (with B-G7) | B-G7 | Test with out-of-support target | Open | S |
| B-G9 | No "PDE-only" arm in code; code has a no-PDE comparator instead [O] | Required baseline missing | Define PDE-only exactly, then code it | DEC | B-G3 | Definition frozen | Open | M |
| B-G10 | Pass rule not computable: no copy-last contrast, no max-T, separate bootstrap seeds per contrast [O] | Cannot declare pass/fail | One function: patients × contrasts delta matrix, shared resamples, max-T intervals | FIX | B-G3 (contrast family) | Synthetic test of max-T | Open | M |
| B-G11 | No checkpointing; a ~678 h run lives in one process [O] | One crash loses weeks | Persist per-fold models and forecasts keyed by fingerprint | FIX | — | Resume test | Open | M |
| B-G12 | Draft decisions open: pass-rule reading, fold rule, time mismatch A/B, outliers, minimum 25, preflight read scope [O draft §10] | Freeze cannot happen | User decides | DEC | — | Decisions recorded | Blocked on user | — |
| B-G13 | Fixed-growth spec reads as governing LUMIERE; stale "no timing pilot" blocker [O] | Wrong rules could apply | Scope note added | FIX | — | Note present | Done `08a6bd3` | S |
| B-G14 | Source-gate record said the school "has not been contacted" next to the user-reported OK [O] | Self-contradiction | Reworded; written record still missing; LUMIERE not named | FIX (wording) / EXT (record) | User | Written approval saved | Wording done `08a6bd3`; record open | S |
| B-G15 | Visit 1 → visit 2 alignment is not QC'd; label warp rule (nearest-neighbour, 1→2 mm) not frozen [O] | PDE and U-Net use v1 masks | Add v1 CT1 to preflight scope (new approval) or declare limit; freeze warp rule | DEC | User | Rule in freeze | Open | M |
| B-G16 | Label gate: DeepBraTumIA version and label 0 meaning unstated; no cavity label; prior LUMIERE mask-QC problems (ledger 2ab, 2aj) not in plan limits [O] | Dice measures agreement with one automated tool | Cite 2ab/2aj in limits; keep as stated limit | EXT | Published docs (NITRC page unreachable) | Limit stated | Open | S |
| B-G17 | Licence: non-commercial; derived-output terms unstated [R] | Sharing risk | Summary numbers only | EXT | — | — | Open | — |
| B-G18 | Real-patient registration QC not run; thresholds from phantoms [O] | Cohort unknown | 84-file preflight after written user OK | DEC | User | `qc.json` saved | Blocked | M |
| B-G19 | Registration NCC can drop with real tumor growth, so QC may remove growing patients [O readiness audit] | Selection bias | Report excluded count with week gap and history volume class (outcome-blind) | FIX (in preflight) | B-G18 | Exclusion table | Open | S |
| B-G20 | Inventory JSON `rules` text says "13 weeks ~ 91 days", conflicting with draft L3 [O] | Minor wording | Output JSON: do not hand-edit; fix the generator wording on next header run | FIX (code) | — | — | Open | S |
| B-G21 | Prior-art checks incomplete (Scholar, citation sweeps) [O prior-art follow-up] | Novelty gap unconfirmed | Run searches; log queries | EXT | — | Search log | Open | M |
| B-G22 | No untouched multi-visit cohort without an access request [O external screen] | No confirmatory test possible | Duke stays pending; no download | NO | Duke + school approval | — | Blocked | — |
| B-G23 | Week-000 scans: pre-op or immediate post-op is not documented; "first post-op triplet" label unverified [O] | Eligibility meaning | Find LUMIERE documentation statement | EXT | — | Quote saved | Open | S |

---

## 3. Track A register of all hypotheses, tests and screens (A-G6)

| # | Item | Data | Decided before scoring? | Outcome | Source |
|---|---|---|---|---|---|
| 1 | Random-cell vs patient split (C-GAT) | Old atlas | No | Leak shown; patient-level 0.638 ≈ cells-per-patient 0.646 | START_HERE §3 |
| 2 | Composition audit (script 157) | 7 SNU donors | No (exploratory) | Genes AUC 0.904 vs cell type 0.713; region = library confound | START_HERE line on script 157 |
| 3 | Spatial same-section core/edge | Greenwald, Pai, Ravi | Screen | Paused: no same-section labels | START_HERE "spatial PAUSED" |
| 4 | MES shift, Greenwald distance, Gill GSE59612, BMDM/microglia ratio | various | Screen | Rejected (prior art or blocked) | START_HERE "new-question search" |
| 5 | TMZ-tolerance program (arena winner) | — | Gate 1 | Dropped: no gene list | same |
| 6 | SenMayo senescence (sealed) | Ivy 8 + Angers 7 | Yes (hash) | Adjusted rule failed | `output/trackA_senmayo/result.json` |
| 7 | Hallmark OXPHOS (sealed) | same 15 | Yes | Adjusted rule failed | `output/trackA_oxphos/result.json` |
| 8 | TAM activation (sealed) | E-MTAB-8743, 10 | Yes | Adjusted rule failed | `output/trackA_tam_5ala/result.json` |
| 9 | A5 candidate vs logistic | Yu 2020, 6 donors | Yes (A5) | 5/6 rule failed | `output/trackA_yu2020/result.json` |
| 10 | Ivy GAP malignant-cell preflight | metadata | Screen | Too few donors | `paper/trackA_execution_audit_20261031.md` |
| 11 | External cohort screens (HRA009021, GSE286419, GSE318562/4, SNUH, Harwood, Petti, EGA, others) | metadata | Screen | None passes all gates | `paper/trackA_alternate_dataset_plan_20261101.md` |

No test passed its full rule, so no multiplicity correction is needed for a claim. Tests 6–8 were hypothesis-informed by item 2.

---

## 4. Completed fixes and verification

| Commit | Files | What changed | Verification |
|---|---|---|---|
| `4170221` | `src/fk_mu_patient_runner.py`, `tests/test_fk_mu_pipeline_components.py` | B-G6: two baseline modules added to result code hashes | Scratch copy, Python 3.14.4: `tests/test_fk_mu_pipeline_components.py` 14 passed. No patient data |
| `08a6bd3` | `paper/track_b_no_contact_dataset_plan_2026-10-31.md`, `paper/track_b_source_gate_evidence_2026-10-31.md`, `paper/track_b_fixed_growth_baseline_spec_2026-10-31.md` | B-G1, B-G4, B-G13, B-G14 | Diff read back; each replacement matched exactly once |
| `698c17b` | `paper/trackA_null_result.md`, `paper/trackA_composition_finding.md`, `paper/trackA_alternate_dataset_plan_20261101.md`, `paper/trackA_gse318562_no_contact_plan_20261101.md`, `paper/trackA_execution_audit_20261031.md`, `paper/trackA_dataset_audit_20261031.md`, `paper/trackABC_execution_plan_2026-10-31.md` | A-G5, A-G7, A-G8, A-G9, A-G10, A-G11 | Diff read back; none of these files is in any `SHA256SUMS` |
| this commit | `paper/trackAB_audit_and_gap_closure.md`, `START_HERE.md`, `vault/PAPER_FINDINGS_LEDGER.md`, `paper/trackABC_execution_plan_2026-10-31.md` | A-G1, A-G2, A-G7, A-G12, A-G13, A-G19 label, B-G1/B-G4 in START_HERE, ledger 2by and 2bz | Numbers re-read from result JSONs (2,519 cells; shuffle 0.427/0.522/0.506; AURC 0.1473/0.1479) |

Seal check after all edits: no file listed in `paper/trackA_*/SHA256SUMS` was changed (`PROTOCOL.md` files, gene lists, `donors.tsv`, `src/189_trackA_yu2020_lodo.py`, the appendix).

---

## 5. Blocked items and exact needs

| Item | Exact decision or evidence needed |
|---|---|
| B-G12 draft decisions | Your answers to the six choices in `paper/track_b_lumiere_scoring_amendment_draft_2026-11-01.md` §10, after the other session fills §5a |
| B-G18 preflight | Written OK to read 84 files: visit-2 and visit-3 atlas CT1, visit-2 atlas brain mask, 28 patients (`paper/track_b_readiness_audit_2026-11-01.md` §5). If B-G15 is chosen as "QC v1 too", the scope grows to 112 files and needs its own OK |
| B-G14 school | A saved written approval with date, approver, conditions, and the LUMIERE dataset named |
| B-G16, B-G23 labels | Published DeepBraTumIA version, label-0 meaning, week-000 meaning. If none exist, accept as stated limits |
| B-G2, A-G15 appendix | Your choice: amend the shared appendix (then re-verify the A5 seal) or keep LUMIERE rules only in the LUMIERE amendment |
| A-G20, A-G21, B-G22 cohorts | A named new accession (Track A) or Duke + school approval (Track B). Contact is not allowed under current rules |
| A-G7 | Permission to read matrix row IDs (no values) for gene identity |

---

## 6. Dependency order and parallel lanes

```
Lane A (Track A, independent of B)
  A-G14 runner tests ─┐
  A-G16 Table S1 refs ├─> one-page Track A negative/methods write-up (story clarity)
  A-G17 endpoint note ┘

Lane B (Track B)
  [running] synthetic gap timing ─> draft §5a filled (other session)
        └─> B-G12 user decisions ─> B-G3 governing protocol ─┬─> B-G9 PDE-only definition
                                                             ├─> B-G10 max-T function
                                                             └─> B-G7 LUMIERE loader + B-G8 miss rule
  B-G11 checkpointing (parallel, no dependency)
  B-G18 preflight (needs written OK) ─> B-G19 exclusion table ─> cohort freeze
  freeze (code + cohort + config hashes) ─> any scoring decision
```

Lanes A and B can run in parallel. Within B, B-G11 and the B-G10 function body can start now; their contrast list waits on B-G3.

**Cost if every drafted Track B option is adopted [I]:** 677.6 h (adopted) + about 250 h (history-gap sensitivity, own U-Net) + about 447 h (time-aware baselines) ≈ 1,375 wall-clock hours ≈ 57 days on this PC. With the max step instead of the mean, higher.

---

## 7. Quality scores (0–10)

Scores are the lead's reconciliation of five subagent proposals and direct checks. "Before" = state at commit `e4df62c`. "After" = after the fixes above.

### Track A

| Dimension | Before | After | Evidence | Confidence | What would change it |
|---|---|---|---|---|---|
| Novelty and uniqueness | 2 | 2 | Core/periphery cell differences and the classification task are published (Darmanis 2017, PMC5810554; Lopes & Vinga 2020, PMC7029554; Pai 2025, PMC12603948; Sonpatki 2026, 10.1038/s41467-026-69716-2). Only narrow angle: donor-held-out benchmark | High | Only a new independent ≥ 6-donor malignant-cell test. **8 is not reachable with current data** |
| Strength of evidence | 2 | 2 | Every primary rule failed; cohorts exposed; 6–15 donors | High | Same as above. **8 not reachable** |
| Scientific rigor | 7 | 7 | Hash seals, nested donor CV, frozen rules; untested sealed runners, CRLF-dependent hashes | Medium | A-G14 tests + CRLF hash rule + A-G16 sources → 8 |
| Claim validity | 6 | 8 | Overclaim wording fixed (A-G5, A-G9); weakest points added (A-G1) | Medium | Stays 8 if write-ups follow §1 wording |
| Honesty about limits | 8 | 9 | All failures reported; corrections disclosed; register in §3 | High | — |
| Story clarity | 4 | 6 | Register (§3) and single status; still many strands across files | Medium | One-page negative/methods write-up → 8 |
| Feasibility and reproducibility | 6 | 6 | 2.8 h run; hashes recorded; no environment lock; Python 3.14 only; CRLF hash dependency | Medium | Pinned environment + one-command hash verifier → 8 |

### Track B

| Dimension | Before | After | Evidence | Confidence | What would change it |
|---|---|---|---|---|---|
| Novelty and uniqueness | 2 | 2 | SoC-DT, PREDICT-GBM and patient RD forecasting exist; LUMIERE inside PREDICT-GBM; gap unconfirmed (B-G21) | Medium | Completed prior-art search could move it ±1. **8 not reachable** on exposed public data |
| Strength of evidence | 1 | 1 | No new patient result; MU result inconclusive | High | A scored LUMIERE run gives exploratory evidence only; ceiling about 3–4. **8 not reachable** without an untouched cohort |
| Scientific rigor | 5 | 5 | Strong gate design and synthetic tests; but no LUMIERE gate code, no max-T, no PDE-only arm, phantom thresholds | Medium | B-G3, B-G7, B-G9, B-G10, B-G11 done and freeze hashed → 7–8 |
| Claim validity | 6 | 7 | Exposure list added; units fixed; no overclaim found | Medium | Governing protocol named (B-G3) → 8 |
| Honesty about limits | 7 | 8 | Exposure, cost and contradictions now stated | High | — |
| Story clarity | 4 | 5 | Many overlapping drafts (MU, UCSD, CFB, Duke, LUMIERE); START_HERE now points here | Medium | One adopted LUMIERE protocol and archive banners on MU docs → 7 |
| Feasibility and reproducibility | 3 | 3 | ~28 days continuous minimum; up to ~57 days with drafted options; no checkpointing | High | Reduced scope (§8 option) + checkpointing → 6–7 |

---

## 8. Route toward 8, and honest scope where 8 is not reachable

- **Novelty and evidence, both tracks:** 8 is not achievable under current data and the no-contact rule. Every available cohort is exposed or below the donor gate. Honest scope: Track A = "a sealed negative and leakage/composition methods study"; Track B = "an exploratory, pre-registered re-analysis of an exposed public cohort with strong baselines".
- **Track A rigor, story, reproducibility to 8:** A-G14 tests (M), A-G16 source rows (S), CRLF hash rule (S), a one-page write-up (M), a pinned environment file and hash verifier (M). No new data needed.
- **Track B rigor and claim validity to 8:** user decisions (B-G12), then B-G3, B-G9, B-G10, B-G7/B-G8, B-G11, all synthetic-tested, then a hashed freeze. About 3–5 days of coding.
- **Track B feasibility:** two honest scopes for the user [P]:
  1. **Full adopted design** (677.6 h, about 28 days; more with drafted options). Needs checkpointing first.
  2. **Reduced design:** run QC first; then one analysis (primary prior support only) with copy-last, fixed growth, PDE-only and the U-Net at full recipe; central-gap and sensitivity runs only if time remains. About 281 h (223.3 + 57.8), about 12 days. This changes the adopted four-analysis rule, so it needs a dated amendment before any target read.

---

## 9. Remaining authorized next actions (no data, no contact)

1. A-G14: synthetic tests for sealed Track A runner functions (code untouched).
2. B-G11: checkpointing in the Track B runner, synthetic resume test.
3. B-G10: max-T shared-resample function with synthetic tests (contrast list as a parameter).
4. B-G20: fix the header generator's "~91 days" text (output JSON regenerated only on its next metadata run).
5. B-G16: add ledger 2ab/2aj mask-QC history to the LUMIERE plan limits.
6. B-G21: prior-art searches with logged queries (metadata only).
7. After the other session's timing job finishes: it fills draft §5a; then present the §10 decisions.

**Verification statement.** Read back after writing: every changed file's diff; result JSON values quoted in §1, §3 and §7; seal files unchanged. No restricted data were opened, no contacts were made, no figures changed, and the running job and uncommitted draft were not touched.
