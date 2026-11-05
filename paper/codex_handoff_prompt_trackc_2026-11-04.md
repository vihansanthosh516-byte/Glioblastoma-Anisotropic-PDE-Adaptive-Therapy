# Codex handoff prompt: Track C (copy everything below the line into Codex)

---

You are taking over Track C of a glioblastoma (GBM) computational-biology project in `C:\Users\vihan\20206 science fair` (git repo, branch `main`, Windows). Work at maximum effort. Use as many parallel sub-agents/workers as your tools allow. Do not stop until the goal below is met or you can show it cannot be met.

## 0. Read first (in this order)
1. `git log --oneline -15` and `git status --short` (Claude may have changed files; verify on disk).
2. `START_HERE.md` (durable handoff; update it whenever the decision or next step changes).
3. `CLAUDE.md` (project rules; includes the Uniqueness rule below), `paper/trackc_overview_for_review_2026-11-04.md` (best summary), `paper/trackc_stage2_report_2026-11-04.md` (full log), `paper/trackc_anchored_spec_2026-11-04.md` (design), `paper/trackc_calibration_targets_2026-11-03.md` (published trial numbers), `vault/PAPER_FINDINGS_LEDGER.md` (rows C1 to C40).

## 1. The user's goal (non-negotiable)
"New biology" means a NEW METHOD that is clearly unique, like the user's hybrid method in Track A. It must contain or be benchmarked against ALL existing methods in the field, on the same data and the same score, and must beat them or add something none of them has. Done, verified or well written is NOT enough. Do not stop until such a method exists. If it does not beat the baselines, say so plainly and propose the next attempt. Never fabricate, never tune toward a positive result, never hand-edit output JSONs.

## 2. What is done (simulation only; all interim; development seeds)
Simulator `src/trackc_sim.py` (sensitive/resistant cells, 28-day cycles, 12-course budget, RANO-like and baseline-referenced progression calls). Tests A to M and G are finished and committed (scripts 190 to 205). Findings:
1. Matching the published RTOG 0525 standard-arm median PFS (5.5 mo, 143 to 186 d) plus the two negative controls (dose-dense HR 0.87, 12-vs-6 cycles HR 0.99) selects WEAK drug kill. 129 "controlled" worlds in 5 configurations (`single`, `single_rep`, `mixture`, `single_mix`, `mixture_mix`).
2. In those worlds NO scan-reading rule tested (raw on/off, nadir, response-guided stopping) beats the best fixed schedule at equal drug use by more than about 1 day (median gap; 164 gaps computed). Scan-reading rules are mostly below the fixed-dosing frontier.
3. An unbuildable "informed" controller (knows each patient's parameters) is far above fixed dosing under the RANO-like rule (+28.8 to +47.9 d) but only +0.3 to +8 d under the baseline-referenced rule. Test M: most of the RANO-like gain is not a better true disease course; it is scan timing / the nadir reference. Test G: a finer search does not change this.
4. Test H: in the mixture families, knowing the responder class (a threshold on kill strength m) carries about half of the perfect-information value; in single families nothing is worth much.
5. Not done: toxicity/drop-out/death, particle-filter (belief) rules on the anchored worlds, measurement error in the endpoint, any patient data, a replicate of the mixture families, a full prior-art search (only snippet-level), sealed test (needs the user's decision; do NOT run one).

Honest status: Track C is currently a conditional MAP plus mostly negative findings. It is not yet a unique method.

## 3. What you must do to make it unique (your main job)
A. FULL PRIOR-ART SEARCH (not snippets). Open and read the papers. Cover: adaptive therapy (Zhang/Gatenby, Viossat and Noble, Sargun 2026, Corigliano 2024, Segura-Collar 2022, Browning 2026), model-predictive/optimal-control dosing in GBM and other cancers, RANO/Macdonald response modelling, TMZ schedules (RTOG 0525, GEINO 14-01, Stupp), MGMT-stratified outcomes, digital-twin and reaction-diffusion GBM papers, value-of-information in treatment design. Write `paper/trackc_prior_art_full_2026-11-XX.md` with a table: method, what it assumes, whether it handles resistance/drug cost/scan noise/tail risk, source URL actually opened. Mark anything not opened as unverified.
B. BUILD THE BASELINE SET: implement or reuse every existing dosing/control method as an arm in the same simulator, same patients, same score (restricted mean time to progression at 365/730 d, courses used, paired worst-decile tail). Already present: standard 6 courses, front-loaded, paced, raw/nadir/response-guided rules, informed oracle. Add: Zhang-style adaptive thresholds, Gatenby dose-skipping, MPC with state estimation, particle-filter belief controller, and any method from step A.
C. INVENT THE UNIQUE METHOD. Candidates (pick by evidence, try several, report all): (1) marker-conditioned dosing: a pre-treatment marker (MGMT-like) with a given accuracy picks a fixed schedule per class; find how accurate the marker must be to capture the responder-class value from test H, and compare with every baseline; (2) a hybrid controller that contains the existing methods as special cases and learns which to use per patient from early readings, with a cost for drug use; (3) a belief-state controller that is robust to the endpoint definition (RANO-like vs baseline-referenced) so it does not exploit the scan grid; (4) anything the prior-art search shows is missing. It must be buildable (no true parameters) and must win under BOTH progression definitions or say exactly where it loses.
D. EVALUATE honestly: development seeds first (new block inside the registry's reserved range); paired comparisons; intervals not just sign counts; the worst-decile tail; sensitivity to the assumed ranges. Include the negative results. Only after the method is frozen and the user agrees, a single sealed run (user decision; amend the 2026-10-31 stop rule; seeds reserved 20262300 to 20262499).
E. Write-up: update ledger rows, `START_HERE.md`, the overview, and the research-plan form `paper/research_plan_trackC_v2_2026-11-04.md`. State novelty claims with how they were checked and what is unverified.

## 4. Rules you must follow (from the repo)
- Every number must cite a file in `output/` or `vault/PAPER_FINDINGS_LEDGER.md`; put the weakest point next to each result; label simulation / assumed / post hoc / exploratory / interim.
- No clinical-benefit, patient-specific, or digital-twin claims. Keep Tracks A, B, C separate (no patient-level join exists).
- H-2/H-3 test sets and sealed seeds of scripts 167 to 173 are spent. Use `src/trackc_seed_registry.py`: launch run, wait for "manifest written ... git clean True", then register seeds, run `python src/trackc_seed_registry.py --write`, run registry tests in a scratch copy, commit, push. Code tree (`src`, `tests`, `output/trackc_seed_registry.json`) must be clean at launch. Dev seed block 20261800 to 20262099 (used up to 20261879).
- Derive paths from `__file__`. Run pytest only in a scratch copy, never in the working tree. Never commit `*_smoke` folders. Do not stage Track A files (`output/patient_split/*`) or `paper/Track A Research Plan Science Fair MD Final.md`.
- Git: check `git status --short`, `git add` named files only (never `git add -A`), commit message `Script N: ...` or `docs(track-c): ...`, push immediately. Commit and push often; the user switches between Claude and Codex when a limit hits.
- Check every ChatGPT/reviewer claim against the files before acting; log replies in `paper/chatgpt_review_response.md`.
- No dataset download or outside contact without the user's OK. CPU jobs run on this PC (8 logical cores).

## 5. How to use sub-agents / parallel workers (maximum)
Run these in parallel, each reading files on disk and not editing the same files: (1) literature readers x3 (one per topic in step A), (2) baseline implementers (step B), (3) method inventors x3 (one per candidate in step C), (4) statistician (power, intervals, multiplicity), (5) skeptic and verifier on every result before it is written (they must open the files and recompute numbers), (6) code reviewer for new code and tests, (7) writer. You merge, run the code, and decide. Cheap models are fine for reading, checking and formatting; use the strongest model for method design and for the final skeptic pass. Re-run a verifier after every doc change.

## 6. What the USER must do (tell them clearly when each is needed)
1. Keep the PC awake and on power during long runs.
2. Give ChatGPT these files for an independent rating and paste the reply back: `paper/trackc_overview_for_review_2026-11-04.md`, `paper/trackc_calibration_targets_2026-11-03.md`, `paper/trackc_anchored_spec_2026-11-04.md` (and the new prior-art file once written).
3. Decide on any sealed run (not before the method is frozen and reviewed).
4. Decide whether to request cycle-level patient data (NCTN/dbGaP route; needs an adult or school sponsor; takes months).
5. OK before any dataset download.

## 7. Output style for the user
The user is tired: report in short plain sentences (what you did, did it work, what to do now). Offer at most 2 options with a recommendation when a decision is needed.
