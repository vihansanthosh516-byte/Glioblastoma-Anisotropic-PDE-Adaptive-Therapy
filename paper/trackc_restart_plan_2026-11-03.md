# Track C restart plan (2026-11-03)

**Status (updated 2026-11-04):** approved by the user on 2026-11-03. Stage 0 (clean record) is done and committed. Stage 1 (targets table and anchored spec) is drafted and waits for the outside review; its reply must be logged in `paper/chatgpt_review_response.md`. Stage 2 is in progress on development seeds: script 190 (test A) is running, 191 (test B) is done, 192 (test C, calibration) and 193 (test D, negative controls) are written and 192 is running, 194 (test E, regime map, simple rules only) is written. Scripts 195 to 198 are claimed for tests F to I. Stages 3 to 5 not started. No sealed run is authorized. This file is the durable copy of the approved plan; `START_HERE.md` points here.

Deviations from the plan as approved (recorded 2026-11-04): (1) the parallel harness of script 177 was not moved to `src/` because its execution manifest pins its path and hash; new scripts use their own process pools and the run manifest of `src/trackc_manifest.py`. (2) Test C uses an interpolated Kaplan-Meier median because a scan every 84 days puts the plain median on the grid 84, 168, 252 days. (3) Test E first uses simple rules only; the particle-filter rules are too slow for a first pass over many worlds.

Everything below is simulation work unless it says otherwise. No clinical claim, no patient-specific effect, no digital-twin claim.

## Context (in one minute)
- Goal: a Track C result that survives a hostile outside review. I do not promise a positive result. I will not tune toward one.
- Finding: old Track C is a toy world. It does not match real GBM. Its tail test was biased. The mean gain is real inside the model.
- Plan: 5 stages. Stages 0 to 2 are docs, code and dev-only runs on this PC. Stage 3 (one sealed run) needs the user's OK later.
- Best honest outcome: a map of when feedback dosing helps, checked against real trial numbers, plus the one experiment that would settle it.

## What 14 Haiku 5.5 experts found (read-only; the marked items were re-checked from files)
Words: TTP = time to progression in the model. CVaR10 = mean of the worst 10% of cases. Open loop = the rule ignores readings.

Good
- The code audit found no bug that breaks the math. Budget, noise and score code look right.
- Headline numbers match the files: 167, 168, 169, 170, 171, 173, 177, 178.
- Mean gain is real in the model: +24.7 d (95% CI 19.6 to 29.8) vs the best simple rule. `output/trackc_sealed_containment/sealed_results.json`
- Math: with no resistance cost and no switching, the resistant share depends only on total drug-days, not on timing (Lemma 1). Not yet proved.

Bad (each one is a reviewer attack)
1. Toy scale. Model mean time-to-progression is 596 to 702 d. Real median progression-free survival is about 7 months (about 210 d; to verify). Scans: model every 28 d, clinic about every 3 months. Courses: model 12, standard 6. `src/trackc_sim.py` lines 31-36 (checked).
2. The score ignores drug use. Using all courses early costs nothing. Less drug at equal control is the point of adaptive therapy. It is not scored.
3. Nothing measures the switching rate (mu) or the reversion rate (nu). Two lab studies disagree. No open patient dataset passes the gate. Closest: Ho 2020 (141 patients, TMZ length, white-cell nadir, survival). It is weak and confounded.
4. The tail test is biased against feedback rules. Worst-decile score (CVaR10) of wide containment rises with more repeats: 171.1, 174.5, 176.0 at 1, 2, 4 repeats. Front-loaded stays 182.2 (open loop). The tail cost reported in 171 and 173 (1.6 to 4.2 d) is the same size as this bias. Re-run and checked by `src/trackc_doc_check.py`. `output/trackc_sealed_containment/sealed_ttp_by_scenario_replicate.csv`
5. The paired view is worse: wide containment loses to front-loaded in 21.0% of scenarios. Mean paired loss in the worst 10% of those: -42.1 d. Also checked by `src/trackc_doc_check.py`.
6. Even a perfect policy gains only about 3.4 d on the tail (CI 1.3 to 6.9). It gains 83.5 d on the mean. Both are lower bounds, because the oracle's schedules are achievable, not proven optimal. `output/trackc_headroom/headroom_results.json`
7. Uneven comparison. Filter rules get the true noise level and a prior that matches the generator. Simple rules get neither. Nominal tests are the best case for filters.
8. Record problems. Seed 20261121 tuned four sealed tests (167, 169, 171, 173). Script 173 re-seals the same arm as 171. A sealed H1 also failed in 169. Script 170 was written after 169. The final-status file omitted 169 and 170. In 168 the narrow arm harms with weak drugs (-59.3, -81.7, -19.8 d). Solver mode changed from "grid" to "refine" at 167. CIs for 171 and 173 cannot be rebuilt (no per-scenario files).
9. Prior art. Closest: Viossat and Noble 2021, Sargun 2026, Corigliano 2024, Segura-Collar 2022, Browning 2026. A worst-decile score in adaptive therapy was not found (weak search, snippets only).
10. MGMT. Old test: 130 patients, hazard ratio 1.11 (0.75 to 1.66). Too small to see an effect. `output/mgmt_resistance/results.json`. Use only published MGMT-stratified numbers as targets.

My errors in the brief (so a reader can discount them)
- I wrote that forced holidays improve the worst decile. The files say the opposite: higher mean in 8 of 9 cells, lower worst-decile in 7 of 9. `paper/trackc_final_status_2026-10-31.md` section 2.
- I wrote "about 2 years". The cap is 1,456 d (about 4 years).
- Web facts from agents are snippets. Treated as "to verify" until Stage 1.

## Target result (one sentence)
When does feedback dosing beat fixed dosing, once the model matches real GBM numbers and drug use counts? And which unmeasured quantity decides it?

## Stages
| Stage | Work | Gate to leave |
|---|---|---|
| 0 Clean record | Fix docs and labels. Seed registry. Fix empty tests. Label the oracle "lower bound". | Every doc number re-checked by script. Tests green in a scratch copy. Commit and push. |
| 1 Anchored spec | Verify trial targets from opened papers. Write the spec: clock, cadence, 6 courses, RANO-like progression, endpoints, margins. No simulation. | No "UNVERIFIED" left. The user sends it to ChatGPT. The reply is logged in `paper/chatgpt_review_response.md`. |
| 2 Dev-only runs | Tests A to I below. New dev seed range. Runs on this PC. | Each result read back from disk, with the 4-question report. |
| 3 Gate, then one sealed run | Frozen policy. Fresh seeds. At least 8 repeats. Score once. | Rules at Stage 3 below. |
| 4 Real-data side branch | Check Ho 2020 terms. MGMT as published targets only. | Cheap. Runs in parallel. |
| 5 Write-up | Ledger, `START_HERE.md`, Track C research-plan draft in the A and B format. | Negatives included. |

### Stage 0 details
- Edit: `START_HERE.md` (line 131 "H1 failed twice"), `paper/trackc_final_status_2026-10-31.md`, `paper/trackc_mechanism_review.md` (line 31 "cut"), `vault/PAPER_FINDINGS_LEDGER.md` (line 35), `paper/trackc_prior_art.md`, `src/trackc_oracle.py` (header said "upper"), `tests/test_trackc_sim.py` (always-true asserts near lines 120, 173, 326).
- New: `src/trackc_seed_registry.py` writes the seed registry from the `frozen_protocol.json` files. Move the parallel harness from `output/trackc_dev_plasticity_sensitivity/` to `src/`.
- New protocol hash covers solver constants (`EVENT_MODE`, `DT`, `RANGES`, `PROG_FACTOR`) and checks a clean git tree.
- Delete or relabel claims the referee expert flagged: "First positive result" (`Claude outputs/trackC_results_report_claude.md`), "170 passed" as a sealed result, 168 "no harm" without the narrow-arm caveat, mixture H1 "0 [0, 0]" as a pass, and "robust across tested departures". Add the 169 H1 failure, the 170 note, and the holiday direction. Old results stay as they were; only the wording changes.

### Stage 2 tests (dev seeds only; tune nothing on sealed seeds)
| Test | Question | Kill rule |
|---|---|---|
| A Tail estimator | How does CVaR10 move with 1 to 16 repeats? Which estimator is unbiased? | None. Cheap. Do first. |
| B Lemmas | Prove Lemmas 1 to 4. Add numeric tests. Haiku proof-checker. | Wrong algebra: drop the lemma. |
| C Clock | Tune only the standard-of-care arm (6 courses, every 28 d) to published PFS. Do not look at policy contrasts. | Cannot match inside the anchored ranges: report "model class cannot reproduce GBM". |
| D Negative controls | Does the calibrated model show no benefit for dose-dense TMZ and for 12 vs 6 cycles, as real trials did? | Big model benefit outside published CIs for every parameter set that passed C: model lacks plasticity. Say so. |
| E Regime map | Sweep mu, nu, cost, kill, resistant start. Policy contrasts: restricted mean TTP at 365 and 730 d, courses used, paired-regret tail. | Sign flips: report the flip boundary. That is the result. |
| F Cost-aware frontier | TTP vs courses vs tail. Raw-rule grid below 0.6. Give raw rules the noise level; give filters a wrong prior. Add theory-driven rules. | A rule that wins only under matched information: say so. |
| G Oracle | Finer dynamic-programming grid. Count cases where a policy beats the oracle. | Any violation: oracle is not a bound. |
| H Value of information | How many washout durations and repeats place mu and nu in the right regime with 90% chance? | None. Output: the experiment design. |
| I Numerics | dt 0.25 vs 0.1. Batch size 1 vs 150. Both event modes. Scratch copy only. | Verdict flips with numerics: it was numerical. |

### Stage 3 rules
- Gate: C passes. D passes or is explained. One frozen policy chosen on dev seeds. Margin and power set before seeds. Outside review done.
- The user must also amend the 2026-10-31 stop rule (item 4 asks for a measured washout curve). Proposed text: "published trial constraints plus outside review, results conditional on the declared mu-nu region".
- Run: fresh seed range (for example 20262300 and up), new noise keys, at least 8 repeats. Co-primary: restricted mean TTP at 365 d and 730 d, paired-regret CVaR10, courses used. Save per-scenario arrays. Score once. About 17 h on this PC (agent estimate from the 171 run time, 6,295 s). Keep the PC awake and on power.
- Pass or fail: report unchanged. If negative, an improvement expert goes first. No retuning.

## Experts for later waves (all Haiku 5.5)
Reference verifier (two), proof checker, statistician (power and pre-registration), skeptic and verifier at every gate, code reviewer for new code, neuro-oncology realism reviewer for the spec, improvement expert after any negative, visualizer, writer. The main session writes and runs the code.

## Rules
- Spent forever: H-2 and H-3 sets. Sealed seeds of 167 to 173. Old margins and endpoints stay as they were. They are not changed.
- Every number cites a file in `output/` or `vault/PAPER_FINDINGS_LEDGER.md`. Weakest point beside each result. Say simulation, assumed, post hoc.
- No clinical claim. No patient-specific effect. No digital-twin claim.
- No dataset download and no contact without the user's OK. CPU jobs run on this PC.
- Commit after each verified step. Push at once. Format "Script N: ..." for scripts, "docs(track-c): ..." for docs. Update `START_HERE.md` when the decision changes.
- Interim numbers are labelled interim.

## Checks that prove it works
- Stage 0: run a script that compares each cited doc number with its JSON (`src/trackc_doc_check.py`). Run `pytest tests/test_trackc_sim.py tests/test_trackc_oracle.py -q` in a scratch copy. Check `git log origin/main -1`.
- Stage 1: each target has the URL of a page that was opened.
- Stage 2: each output folder has a manifest (seed range from the registry, code hash, git commit, clean tree). Read each file back: exists, size over 0, values sane, same budget and same seeds in every arm. One independent Haiku verifier re-runs one test.
- Stage 3: protocol committed before scoring. Seal marker present. One scoring pass.

## Decisions for the user (defaults apply if only approved)
1. Sealed run: decide at the Stage 3 gate. Default: ask then. Other option: never run another sealed test; finish with dev-only work and the lemmas.
2. Patient data: default open-download only. Other option: also start an NCTN/dbGaP request now (only route with cycle-level data; takes months; needs an adult or school sponsor).
