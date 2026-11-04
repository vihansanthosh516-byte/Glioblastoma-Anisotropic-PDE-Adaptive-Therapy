# Track C anchored specification (draft for review, 2026-11-04)

Status: Stage 1 of `paper/trackc_restart_plan_2026-11-03.md`. Design only. No simulation was run for this file. It is not a protocol, it freezes no hypothesis, and it authorises no sealed run. Published numbers and their verification levels are in `paper/trackc_calibration_targets_2026-11-03.md`. Everything here is a computer simulation with assumed parameters. No clinical claim, no patient-specific claim, no digital-twin claim.

## 1. The question

When does feedback dosing (a rule that reads noisy scans) beat fixed dosing, once the simulated tumour world is made to match published standard-of-care numbers and drug use counts? And which unmeasured quantity decides it?

## 2. What the Stage 2 development tests have found so far

All SIMULATION; development seeds or saved development data; exploratory. Each is a lower-level fact that changes the design below.

1. **The tail estimate was biased against feedback rules.** In the sealed 167 per-scenario file, the worst-decile score (CVaR10) of wide containment rises with the number of noise repeats (171.1, 174.5, 176.0 at 1, 2, 4 repeats) and stays at 182.2 for front-loaded dosing. Scripts 171 and 173 used 2 repeats. `output/trackc_sealed_containment/sealed_ttp_by_scenario_replicate.csv`; checked by `src/trackc_doc_check.py`. Test A (script 190, running on fresh development seeds) re-measures this with up to 16 repeats; result pending.
2. **With no crowding and no cost of resistance, front-loaded dosing is optimal** for time to progression, in the no-crowding variant of the equations (`paper/trackc_lemmas_2026-11-04.md`, Proposition; algebra independently checked by a read-only Haiku proof-checker). By the contrapositive, any gain of feedback in the full model must come from crowding (starting burden B0 close to the carrying capacity), the cost c of resistance, or switching. Brute force on short horizons agrees: the best schedule beat front-loaded dosing by 1.3 d on average at B0 0.05-0.10 and by 21.9 d at B0 0.40-0.50 (`tests/test_trackc_lemmas.py`).
3. **In the saved oracle headroom, B0 is the largest single driver.** Mean gain from known-parameter dosing over front-loaded dosing: 29 d in the lowest B0 quarter, 132 d in the highest (`output/trackc_headroom/headroom_per_scenario.csv`, 200 development scenarios, exploratory linear model R-squared 0.33). Cost c matters mainly when B0 is high (105 d vs 156 d at high B0; 31 d vs 38 d at low B0).
4. **The old world is far from the published standard arm.** The model's mean time to progression is 596 to 702 days. RTOG 0525's standard arm has a median PFS of 5.5 months (143 to 186 days at the 95% interval) from the start of adjuvant treatment, with scans about every 3 months and a median of 3 cycles given. The old model scans every 28 days, gives up to 12 courses, and calls progression at 1.2 times the baseline burden, not at the RANO nadir rule.

Weakest point of this section: items 2 and 3 are derived inside the assumed model. They say what the equations allow, not that real tumours behave this way.

## 3. Changes from the old world (what is new, what stays frozen)

Stays frozen: the equations, parameter symbols, particle filter, and all code in `src/trackc_sim.py` unless a change below needs a code addition (each addition is off by default, so old outputs reproduce). Old sealed seeds and old margins stay spent.

| # | Change | Why | Code needed |
|---|---|---|---|
| 1 | Primary scan interval 84 days (every third cycle); sensitivity 28, 56, 112 days | Both reference trials scanned about every 3 months; EANO calls 2 to 6 months common practice | New option `scen["scan_every"]`: scans visible only every k-th cycle, controller told (the 170 "told" mode, but regular). Implemented, off by default, tested in `tests/test_trackc_scan_endpoint.py` |
| 2 | Standard arm = 6 courses (label), control arm = 12 courses | Label: 6 maintenance cycles; 12-vs-6 trials show no benefit | None (budget is a parameter) |
| 3 | Two progression endpoints, both reported: (a) old rule, true burden reaches 1.2 times baseline; (b) RANO-like, at a scan the burden is at least 1.4 times the smallest earlier scan value (volume equivalent of +25% in the product of diameters) | The old rule follows the adaptive-therapy convention; RANO measures from the nadir. The two can differ by hundreds of days. | New option `simulate(..., endpoint="rano")`, off by default. A floor of 0.1 x baseline stands in for the smallest measurable lesion (assumed). Also added: `simulate(..., budget=k)` for the course cap and `scen["dose"]` for a drug-exposure multiplier. All tested in `tests/test_trackc_scan_endpoint.py` |
| 4 | Score includes drug use: courses used, and TTP-versus-courses frontier | Less drug at equal control is the point of adaptive therapy; the old score ignored it | Saved already in script 190 output |
| 5 | Tail endpoints: joint CVaR10 and paired-regret CVaR10 against the standard arm; plug-in CVaR10 of scenario means only with 8 or more repeats | Test A shows the plug-in is biased low at few repeats | `src/trackc_tail_estimators.py` (done) |
| 6 | Mean endpoint: restricted mean TTP at 365 and 730 days (not the 1,456-day cap) | Cap is about 4 years; trials report 2-year horizons | Analysis only |
| 7 | Information parity: raw rules get the true noise level; filter rules also run with a wrong prior | Filters had the true noise sd and a matching prior, raw rules had neither | Small policy option |
| 8 | Starting burden B0, cost c, switching mu and reversion nu become named regime axes, with growth rate r and kill m as calibration axes | Section 2 items 2 and 3 | None |

## 4. Worlds, calibration and the regime map

A "world" is one setting of the parameter ranges. Calibration means the standard arm (6 courses, one per 28-day cycle, scans every 84 days) must reproduce **median PFS 143 to 186 days** on the RANO-like endpoint (3b), measured from the first adjuvant course. This uses only the standard arm. Policy contrasts are not looked at while choosing worlds.

1. Draw worlds by Latin hypercube over seven axes (`src/trackc_calibration.py`, all ranges assumed): a shift of the growth-rate range r (x0.5 to x2.5), a scale of the kill-to-growth range m (x0.05 to x1; kill strength is not anchored), a shift of the initial resistant fraction f0 (x0.5 to x10, capped at 0.9), the centre of the starting burden B0 (0.12 to 0.55, half width 0.10), the cost c (0 to 0.25), and the switching and reversion rates mu and nu (levels 0, 0.01, 0.04 per day, as in script 177). Each world gets its own patients. The standard arm is open loop, so no particle filter is needed and this step is cheap. The first run (script 192) uses 500 worlds with 400 patients each.
2. Keep worlds whose standard-arm median PFS falls in 143 to 186 days. The median is the Kaplan-Meier median with linear interpolation between scan days, because a scan every 84 days puts the plain median on the grid 84, 168, 252 days. If none do inside the ranges allowed by the targets file, the finding is "this model class cannot reproduce the standard arm" and Stage 2 stops there for the old equations. The same worlds are re-run with the RANO floor and factor changed, with scans every 56 and 112 days, and with the old 1.2 x baseline rule, to show how much the answer depends on these assumptions.
3. Negative controls on the kept worlds (script 193, fresh patients): (i) dose-dense versus standard dosing should show a PFS hazard ratio inside the RTOG 0525 interval, HR 0.87 (0.75 to 1.00); the dose-dense arm is the standard arm with the drug effect multiplied by 2.0 (an assumed ratio of drug exposure; 1.5 and 3.0 are sensitivity runs); (ii) a landmark comparison at the day-168 scan, 12 versus 6 courses among patients without a progression call, should fall inside GEINO 14-01, PFS HR 0.99 (0.71 to 1.38). A world with a model hazard ratio outside these intervals is marked "fails negative control". Hazard ratios use a Cox model with Efron ties (`src/trackc_survival.py`).
4. In the controlled worlds (script 194 runs all of them, because the rules are simple; the first plan was a spanning subset of 6 to 12), run the policy contrasts: feedback rules against front-loaded dosing and against the best fixed schedule, on mean RMST, joint and paired tail, and courses used.
5. The result is the map: where in (B0, c, mu, nu) each contrast is positive, zero or negative, and the boundary. A sign flip is the finding.

## 5. Stage 2 tests (development seeds 20261800 to 20262099 and noise keys 1,800,000 to 1,899,999; nothing is tuned on sealed seeds)

| Script | Test | Status |
|---|---|---|
| 190 | A: tail estimator versus number of repeats | done 2026-11-04 (interim): the bias from few repeats is about 1 to 4 d and saturates; the tail deficit of feedback rules remains at 16 repeats |
| 191 | B: exact identities and gain-versus-burden checks | done: 13 tests pass (scratch copy); algebra independently checked; wording corrected |
| 192, 193 | C: calibrate the standard arm; D: negative controls | done for config `single` (81 of 500 worlds pass, 39 controlled). `--config single_rep` (fresh-seed replicate, adds 3- and 4-course variants) and `--config mixture` (responder/non-responder family) are defined in `src/trackc_calibration.py`; their status is in `START_HERE.md` |
| 194 | E: regime map on all controlled worlds (simple rules, both endpoints, cost-aware arms, informed family maximum) | done for `single` (interim); `mixture` pending |
| to write | F: cost-aware frontier and information-matched rules; G: oracle on a finer grid | F is partly covered by the cost-aware arms of 194; G by the share of patients who beat the informed family maximum (194) |
| to write | H: how many washout durations and repeats place mu and nu in the right region with 90% probability | after the mixture result of 194 |
| 198 | I: numerics (step 0.25 vs 0.1) | done for `single`: no changed call, no flipped flag. Batch size and the old event mode were not tested |
| 199 | J (added): how deep a response the controlled worlds predict | done for `single`: no deep responders predicted; `mixture` pending |
| 200 | post hoc analysis of C and D from the saved files | done for `single` |

Kill rules are in the restart plan. Every output folder gets a manifest (seed block from `src/trackc_seed_registry.py`, protocol hash, code hashes, git commit, clean-tree flag).

## 6. What is still open before any sealed run (Stage 3)

1. The 2026-10-31 stop rule asks for a measured washout curve. The restart plan proposes "published trial constraints plus outside review, results conditional on the declared mu-nu region". This needs the user's amendment.
2. Margins and power: set before seeds, by the statistician checker, from Test A and the regime map.
3. Outside review of this spec and the targets file, and the reply logged in `paper/chatgpt_review_response.md`.
4. A decision on the nadir-referenced endpoint: is the primary endpoint 3a (adaptive-therapy convention) or 3b (RANO-like)? Both are reported either way.

## 7. What this cannot show

- It cannot show that any dosing rule helps patients. The tumour is a two-population equation set with assumed ranges.
- B0, c, mu and nu are not measured here. Nothing in the open literature checked so far gives a rate for mu or nu (abstract-level only).
- One calibration number (a median) leaves many parameter combinations that fit equally well; the map shows how much the answer depends on that freedom.
- The volume-change-per-cycle data needed to anchor the kill strength were not found.
