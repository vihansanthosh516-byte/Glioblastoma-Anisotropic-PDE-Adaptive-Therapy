# Response to the outside review (ChatGPT)

Written 2026-10-04. The reviewer saw only `overview_for_review.md`. It could not open any file. Below are the points I think are wrong or overstated, with the evidence. Points I agree with are listed at the end.

## 1. WRONG: "The model also loses on LUMIERE"
The review says the forecast "also loses on LUMIERE" and builds advice on it ("did not maintain the advantage on independent data").
- What LUMIERE tested (script 91, `output/lumiere_volume_forecast.json`) was a **growth-only volume guess**: extend the last growth rate to the next scan. It was compared with "no change" on tumour volume.
- It was **not** the PDE model, **not** Dice, and used **no** masks (masks are in a 30 GB zip that was not downloaded).
- The PDE forecast (Dice 0.259 vs 0.233, `output/forecast_labels_core/results.json`) has **never been tested on outside data**.
- Correct statement: "The PDE forecast has no independent test yet. A simpler growth-only volume guess also loses to no-change on LUMIERE, which matches the MU-Glioma finding that growth-only forecasts fail on shrinking tumours."
- Do not write "failed to generalise". Nothing was tested for generalisation.

## 2. OVERSTATED: "Most publishable: the MRI forecast ... an explicit prediction task ... suitable for testing"
The review lists strengths and leaves out the weaknesses that are in the vault:
- The target (cellular core, labels 1 and 3) was chosen **after** a volumes-only probe. It is disclosed in the script header. The earlier `mask>0` target was negative. (`vault/PAPER_FINDINGS_LEDGER.md` B2, L17)
- The median gain is +0.005. The mean gain is +0.026 (CI 0.008-0.044). The model is better in only 55% of patients.
- Absolute Dice is low (0.26).
- 19 patients were skipped (empty core).
- It has no outside test.
It is still the best Track B result. It is a moderate result, not a strong one.

## 3. ALREADY DONE: "Choose the threshold using training data only; test it on held-out patients"
Script 87 (`output/limit_checks.json`) did this for the 0.02/day cutoff:
- Leave-one-out: the cutoff is chosen without patient i, then used on patient i. 60 of 61 correct (98.4%).
- Permutation test: a perfect split by chance has p = 5e-5 (20,000 shuffles).
- Exact CIs: 10/10 above the cutoff (CI 0.69-1.00), 0/51 below (CI 0-0.070).
What is still true: n = 10, the cutoff was first picked after seeing the data, and no other cohort has enough fast-growing patients (LUMIERE has only 3-4 above 0.02). The result stays exploratory.

## 4. INCOMPLETE: "CGGA gives HR 1.009, p = .08, so age should be cut"
That is the pooled all-GBM row (n = 374). It includes IDH-mutant secondary tumours, which behave differently. In the cohort closest to TCGA (IDH-wild-type, primary, pooled n = 179) the age HR is 1.018 (1.005-1.032), Holm p = 0.016 (`output/cgga_validation.json`). So age replicates in direction, with a weaker size. Cutting it from the headline is still right, because age is a known factor. But it did not fail to replicate.

## 5. UNVERIFIED: "If patient-level C-GAT ends up around 68-74%, that is believable"
The reviewer took this from my interim numbers (2 of 5 folds). The final result is not in yet (script 88). Do not use this range. The first two folds already vary from 58% to 78% for the clean variants.

## 6. MISSING CONTEXT: "The simple paced rule beating the standard schedule is more interesting if it survives validation"
It leaves out what the data say:
- The rule is evaluated only in simulation, with **assumed** kill rates.
- It extends time to progression (+8 to +9 days on 91 MU patients, +16 to +19 days on LUMIERE growth rates, all CIs above 0, random kill rates, script 93) but it **loses** on day-90 tumour volume.
- There is nothing to "validate" against a real patient outcome. No real treatment-response data are in the project.
Say "simulation evidence" every time.

## 7. MINOR
- The review rates against "ISEF" and "high school". That came from my overview, which assumed a level. It is the author's call, not a fact.
- "Honesty 9/10" and the ratings are opinions from a text summary. The reviewer could not check any number.
- "97 scripts producing a huge number of claims" is fair as style advice. The paper plan already limits the main text to the top 3 findings per track (`PAPER_FINDINGS_LEDGER.md`).

## Where the review is right (keep)
- The title overpromises. Use the conservative title.
- Track A does not feed Track B, so "multi-scale" is not earned.
- Track A's 78.7% must not be a headline. Report the drop under patient-level splitting.
- RL adds nothing over a one-line rule. Say so.
- Cut or downgrade: age to survival, simulated elongation, DTI, the dual-drug rescue.
- Make Track B the centre. Make Track A exploratory and Track C a simulation study.
- Every main finding needs: baseline, effect size, uncertainty, limit.


---

# Review of analysis_plan_v1 (2026-10-04, pre-commit checklist from ChatGPT)

Each item was checked against `analysis_plan_v1.md` and the code. The "9/10" rating was not used. All fixes are in Amendment 1 of the plan.

| # | Review claim | Check | Verdict |
|---|---|---|---|
| 1 | Patient-level aggregation undefined | Plan said "patient is the unit" but not how several forecasts become one value | Right. A1.1 |
| 2 | Eligibility undefined | No criteria in v1 | Right. A1.2 (from `forecast_pairs()` and script 81) |
| 3 | Missing / failed-fit rules | Absent | Right. A1.3 |
| 4 | Baselines not defined | Only named. Also: script 81's `uniform_dilation` uses the true scan-2 volume (oracle) | Right, and worse than stated. A1.4 |
| 5 | PDE spec not frozen | Absent; no optimizer exists (grid search) so "optimizer / regularization" are "none" | Right. A1.5, written from `run_improved_aniso.py` |
| 6 | H2 "same sign" too weak | v1 said "same sign, CI reported" | Partly right (CI was already reported). Tiers added. A1.7 |
| 7 | H3 not stated frozen | v1 had a general freeze rule | Partly right. A1.8. Also LUMIERE may have no tensors, so H3/H4 may be untestable externally |
| 8 | H4 should be exploratory | Result already seen (-0.001, p 0.97) | Right. A1.9 |
| 9 | H8 needs patient linkage | v1 wording implied a patient-level gain | Right. A1.10 |
| 10 | Utility function undefined | Absent | Right. A1.11. Weights and progression rule are our choices, flagged as such |
| 11 | H10 should be neutral | v1: "RL beats it" | Right. A1.12 |
| 12 | Controller information | Absent | Right. A1.12 |
| 13 | Subgroups mathematical | v1 had cutoffs +/-10% but "previous two scans" | Right, and a real flaw: first forecasts have no previous scan, so they cannot be subgrouped. A1.13 |
| 14 | Early post-RT is not RANO | v1 said "(RANO 2.0)" | Right. A1.14 |
| 15 | Catastrophic threshold wording | v1 had threshold, no label | Right. A1.15 |
| 16 | H1 needs mean AND median | v1 already required both | Already in the plan. Wording made explicit. A1.6 |
| 17 | Multiple testing | v1 already had Holm / BH, primary separate | Already in the plan. A1.16 restates |
| 18 | Freeze list | v1 had tag only | Partly right. A1.17 |
| 19 | External reruns | v1 said "re-runs only for crashes, logged" | Already in the plan. A1.18 restates |
| 20 | Solver verification gate | v1 said no biological result until fixed | Already in the plan. A1.19 restates |

**What the outside review missed (found in the files, not by the reviewer)**
1. Four of five CV folds in script 81 chose the grid edge (rho 0.1, d 0.01). The fit may lie outside the grid. A1.5 adds a development-only grid-extension check before the freeze.
2. The isotropic homogeneous model (+0.028) is as good as the anisotropic one (+0.026). The gain comes from growth, not from direction (`output/forecast_labels_core/results.json`).
3. Script 81 forecasts only scan 1 to 2. Rolling-origin needs the count of patients with 3 or more scans, which is not known yet. A1.1 sets a rule.
4. The "grew / shrank-or-same" split uses the target scan, so it is outcome-defined.

## Round 2 (2026-10-05): review of `paper/review_packet_phase0-3.md`
Each claim checked against files. No agreement by deference.

| # | ChatGPT claim | Checked against | Verdict | Action |
|---|---|---|---|---|
| 1 | PDE must beat persistence AND the geometric rule; PDE vs geometric is the key secondary | `analysis_plan_v1.md` A3.2 | Right, and already registered (A3.2). | None. Report both, paired per patient. |
| 2 | Geometric rule needs a leakage audit | `src/98_baseline_ladder.py` lines 133-156 | Audited. Inputs: input mask, input-scan brain mask, input volume, interval dt, and g = median log growth rate of pairs from patients NOT in the test fold (`prim["fold"] != f`). No future mask, no future volume, no per-patient fit. dt uses the date of the target scan, as the PDE does. The oracle uses the true volume and is labelled. No leak found. | None. |
| 3 | The geometric rule is "potentially your biggest discovery" | `baseline_ladder.json` | Overstated. Gain is +0.0136 on a persistence Dice of 0.277. It is a population mean growth rate applied to the mask. It shows that per-patient extrapolation (last_rate -0.0445, linear -0.028) is worse than a population rate. It is a useful baseline, not a discovery. | Keep wording modest. |
| 4 | Cell selection must use training patients only; state it in the plan | `src/100_pde_manifest.py` `select_cells()` lines 244-259 | Right and already true: the cell for fold f is chosen on first pairs of patients with `fold != f`. The plan did not say this in one sentence. | Add the sentence in Amendment 4. |
| 5 | Do not call training-patient results "held-out" | `analyze()` | Right. Reported Dice are out-of-fold only. Terms to use: training, development (MU folds), external test (LUMIERE). | Wording rule in Amendment 4. |
| 6 | 2 mm vs 1 mm sensitivity; connect numerical error to Dice | W22, `PILOT_N` | Right and queued. But the pilot takes the first 20 fold-0 patients by ID, not a spread of cases. That is a real gap. | Change the pilot to a pre-declared stratified subset (by input volume and interval, both known before forecast), after the main run ends. Not done yet. |
| 7 | Segmentation perturbation is more valuable than another model | Masterplan 3.3 | Agree. Planned. | Phase 3.3. |
| 8 | Do not call 57% of pairs "unforecastable" | ledger 2s | We do not. The ledger says "below the assumed noise floor". | None. |
| 9 | DTI wording: no evidence that the atlas direction explains centroid displacement | ledger 2t | Agree. Ledger already says it does not prove tracts are irrelevant. | None. |
| 10 | Cohort claim: say "eligible longitudinal cohort" | A2.2, ledger 2p | Right and already in A2.2. | None. |
| 11 | Make GBM-only the primary analysis | A1.2, A2.1 | Fair point. Registered primary is all 134 eligible; GBM-only (n=104) is a sensitivity run. No PDE result exists yet, so a change now would be clean. Baseline GBM-only values (geometric +0.0149 vs +0.0136) are already seen and are close. | User decision (2 options). |
| 12 | LUMIERE masks are machine-made; do not call it "confirmatory" | A1.18, plan lines 10 and 32 | Right on the masks (already stated). Plan text does say "confirmatory". | Replace with "prespecified external evaluation" in Amendment 4; add external segmentation shift as a failure mode. |
| 13 | Treatment assumptions must not contaminate the MRI evidence | A1.11 | Agree; treatment is simulation only. | Keep separate in the paper. |
| 14 | 93 tests are engineering evidence, not validation | ledger 2u | Agree. We say so. | None. |
| 15 | Priority order: PDE rerun, 1 mm, segmentation, GBM-only, negative controls, freeze, LUMIERE | masterplan Phase 3 | Same order as the plan. | None. |

What ChatGPT did not say: the packet's own correction (+0.031 was wrong, +0.0136 is right) came from checking the ledger against the JSON. Other ledger entries have not all been re-checked against their JSON yet (open task).

