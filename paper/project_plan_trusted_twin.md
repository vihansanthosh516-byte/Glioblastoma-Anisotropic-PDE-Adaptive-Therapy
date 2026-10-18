# Project plan: a glioblastoma twin that knows when to trust itself (DRAFT, 2026-10-18)

Status: plan only. Nothing here has been run. Not yet declared as an amendment. Written for the user's goal: keep the GBM project and its parts, change a lot, add a new mathematical core, for GARSEF (February 2027).

## 1. Honest novelty check (searches, snippet level; unverified)
- Uncertainty-aware glioma digital twins already exist: Bayesian calibration with probabilistic forecasts and risk-aware radiotherapy (https://arxiv.org/pdf/2308.12429; https://arxiv.org/pdf/2505.08927). Strong groups work here. A plain "digital twin" is NOT new.
- NASEM 2023 stresses VVUQ and "fit-for-purpose" validation for digital twins (secondary sources only; exact wording not read). A cancer twin is a named example. No oncology VVUQ benchmark was found (https://www.frontiersin.org/journals/physiology/articles/10.3389/fphys.2024.1473125/full).
- NOT found in about 9 searches: (a) a per-patient rule that withholds a growth forecast when it cannot beat "no change" given measurement noise; (b) a forecast skill score against persistence for glioma growth models; (c) a closed-form forecast horizon from noise and growth speed; (d) conformal coverage guarantees for glioblastoma radiotherapy margins.
- Closest work: a rat-simulation study of imaging quality (https://pmc.ncbi.nlm.nih.gov/articles/PMC12564088/, no theorem, no baseline); Bayesian glioma growth that folds segmentation error into model error (https://pmc.ncbi.nlm.nih.gov/articles/PMC10599765).
- Conclusion: the new part is narrow: a trust criterion with a theorem behind it, tested against a no-change baseline on real patients. The twin around it is not new. Claim only that.

## 2. What stays, changes, goes
| Item | Decision | Reason |
|---|---|---|
| Recurrence model and two-centre all-recurrence gain (+0.081, +0.085) | Keep as module 1 ("where") | Verified in files; not new alone |
| Exact monotone solver (Selling) | Keep as the verified simulator | Needed to check the theorems |
| Measured segmentation noise floor | Keep, extend to all longitudinal sets | Core input to the theory |
| 4D viewer | Keep, add a trust indicator and a confidence band | Demo for the poster |
| Persistence tests (H1, H3, volume forecast) | Keep as evidence for the theory | Real data, pre-declared |
| Drug what-if | Keep as labelled simulation only | Assumed rates |
| Triage, H-4 test-set hunt, enhancing-recurrence chase | Drop | No headroom, no data |
| Headline "digital twin predicts growth" | Drop | Not supported |

## 3. The project in one sentence
A glioblastoma twin with three honest parts: it shows WHERE the tumour returns (validated), says WHEN a growth forecast is meaningful at all (new theory), and shows WHAT-IF drug scenarios labelled as simulation.

## 4. Mathematics to prove (the new core)
T1. Exact identity: Dice(A,B) = 1 - |A sym-diff B| / (|A|+|B|). First-order law: 1 - Dice ~ (mean |u|) S / (2V) for a boundary displaced by a normal field u (S surface area, V volume). State the error term and the smoothness assumptions.
T2. Noise ceiling: independent boundary noise of standard deviation sigma on both scans bounds the expected Dice of every forecast, a closed form with folded-normal moments.
T3. Skill score SS = (Dice_model - Dice_persistence) / (1 - Dice_persistence). Closed form for the best attainable SS (SS_max) as a function of displacement over noise.
T4. Fisher-KPP front speed c = 2 sqrt(D rho) (known) turns displacement into days: forecast horizon T*(sigma, D, rho) at which SS_max reaches a chosen level. Per-patient Forecastability Index FI = expected displacement / noise.
T5. Error budget: forecast error = noise + speed error + direction error. State when persistence provably wins.
Each theorem gets: statement, assumptions, proof, numerical check.

## 5. Simulations (verify the theory before touching patients)
- S1: random smooth tumours (spheres, ellipsoids, brain-shaped domains). Displace boundaries by known fields, add noise, compare measured Dice with T1-T3. Pass rule fixed first: relative error under 5% in the regime where the assumptions hold; show where they fail.
- S2: Fisher-KPP with the exact monotone solver on a brain-shaped domain; check that front displacement matches c Delta t and that SS_max matches T3/T4.

## 6. Real-data tests (each declared in an amendment BEFORE it is run)
- R1 (parameter-free prediction): the theory predicts, from measured noise and measured displacement alone, the share of real scan pairs where any forecast can beat persistence. Compare with the observed share of pairs where the FK forecast beats persistence (data: MU 334 pairs and the longitudinal sets already on disk). No fitted parameter. Success rule to be written in the amendment.
- R2: does FI rank pairs by forecast gain (AUC with bootstrap CI)? Leave-one-dataset-out.
- R3: guaranteed-coverage plan (conformal) for the recurrence map: smallest volume that reaches a stated coverage with a finite-sample guarantee; test coverage leave-one-centre-out; report volume against the 15 mm plan.
All data are spent sets: label results "exploratory, leave-one-dataset-out", never "locked".

## 7. Deliverables
Paper (theorems with proofs, simulations, data), poster (one clear question, three figures), 4D viewer with trust indicator, interview notes (limits answered fully when asked).

## 8. Timeline to February (fair date and forms UNVERIFIED: check the GARSEF handbook)
- Weeks 1-2: theorems T1-T4 drafted, S1 simulations, amendment 15 declaring S1/S2 pass rules.
- Weeks 3-4: S2 with the exact solver; measure noise and displacement on all longitudinal pairs.
- Weeks 5-6: declare and run R1, R2 once each.
- Weeks 7-8: declare and run R3.
- Weeks 9-10: viewer trust indicator and confidence band.
- Weeks 11-13: paper, poster, mock interviews.
- Weeks 14-15: buffer, fair forms, fixes.

## 9. Risks
- The theory may already exist in forecasting or image-analysis literature (not found, not proven absent). Mitigation: deeper search in week 1 before proving anything; ask a faculty member.
- R1 could fail (theory and data disagree). Then report the disagreement and what it says about the model; the theory and simulations still stand.
- Proofs must be correct. Check each numerically before claiming.
- Novelty is narrow. A win is not predictable.
