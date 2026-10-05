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
