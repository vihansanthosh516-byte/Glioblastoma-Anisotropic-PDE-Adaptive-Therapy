# Track A review follow-up (2026-11-10)

Checks of the 2026-11-10 read-only review findings against saved files. A5 code, result, and its failed 5-of-6 rule are unchanged. No expression values were read.

## 1. AURC scale (finding A-4, verified)

`class_balanced_aurc` in `output/trackA_yu2020/result.json` integrates class-balanced selective risk over coverage from ceil(N/2)/N to 1 with the trapezoid rule and no normalization (`src/189_trackA_yu2020_lodo.py:486-493`; PROTOCOL section D). The coverage span is about 0.5, so values are **not** on the usual 0-1 AURC scale. A classifier with constant risk 0.5 (chance for class-balanced risk) scores about 0.25. Lower is better. Every A5 AURC number (for example candidate 0.1479, logistic 0.1764) must be read on this scale.

## 2. Within-cell-type check (finding A-1, verified; previously unreported)

PROTOCOL sections 8-9 require it. Only the Myeloid stratum was eligible (`result.json` `cell_type_stratum_eligibility`; donors GS1, GS2, GS3, GS5, GS13; GS6 lacked periphery cells). Per-donor Myeloid AURC (same scale as above):

| Donor | Candidate | Logistic | Gradient boosting |
|---|---|---|---|
| GS1 | 0.2498 | 0.1432 | 0.1977 |
| GS2 | 0.2507 | 0.3270 | 0.2855 |
| GS3 | 0.2465 | 0.1914 | 0.2166 |
| GS5 | 0.2173 | 0.0332 | 0.0204 |
| GS13 | 0.2494 | 0.0403 | 0.0189 |

Within Myeloid cells, the candidate scored near the chance level of about 0.25 in every donor. Logistic regression was lower than the candidate in 4 of 5 donors. Reading: the candidate's mixed-cell advantage did not appear inside this cell type. This is a descriptive check with small per-donor counts (for example GS5 has 23 core and 11 periphery Myeloid cells), not a new test.

## 3. Per-donor drivers (finding A-2, corrected)

Candidate minus logistic AURC: GS1 −0.104, GS2 −0.093, GS5 −0.054, GS6 −0.005, GS3 +0.024, GS13 +0.062. Three donors (GS1, GS2, GS5) drive the lower mean; the reviewer's statement that GS2 was close is wrong. GS6 was nearly tied.

## 4. Ablation name (finding A-6, verified; no code change)

`ablation_no_gate` sets the gate's L1 weight lambda_1 to 0 but still trains the gate (`189_...:643-644`). It matches protocol ablation 1. Renaming the key would change the sealed runner hash and the saved result keys, so the name stays. Read it as "no gate sparsity penalty", not "gate removed".

## 5. GBmap version (finding A-3, checked; no code change)

The runner hard-codes `CELLxGENE dataset c888b684-6c51-431f-972a-6c963044cef0; version 861acfd8-25f0-418b-a445-aa96da232827` (`189_...:296`). The public CELLxGENE curation API (queried 2026-11-10) lists this dataset as "Core GBmap" in the pinned collection `999f2a15-3d7e-440b-96ae-2c806799c08c`, with the same current version ID. The value is consistent with the protocol pin. The remaining weakness is provenance form: the string is typed, not read from the saved snapshot. Not changed, because the runner is sealed.
