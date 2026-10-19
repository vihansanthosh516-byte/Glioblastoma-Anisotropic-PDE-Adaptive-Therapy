# Probabilistic RANO audit (A–D, G scope)

## Scope and checks

Read amendments 14–21, scripts 140–145, `src/rano_calc.py`, its seven declared tests, ledger 2br–2bx, and cited RANO outputs. Did not execute scripts 140–145 or pytest. Independently recomputed A15 flip rate from its saved pair CSV and A21 conversion arithmetic from the published inputs in the amendment.

## Findings

| Finding | File:line | Evidence | Severity |
|---|---|---|---|
| A21 simulation grid is not the union of the conversion values. The plan says “constant change-SD s in the union of the audit values” and PU1/PU2 are evaluated over that list. Script 145 simulates only 7 rounded values (`0.08, 0.111, 0.126, 0.139, 0.168, 0.193, 0.262`). The saved conversion table contains 18 distinct values across both versions and all three readings (including 0.099, 0.143, 0.146, 0.160, 0.199, 0.202, 0.230, 0.235, 0.291, 0.342; 0.079 is rounded to 0.08). Thus the confidence-range endpoints for multiple readings are not simulated. PU1/PU2 pass for the declared seven-point grid, but do not establish the declared “union” coverage or each reading’s full range. | `analysis_plan_v1.md:624`; `src/145_conversion_audit.py:43`; `output/conversion_audit/conversion_table.csv:2-7` | Script grid compared with all six conversion rows; independently recalculated R1/R2/R3 using the printed RC%, RC mL, and mean volumes. | **High** |
| A15’s original headline calls the consecutive-scan proxy a “RANO 2.0 volumetric class.” Script 140 classifies `V1/V0` for each consecutive pair. That does not use the RANO progression reference (nadir) and response reference (baseline). Amendment 18 later states this erratum and the ledger 2bt corrects the framing, but ledger 2br remains overbroad when read alone. | `src/140_rano_tool_dependence.py:10`; `vault/PAPER_FINDINGS_LEDGER.md:124,126`; `analysis_plan_v1.md:550` | Compared declared implementation and data logic with A18 erratum: A15 tests consecutive ratios; A18 says progression uses nadir, response baseline. | **Medium** |
| A16/A18/A19/A21 trial simulations model only a simplified volumetric progression process. They omit clinical confirmation timing, treatment-era rules, human approval of segmentations, and non-volumetric progression criteria. The ledger identifies simplified confirmation and volumetric-only scope; results must remain labelled simulation under assumed parameters. This is a scope limit, not a code-versus-amendment deviation. | `src/141_trial_sim.py:42-60,94-115`; `analysis_plan_v1.md:509-516,578-585,624-630`; `vault/PAPER_FINDINGS_LEDGER.md:125,127,130` | Traced event rule, confirmation rule, and estimands in code; compared with declared limits and current RANO 2.0 workflow description. | **Medium** |

## Reconciliation and method checks

- Headline output figures in ledger entries 2br–2bx reconcile with their cited summaries/CSV rollups at reported precision. A15 saved pairs independently give 121 eligible rows and flip rate 0.5454545; its saved bootstrap interval is 0.4559–0.6381. A15 MDC outputs give the ledger’s bands and readings. A16 summary gives H-S1/H-S2 supported and H-S3/H-S4 unsupported; the four out-of-band type-I rates match 2bs. A18 summary gives PA1, PA2, PB1, PD1 supported and exact A16 reproduction; A19 gives PX1–PX4 supported; A21 gives PU1–PU4 supported for its seven-point grid. Calculator and viewer saved check fields agree with 2bv/2bw. No headline-number mismatch found in these outputs.
- Commit order passes for Amendments 15–21: each declaration commit precedes code commit, which precedes result commit, per `git log`. Examples: 51fd053 → ad9c55c → 945e2ab; 2d3e6a4 → 48c5933 → 131495f.
- `src/140_rano_tool_dependence.py:57-69`: bootstrap resamples whole patients with replacement and computes the pair-weighted statistic, matching the declared patient-cluster bootstrap. Saved pair CSV rollup matches the JSON rate.
- `src/141_trial_sim.py:42-60,63-76,94-115`: event detection uses prior-scan nadir; two-scan rule confirms against same prior nadir. KM scan-grid risk sets and log-rank score/variance are internally consistent with the simulated discrete scan times. `src/142_rano_sensitivity.py:94-101` implements stationary marginal-SD AR(1) recurrence correctly. `src/145_conversion_audit.py:47-53` matches the three stated conversions; independent values agree with `conversion_table.csv` to displayed precision.
- `src/rano_calc.py:28-78` implements the stated normal-error forward and flat-prior backward probabilities and MDC. `tests/test_rano_calc.py` contains the seven declared tests; tests were not run in this audit.
- The saved predictions are reproduced by summary booleans: H-S1 true, H-S2 true, H-S3 false, H-S4 false; PA1/PA2/PB1/PD1 true; PX1–PX4 true; PU1–PU4 true. Ledger reports these statuses unchanged. PU statuses are qualified by the A21 grid-coverage finding above.

## Cannot verify

- Exact run reproducibility to three decimals: scripts 140–145 were not run, per instruction.
- Calculator test pass claim: tests were inspected, not run.
- Clinical validity or applicability of the assumed segmentation-noise distributions and the one-study RC% conversion: no patient validation was performed here.
