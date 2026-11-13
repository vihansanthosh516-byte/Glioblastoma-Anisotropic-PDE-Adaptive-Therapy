# Integration (BAMD) Gate 1b result: biopsy-site MRI does not predict biopsy biology in held-out patients (2026-11-12)

Protocol: `paper/integration_protocol_bamd_gate1b_2026-11-12.md` (committed f451300 before analysis). Result: `output/integration_bamd/gate1b.json`. Data: Hu et al. 2023 public table, 95 MRI-localized biopsies, 25 patients, 0 rows dropped. No survival data.

## Result (leave-one-patient-out ridge; 2,000 CE/NE-stratified permutations)

| Pathway score | rho, with CE flag | Holm p | rho, without CE flag | Holm p |
|---|---|---|---|---|
| GPM (glycolytic / plurimetabolic) | -0.386 | 1.0 | -0.356 | 1.0 |
| MTC (mitochondrial) | -0.054 | 1.0 | -0.037 | 0.88 |
| NEU (neuronal) | 0.238 | 1.0 | 0.210 | 0.43 |
| PPR (proliferative / progenitor) | 0.178 | 1.0 | 0.051 | 1.0 |

**No target passed (rule: rho >= 0.25 and Holm p < 0.05).** Stage 2 (transfer to UPENN) is not run. UPENN survival stays unopened.

## Reading the numbers

- The negative GPM rho is a known leave-one-group-out effect: when a held-out patient's mean differs from the training mean, out-of-patient predictions can move against the truth. It shows there is no stable patient-transferable signal, not a real inverse link.
- NEU (0.24) was the strongest positive value but did not beat the permutation null. The null is wide because the permutation keeps the CE/NE split, which already carries much of the variance.
- Hu et al. 2023 reported within-cohort associations with patient random effects. This test asks a harder question (prediction for a new patient) and the answer, at n = 25 patients, is no.

## Weakest points

25 patients; only 7 local MRI values; ridge only; pathway scores as published. A within-patient (not across-patient) link may still exist; this test cannot see it.

## Combined integration status

Two independent, pre-registered tests of "MRI carries tumor biology" both failed:
1. Whole-tumor MRI vs bulk RNA cell states, TCGA, n = 91 (Gate 1).
2. Biopsy-site MRI vs biopsy RNA pathways, Hu 2023, 25 patients / 95 biopsies (Gate 1b).
