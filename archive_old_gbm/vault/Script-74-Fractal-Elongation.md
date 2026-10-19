> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Script 74 — Matched-isotropic paired test (D_f and elongation)

**Source:** `output/fractal_aniso_vs_iso.json` (script `src/74_fractal_aniso_vs_iso.py`). n=103 paired simulations using MU-Glioma-Post masks; this is simulation-vs-simulation, not measured longitudinal tumor growth. Design: paired, same patient, rho, seed, dt and steps. Only the tensor differs. Iso tensor = trace(D)/2 * I. Values below were read from the JSON on 2026-10-04.

| Endpoint (150 d) | Aniso | Iso | Diff | Aniso > iso | Cohen dz | Wilcoxon p |
|---|---|---|---|---|---|---|
| Fractal dimension D_f | 0.9248 | 0.9429 | -0.0181 | 8/103 | -1.47 | 2.9e-15 |
| **Elongation** | **1.2119** | **1.0337** | **+0.1782** | **98/103** | **2.89** | **4.1e-19** |
| Tract alignment | 0.4535 | 0.5011 | -0.0475 | 1/103 | -4.39 | 2.4e-19 |

## What this means
- **Elongation is the clean positive result.** Anisotropy stretches the tumour shape in 98 of 103 patients.
- **D_f is not a valid anisotropy metric here.** It goes the wrong way. Tumours are about 120-140 px, so box counting measures size. Drop D_f.
- **Tract alignment goes the wrong way too** (0.45 vs 0.50, 1/103 higher). This was not in the vault before. Do not claim that anisotropic growth "aligns better with tracts" without explaining this. Likely cause: not checked. State it as an open question.
- At 600 d, elongation is 1.56 vs 1.13 (checked2026-10-19 in output/fractal_aniso_vs_iso.json:secondary.elongation_600d, mean1.5592 vs1.1284,99/103 higher; paired simulation).
- Elongation here is on the **patient-tensor model**. Script 78 shows elongation changes sign across horizons in a different setup. Keep the two separate.

## README claim that is retracted
README says "D_f = 1.20-1.55 vs isotropic ~0.0, Cohen's d = 10.43". Script 74 replaces it. Do not cite the README value.

## Script 42 note
Historical script42 median0.909/0-of103-above1.2 claim is unverified here; exact source output and keys are missing. It does not affect the supported script74 retraction above.
