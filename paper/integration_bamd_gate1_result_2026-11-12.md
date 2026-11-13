# Integration (BAMD) Gate 1 result: MRI does not recover tumor cell states (2026-11-12)

Protocol: `paper/integration_protocol_bamd_2026-11-12.md` (SHA-256 0dc8d720...) and note 1. Result file: `output/integration_bamd/gate1.json`. Features: `output/integration_bamd/features.csv` (505 cases: 102 TCGA, 403 UPENN; SHA-256 01d6d525...). No survival data was read.

## Result

TCGA anchor set: n = 91 patients with BraTS 2021 MRI + expert masks + U133A RNA.

| Cell state (Neftel) | CV Spearman rho (MRI-predicted vs RNA score) | Holm p (1,000 permutations) | Pass (rho >= 0.25 and p < 0.05) |
|---|---|---|---|
| MES | 0.072 | 0.56 | no |
| AC | -0.071 | 0.96 | no |
| OPC | -0.034 | 0.96 | no |
| NPC | -0.054 | 0.96 | no |

**No state passed. Per the protocol, BAMD stops here.** H1 (predicted states add survival information in UPENN) is not tested. The residual split (H2) needs the frozen map, so it is not tested either. UPENN survival was not opened and stays a fresh cohort.

## What this means

With 59 hand-made MRI features and 91 patients, routine pre-operative MRI showed no detectable linear link to bulk-RNA Neftel cell-state scores. This is a clear negative result for the first step of the new method.

## Weakest points (why the negative may not be final)

- Small anchor set (91). A weak true link (rho about 0.2) would be hard to detect.
- Bulk RNA scores mix tumor and normal cells; the RNA sample site is not matched to MRI regions.
- U133A microarray from about 2005-2010; MRI from mixed scanners.
- Features are simple (volumes, shape, intensity statistics), not texture or deep features.
- Linear ridge map only.

## Honest next options (each would need a new dated protocol; none is a re-test of this one)

1. A larger anchor: add CPTAC-GBM (33 with MRI; RNA open at GDC, not yet checked) and IvyGAP (30). Even then n is about 150.
2. Region-matched biology: IvyGAP has RNA from anatomic regions (leading edge, cellular tumor, necrosis) that could be matched to MRI habitats. This is a different, more direct question.
3. Accept the negative and use UPENN for a different pre-registered question.

## Verdict for this method

BAMD as specified does not work at this data size: the MRI-to-biology anchor failed its outcome-free gate.
