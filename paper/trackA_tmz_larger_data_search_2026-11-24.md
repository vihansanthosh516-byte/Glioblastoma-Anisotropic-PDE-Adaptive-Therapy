# Track A TMZ: search for a larger dataset (2026-11-24)

Goal: more patients with untreated-tumour expression and a TMZ response, with a no-TMZ comparison group. Counts below come from metadata only; no survival analysis was run.

| Source | Expression | TMZ yes / no | Access | Main problem |
|---|---|---|---|---|
| TCGA-GBM (GDC) | U133A array about 500, RNA-seq about 170 (from project knowledge, verify) | 617 cases; 356 name temozolomide in GDC treatment records; 574 have treatment yes/no recorded (GDC API, 2026-11-24) | Open | Not randomized (era, fitness, choice). Exact no-TMZ count with expression still to count. |
| CGGA mRNAseq 693 + 325 | RNA-seq | Primary GBM 225: TMZ 170, untreated 41, NA 14 (clinical files, 2026-11-24) | Open | Not randomized; small untreated group; 35 IDH-mutant included. |
| GSE7696 | U133 Plus 2 | 42 / 26 | Open | Spent for the repair-axis score. |
| RTOG 0525 / NRG | Subset | All got TMZ (dose comparison) | Controlled (NCTN archive) | No no-TMZ group. |
| EORTC 26981 / NCIC CE.3 full | Subset | Randomized TMZ+RT vs RT | Controlled, application | Access by application only. |

Assessment: TCGA-GBM is the largest open set with TMZ yes/no and expression. Combined with CGGA it may give roughly 150-250 no-TMZ patients, but treatment was not randomized, so a treatment-by-score interaction is confounded by why patients did or did not get TMZ. The only clean (randomized) large set is the EORTC/NCIC trial, which needs an application.

Next: count TCGA-GBM patients with array expression and a clear TMZ yes/no (metadata only), then decide with the user between (1) open but confounded TCGA+CGGA, and (2) applying for EORTC trial data.
