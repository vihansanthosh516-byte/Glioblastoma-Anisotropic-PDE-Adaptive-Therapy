# Track A: Core-vs-edge programs in labeled GBM tissue are explained by known confounds

Status: final for this cohort (2026-10-29). Two pre-registered tests; both fail their own confound rule. No third test on this cohort.

## Question
In patient-matched tumor core and edge samples, do fixed external gene programs differ beyond (a) expression-matched random gene sets, (b) the hypoxia gradient, and (c) normal-brain neuron content?

## Cohort (labels not derived from expression)
- Ivy GAP, 8 donors: LE-reference-histology (edge) vs CT-reference-histology (core). Histology labels. RNA-seq (FPKM).
- Angers E-MTAB-5092, 7 patients: interface (edge) vs tumor zone (core). Preoperative MRI, neuronavigation, pathologist review (Aubry 2015). Agilent arrays.
- Total 15 patients. Each patient gives one paired value. Cohorts analysed on their own platform, then pooled by within-cohort scaling.

## Method (sealed before scoring)
- Score = mean within-cohort z-score of the list's genes.
- Per patient: E = (edge - core score) - median (edge - core) over 1,000 random gene sets matched on mean-expression decile (seed 20261029).
- Exact two-sided sign-flip test over all 2^15 patterns (p floor 6.1e-5).
- Success also needed: same sign in both cohorts, and significance after regressing out confound scores computed the same way.
- Protocols and gene lists hashed and committed before scores were computed: `paper/trackA_senmayo/` (commit c02217b), `paper/trackA_oxphos/` (commit 162e55a).

## Results

| Test | Direction | p (unadjusted) | Cohorts agree | Adjusted for | p (adjusted) | Source |
|---|---|---|---|---|---|---|
| 1. SenMayo senescence/SASP (125 genes; alpha 0.05) | core > edge | 0.0035 | yes | hypoxia | 0.37 | `output/trackA_senmayo/result.json` |
| 2. Hallmark OXPHOS (200 genes; alpha 0.025) | edge > core | 0.0024 | yes | hypoxia + neuron markers | 0.67 | `output/trackA_oxphos/result.json` |

Per-cohort, test 1: Ivy GAP p = 0.063, Angers p = 0.047. Test 2: Ivy GAP p = 0.0078, Angers p = 0.23. (Same files, `secondaries`.)

## Interpretation
- Both programs show strong, reproducible core-vs-edge shifts in raw tissue.
- Both shifts vanish once known tissue gradients are removed: the senescence/SASP signal tracks hypoxia in the core; the OXPHOS signal tracks neuron-rich normal brain at the edge.
- So, in bulk labeled tissue, these two programs give no evidence of tumor-cell-intrinsic core-vs-edge biology beyond those gradients.
- Value for the framework: it shows that bulk core-vs-edge contrasts must be adjusted for hypoxia and normal-brain content before any claim of new biology. Tracks B and C do not use Track A outputs; there is no patient-level join.

## Limits (stated in full)
- Bulk samples mix tumor and normal cells; adjustment by marker scores is partial and can remove real biology that co-varies with the confounds (a null after adjustment is "not shown", not "absent").
- 15 patients; two platforms; Angers interface samples are often mostly brain or mostly tumor (Aubry 2015).
- Only two programs tested; other programs were not tested and the cohort is now closed to avoid multiple testing.
- Neuron covariate (LEIN_NEURON_MARKERS) does not capture astrocytes or oligodendrocytes.
- Prior-art searches were limited web searches, not systematic reviews.
- Angers preprocessing (gProcessedSignal, log2, quantile normalization, probe mean) was fixed in code before scoring but was not written into the first sealed protocol.
- Earlier rejected or dropped Track A candidates (MES shift, myeloid ratio, TMZ tolerance, within-slide distance, Gill cohort) are recorded in `START_HERE.md`.
