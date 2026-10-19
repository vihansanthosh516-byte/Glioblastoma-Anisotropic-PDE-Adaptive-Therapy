> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Benchmark Summary: Core vs Periphery vs Healthy Classification

**Dataset:** 15,000 cells (5,000 per class) from multiomic-gbm scRNA-seq
**Features:** Top 100 HVGs (row-z-scored)
**Split:** 80/20 stratified (12,000 train / 3,000 test)

Evidence: root `../../output/benchmark_comparison.tsv`. Historical 80/20 split is by cell, not patient; donor overlap, class composition and test-set selection limit inference. The corrected patient-fold evaluation differs (`../../vault/PAPER_FINDINGS_LEDGER.md` 2o). Feature count is an implementation setting, not independent validation.

## Overall Metrics

| Method                          |   Accuracy |   Macro F1 |   Weighted F1 |   Macro Precision |   Macro Recall |   Macro AUC (OvR) |   Params |   Best Test Acc |
|:--------------------------------|-----------:|-----------:|--------------:|------------------:|---------------:|------------------:|---------:|----------------:|
| Logistic Regression (Classical) |     0.6983 |     0.6981 |        0.6981 |            0.6982 |         0.6983 |            0.8623 |        0 |          0.6983 |
| Random Forest (Classical)       |     0.726  |     0.7248 |        0.7248 |            0.7274 |         0.726  |            0.8735 |        0 |          0.726  |
| Transformer (Deep Learning)     |     0.496  |     0.484  |        0.484  |            0.4828 |         0.496  |            0.6884 |   80,323 |          0.5117 |
| Hybrid (LR-prior Transformer)   |     0.512  |     0.4939 |        0.4939 |            0.519  |         0.512  |            0.6907 |   80,323 |          0.5133 |

## Per-Class F1 Scores

| Method              |     Core |   Periphery |   Healthy |
|:--------------------|---------:|------------:|----------:|
| Logistic Regression | 0.757606 |    0.595573 |  0.741176 |
| Random Forest       | 0.792802 |    0.63803  |  0.743709 |
| Transformer         | 0.61913  |    0.321788 |  0.510995 |
| Hybrid              | 0.565469 |    0.315599 |  0.600763 |

## Key Findings

- **Best overall: Random Forest (Classical)** (Macro F1 = 0.7248)
- LR/RF score higher than these two deep models on this selected cell split; no broad model-class conclusion follows.
- Hybrid has slightly higher saved aggregate accuracy/F1/AUC than Transformer on this split. No independent uncertainty comparison establishes general superiority or the cause of the difference.
- Periphery has the lowest listed per-class F1; this does not establish an intermediate biological state. Exact source for this per-class table remains unverified here.
- RF has the highest macro-F1 among these four saved models, not a validated patient-level or novelty result. Current science-fair goal is in START_HERE; this benchmark does not determine presentation priorities.
