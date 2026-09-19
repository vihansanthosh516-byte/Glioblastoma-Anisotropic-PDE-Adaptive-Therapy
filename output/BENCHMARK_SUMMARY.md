# Benchmark Summary: Core vs Periphery vs Healthy Classification

**Dataset:** 15,000 cells (5,000 per class) from multiomic-gbm scRNA-seq
**Features:** Top 100 HVGs (row-z-scored)
**Split:** 80/20 stratified (12,000 train / 3,000 test)

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
- Classical methods (LR, RF) outperform deep learning on this small-feature, tabular-style data.
- The Hybrid method (injecting LR coefficients as attention bias) did not improve over the Transformer baseline, likely due to the limited sequence length (100 genes) and shallow model.
- All methods struggle most with **Periphery** class (intermediate biology between Core and Healthy).
- Recommendation for science fair: Lead with **Random Forest** as the primary result, show Transformer as 'deep learning attempt', Hybrid as 'novel architecture exploration'.
