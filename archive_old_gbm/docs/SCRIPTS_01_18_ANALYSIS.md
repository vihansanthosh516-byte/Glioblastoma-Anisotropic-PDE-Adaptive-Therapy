> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Scripts 01–18: historical implementation and output notes

Corrected 2026-10-19. These are legacy exploratory scripts, not the current
probabilistic-RANO pipeline. Use exact filenames because numbers 15 and 16 are
reused. Paths below refer to the repository root. No training, preprocessing,
large H5AD recount or scientific analysis was run for these corrections.

## Data and representation

`src/01_load_and_filter.py` defines UCSC multiomic-gbm singlet/Core/Periphery/
NormalBrain filtering. Local data paths are machine-dependent; current DATA_DIR
is a fixed Windows path, not the old WSL example. Both 01 and 18 still contain
working-directory-relative paths and are not certified as portable current runs.
The source comments report 140355 cells by 29661 genes; exact full-H5AD dimensions
and class counts need a saved input/count manifest before being quoted as checked
outputs. Memory/RSS assertions lack a run log and are withdrawn.

`src/02_preprocess_umap_de.py` implements 5000/class subsampling for visualization
and full-data cell-level DE. HVG selection uses Seurat-v3, 2500 genes and class as
batch key. This exploratory cell-level comparison does not establish independent
patient inference. Script 03 compares marker overlap with a chosen list; overlap
is not validation against a publication. The old claimed 300 combined rows and
published-signature agreement need exact table counts and a paper reference.

`output/nn_X.npy` header is (15000,2500), float32. Its numeric payload is
150000000 bytes (150 MB decimal, about 143 MiB), not a measured process peak.
Script 04 implements row standardization and a seed-42 80/20 stratified random
cell split. That split is not patient-held-out validation.

Scripts 10–11 implement a chosen cVAE objective with patient/region-defined
training-positive pairs. Those labels do not establish biological positive pairs
or meaningful latent causality. Array headers read without model execution:
`output/cgat/cvae_latent.npy` (15000,32), `cvae_latent_full.npy` (140355,32).
The latter float32 payload is 17965440 bytes (18.0 MB decimal, 17.1 MiB).
Checkpoint presence does not verify every training loss/epoch or metadata claim.

## Saved model metrics, not patient-independent performance

| Method / artifact | Accuracy | Macro F1 | Macro AUC | Parameters |
|---|---:|---:|---:|---:|
| LR / `output/benchmark_comparison.tsv` | 0.6983 | 0.6981 | 0.8623 | 0 |
| RF / same TSV | 0.7260 | 0.7248 | 0.8735 | 0 |
| Transformer / same TSV | 0.4960 | 0.4840 | 0.6884 | 80323 |
| Hybrid / same TSV | 0.5120 | 0.4939 | 0.6907 | 80323 |
| Custom VAE / `output/scvi_metrics.json` | 0.729667 | 0.728794 | 0.880037 | not recorded here |
| NMF / `output/nmf_metrics.json` | 0.610333 | 0.611204 | 0.789947 | not recorded here |
| C-GAT / `output/cgat/gat_metrics.json` | 0.787333 | 0.783116 | 0.921789 | 41499 |

`output/nn_metrics.json` separately gives attention NN best_test_accuracy=0.79,
final_test_accuracy=0.7843333, n_params=488899. It is not a 3.2-million-parameter
model. Attention weights are model attribution, not a causal biological explanation.
Hybrid has slightly higher saved aggregate scores than Transformer; one split
with no fair independent uncertainty comparison does not establish superiority.
The old 94–96% LR/RF and 92–95% GAT ranges are unsupported by these artifacts.

The source named `src/15_baseline_scvi.py` is a custom VAE, not scvi-tools. Its
filename/artifact label does not identify the external package. NMF components
are 0–3; positional AC/MES/NPC/OPC names are not independently validated cell
identities. Former spatial/biological representation claims are withdrawn.

Script 12 specifies K=15 excluding self in the subsample graph. Expected
15000x15 directed edges are an implementation expectation, not recovered graph
metadata. `output/cgat/gat_edge_index.npy`, `gat_edge_attr.npy` and `gat_meta.json`
are absent. The graph uses patient/region metadata, not demonstrated spatial
coordinates. Scripts 13 and 16 use a transductive full graph and random cell
splits with patient overlap. Test accuracy is consulted repeatedly for checkpoint
selection, so it is a selected-set score. The full-graph branch also includes
label-derived `class_diff`, a direct leakage source. These are historical defects,
not independent validation.

Ledger 2o (`output/cgat_leak_free.json`, `output/trackA_audit.json`) gives the
leak-free patient-fold C-GAT mean 63.8% with large fold variation and donor,
diagnosis and fresh/frozen confounding. Patient-ID lookup alone reaches 82.8% on
the original split. Cell-level accuracy does not establish strong biological
separability or a graph advantage. Script 14 was skipped; the old final
leaderboard was assembled manually, not emitted as a verified final comparison.
`src/09_benchmark_comparison.py` is excluded from the runnable CLI inventory;
no complete 01–18 reproduction chain is certified.

## Exploratory state-score summaries

Script 17's entropy and NMF-gradient claims have a figure but no cited numerical
table/test here; they are exploratory visual impressions, not a confirmed
intermediate biological state. Script 18 defines a researcher-chosen composite
score combining NMF entropy and distance to the healthy centroid. Weights
alpha=beta=0.5 and sigma=1 are assumptions (`output/csgt_metrics.json`).

That output's Kruskal-Wallis H=141.717101, p=1.6846998e-31, with 5000 cells per
zone, detects distributional differences. Means are Healthy 0.586675, Periphery
0.649461, Core 0.646952. They are not monotonically Healthy<Periphery<Core.
The test proves neither continuity, ordering, a trajectory, nor a biological
intermediate state. “CSGT proof,” “true intermediate” and validated-transition
claims are withdrawn.

## Missing and moved artifacts

The former final-benchmark flow and output list is not a record that all outputs
were produced. At the named root paths these are unavailable:
`output/cgat/cvae_labels_full.npy`, `output/cgat/gat_edge_index_full.npy`,
`output/final_benchmark_bar.png`, `output/final_benchmark_comparison.tsv`,
`output/final_benchmark_confusion.png`, `output/final_benchmark_perclass_heatmap.png`,
`output/nn_model.pt`, and `output/scvi_model.pt`.
`output/BENCHMARK_SUMMARY.md` was moved to `../output/BENCHMARK_SUMMARY.md` inside
this archive. The extra missing graph metadata named above also needs recovery
before exact graph certification. No missing artifact was generated by this audit.

Current trust hierarchy: `../../START_HERE.md`; corrected ledger:
`../../vault/PAPER_FINDINGS_LEDGER.md`; correction tracking:
`../../paper/codex_audit2/fix_progress.md`. Historical setup/outputs are background,
not current safe-runnable instructions or clinically validated findings.
