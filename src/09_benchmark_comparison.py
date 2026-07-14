"""
09_benchmark_comparison.py
==========================
Head-to-head comparison: Method 1 (Classical) vs Method 2 (Transformer) vs Method 3 (Hybrid)
Produces:
  - Combined metrics TSV
  - Bar chart of accuracy / macro F1 / AUC
  - Per-class heatmap
  - Summary markdown for science fair poster
"""
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

ROOT    = "/mnt/c/Users/vihan/20206 science fair"
OUT_DIR = os.path.join(ROOT, "output")

# Input metric files
M1_FILE = os.path.join(OUT_DIR, "method1_metrics.json")
M2_FILE = os.path.join(OUT_DIR, "method2_metrics.json")
M3_FILE = os.path.join(OUT_DIR, "method3_metrics.json")

# Outputs
TSV_OUT  = os.path.join(OUT_DIR, "benchmark_comparison.tsv")
PNG_BAR  = os.path.join(OUT_DIR, "benchmark_bar_chart.png")
PNG_HEAT = os.path.join(OUT_DIR, "benchmark_per_class_heatmap.png")
MD_OUT   = os.path.join(OUT_DIR, "BENCHMARK_SUMMARY.md")

def load_metrics(path):
    with open(path) as f:
        return json.load(f)

# Load all
m1 = load_metrics(M1_FILE)
m2 = load_metrics(M2_FILE)
m3 = load_metrics(M3_FILE)

# We'll use LogisticRegression as the representative of Method 1
lr_metrics = m1["LogisticRegression"]
rf_metrics = m1["RandomForest"]

methods = [
    ("Logistic Regression (Classical)", lr_metrics),
    ("Random Forest (Classical)", rf_metrics),
    ("Transformer (Deep Learning)", m2),
    ("Hybrid (LR-prior Transformer)", m3),
]

# Build comparison table
rows = []
for name, m in methods:
    rows.append({
        "Method": name,
        "Accuracy": f"{m['accuracy']:.4f}",
        "Macro F1": f"{m['macro_f1']:.4f}",
        "Weighted F1": f"{m['weighted_f1']:.4f}",
        "Macro Precision": f"{m['macro_precision']:.4f}",
        "Macro Recall": f"{m['macro_recall']:.4f}",
        "Macro AUC (OvR)": f"{m['macro_auc_ovr']:.4f}",
        "Params": f"{int(m.get('n_params', 0)):,}" if isinstance(m.get('n_params', 0), (int, float)) else "N/A",
        "Best Test Acc": f"{m.get('best_test_acc', m['accuracy']):.4f}",
    })

df = pd.DataFrame(rows)
df.to_csv(TSV_OUT, sep="\t", index=False)
print(f"Saved {TSV_OUT}")
print(df.to_string(index=False))

# ---- 1. Bar chart: Accuracy / Macro F1 / AUC ----
fig, axes = plt.subplots(1, 3, figsize=(14, 4))
colors = ["#4c72b0", "#55a868", "#c44e52", "#8172b3"]
metrics_to_plot = [
    ("Accuracy", "accuracy"),
    ("Macro F1", "macro_f1"),
    ("Macro AUC (OvR)", "macro_auc_ovr"),
]
for ax, (title, key) in zip(axes, metrics_to_plot):
    vals = [m[key] for _, m in methods]
    bars = ax.bar(range(len(methods)), vals, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels([n.split("(")[0].strip() for n, _ in methods], rotation=15, ha="right")
    ax.set_ylim(0, 1.0)
    ax.set_title(title)
    ax.set_ylabel("Score")
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f"{val:.3f}", ha="center", va="bottom", fontsize=9)
fig.suptitle("Method Comparison: Classical vs Deep vs Hybrid", fontsize=14, fontweight="bold")
fig.tight_layout()
fig.savefig(PNG_BAR, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved {PNG_BAR}")

# ---- 2. Per-class heatmap (Precision, Recall, F1) ----
class_names = ["Core", "Periphery", "Healthy"]
heat_data = []
for name, m in methods:
    per = m["per_class"]
    for cn in class_names:
        heat_data.append({
            "Method": name.split("(")[0].strip(),
            "Class": cn,
            "Precision": per[cn]["precision"],
            "Recall": per[cn]["recall"],
            "F1": per[cn]["f1-score"],
        })
heat_df = pd.DataFrame(heat_data)

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
for ax, metric in zip(axes, ["Precision", "Recall", "F1"]):
    pivot = heat_df.pivot(index="Method", columns="Class", values=metric)
    pivot = pivot.loc[[m.split("(")[0].strip() for m, _ in methods], class_names]
    sns.heatmap(pivot, annot=True, fmt=".3f", cmap="RdYlGn", vmin=0, vmax=1,
                cbar_kws={"label": metric}, ax=ax)
    ax.set_title(metric)
fig.suptitle("Per-Class Performance", fontsize=14, fontweight="bold")
fig.tight_layout()
fig.savefig(PNG_HEAT, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved {PNG_HEAT}")

# ---- 3. Confusion matrix comparison (side by side) ----
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
for ax, (name, m) in zip(axes.flatten(), methods):
    cm = np.array(m["confusion_matrix"])
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues",
                xticklabels=class_names, yticklabels=class_names, ax=ax, cbar=False)
    ax.set_title(name.split("(")[0].strip())
    ax.set_xlabel("Predicted")
    ax.set_ylabel("True")
fig.suptitle("Confusion Matrices", fontsize=14, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "benchmark_confusion_matrices.png"), dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved benchmark_confusion_matrices.png")

# ---- 4. Markdown summary for poster ----
with open(MD_OUT, "w") as f:
    f.write("# Benchmark Summary: Core vs Periphery vs Healthy Classification\n\n")
    f.write("**Dataset:** 15,000 cells (5,000 per class) from multiomic-gbm scRNA-seq\n")
    f.write("**Features:** Top 100 HVGs (row-z-scored)\n")
    f.write("**Split:** 80/20 stratified (12,000 train / 3,000 test)\n\n")
    f.write("## Overall Metrics\n\n")
    f.write(df.to_markdown(index=False))
    f.write("\n\n## Per-Class F1 Scores\n\n")
    f1_pivot = heat_df.pivot(index="Method", columns="Class", values="F1")
    f1_pivot = f1_pivot.loc[[m.split("(")[0].strip() for m, _ in methods], class_names]
    f.write(f1_pivot.to_markdown())
    f.write("\n\n## Key Findings\n\n")
    best = df.loc[df["Macro F1"].astype(float).idxmax(), "Method"]
    f.write(f"- **Best overall: {best}** (Macro F1 = {df['Macro F1'].max()})\n")
    f.write(f"- Classical methods (LR, RF) outperform deep learning on this small-feature, tabular-style data.\n")
    f.write(f"- The Hybrid method (injecting LR coefficients as attention bias) did not improve over the Transformer baseline, likely due to the limited sequence length (100 genes) and shallow model.\n")
    f.write(f"- All methods struggle most with **Periphery** class (intermediate biology between Core and Healthy).\n")
    f.write(f"- Recommendation for science fair: Lead with **Random Forest** as the primary result, show Transformer as 'deep learning attempt', Hybrid as 'novel architecture exploration'.\n")

print(f"Saved {MD_OUT}")
print("\nDONE.")