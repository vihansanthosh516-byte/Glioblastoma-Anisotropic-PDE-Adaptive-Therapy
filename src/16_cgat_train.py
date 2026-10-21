"""
16_cgat_train.py
================
Train Edge-Aware GAT on FULL 140k-cell graph using full-batch training.
No neighbor sampling needed - full graph fits in memory.

NOTE (2026-10-19 audit): historical legacy Track A script. Test accuracy is evaluated every epoch and used for checkpoint selection ((UPDATED 2026-10-19: the epoch is now chosen on a validation split carved from train and test is scored once; the text in this sentence describes the older behaviour: the test set was part of model selection)); no torch seed is set; the split is per cell. Needs torch and torch_geometric (not in Requirements.txt). The root was a hard-coded WSL path and now derives from __file__.
"""
import os, time, resource, gc, json
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.nn import GATv2Conv
from sklearn.metrics import (accuracy_score, f1_score, precision_recall_fscore_support,
                             confusion_matrix, classification_report, roc_auc_score)
from sklearn.preprocessing import label_binarize
from sklearn.model_selection import train_test_split
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
# GLOBAL_SEED_2026_10_19: fixed seed added by the audit corrections (legacy script was unseeded); outputs from earlier runs came from unseeded runs.
_SEED = int(os.environ.get("GBM_SEED", "42"))  # seed sweep; default 42
import random as _random_seed
_random_seed.seed(_SEED)
try:
    import numpy as _np_seed
    _np_seed.random.seed(_SEED)
except ImportError:
    pass
try:
    import torch as _torch_seed
    _torch_seed.manual_seed(_SEED)
except ImportError:
    pass

ROOT = __import__("os").path.dirname(__import__("os").path.dirname(__import__("os").path.abspath(__file__)))  # was a hard-coded WSL path
OUT_DIR = os.path.join(ROOT, "output", "cgat")
# GBM_SPLIT_MODE=patient: patient-disjoint split, results under output/cgat/patient_split/ (default "cell" keeps the original behaviour and names).
SPLIT_MODE = os.environ.get("GBM_SPLIT_MODE", "cell")
RES_DIR = OUT_DIR if SPLIT_MODE == "cell" else os.path.join(OUT_DIR, "patient_split_seeds", f"seed{_SEED}")  # one folder per seed: no overwrite
os.makedirs(RES_DIR, exist_ok=True)

# upstream inputs (scripts 10-15) stay in the shared folder; only outputs are per-seed
INPUT_DIR = OUT_DIR if SPLIT_MODE == "cell" else os.path.join(OUT_DIR, "patient_split")
EDGE_INDEX_PATH = os.path.join(INPUT_DIR, "gat_edge_index_full.npy")
EDGE_ATTR_PATH = os.path.join(INPUT_DIR, "gat_edge_attr_full.npy")
LATENT_PATH = os.path.join(INPUT_DIR, "cvae_latent_full.npy")
LABELS_PATH = os.path.join(INPUT_DIR, "cvae_labels_full.npy")

MODEL_PATH = os.path.join(RES_DIR, "gat_full_model.pt")
METRICS_PATH = os.path.join(RES_DIR, "gat_full_metrics.json")
PRED_PATH = os.path.join(RES_DIR, "gat_full_predictions.tsv")
PNG_CM = os.path.join(RES_DIR, "gat_full_confusion.png")
PNG_LOSS = os.path.join(RES_DIR, "gat_full_training_loss.png")

# Hyperparameters
HIDDEN_DIM = 64
HEADS = 8
N_LAYERS = 2
DROPOUT = 0.3
LR = 1e-3
WEIGHT_DECAY = 1e-4
EPOCHS = 50
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

def mem_mb(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}  RSS={mem_mb():.0f}MB  DEV={DEVICE}", flush=True)

# ---- 1. Load Data ----
log("Loading full graph...")
edge_index = np.load(EDGE_INDEX_PATH).astype(np.int64)
edge_attr = np.load(EDGE_ATTR_PATH).astype(np.float32)
latent = np.load(LATENT_PATH).astype(np.float32)
labels = np.load(LABELS_PATH, allow_pickle=True)

class_names = ["Healthy", "Periphery", "Core"]
class_to_idx = {c: i for i, c in enumerate(class_names)}
y = np.array([class_to_idx[c] for c in labels], dtype=np.int64)

n_cells, in_dim = latent.shape
log(f"Data: {n_cells} cells, {in_dim} features, {edge_index.shape[1]} edges")
log(f"Class dist: {np.bincount(y)}")

# Convert to tensors
x = torch.from_numpy(latent).float()
edge_index_t = torch.from_numpy(edge_index).long()
edge_attr_t = torch.from_numpy(edge_attr).float()
y_t = torch.from_numpy(y).long()

# Train/test split. "cell" = original stratified random cell split (leaky: cells of one patient on both sides).
# "patient" = StratifiedGroupKFold fold 0 with groups = cvae_patient_full.npy (no patient in both sets).
if SPLIT_MODE == "cell":
    train_idx, test_idx = train_test_split(np.arange(n_cells), test_size=0.2, random_state=42, stratify=y)
else:
    from sklearn.model_selection import StratifiedGroupKFold
    _pf = os.path.join(INPUT_DIR, "cvae_patient_full.npy")
    if not os.path.exists(_pf):  # script 11 never ran at 140k; use the same 15k-cell patient IDs as scripts 05/07/08/13
        _pf = os.path.join(OUT_DIR, "nn_patient_ids.npy")
    if not os.path.exists(_pf):
        raise SystemExit(f"patient split needs {_pf} (written by script 11); not found")
    groups = np.load(_pf, allow_pickle=True).astype(str)
    train_idx, test_idx = next(StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42).split(np.arange(n_cells), y, groups))
    assert not (set(groups[train_idx]) & set(groups[test_idx])), "patient overlap between train and test"
# Validation cells carved from TRAIN (10%, stratified, seed 42); epoch chosen on validation, test scored once.
if SPLIT_MODE == "patient":
    # Patient mode: validation grouped by patient, so no donor is in both train and validation.
    _a, _b = next(StratifiedGroupKFold(n_splits=10, shuffle=True, random_state=42).split(np.zeros(len(train_idx)), y[train_idx], groups[train_idx]))
    assert not (set(groups[train_idx[_a]]) & set(groups[train_idx[_b]])), "patient overlap between train and validation"
    train_idx, val_idx = train_idx[_a], train_idx[_b]
else:
    train_idx, val_idx = train_test_split(train_idx, test_size=0.10, random_state=42, stratify=y[train_idx])
train_mask = torch.zeros(n_cells, dtype=torch.bool)
train_mask[train_idx] = True
val_mask = torch.zeros(n_cells, dtype=torch.bool)
val_mask[val_idx] = True
test_mask = torch.zeros(n_cells, dtype=torch.bool)
test_mask[test_idx] = True

log(f"Train: {train_mask.sum().item()}, Test: {test_mask.sum().item()}")

data = Data(x=x, edge_index=edge_index_t, edge_attr=edge_attr_t, y=y_t)
data = data.to(DEVICE)

# ---- 2. Edge-Aware GAT ----
class EdgeAwareGAT(nn.Module):
    def __init__(self, in_dim, hidden_dim, out_dim, heads=8, n_layers=2, dropout=0.3, edge_dim=6):
        super().__init__()
        self.n_layers = n_layers
        self.dropout = dropout
        
        self.convs = nn.ModuleList()
        self.bns = nn.ModuleList()
        
        # Layer 1
        self.convs.append(GATv2Conv(in_dim, hidden_dim, heads=heads, dropout=dropout,
                                    edge_dim=edge_dim, add_self_loops=True, concat=True))
        self.bns.append(nn.BatchNorm1d(hidden_dim * heads))
        
        # Hidden layers
        for _ in range(n_layers - 2):
            self.convs.append(GATv2Conv(hidden_dim * heads, hidden_dim, heads=heads, dropout=dropout,
                                        edge_dim=edge_dim, add_self_loops=True, concat=True))
            self.bns.append(nn.BatchNorm1d(hidden_dim * heads))
        
        # Output layer (single head, no concat)
        self.convs.append(GATv2Conv(hidden_dim * heads, out_dim, heads=1, dropout=dropout,
                                    edge_dim=edge_dim, add_self_loops=True, concat=False))
    
    def forward(self, x, edge_index, edge_attr):
        for i, conv in enumerate(self.convs):
            x = conv(x, edge_index, edge_attr=edge_attr)
            if i < len(self.convs) - 1:
                x = F.elu(x)
                x = self.bns[i](x)
                x = F.dropout(x, p=self.dropout, training=self.training)
        return x

model = EdgeAwareGAT(
    in_dim=in_dim,
    hidden_dim=HIDDEN_DIM,
    out_dim=len(class_names),
    heads=HEADS,
    n_layers=N_LAYERS,
    dropout=DROPOUT,
    edge_dim=edge_attr.shape[1]
).to(DEVICE)

log(f"Model params: {sum(p.numel() for p in model.parameters()):,}")

optimizer = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=EPOCHS)
criterion = nn.CrossEntropyLoss()

# ---- 3. Training (Full-Batch) ----
log("Training (full-batch)...")
train_losses = []
test_losses = []
test_accs = []
best_acc = 0

for epoch in range(1, EPOCHS + 1):
    model.train()
    t0 = time.time()
    
    optimizer.zero_grad()
    out = model(data.x, data.edge_index, data.edge_attr)
    loss = criterion(out[train_mask], data.y[train_mask])
    loss.backward()
    optimizer.step()
    scheduler.step()
    train_loss = loss.item()
    train_losses.append(train_loss)
    
    # Validation
    model.eval()
    with torch.no_grad():
        out = model(data.x, data.edge_index, data.edge_attr)
        test_loss = criterion(out[val_mask], data.y[val_mask]).item()  # validation loss (name kept)
        test_losses.append(test_loss)
        
        pred = out[val_mask].argmax(dim=-1)
        acc = (pred == data.y[val_mask]).float().mean().item()  # validation accuracy
        test_accs.append(acc)
        
        if acc > best_acc:
            best_acc = acc
            torch.save(model.state_dict(), MODEL_PATH)
    
    log(f"Epoch {epoch:3d}/{EPOCHS}  TrainLoss={train_loss:.4f}  ValLoss={test_loss:.4f}  ValAcc={acc:.4f}  Best={best_acc:.4f}  ({time.time()-t0:.1f}s)")

# ---- 4. Final Evaluation ----
log("Final evaluation...")
model.load_state_dict(torch.load(MODEL_PATH))
model.eval()
with torch.no_grad():
    out = model(data.x, data.edge_index, data.edge_attr)
    probs = F.softmax(out[test_mask], dim=-1).cpu().numpy()
    pred = probs.argmax(axis=1)
    true = data.y[test_mask].cpu().numpy()

# Metrics
acc = accuracy_score(true, pred)
np.savez(os.path.join(RES_DIR, f"test_predictions_cgat16_seed{_SEED}.npz"), y_true=true, y_pred=pred, test_index=test_idx)
macro_f1 = f1_score(true, pred, average="macro")
weighted_f1 = f1_score(true, pred, average="weighted")
prec, rec, _, _ = precision_recall_fscore_support(true, pred, average="macro")
y_true_bin = label_binarize(true, classes=[0,1,2])
auc = roc_auc_score(y_true_bin, probs, average="macro", multi_class="ovr")
cm = confusion_matrix(true, pred, labels=[0,1,2])
report = classification_report(true, pred, target_names=class_names, output_dict=True)

metrics = {
    "method": "C-GAT (Full Graph, Full-Batch)",
    "accuracy": float(acc),
    "macro_f1": float(macro_f1),
    "weighted_f1": float(weighted_f1),
    "macro_precision": float(prec),
    "macro_recall": float(rec),
    "macro_auc_ovr": float(auc),
    "confusion_matrix": cm.tolist(),
    "per_class": report,
    "best_val_acc": best_acc,
    "selection": "epoch chosen on a validation split carved from train; test scored once",
    "split_mode": SPLIT_MODE,
    "final_test_acc": float(acc),
    "n_params": sum(p.numel() for p in model.parameters()),
    "n_cells": int(n_cells),
    "n_edges": int(edge_index.shape[1]),
    "epochs": EPOCHS
}

log(f"C-GAT Full: acc={acc:.4f}, macro_f1={macro_f1:.4f}, AUC={auc:.4f}")

# Save metrics
with open(METRICS_PATH, "w") as f:
    json.dump(metrics, f, indent=2)

# Save predictions
import pandas as pd
pred_df = pd.DataFrame({
    "true": [class_names[i] for i in true],
    "pred": [class_names[i] for i in pred],
    **{f"prob_{cn}": probs[:, i] for i, cn in enumerate(class_names)}
})
pred_df.to_csv(PRED_PATH, sep="\t", index=False)

# Plots
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].plot(train_losses, label="train")
axes[0].plot(test_losses, label="val")
axes[0].set_title("Loss"); axes[0].legend()
axes[1].plot(test_accs, label="val acc", color="green")
axes[1].set_title("Val Accuracy"); axes[1].legend()
fig.tight_layout()
fig.savefig(PNG_LOSS, dpi=180); plt.close(fig)

fig, ax = plt.subplots(figsize=(5, 4))
sns.heatmap(cm, annot=True, fmt="d", cmap="Blues",
            xticklabels=class_names, yticklabels=class_names, ax=ax)
ax.set_title("C-GAT Full Graph"); ax.set_xlabel("Predicted"); ax.set_ylabel("True")
fig.tight_layout()
fig.savefig(PNG_CM, dpi=180); plt.close(fig)

log(f"Saved: {METRICS_PATH}, {PRED_PATH}, {PNG_CM}, {PNG_LOSS}")
log("DONE.")

if __name__ == "__main__":
    import sys
    sys.exit(0)