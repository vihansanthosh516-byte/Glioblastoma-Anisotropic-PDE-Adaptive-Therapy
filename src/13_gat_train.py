"""
13_gat_train.py
================
Train Graph Attention Network (GAT) with edge-aware attention on the cVAE latent graph.

Architecture:
- 2-layer GAT with 8 heads + 1 output head
- Edge features injected via attention bias (edge-aware attention)
- Trained on 15k balanced cells (5k per class)
- 80/20 train/test split (same as other methods for fair comparison)

NOTE (2026-10-19 audit): historical Track A script. The split is per cell, not patient-held-out, and test accuracy is
evaluated every epoch with the checkpoint saved when test accuracy improves, so the reported accuracy is selected on
the test set (optimistic). No random seed is set. Needs torch, torch_geometric and anndata (not in Requirements.txt).
The edge attributes come from script 12 and encode tissue-label information (see its note). Paths now derive from __file__.
The (15000, 32) latent shape is an assumption, not checked.
"""
import os, time, gc, json
try:
    import resource
except ImportError:  # native Windows Python has no resource module; RSS then reports 0
    resource = None
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.nn import GATv2Conv
from sklearn.metrics import (accuracy_score, f1_score, precision_recall_fscore_support,
                             confusion_matrix, classification_report, roc_auc_score)
from sklearn.preprocessing import label_binarize
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

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(ROOT, "output", "cgat")
# GBM_SPLIT_MODE=patient: patient-disjoint split, results under output/cgat/patient_split/ (default "cell" keeps the original behaviour and names).
SPLIT_MODE = os.environ.get("GBM_SPLIT_MODE", "cell")
RES_DIR = OUT_DIR if SPLIT_MODE == "cell" else os.path.join(OUT_DIR, "patient_split")
os.makedirs(RES_DIR, exist_ok=True)

EDGE_INDEX_PATH = os.path.join(RES_DIR, "gat_edge_index.npy")
EDGE_ATTR_PATH = os.path.join(RES_DIR, "gat_edge_attr.npy")
LATENT_PATH = os.path.join(RES_DIR, "cvae_latent.npy")
ADATA_SUB_PATH = os.path.join(ROOT, "output", "02_adata_subsampled.h5ad")
SPLIT_PATH = os.path.join(ROOT, "output", "nn_train_test_split.tsv" if SPLIT_MODE == "cell" else "nn_train_test_split_patient.tsv")

MODEL_PATH = os.path.join(RES_DIR, "gat_model.pt")
METRICS_PATH = os.path.join(RES_DIR, "gat_metrics.json")
PRED_PATH = os.path.join(RES_DIR, "gat_predictions.tsv")
PNG_CM = os.path.join(RES_DIR, "gat_confusion.png")
PNG_LOSS = os.path.join(RES_DIR, "gat_training_loss.png")

# Hyperparameters
HIDDEN_DIM = 64
HEADS = 8
N_LAYERS = 2
DROPOUT = 0.3
LR = 1e-3
WEIGHT_DECAY = 1e-4
EPOCHS = 100
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

def mem_mb(): return (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss if resource else 0) / 1024
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}  RSS={mem_mb():.0f}MB  DEV={DEVICE}", flush=True)

# ---- 1. Load data ----
log("Loading data...")
latent = np.load(LATENT_PATH).astype(np.float32)  # (15000, 32)
edge_index = np.load(EDGE_INDEX_PATH).astype(np.int64)  # (2, E)
edge_attr = np.load(EDGE_ATTR_PATH).astype(np.float32)  # (E, 5)

import anndata as ad
adata = ad.read_h5ad(ADATA_SUB_PATH)
labels = adata.obs["class"].values  # Core, Periphery, Healthy
class_names = ["Core", "Periphery", "Healthy"]
class_to_idx = {c: i for i, c in enumerate(class_names)}
y = np.array([class_to_idx[c] for c in labels], dtype=np.int64)

# Train/test split (same as other methods)
split_df = pd.read_csv(SPLIT_PATH, sep="\t")
split_map = dict(zip(split_df["index"], split_df["split"]))
train_mask = np.array([split_map.get(i, "train") == "train" for i in range(len(y))])
test_mask = ~train_mask
# Validation cells carved from TRAIN (10%, stratified, seed 42); epoch chosen on validation, test scored once.
from sklearn.model_selection import train_test_split as _tts, StratifiedGroupKFold as _sgkf
if SPLIT_MODE == "patient":
    # Patient mode: validation grouped by patient, so no donor is in both train and validation.
    _tri = np.where(train_mask)[0]
    _pid = np.load(os.path.join(ROOT, "output", "nn_patient_ids.npy"), allow_pickle=True)[_tri]
    _a, _b = next(_sgkf(n_splits=10, shuffle=True, random_state=42).split(np.zeros(len(_tri)), y[_tri], _pid))
    assert not (set(_pid[_a]) & set(_pid[_b])), "patient overlap between train and validation"
    _tr_i, _va_i = _tri[_a], _tri[_b]
else:
    _tr_i, _va_i = _tts(np.where(train_mask)[0], test_size=0.10, random_state=42, stratify=y[train_mask])
val_mask = np.zeros_like(train_mask); val_mask[_va_i] = True
train_mask = np.zeros_like(train_mask); train_mask[_tr_i] = True

log(f"Data: {len(y)} cells, {latent.shape[1]} features")
log(f"Train: {train_mask.sum()}, Test: {test_mask.sum()}")
log(f"Class dist (train): {np.bincount(y[train_mask])}")
log(f"Class dist (test): {np.bincount(y[test_mask])}")

# ---- 2. PyG Data object ----
x = torch.from_numpy(latent).float()
edge_index_t = torch.from_numpy(edge_index).long()
edge_attr_t = torch.from_numpy(edge_attr).float()
y_t = torch.from_numpy(y).long()

data = Data(x=x, edge_index=edge_index_t, edge_attr=edge_attr_t, y=y_t)
data = data.to(DEVICE)

# Train/test masks
train_mask_t = torch.from_numpy(train_mask).bool().to(DEVICE)
test_mask_t = torch.from_numpy(test_mask).bool().to(DEVICE)
val_mask_t = torch.from_numpy(val_mask).bool().to(DEVICE)

# ---- 3. Edge-aware GAT ----
class EdgeAwareGAT(torch.nn.Module):
    """
    GATv2 with edge-feature bias in attention.
    Attention score: alpha_ij = softmax( LeakyReLU( a^T [W h_i || W h_j] + b_edge(e_ij) ) )
    """
    def __init__(self, in_dim, hidden_dim, out_dim, heads=8, n_layers=2, dropout=0.3, edge_dim=5):
        super().__init__()
        self.n_layers = n_layers
        self.dropout = dropout
        
        self.convs = nn.ModuleList()
        self.batch_norms = nn.ModuleList()
        
        # First layer
        self.convs.append(GATv2Conv(in_dim, hidden_dim, heads=heads, dropout=dropout,
                                    edge_dim=5, add_self_loops=True, concat=True))
        self.batch_norms.append(nn.BatchNorm1d(hidden_dim * heads))
        
        # Hidden layers
        for _ in range(n_layers - 2):
            self.convs.append(GATv2Conv(hidden_dim * heads, hidden_dim, heads=heads, dropout=dropout,
                                        edge_dim=5, add_self_loops=True, concat=True))
            self.batch_norms.append(nn.BatchNorm1d(hidden_dim * heads))
        
        # Output layer (single head, no concat)
        self.convs.append(GATv2Conv(hidden_dim * heads, out_dim, heads=1, dropout=dropout,
                                    edge_dim=5, add_self_loops=True, concat=False))
        
    def forward(self, x, edge_index, edge_attr):
        for i, conv in enumerate(self.convs):
            x = conv(x, edge_index, edge_attr=edge_attr)
            if i < len(self.convs) - 1:
                x = F.elu(x)
                x = self.batch_norms[i](x)
                x = F.dropout(x, p=self.dropout, training=self.training)
        return x

model = EdgeAwareGAT(
    in_dim=latent.shape[1],
    hidden_dim=HIDDEN_DIM,
    out_dim=len(class_names),
    heads=HEADS,
    n_layers=N_LAYERS,
    dropout=DROPOUT,
    edge_dim=5
).to(DEVICE)

log(f"Model params: {sum(p.numel() for p in model.parameters()):,}")

optimizer = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=EPOCHS)
criterion = nn.CrossEntropyLoss()

# ---- 4. Training ----
log("Training...")
train_losses = []
test_losses = []
test_accs = []
best_acc = 0

for epoch in range(1, EPOCHS + 1):
    model.train()
    t0 = time.time()
    
    optimizer.zero_grad()
    out = model(data.x, data.edge_index, data.edge_attr)
    loss = criterion(out[train_mask_t], data.y[train_mask_t])
    loss.backward()
    optimizer.step()
    scheduler.step()
    
    # Validation
    model.eval()
    with torch.no_grad():
        out = model(data.x, data.edge_index, data.edge_attr)
        test_loss = criterion(out[val_mask_t], data.y[val_mask_t])  # validation loss (name kept)
        pred = out[val_mask_t].argmax(dim=-1)
        acc = (pred == data.y[val_mask_t]).float().mean().item()  # validation accuracy
    
    train_losses.append(loss.item())
    test_losses.append(test_loss.item())
    test_accs.append(acc)
    
    if acc > best_acc:
        best_acc = acc
        torch.save(model.state_dict(), MODEL_PATH)
    
    log(f"Epoch {epoch:3d}/{EPOCHS}  TrainLoss={loss.item():.4f}  ValLoss={test_loss.item():.4f}  ValAcc={acc:.4f}  ({time.time()-t0:.1f}s)")

# ---- 5. Final evaluation ----
log("Final evaluation...")
model.load_state_dict(torch.load(MODEL_PATH))
model.eval()
with torch.no_grad():
    out = model(data.x, data.edge_index, data.edge_attr)
    probs = F.softmax(out[test_mask_t], dim=-1).cpu().numpy()
    pred = probs.argmax(axis=1)
    true = data.y[test_mask_t].cpu().numpy()

# Metrics
acc = accuracy_score(true, pred)
np.savez(os.path.join(RES_DIR, f"test_predictions_gat13_seed{_SEED}.npz"), y_true=true, y_pred=pred, test_index=np.where(test_mask)[0])
macro_f1 = f1_score(true, pred, average="macro")
weighted_f1 = f1_score(true, pred, average="weighted")
prec, rec, _, _ = precision_recall_fscore_support(true, pred, average="macro")
y_true_bin = label_binarize(true, classes=[0,1,2])
auc = roc_auc_score(y_true_bin, probs, average="macro", multi_class="ovr")
cm = confusion_matrix(true, pred, labels=[0,1,2])
report = classification_report(true, pred, target_names=class_names, output_dict=True)

metrics = {
    "method": "C-GAT",
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
    "final_test_acc": float(acc),
    "split_mode": SPLIT_MODE,
    "n_params": sum(p.numel() for p in model.parameters())
}

log(f"C-GAT: acc={acc:.4f}, macro_f1={macro_f1:.4f}, AUC={auc:.4f}")

# Save metrics
with open(METRICS_PATH, "w") as f:
    json.dump(metrics, f, indent=2)

# Save predictions
pred_df = pd.DataFrame({
    "true": [class_names[i] for i in true],
    "pred": [class_names[i] for i in pred],
    **{f"prob_{cn}": probs[:, i] for i, cn in enumerate(class_names)}
})
pred_df.to_csv(PRED_PATH, sep="\t", index=False)

# Plot training curves
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].plot(train_losses, label="train")
axes[0].plot(test_losses, label="validation")
axes[0].set_title("Loss"); axes[0].legend()
axes[1].plot(test_accs, label="validation acc", color="green")
axes[1].set_title("Validation Accuracy"); axes[1].legend()
fig.tight_layout()
fig.savefig(PNG_LOSS, dpi=180); plt.close(fig)

# Confusion matrix
fig, ax = plt.subplots(figsize=(5, 4))
sns.heatmap(cm, annot=True, fmt="d", cmap="Blues",
            xticklabels=class_names, yticklabels=class_names, ax=ax)
ax.set_title("C-GAT Confusion Matrix")
ax.set_xlabel("Predicted"); ax.set_ylabel("True")
fig.tight_layout()
fig.savefig(PNG_CM, dpi=180); plt.close(fig)

log(f"Saved {METRICS_PATH}, {PRED_PATH}, {PNG_CM}, {PNG_LOSS}")
log("DONE.")