"""
05_attention_gated_network.py
=============================
A PyTorch implementation of an Attention-Gated Neural Network that:
  - Classifies a 2,500-dim gene-expression vector into the 3 classes
    (Core / Periphery / Healthy).
  - Uses a *gene-level attention head* that learns which genes matter
    for each prediction. The attention weights are saved per-cell,
    so we can compare them against the paper's marker genes.

Architecture (attention-gated):
  x  (N, 2500)         ->  gene embedding (N, 2500, D)
  gene_attn            ->  scores over 2500 genes  (logistic / softmax)
  gated_x  = a * v     ->  weighted per-gene features (N, D)
  dropout + 2 dense    ->  logits (N, 3)

Intent (original author notes; NOT demonstrated by this script, 2026-10-19 audit):
  - Attention weights are saved per cell, but nothing here shows they explain
    predictions (no ablation or comparison code).
  - No comparison against a 10,000-gene model is run here.
  - No test of agreement with DE markers beyond listing top weights.
  - Test accuracy is read every epoch and the best epoch is reported, on a
    random cell-level split (not patient-held-out), so it is optimistic.

Size: output/nn_metrics.json reports n_params = 488899 (the earlier "~3.2M"
was wrong). No torch.save is called, so nn_model.pt is not produced here.

Run:
    cd "."
    python src/05_attention_gated_network.py
"""
import os, time, json, gc
try:
    import resource
except ImportError:  # native Windows Python has no resource module; RSS then reports 0
    resource = None
os.environ["OMP_NUM_THREADS"]      = "2"
os.environ["MKL_NUM_THREADS"]      = "2"

import numpy as np
import pandas as pd
import torch, torch.nn as nn, torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
# GLOBAL_SEED_2026_10_19: fixed seed added by the audit corrections (legacy script was unseeded); outputs from earlier runs came from unseeded runs.
import random as _random_seed
_random_seed.seed(42)
try:
    import numpy as _np_seed
    _np_seed.random.seed(42)
except ImportError:
    pass
try:
    import torch as _torch_seed
    _torch_seed.manual_seed(42)
except ImportError:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # project root (was "."; cwd-dependent)
OUT_DIR = os.path.join(ROOT, "output")
# GBM_SPLIT_MODE=patient: use the patient-disjoint split from script 04 and write all results under output/patient_split/
# (the default "cell" mode is the original leaky random cell split and keeps the original output names).
SPLIT_MODE = os.environ.get("GBM_SPLIT_MODE", "cell")
RES_DIR = OUT_DIR if SPLIT_MODE == "cell" else os.path.join(OUT_DIR, "patient_split")
os.makedirs(RES_DIR, exist_ok=True)
NPY_X   = os.path.join(OUT_DIR, "nn_X.npy")
NPY_Y   = os.path.join(OUT_DIR, "nn_y.npy")
TSV_GENE_NM  = os.path.join(OUT_DIR, "nn_gene_names.tsv")
TSV_CLS_NM   = os.path.join(OUT_DIR, "nn_class_names.tsv")
TSV_SPLIT    = os.path.join(OUT_DIR, "nn_train_test_split.tsv" if SPLIT_MODE == "cell" else "nn_train_test_split_patient.tsv")

PT_MODEL     = os.path.join(RES_DIR, "nn_model.pt")
TSV_ATTEN    = os.path.join(RES_DIR, "nn_attention_weights.tsv")
TSV_CLASSIF_REQ  = os.path.join(RES_DIR, "nn_classification_report.txt")
TSV_METRICS  = os.path.join(RES_DIR, "nn_metrics.json")
PNG_LOSS     = os.path.join(RES_DIR, "nn_training_loss.png")
PNG_CONFUSN  = os.path.join(RES_DIR, "nn_confusion_matrix.png")
TOP_N_DE     = 100
EMBED_DIM    = 64
DROPOUT      = 0.3
HIDDEN       = 128
BATCH        = 128
LR           = 1e-3
EPOCHS       = 25
DEVICE       = "cuda" if torch.cuda.is_available() else "cpu"

def mem_mb(): return (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss if resource else 0)/1024
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}  RSS={mem_mb():.0f}MB  DEV={DEVICE}", flush=True)

# ---- 1. Load -------------------------------------------------------------
log("Loading tensors...")
X = np.load(NPY_X); y = np.load(NPY_Y)
gene_names = pd.read_csv(TSV_GENE_NM, sep="\t")["gene"].to_list()
class_df   = pd.read_csv(TSV_CLS_NM,  sep="\t")
class_names = class_df["class_name"].tolist()
split_df   = pd.read_csv(TSV_SPLIT, sep="\t")
split_map  = dict(zip(split_df["index"], split_df["split"]))
log(f"X {X.shape}, y {y.shape}; class_names={class_names}")
log(f"train n={sum(v=='train' for v in split_map.values())}, test n={sum(v=='test' for v in split_map.values())}")

# ---- 2. PyTorch tensors / split ------------------------------------------
train_mask = np.array([split_map[i] == "train" for i in range(len(y))])
test_mask  = ~train_mask
X_train_t  = torch.from_numpy(X[train_mask]).float()
y_train_t  = torch.from_numpy(y[train_mask]).long()
X_test_t   = torch.from_numpy(X[test_mask]).float()
y_test_t   = torch.from_numpy(y[test_mask]).long()
log(f"X_train_t {X_train_t.shape}, X_test_t {X_test_t.shape}")
# Validation split carved from TRAIN (10%, stratified, seed 42): the epoch is chosen on validation; the test set is scored once.
from sklearn.model_selection import train_test_split as _tts, StratifiedGroupKFold as _sgkf
if SPLIT_MODE == "patient":
    # Patient mode: validation is grouped by patient, so no donor is in both train and validation.
    _pid = np.load(os.path.join(OUT_DIR, "nn_patient_ids.npy"), allow_pickle=True)[train_mask]
    _tr_i, _va_i = next(_sgkf(n_splits=10, shuffle=True, random_state=42).split(np.zeros(len(y_train_t)), y_train_t.numpy(), _pid))
    assert not (set(_pid[_tr_i]) & set(_pid[_va_i])), "patient overlap between train and validation"
    print(f"patient-grouped validation: {len(set(_pid[_va_i]))} patients, {len(_va_i)} cells", flush=True)
else:
    _tr_i, _va_i = _tts(np.arange(len(y_train_t)), test_size=0.10, random_state=42, stratify=y_train_t.numpy())
X_val_t, y_val_t = X_train_t[_va_i], y_train_t[_va_i]
X_train_t, y_train_t = X_train_t[_tr_i], y_train_t[_tr_i]

# ---- 3. Model: attention-gated ---------------------------------------------
class AttnGated(nn.Module):
    """
    Concept:
      Each gene gets a small embedding (D).
      The "gate" produces a scalar per gene from that embedding.
      The "value" is the embedding itself. We multiply value by gate.
    """
    def __init__(self, n_genes, n_classes, emb=EMBED_DIM, h=HIDDEN, p=DROPOUT):
        super().__init__()
        # shared gene-embedding (treating each gene index as a token id 0..n_genes-1)
        self.gene_emb = nn.Embedding(n_genes, emb)
        # attention head per gene
        self.gate    = nn.Linear(emb, 1, bias=False)
        self.norm    = nn.LayerNorm(emb)
        self.dropout = nn.Dropout(p)
        # input to MLP = pooled (emb) + raw x (n_genes)
        self.fc1     = nn.Linear(emb + n_genes, h)
        self.fc2     = nn.Linear(h, n_classes)
        # init
        nn.init.normal_(self.gene_emb.weight, std=0.02)
        nn.init.normal_(self.gate.weight, std=0.02)

    def forward(self, x):
        # x: (B, G) where x[i,j] = z-scored expression of gene j for cell i
        B, G = x.shape
        ids   = torch.arange(G, device=x.device)
        emb   = self.gene_emb(ids)               # (G, D)
        emb   = self.norm(emb)
        gate  = self.gate(emb).squeeze(-1)        # (G,)
        # mask: ignore genes with very low |z|-score in this cell (paper-style noise filter)
        mask  = (x.abs() > 0.1).float()
        # per-cell attention logits: a_g += gate_g + cell-specific score (linear add)
        # we make cell-conditioned attention by adding x @ W_cell, where W_cell is a vector
        # simple option: score_ij = gate_g + log(mask_ij + eps)  (truly attention-gated)
        score = gate.unsqueeze(0) + torch.log(mask + 1e-6)   # (B, G)
        attn  = F.softmax(score, dim=-1)                     # (B, G)
        # weighted sum:
        # treat each gene's contribution as emb[g] * x[i,g] (so we let real values flow)
        # then average by attention weight
        v = x.unsqueeze(-1) * emb.unsqueeze(0)              # (B, G, D)
        gated = attn.unsqueeze(-1) * v                      # (B, G, D)
        pooled = gated.sum(dim=1) / (attn.sum(dim=1, keepdim=True) + 1e-6)  # (B, D)
        # concat with global stats to give the MLP richer features
        feat = torch.cat([pooled, x], dim=-1)                # (B, D+G)
        feat = self.dropout(feat)
        h    = F.relu(self.fc1(feat))
        h    = self.dropout(h)
        return self.fc2(h), attn

    def attention(self, x):
        with torch.no_grad():
            B, G = x.shape
            ids  = torch.arange(G, device=x.device)
            emb  = self.gene_emb(ids)
            emb  = self.norm(emb)
            gate = self.gate(emb).squeeze(-1)
            mask = (x.abs() > 0.1).float()
            score = gate.unsqueeze(0) + torch.log(mask + 1e-6)
            return F.softmax(score, dim=-1)

model = AttnGated(n_genes=X.shape[1], n_classes=len(class_names)).to(DEVICE)
opt   = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=1e-5)
crit  = nn.CrossEntropyLoss()
nparams = sum(p.numel() for p in model.parameters())
log(f"model built, {nparams:,} params")

train_ds = TensorDataset(X_train_t, y_train_t)
val_ds   = TensorDataset(X_val_t, y_val_t)
test_ds  = TensorDataset(X_test_t,  y_test_t)
train_dl = DataLoader(train_ds, batch_size=BATCH, shuffle=True,  num_workers=0)
test_dl  = DataLoader(test_ds,  batch_size=BATCH, shuffle=False, num_workers=0)
val_dl   = DataLoader(val_ds,   batch_size=BATCH, shuffle=False, num_workers=0)

# ---- 4. Train -------------------------------------------------------------
log("Begin training...")
import copy as _copy
best_val, best_epoch, best_state = -1.0, 0, None
tr_loss, te_loss, te_acc = [], [], []  # NOTE: te_* now hold VALIDATION values per epoch
for epoch in range(1, EPOCHS + 1):
    model.train(); t0=time.time(); loss_acc=0; n_acc=0
    for xb, yb in train_dl:
        xb, yb = xb.to(DEVICE), yb.to(DEVICE)
        opt.zero_grad()
        out, _ = model(xb)
        loss = crit(out, yb)
        loss.backward()
        opt.step()
        loss_acc += loss.item() * xb.size(0)
        n_acc    += xb.size(0)
    tr_loss.append(loss_acc / max(n_acc,1))

    model.eval()
    with torch.no_grad():
        loss_acc=0; n_acc=0; correct=0
        for xb, yb in val_dl:
            xb, yb = xb.to(DEVICE), yb.to(DEVICE)
            out, _ = model(xb)
            loss = crit(out, yb)
            loss_acc += loss.item() * xb.size(0)
            n_acc    += xb.size(0)
            correct  += (out.argmax(-1) == yb).sum().item()
    te_loss.append(loss_acc / max(n_acc,1))
    te_acc.append(correct / max(n_acc,1))
    if te_acc[-1] > best_val:
        best_val, best_epoch, best_state = te_acc[-1], epoch, _copy.deepcopy(model.state_dict())
    log(f"epoch {epoch:02d}/{EPOCHS} trLoss={tr_loss[-1]:.4f} valLoss={te_loss[-1]:.4f} valAcc={te_acc[-1]:.4f} ({time.time()-t0:.1f}s)")

model.load_state_dict(best_state)  # epoch selected on validation
log(f"selected epoch {best_epoch} (val acc {best_val:.4f})")

# ---- 5. Plot training -----------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(tr_loss, label="train", color="tab:blue")
ax[0].plot(te_loss, label="validation",  color="tab:orange")
ax[0].set_title("Loss")
ax[0].set_xlabel("epoch"); ax[0].legend()
ax[1].plot(te_acc, color="tab:green")
ax[1].set_title("Validation accuracy")
ax[1].set_xlabel("epoch")
fig.tight_layout()
fig.savefig(PNG_LOSS, dpi=180); plt.close(fig)

# ---- 6. Confusion matrix --------------------------------------------------
from sklearn.metrics import confusion_matrix, classification_report
model.eval()
preds=[]; true=[]
with torch.no_grad():
    for xb, yb in test_dl:
        out, _ = model(xb.to(DEVICE))
        preds.append(out.argmax(-1).cpu().numpy())
        true .append(yb.numpy())
preds = np.concatenate(preds); true = np.concatenate(true)

cm = confusion_matrix(true, preds, labels=list(range(len(class_names))))
import seaborn as sns
fig, ax = plt.subplots(figsize=(5,4))
sns.heatmap(cm, annot=True, fmt="d", xticklabels=class_names,
            yticklabels=class_names, cmap="Blues", ax=ax)
ax.set_xlabel("Predicted"); ax.set_ylabel("True")
fig.tight_layout(); fig.savefig(PNG_CONFUSN, dpi=180); plt.close(fig)
log(f"Wrote {PNG_CONFUSN}")

report = classification_report(true, preds, target_names=class_names)
with open(TSV_CLASSIF_REQ, "w") as f:
    f.write(report)
log(f"Wrote {TSV_CLASSIF_REQ}\n{report}")

metrics = {
    "best_val_accuracy": best_val,
    "selected_epoch": best_epoch,
    "selection": "epoch chosen on a validation split carved from train; test scored once at that epoch",
    "final_epoch_val_accuracy": te_acc[-1],
    "final_train_loss":   tr_loss[-1],
    "final_epoch_val_loss": te_loss[-1],
    "n_params":            nparams,
    "epochs":              EPOCHS,
    "batch_size":          BATCH,
    "embedding_dim":       EMBED_DIM,
    "hidden_dim":          HIDDEN,
    "dropout":             DROPOUT,
    "lr":                  LR,
    "device":              DEVICE,
}
with open(TSV_METRICS, "w") as f:
    json.dump(metrics, f, indent=2)

# ---- 7. Extract attention weights and save -------------------------------
log("Extracting attention weights...")
model.eval()
with torch.no_grad():
    A_test = model.attention(X_test_t.to(DEVICE)).cpu().numpy()  # (N_test, G)
A_test_df = pd.DataFrame(A_test, columns=gene_names)
A_test_df.insert(0, "true_class",      [class_names[i] for i in true])
A_test_df.insert(1, "predicted_class", [class_names[i] for i in preds])
A_test_df.insert(2, "correct",         (true == preds).astype(int))
A_test_df.to_csv(TSV_ATTEN, sep="\t", index=False)
log(f"Wrote {TSV_ATTEN}: {A_test_df.shape}")

# ---- 8. Identify top attention-weighted genes (averaged per class) --------
log("Computing class-averaged attention...")
A_full = model.attention(
    torch.from_numpy(X).float().to(DEVICE)
).cpu().numpy()  # (N, G)
y_all = y
records = []
for ci, name in enumerate(class_names):
    sel = (y_all == ci)
    mean_a = A_full[sel].mean(axis=0)
    top_idx = np.argsort(mean_a)[::-1][:TOP_N_DE]
    for rank, gi in enumerate(top_idx, 1):
        records.append({"class": name, "rank": rank,
                        "gene": gene_names[gi],
                        "mean_attention": float(mean_a[gi])})
attn_top = pd.DataFrame(records)
attn_top.to_csv(os.path.join(RES_DIR, "nn_attention_top100_per_class.tsv"),
                sep="\t", index=False)
log(f"Wrote top-100 attention per class")
log("DONE.")
