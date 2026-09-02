"""
08_method3_hybrid.py
====================
Method 3 (Hybrid): Combine classical + deep learning.

Idea: Use Logistic Regression coefficient magnitudes as a *prior bias* for 
the Transformer's self-attention mechanism. Genes with higher |coef| in LR
get boosted attention scores.

Architecture:
  - Same Transformer encoder as Method 2
  - Attention bias matrix: for each head, add (scaled) LR coefficient magnitude
    as a static bias to the attention logits before softmax
  - This injects the "statistical knowledge" from Method 1 into Method 2
  - Final classification head same as before
"""
import os, time, resource, json, gc
os.environ["OMP_NUM_THREADS"] = "2"
os.environ["MKL_NUM_THREADS"] = "2"

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
from sklearn.metrics import (accuracy_score, f1_score, precision_recall_fscore_support,
                             confusion_matrix, classification_report, roc_auc_score)
from sklearn.preprocessing import label_binarize

ROOT    = "/mnt/c/Users/vihan/20206 science fair"
OUT_DIR = os.path.join(ROOT, "output")
NPY_X   = os.path.join(OUT_DIR, "nn_X.npy")
NPY_Y   = os.path.join(OUT_DIR, "nn_y.npy")
TSV_GENE_NM = os.path.join(OUT_DIR, "nn_gene_names.tsv")
TSV_CLS_NM  = os.path.join(OUT_DIR, "nn_class_names.tsv")
TSV_SPLIT   = os.path.join(OUT_DIR, "nn_train_test_split.tsv")

# Load LR coefficients from Method 1
TSV_LR_COEF = os.path.join(OUT_DIR, "method1_lr_coefficients.tsv")

PT_MODEL     = os.path.join(OUT_DIR, "method3_hybrid.pt")
TSV_METRICS  = os.path.join(OUT_DIR, "method3_metrics.json")
TSV_PRED     = os.path.join(OUT_DIR, "method3_predictions.tsv")
PNG_CM       = os.path.join(OUT_DIR, "method3_confusion.png")
PNG_LOSS     = os.path.join(OUT_DIR, "method3_training_loss.png")

# Hyperparameters - CPU friendly (match Method 2)
N_GENES      = 100
N_CLASSES    = 3
EMBED_DIM    = 64
N_HEADS      = 4
N_LAYERS     = 2
FF_DIM       = 128
DROPOUT      = 0.2
BATCH_SIZE   = 256
LR           = 1e-3
EPOCHS       = 5
WEIGHT_DECAY = 1e-4
DEVICE       = "cuda" if torch.cuda.is_available() else "cpu"

# Hybrid-specific: scale factor for LR coefficient bias
LR_BIAS_SCALE = 0.5  # how strongly to inject LR prior

def mem_mb(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}  RSS={mem_mb():.0f}MB  DEV={DEVICE}", flush=True)

# ---- 1. Load data + LR coefficients ----------------------------------------
log("Loading dense tensors + LR coefficients...")
X = np.load(NPY_X)[:, :N_GENES]
y = np.load(NPY_Y)
gene_names = pd.read_csv(TSV_GENE_NM, sep="\t")["gene"].tolist()[:N_GENES]
class_df   = pd.read_csv(TSV_CLS_NM, sep="\t")
class_names = class_df["class_name"].tolist()
split_df   = pd.read_csv(TSV_SPLIT, sep="\t")
split_map = dict(zip(split_df["index"], split_df["split"]))

# Load LR coefficients (3 classes x 2500 genes) -> subset to N_GENES
lr_coef_df = pd.read_csv(TSV_LR_COEF, sep="\t", index_col=0)
lr_coef_df = lr_coef_df[gene_names]  # reorder + subset
# Take mean absolute coefficient across the 3 classes as "importance prior"
lr_prior = lr_coef_df.abs().mean(axis=0).values  # (N_GENES,)
# Normalize to [0, 1] range
lr_prior = (lr_prior - lr_prior.min()) / (lr_prior.max() - lr_prior.min() + 1e-8)
lr_prior_t = torch.from_numpy(lr_prior).float().to(DEVICE)  # (G,)
log(f"LR prior loaded: shape={lr_prior_t.shape}, range=[{lr_prior_t.min():.3f}, {lr_prior_t.max():.3f}]")

train_mask = np.array([split_map[i] == "train" for i in range(len(y))])
test_mask  = ~train_mask
X_train, X_test = X[train_mask], X[test_mask]
y_train, y_test = y[train_mask], y[test_mask]
log(f"train {X_train.shape}, test {X_test.shape}")

# ---- 2. PyTorch tensors / DataLoaders --------------------------------------
X_train_t = torch.from_numpy(X_train).float()
y_train_t = torch.from_numpy(y_train).long()
X_test_t  = torch.from_numpy(X_test).float()
y_test_t  = torch.from_numpy(y_test).long()

train_ds = TensorDataset(X_train_t, y_train_t)
test_ds  = TensorDataset(X_test_t,  y_test_t)
train_dl = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True,  num_workers=0)
test_dl  = DataLoader(test_ds,  batch_size=BATCH_SIZE, shuffle=False, num_workers=0)

# ---- 3. Hybrid Transformer with LR-biased attention ------------------------
class HybridGeneTransformer(nn.Module):
    def __init__(self, n_genes, n_classes, embed_dim, n_heads, n_layers, ff_dim, dropout, lr_prior):
        super().__init__()
        self.n_genes = n_genes
        self.embed_dim = embed_dim
        self.n_heads = n_heads
        self.register_buffer("lr_prior", lr_prior)  # (G,)
        
        self.gene_emb = nn.Embedding(n_genes, embed_dim)
        self.val_proj = nn.Linear(1, embed_dim)
        self.cls_token = nn.Parameter(torch.randn(1, 1, embed_dim))
        self.pos_emb = nn.Parameter(torch.randn(1, n_genes + 1, embed_dim))
        
        # Custom encoder layers that accept attention bias
        self.layers = nn.ModuleList([
            nn.TransformerEncoderLayer(
                d_model=embed_dim, nhead=n_heads, dim_feedforward=ff_dim,
                dropout=dropout, batch_first=True, activation="gelu"
            ) for _ in range(n_layers)
        ])
        self.dropout = nn.Dropout(dropout)
        self.norm = nn.LayerNorm(embed_dim)
        self.head = nn.Linear(embed_dim, n_classes)
        nn.init.normal_(self.gene_emb.weight, std=0.02)
        nn.init.normal_(self.cls_token, std=0.02)
        nn.init.normal_(self.pos_emb, std=0.02)

    def forward(self, x):
        B, G = x.shape
        gene_ids = torch.arange(G, device=x.device).unsqueeze(0)
        gene_e = self.gene_emb(gene_ids)
        val_e = self.val_proj(x.unsqueeze(-1))
        tokens = gene_e + val_e
        cls = self.cls_token.expand(B, -1, -1)
        tokens = torch.cat([cls, tokens], dim=1)
        tokens = tokens + self.pos_emb
        tokens = self.dropout(tokens)
        
        # Create attention bias from LR prior: shape (G+1, G+1)
        # cls token has no prior (0), genes have lr_prior
        # We want: bias for cls->gene = lr_prior, gene->cls = 0, gene->gene = lr_prior
        # But standard MHA expects (B*n_heads, T, T) or (T, T)
        # We'll add to the attention logits inside the encoder layers
        
        # Build bias matrix: (1, 1, T, T) where T = G+1
        T = G + 1
        bias = torch.zeros(1, 1, T, T, device=x.device)
        # cls attends to genes with LR prior
        bias[:, :, 0, 1:] = self.lr_prior.unsqueeze(0) * LR_BIAS_SCALE
        # genes attend to cls with no bias (or small)
        # genes attend to other genes with LR prior of target gene
        bias[:, :, 1:, 1:] = self.lr_prior.unsqueeze(0).unsqueeze(0) * LR_BIAS_SCALE
        
        for layer in self.layers:
            # We need to manually call self_attn with attn_mask
            # TransformerEncoderLayer doesn't easily expose this, so we manually iterate
            # Standard path: layer.self_attn(q, k, v, attn_mask=bias)
            tokens = self._encoder_layer_with_bias(layer, tokens, bias)
        
        cls_out = tokens[:, 0]
        cls_out = self.norm(cls_out)
        logits = self.head(cls_out)
        return logits

    def _encoder_layer_with_bias(self, layer, x, bias):
        # LayerNorm -> SelfAttn -> Residual -> LayerNorm -> FFN -> Residual
        residual = x
        x = layer.norm1(x)
        # Self-attention with bias
        # layer.self_attn is nn.MultiheadAttention
        # We need to pass attn_mask of shape (T, T) or (B*n_heads, T, T)
        # Our bias is (1, 1, T, T) - expand to (B, n_heads, T, T)
        B, T, E = x.shape
        attn_mask = bias.expand(B, self.n_heads, T, T).reshape(B * self.n_heads, T, T)
        x, _ = layer.self_attn(x, x, x, attn_mask=attn_mask, need_weights=False)
        x = layer.dropout1(x)
        x = residual + x
        
        residual = x
        x = layer.norm2(x)
        x = layer.linear2(layer.dropout(F.gelu(layer.linear1(x))))
        x = layer.dropout2(x)
        x = residual + x
        return x

model = HybridGeneTransformer(
    n_genes=N_GENES, n_classes=N_CLASSES,
    embed_dim=EMBED_DIM, n_heads=N_HEADS,
    n_layers=N_LAYERS, ff_dim=FF_DIM, dropout=DROPOUT,
    lr_prior=lr_prior_t
).to(DEVICE)

n_params = sum(p.numel() for p in model.parameters())
log(f"Hybrid Model: {n_params:,} params (with LR prior injection)")

criterion = nn.CrossEntropyLoss()
optimizer = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=EPOCHS)

# ---- 4. Training loop -------------------------------------------------------
log("Training Hybrid...")
train_losses, test_losses, test_accs = [], [], []
for epoch in range(1, EPOCHS + 1):
    model.train()
    t0 = time.time()
    tr_loss_acc = 0.0
    for xb, yb in train_dl:
        xb, yb = xb.to(DEVICE), yb.to(DEVICE)
        optimizer.zero_grad()
        out = model(xb)
        loss = criterion(out, yb)
        loss.backward()
        optimizer.step()
        tr_loss_acc += loss.item() * xb.size(0)
    tr_loss = tr_loss_acc / len(train_ds)
    train_losses.append(tr_loss)

    model.eval()
    with torch.no_grad():
        te_loss_acc = 0.0
        correct = 0
        for xb, yb in test_dl:
            xb, yb = xb.to(DEVICE), yb.to(DEVICE)
            out = model(xb)
            loss = criterion(out, yb)
            te_loss_acc += loss.item() * xb.size(0)
            correct += (out.argmax(-1) == yb).sum().item()
    te_loss = te_loss_acc / len(test_ds)
    te_acc = correct / len(test_ds)
    test_losses.append(te_loss)
    test_accs.append(te_acc)

    scheduler.step()
    log(f"epoch {epoch:02d}/{EPOCHS}  trLoss={tr_loss:.4f}  teLoss={te_loss:.4f}  teAcc={te_acc:.4f}  ({time.time()-t0:.1f}s)")

# ---- 5. Evaluation ----------------------------------------------------------
log("Evaluating on test set...")
model.eval()
all_preds, all_true, all_proba = [], [], []
with torch.no_grad():
    for xb, yb in test_dl:
        xb = xb.to(DEVICE)
        out = model(xb)
        proba = F.softmax(out, dim=-1).cpu().numpy()
        pred = out.argmax(-1).cpu().numpy()
        all_proba.append(proba)
        all_preds.append(pred)
        all_true.append(yb.numpy())
y_pred = np.concatenate(all_preds)
y_true = np.concatenate(all_true)
y_proba = np.concatenate(all_proba)

# Metrics
acc = accuracy_score(y_true, y_pred)
macro_f1 = f1_score(y_true, y_pred, average="macro")
weighted_f1 = f1_score(y_true, y_pred, average="weighted")
prec, rec, _, _ = precision_recall_fscore_support(y_true, y_pred, average="macro")
y_true_bin = label_binarize(y_true, classes=[0,1,2])
auc = roc_auc_score(y_true_bin, y_proba, average="macro", multi_class="ovr")
cm = confusion_matrix(y_true, y_pred, labels=[0,1,2])
report = classification_report(y_true, y_pred, target_names=class_names, output_dict=True)

metrics = {
    "method": "Hybrid_Transformer_LR_Prior",
    "accuracy": float(acc),
    "macro_f1": float(macro_f1),
    "weighted_f1": float(weighted_f1),
    "macro_precision": float(prec),
    "macro_recall": float(rec),
    "macro_auc_ovr": float(auc),
    "confusion_matrix": cm.tolist(),
    "per_class": report,
    "n_params": n_params,
    "epochs": EPOCHS,
    "best_test_acc": max(test_accs),
    "final_test_acc": test_accs[-1],
    "lr_bias_scale": LR_BIAS_SCALE
}
log(f"Hybrid: acc={acc:.4f}, macro_f1={macro_f1:.4f}, AUC={auc:.4f}")

# ---- 6. Save artifacts ------------------------------------------------------
with open(TSV_METRICS, "w") as f:
    json.dump(metrics, f, indent=2)
log(f"Saved {TSV_METRICS}")

torch.save(model.state_dict(), PT_MODEL)
log(f"Saved model -> {PT_MODEL}")

pred_df = pd.DataFrame({
    "true": [class_names[i] for i in y_true],
    "pred": [class_names[i] for i in y_pred],
    **{f"prob_{cn}": y_proba[:,i] for i, cn in enumerate(class_names)}
})
pred_df.to_csv(TSV_PRED, sep="\t", index=False)
log(f"Saved {TSV_PRED}")

# Plot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].plot(train_losses, label="train")
axes[0].plot(test_losses, label="test")
axes[0].set_title("Loss"); axes[0].legend()
axes[1].plot(test_accs, label="test acc", color="green")
axes[1].set_title("Test Accuracy"); axes[1].legend()
fig.tight_layout(); fig.savefig(PNG_LOSS, dpi=180); plt.close(fig)

fig, ax = plt.subplots(figsize=(5, 4))
sns.heatmap(cm, annot=True, fmt="d", xticklabels=class_names,
            yticklabels=class_names, cmap="Blues", ax=ax)
ax.set_title("Hybrid (LR-prior Transformer)"); ax.set_xlabel("Predicted"); ax.set_ylabel("True")
fig.tight_layout(); fig.savefig(PNG_CM, dpi=180); plt.close(fig)
log(f"Saved {PNG_CM}, {PNG_LOSS}")

log("DONE.")