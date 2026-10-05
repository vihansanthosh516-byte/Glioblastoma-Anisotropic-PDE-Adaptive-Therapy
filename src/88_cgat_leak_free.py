"""Script 88: C-GAT re-evaluated without label leakage (limit L0).

Problems found in the original C-GAT evaluation (scripts 12-13):
  1. Edge attributes (scripts 12) hold same_region, transition_type and region_diff, computed from the TRUE
     tissue region of both endpoints of every edge, including test cells. The attention bias can read labels.
  2. The train/test split is a random cell-level split, so patients appear in both sets.
  3. Script 13 keeps the epoch with the best TEST accuracy (model selection on the test set).
  4. The cVAE (script 10) was trained on random-split training cells with class/patient contrastive pairs.

This script uses the same patient-level folds as script 86 and always evaluates the FINAL epoch (100 epochs).
Variants:
  gat_pca_clean        node features: PCA(32) fitted on training cells only; graph: kNN(15) on those features;
                       edge attributes: [distance, same_patient] only. No label-derived input anywhere.
  gat_cvae_no_label_edges  cVAE latent + original graph, edge attributes [distance, same_patient] only.
                       (cVAE residual leak remains: trained on labelled random-split cells.)
  cgat_original_edges  cVAE latent + original graph + all 5 edge attributes (the published configuration),
                       evaluated under patient folds, final epoch. Shows the size of the label leak.

Input:  output/nn_X.npy, nn_y.npy, cgat/cvae_latent.npy, cgat/gat_edge_index.npy, cgat/gat_edge_attr.npy,
        02_adata_subsampled.h5ad
Output: output/cgat_leak_free.json
"""
import json
import time
from pathlib import Path as _Path

import numpy as np
import scanpy as sc
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.decomposition import PCA
from sklearn.metrics import accuracy_score, f1_score, recall_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.neighbors import NearestNeighbors
from torch_geometric.nn import GATv2Conv

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
OUT_JSON = OUTPUT_DIR / "cgat_leak_free.json"
SEED, N_SPLITS, EPOCHS, K = 42, 5, 100, 15
HIDDEN, HEADS, DROPOUT, LR, WD = 64, 8, 0.3, 1e-3, 1e-4
DEVICE = "cpu"


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


class EdgeAwareGAT(nn.Module):
    """Same architecture as script 13 (2 GATv2 layers, batch norm, ELU), with configurable edge_dim."""

    def __init__(self, in_dim, out_dim, edge_dim):
        super().__init__()
        self.c1 = GATv2Conv(in_dim, HIDDEN, heads=HEADS, dropout=DROPOUT, edge_dim=edge_dim, add_self_loops=True, concat=True)
        self.bn = nn.BatchNorm1d(HIDDEN * HEADS)
        self.c2 = GATv2Conv(HIDDEN * HEADS, out_dim, heads=1, dropout=DROPOUT, edge_dim=edge_dim, add_self_loops=True, concat=False)

    def forward(self, x, ei, ea):
        x = self.bn(F.elu(self.c1(x, ei, edge_attr=ea)))
        x = F.dropout(x, p=DROPOUT, training=self.training)
        return self.c2(x, ei, edge_attr=ea)


def train_eval(x, ei, ea, y, tr, te, seed):
    torch.manual_seed(seed)
    model = EdgeAwareGAT(x.shape[1], 3, ea.shape[1]).to(DEVICE)
    opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WD)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=EPOCHS)
    trm, tem = torch.from_numpy(tr), torch.from_numpy(te)
    for _ in range(EPOCHS):
        model.train()
        opt.zero_grad()
        F.cross_entropy(model(x, ei, ea)[trm], y[trm]).backward()
        opt.step()
        sched.step()
    model.eval()
    with torch.no_grad():
        pred = model(x, ei, ea)[tem].argmax(-1).numpy()
    return pred


def metrics(yt, yp):
    rec = recall_score(yt, yp, labels=[0, 1, 2], average=None, zero_division=0)
    return {"acc": float(accuracy_score(yt, yp)), "macro_f1": float(f1_score(yt, yp, average="macro", zero_division=0)),
            "recall_core": float(rec[0]), "recall_periphery": float(rec[1]), "recall_healthy": float(rec[2])}


def knn_edges(feat, patients, k=K):
    nb = NearestNeighbors(n_neighbors=k + 1, n_jobs=4).fit(feat)
    d, idx = nb.kneighbors(feat)
    rows = np.repeat(np.arange(len(feat)), k)
    cols = idx[:, 1:].reshape(-1)
    dist = d[:, 1:].reshape(-1)
    same = (patients[rows] == patients[cols]).astype(np.float32)
    return np.stack([rows, cols]).astype(np.int64), np.stack([dist, same], axis=1).astype(np.float32)


def rebuild_original_graph(lat, patients, regions):
    """Rebuild script 12's graph exactly (the .npy edge files are git-ignored and absent).
    Edge attrs: [distance, same_patient, same_region, transition_type, region_diff]."""
    order = {"NormalBrain": 0, "peri:GBM": 1, "core:GBM": 2}
    ridx = np.array([order[r] for r in regions])
    nb = NearestNeighbors(n_neighbors=K + 1, metric="euclidean", n_jobs=4).fit(lat)
    d, idx = nb.kneighbors(lat)
    rows = np.repeat(np.arange(len(lat)), K)
    cols = idx[:, 1:].reshape(-1)
    dist = d[:, 1:].reshape(-1)
    pair_type = {("core:GBM", "peri:GBM"): 1, ("NormalBrain", "peri:GBM"): 2, ("NormalBrain", "core:GBM"): 3}
    trans = np.array([0 if regions[i] == regions[j] else pair_type.get(tuple(sorted([regions[i], regions[j]])), 4)
                      for i, j in zip(rows, cols)], dtype=np.float32)
    ea = np.stack([dist, (patients[rows] == patients[cols]).astype(np.float32),
                   (regions[rows] == regions[cols]).astype(np.float32), trans,
                   np.abs(ridx[rows] - ridx[cols]).astype(np.float32)], axis=1).astype(np.float32)
    return np.stack([rows, cols]).astype(np.int64), ea


def summarise(folds):
    out = {}
    for k in folds[0]:
        v = np.array([f[k] for f in folds])
        out[k] = {"mean": float(v.mean()), "sd": float(v.std(ddof=1)), "per_fold": [float(a) for a in v]}
    return out


def main():
    Xg = np.load(OUTPUT_DIR / "nn_X.npy")
    y_np = np.load(OUTPUT_DIR / "nn_y.npy")
    lat = np.load(OUTPUT_DIR / "cgat" / "cvae_latent.npy").astype(np.float32)
    obs = sc.read_h5ad(OUTPUT_DIR / "02_adata_subsampled.h5ad").obs
    patients = obs["ID1"].astype(str).to_numpy()
    ei0, ea0 = rebuild_original_graph(lat, patients, obs["tissue_histology"].astype(str).to_numpy())
    y = torch.from_numpy(y_np).long()
    folds = list(StratifiedGroupKFold(N_SPLITS, shuffle=True, random_state=SEED).split(Xg, y_np, patients))
    res = {"script": "88_cgat_leak_free", "seed": SEED, "epochs": EPOCHS, "n_splits": N_SPLITS,
           "evaluated_epoch": "final (no test-set selection)", "variants": {}}
    ea_cols = ea0[:, :2]  # distance, same_patient
    per = {"gat_pca_clean": [], "gat_cvae_no_label_edges": [], "cgat_original_edges": []}
    for i, (tr, te) in enumerate(folds):
        t0 = time.time()
        pca = PCA(32, random_state=SEED).fit(Xg[tr])
        Zp = pca.transform(Xg).astype(np.float32)
        eip, eap = knn_edges(Zp, patients)
        pred = train_eval(torch.from_numpy(Zp), torch.from_numpy(eip), torch.from_numpy(eap), y, tr, te, SEED)
        per["gat_pca_clean"].append(metrics(y_np[te], pred))
        pred = train_eval(torch.from_numpy(lat), torch.from_numpy(ei0), torch.from_numpy(ea_cols), y, tr, te, SEED)
        per["gat_cvae_no_label_edges"].append(metrics(y_np[te], pred))
        pred = train_eval(torch.from_numpy(lat), torch.from_numpy(ei0), torch.from_numpy(ea0), y, tr, te, SEED)
        per["cgat_original_edges"].append(metrics(y_np[te], pred))
        log(f"fold {i}: " + ", ".join(f"{k} {v[-1]['acc']:.4f}" for k, v in per.items()) + f" ({time.time() - t0:.0f}s)")
        res["variants"] = {k: summarise(v) for k, v in per.items()}
        OUT_JSON.write_text(json.dumps(res, indent=2))
    log("done")


if __name__ == "__main__":
    main()
