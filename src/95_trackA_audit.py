#!/usr/bin/env python3
"""Script 95: Track A leakage and confounding audit (Phase 1.4; masterplan §29-30, §66, §128).

Questions, answered on the same 15,000 cells and the same patient folds as script 88
(StratifiedGroupKFold(5, shuffle, seed 42), groups = ID1):
  A  Who contributes each class? (donor / diagnosis / Fresh-Frozen composition)
  B  What is in each held-out fold?
  C  How far do non-biological features go? Majority, Fresh/Frozen only, sequencing depth only,
     cells-per-patient only, diagnosis only, and PCA(32) + logistic regression as reference.
     Everything is fitted on training patients only.
  D  Within-patient label permutation (keeps each patient's class mix, destroys cell-level signal).
  E  Can the representations predict patient or Fresh/Frozen? (random cell split, so this measures
     identifiability, not generalisation)
  F  Graph audit: the original kNN graph over all cells. Share of edges that join test to train
     nodes, same-patient edges, same-class edges, and a neighbour-vote baseline that uses only
     graph edges and training labels.

Input : output/nn_X.npy, output/nn_y.npy, output/cgat/cvae_latent.npy, output/02_adata_subsampled.h5ad
Output: output/trackA_audit.json, output/trackA_audit.manifest.json
"""
from __future__ import annotations

import json
import time
from pathlib import Path as _Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score, recall_score
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
from sklearn.preprocessing import OneHotEncoder, StandardScaler

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output"
SEED, N_SPLITS, K_PERM = 42, 5, 100
CLASSES = ["Core", "Periphery", "Healthy"]

import sys  # noqa: E402
sys.path.insert(0, str(PROJECT_ROOT))
from src.run_manifest import write_run_manifest  # noqa: E402

K = 15


def rebuild_original_graph(lat, patients, regions):
    """Copy of script 88 rebuild_original_graph (script 12's graph; the .npy edge files are git-ignored).
    Copied, not imported, because script 88 needs torch_geometric at import time.
    Edge attrs: [distance, same_patient, same_region, transition_type, region_diff]."""
    from sklearn.neighbors import NearestNeighbors
    order = {"NormalBrain": 0, "peri:GBM": 1, "core:GBM": 2}
    ridx = np.array([order[r] for r in regions])
    nb = NearestNeighbors(n_neighbors=K + 1, metric="euclidean", n_jobs=4).fit(lat)
    d, idx = nb.kneighbors(lat)
    rows = np.repeat(np.arange(len(lat)), K)
    cols = idx[:, 1:].reshape(-1)
    dist = d[:, 1:].reshape(-1)
    ea = np.stack([dist, (patients[rows] == patients[cols]).astype(np.float32),
                   (regions[rows] == regions[cols]).astype(np.float32),
                   np.zeros(len(rows), np.float32), np.abs(ridx[rows] - ridx[cols]).astype(np.float32)], axis=1)
    return np.stack([rows, cols]).astype(np.int64), ea.astype(np.float32)


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def metrics(yt, yp):
    rec = recall_score(yt, yp, labels=[0, 1, 2], average=None, zero_division=0)
    return {"acc": float(accuracy_score(yt, yp)), "macro_f1": float(f1_score(yt, yp, average="macro", zero_division=0)),
            "recall_core": float(rec[0]), "recall_periphery": float(rec[1]), "recall_healthy": float(rec[2])}


def summarise(folds):
    out = {}
    for k in folds[0]:
        v = np.array([f[k] for f in folds])
        out[k] = {"mean": float(v.mean()), "sd": float(v.std(ddof=1)), "per_fold": [float(a) for a in v]}
    return out


def main():
    t0 = time.time()
    X = np.load(OUT / "nn_X.npy")
    y = np.load(OUT / "nn_y.npy")
    lat = np.load(OUT / "cgat" / "cvae_latent.npy").astype(np.float32)
    obs = sc.read_h5ad(OUT / "02_adata_subsampled.h5ad").obs
    pid = obs["ID1"].astype(str).to_numpy()
    cls = obs["class"].astype(str).to_numpy()
    assert [CLASSES[i] for i in y[:50]] == list(cls[:50]), "nn_y order does not match obs['class']"
    folds = list(StratifiedGroupKFold(N_SPLITS, shuffle=True, random_state=SEED).split(X, y, pid))
    res = {"script": "95_trackA_audit", "seed": SEED, "n_cells": int(len(y)), "n_patients": int(len(set(pid)))}

    # A composition
    g = obs.assign(pid=pid).groupby("pid", observed=True).agg(
        n_cells=("class", "size"), classes=("class", lambda s: ",".join(sorted(set(s)))),
        diagnosis=("Diagnosis", "first"), source=("source", "first"))
    res["A_composition"] = {
        "patients": g.reset_index().to_dict(orient="records"),
        "class_by_diagnosis": pd.crosstab(obs["class"], obs["Diagnosis"]).to_dict(),
        "class_by_source": pd.crosstab(obs["class"], obs["source"]).to_dict(),
        "single_class_patients": int((~g["classes"].str.contains(",")).sum()),
        "healthy_donors": int((g["classes"] == "Healthy").sum()),
        "healthy_donors_diagnoses": g.loc[g["classes"] == "Healthy", "diagnosis"].tolist(),
        "frozen_cells_that_are_healthy": int(((obs["source"] == "Frozen") & (obs["class"] == "Healthy")).sum()),
    }

    # B fold composition
    res["B_folds"] = [{"fold": i, "test_patients": sorted(set(pid[te])), "n_test_cells": int(len(te)),
                       "test_class_counts": {c: int((cls[te] == c).sum()) for c in CLASSES},
                       "train_healthy_donors": int(len(set(pid[tr][cls[tr] == "Healthy"]))),
                       "test_healthy_donors": int(len(set(pid[te][cls[te] == "Healthy"])))}
                      for i, (tr, te) in enumerate(folds)]

    # C baselines (train patients only)
    n_per_patient = pd.Series(pid).map(pd.Series(pid).value_counts()).to_numpy().astype(float)
    qc = np.log1p(obs[["nCount_RNA", "nFeature_RNA"]].to_numpy(dtype=float))
    src = obs[["source"]].astype(str).to_numpy()
    dx = obs[["Diagnosis"]].astype(str).to_numpy()
    base = {k: [] for k in ["majority", "fresh_frozen_only", "sequencing_depth_only", "cells_per_patient_only",
                            "diagnosis_only", "pca32_logreg_reference"]}
    pca_cache = []
    for i, (tr, te) in enumerate(folds):
        maj = np.bincount(y[tr]).argmax()
        base["majority"].append(metrics(y[te], np.full(len(te), maj)))
        enc = OneHotEncoder(handle_unknown="ignore").fit(src[tr])
        m = LogisticRegression(max_iter=1000).fit(enc.transform(src[tr]), y[tr])
        base["fresh_frozen_only"].append(metrics(y[te], m.predict(enc.transform(src[te]))))
        m = HistGradientBoostingClassifier(random_state=SEED, max_iter=100).fit(qc[tr], y[tr])
        base["sequencing_depth_only"].append(metrics(y[te], m.predict(qc[te])))
        m = HistGradientBoostingClassifier(random_state=SEED, max_iter=100).fit(n_per_patient[tr, None], y[tr])
        base["cells_per_patient_only"].append(metrics(y[te], m.predict(n_per_patient[te, None])))
        enc = OneHotEncoder(handle_unknown="ignore").fit(dx[tr])
        m = LogisticRegression(max_iter=1000).fit(enc.transform(dx[tr]), y[tr])
        base["diagnosis_only"].append(metrics(y[te], m.predict(enc.transform(dx[te]))))
        sc_ = StandardScaler().fit(X[tr])
        pca = PCA(32, random_state=SEED).fit(sc_.transform(X[tr]))
        Ztr, Zte = pca.transform(sc_.transform(X[tr])), pca.transform(sc_.transform(X[te]))
        pca_cache.append((Ztr, Zte))
        m = LogisticRegression(max_iter=2000).fit(Ztr, y[tr])
        base["pca32_logreg_reference"].append(metrics(y[te], m.predict(Zte)))
        log(f"C fold {i} done ({time.time() - t0:.0f}s)")
    res["C_baselines_patient_folds"] = {k: summarise(v) for k, v in base.items()}

    # D within-patient permutation null
    rng = np.random.default_rng(SEED)
    null_acc = []
    for r in range(K_PERM):
        yp = y.copy()
        for p in np.unique(pid):
            idx = np.where(pid == p)[0]
            yp[idx] = y[rng.permutation(idx)]
        accs = []
        for i, (tr, te) in enumerate(folds):
            Ztr, Zte = pca_cache[i]
            m = LogisticRegression(max_iter=500).fit(Ztr, yp[tr])
            accs.append(accuracy_score(yp[te], m.predict(Zte)))
        null_acc.append(float(np.mean(accs)))
        if r % 20 == 0:
            log(f"D perm {r}/{K_PERM} ({time.time() - t0:.0f}s)")
    real = res["C_baselines_patient_folds"]["pca32_logreg_reference"]["acc"]["mean"]
    res["D_within_patient_permutation"] = {
        "n_perm": K_PERM, "null_mean_acc": float(np.mean(null_acc)),
        "null_p95": float(np.percentile(null_acc, 95)), "null_max": float(np.max(null_acc)),
        "real_pca32_logreg_acc": real,
        "real_above_null_max": bool(real > np.max(null_acc)),
        "note": "Null keeps each patient's class mix and destroys cell-level signal. 21 patients, so it is coarse."}

    # E identifiability of patient / Fresh-Frozen from representations (random cell split)
    def ident(Z, target, name):
        accs = []
        for tr, te in StratifiedKFold(N_SPLITS, shuffle=True, random_state=SEED).split(Z, target):
            m = LogisticRegression(max_iter=1000).fit(StandardScaler().fit(Z[tr]).transform(Z[tr]), target[tr])
            accs.append(accuracy_score(target[te], m.predict(StandardScaler().fit(Z[tr]).transform(Z[te]))))
        chance = float(np.bincount(pd.factorize(target)[0]).max() / len(target))
        return {"representation": name, "acc": float(np.mean(accs)), "majority_chance": chance}

    pca_all = PCA(32, random_state=SEED).fit_transform(StandardScaler().fit_transform(X))
    pcodes = pd.factorize(pid)[0]
    scodes = pd.factorize(obs["source"].astype(str))[0]
    res["E_identifiability"] = {
        "patient_from_cvae_latent": ident(lat, pcodes, "cvae_latent"),
        "patient_from_pca32": ident(pca_all, pcodes, "pca32"),
        "fresh_frozen_from_cvae_latent": ident(lat, scodes, "cvae_latent"),
        "fresh_frozen_from_pca32": ident(pca_all, scodes, "pca32")}
    log(f"E done ({time.time() - t0:.0f}s)")

    # F graph audit on the published (script 12) graph, rebuilt as in script 88
    regions = obs["tissue_histology"].astype(str).to_numpy()
    ei, ea = rebuild_original_graph(lat, pid, regions)
    rows, cols = ei
    same_pat = ea[:, 1] == 1
    same_cls = y[rows] == y[cols]
    gf = []
    for i, (tr, te) in enumerate(folds):
        is_te = np.zeros(len(y), bool)
        is_te[te] = True
        from_te = is_te[rows]
        to_tr = ~is_te[cols]
        e = from_te & to_tr
        # neighbour vote from TRAIN neighbours only, using no features
        votes = np.zeros((len(y), 3))
        np.add.at(votes, (rows[e], y[cols[e]]), 1)
        has = votes[te].sum(1) > 0
        pred = votes[te].argmax(1)
        gf.append({"fold": i,
                   "share_edges_from_test_nodes_to_train_nodes": float(e.sum() / max(from_te.sum(), 1)),
                   "test_nodes_with_any_train_neighbour": float(has.mean()),
                   "neighbour_vote_acc_on_nodes_with_train_neighbours": float(accuracy_score(y[te][has], pred[has])) if has.any() else None,
                   "neighbour_vote_acc_all_test_nodes_no_neighbour_counts_wrong": float((pred[has] == y[te][has]).sum() / len(te))})
    res["F_graph_audit"] = {
        "n_edges": int(len(rows)), "share_edges_same_patient": float(same_pat.mean()),
        "share_edges_same_class": float(same_cls.mean()),
        "share_edges_same_class_given_diff_patient": float(same_cls[~same_pat].mean()),
        "per_fold": gf,
        "mean_neighbour_vote_acc_on_covered_nodes": float(np.mean([f["neighbour_vote_acc_on_nodes_with_train_neighbours"] for f in gf if f["neighbour_vote_acc_on_nodes_with_train_neighbours"] is not None])),
        "reading": ("The graph is built over ALL cells with the cVAE latent. Test cells have train-labelled neighbours "
                    "(first column). The cVAE (script 10) was trained on labelled random-split cells, so its latent "
                    "space already encodes the labels of test patients' cells. Only gat_pca_clean in script 88 is free "
                    "of both leaks.")}

    ok = OUT / "trackA_audit.json"
    ok.write_text(json.dumps(res, indent=1, default=float))
    write_run_manifest("trackA_audit_95", OUT / "trackA_audit.manifest.json", script="src/95_trackA_audit.py", seed=SEED,
                       config={"n_splits": N_SPLITS, "k_perm": K_PERM, "pca_dim": 32},
                       inputs=[OUT / "nn_X.npy", OUT / "nn_y.npy", OUT / "cgat" / "cvae_latent.npy",
                               OUT / "02_adata_subsampled.h5ad"],
                       dataset="multiomic-gbm 15,000-cell subsample", patient_split="StratifiedGroupKFold(5), groups=ID1",
                       primary_endpoint="n/a (audit)")
    log(f"done {time.time() - t0:.0f}s -> {ok}")


if __name__ == "__main__":
    main()
