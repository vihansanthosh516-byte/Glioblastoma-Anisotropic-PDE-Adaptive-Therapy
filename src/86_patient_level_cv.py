"""Script 86: patient-level cross-validation for the Track A zone classifiers (limit L0).

Scripts 04-16 split the 15,000 cells at random, so cells from one patient sit in both train and
test, and class is confounded with patient (Healthy comes from 3 donors; many patients give one
class). This script compares, on identical data:
  random   : StratifiedKFold(5) over cells
  patient  : StratifiedGroupKFold(5), groups = patient (obs["ID1"]), no patient in both sets
for logistic regression and random forest on the 2,500-gene matrix, logistic regression on the
scVI latent, and a "patient lookup" baseline that predicts a patient's majority training class.

Input:  output/nn_X.npy, nn_y.npy, scvi_latent.npy, 02_adata_subsampled.h5ad (obs["ID1"])
Output: output/patient_level_cv.json
"""
import json
import time
from pathlib import Path as _Path

import numpy as np
import scanpy as sc
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score, recall_score
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
OUT_JSON = OUTPUT_DIR / "patient_level_cv.json"
SEED = 42
N_SPLITS = 5
CLASSES = ["Core", "Periphery", "Healthy"]


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def lookup_predict(train_idx, test_idx, y, groups):
    """Predict each test cell's patient majority class from training cells; fall back to global majority."""
    glob = np.bincount(y[train_idx], minlength=3).argmax()
    maj = {}
    for g in np.unique(groups[train_idx]):
        m = train_idx[groups[train_idx] == g]
        maj[g] = np.bincount(y[m], minlength=3).argmax()
    return np.array([maj.get(groups[i], glob) for i in test_idx])


def metrics(y_true, y_pred):
    rec = recall_score(y_true, y_pred, labels=[0, 1, 2], average=None, zero_division=0)
    return {"acc": float(accuracy_score(y_true, y_pred)),
            "macro_f1": float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
            "recall_core": float(rec[0]), "recall_periphery": float(rec[1]), "recall_healthy": float(rec[2])}


def summarise(folds):
    out = {}
    for k in folds[0]:
        v = np.array([f[k] for f in folds])
        out[k] = {"mean": float(v.mean()), "sd": float(v.std(ddof=1)), "per_fold": [float(x) for x in v]}
    return out


def make_models():
    return {
        "logreg_genes": ("X", lambda: LogisticRegression(C=1.0, max_iter=300, class_weight="balanced")),
        "rf_genes": ("X", lambda: RandomForestClassifier(n_estimators=150, n_jobs=-1, random_state=SEED, class_weight="balanced")),
        "logreg_scvi": ("Z", lambda: LogisticRegression(C=1.0, max_iter=1000, class_weight="balanced")),
    }


def main():
    X = np.load(OUTPUT_DIR / "nn_X.npy")
    Z = np.load(OUTPUT_DIR / "scvi_latent.npy")
    y = np.load(OUTPUT_DIR / "nn_y.npy")
    groups = sc.read_h5ad(OUTPUT_DIR / "02_adata_subsampled.h5ad").obs["ID1"].astype(str).to_numpy()
    assert len(groups) == len(y) == X.shape[0] == Z.shape[0]
    log(f"cells {len(y)}, patients {len(np.unique(groups))}, class counts {np.bincount(y).tolist()}")

    splits = {
        "random": list(StratifiedKFold(N_SPLITS, shuffle=True, random_state=SEED).split(X, y)),
        "patient": list(StratifiedGroupKFold(N_SPLITS, shuffle=True, random_state=SEED).split(X, y, groups)),
    }
    for name, sp in splits.items():
        for i, (tr, te) in enumerate(sp):
            assert name == "random" or not (set(groups[tr]) & set(groups[te])), "patient overlap"
            log(f"{name} fold {i}: train {len(tr)}, test {len(te)}, test classes {np.bincount(y[te], minlength=3).tolist()}, "
                f"test patients {len(np.unique(groups[te]))}")

    res = {"script": "86_patient_level_cv", "seed": SEED, "n_splits": N_SPLITS, "n_cells": int(len(y)),
           "n_patients": int(len(np.unique(groups))), "classes": CLASSES,
           "note_rf": "n_estimators=150 (script 06 used 500) to keep runtime manageable",
           "patients_per_class": {c: sorted(set(groups[y == i])) for i, c in enumerate(CLASSES)},
           "results": {}}
    for sname, sp in splits.items():
        res["results"][sname] = {}
        folds = [metrics(y[te], lookup_predict(tr, te, y, groups)) for tr, te in sp]
        res["results"][sname]["patient_lookup_baseline"] = summarise(folds)
        log(f"{sname} patient_lookup acc {res['results'][sname]['patient_lookup_baseline']['acc']['mean']:.4f}")
        for mname, (src, mk) in make_models().items():
            M = X if src == "X" else Z
            folds = []
            for i, (tr, te) in enumerate(sp):
                t0 = time.time()
                clf = mk().fit(M[tr], y[tr])
                folds.append(metrics(y[te], clf.predict(M[te])))
                log(f"{sname} {mname} fold {i}: acc {folds[-1]['acc']:.4f} ({time.time() - t0:.0f}s)")
            res["results"][sname][mname] = summarise(folds)
            OUT_JSON.write_text(json.dumps(res, indent=2))
    OUT_JSON.write_text(json.dumps(res, indent=2))
    log("done")


if __name__ == "__main__":
    main()
