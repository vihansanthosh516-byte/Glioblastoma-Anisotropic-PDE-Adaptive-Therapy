"""
20_patient_split_methods.py
===========================
Donor-disjoint (patient-split) substitute methods. See
research_notes/plan_patient_split_methods_2026-10-19.md (amendments apply).

Methods (substitutes, not the published tools):
  M1p  PCA latent (50 comps, fit on TRAIN cells only) + LogisticRegression / RandomForest
  M3p  Sparse gene panel: L1-logistic selection on TRAIN cells only, top-k genes, + LogisticRegression

Protocol (fixed in the plan):
  - Train = 16 patients (patient split from script 04). Test = 5 patients, scored once per method.
  - Validation = first fold of StratifiedGroupKFold(10) inside train (3 patients). Hyperparameters chosen on validation macro-F1.
  - Refit on all train patients with the chosen hyperparameters; score test once.
  - Seeds 42-46 (used by RandomForest and L1 solver ordering); mean/min/max reported.
  - Patient-level bootstrap (2000 resamples over the 5 test patients) of test macro-F1.

Output: output/patient_split_v2_methods/ (one JSON per method, plus summary.json).
Run: python src/20_patient_split_methods.py
"""
import os, json
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import accuracy_score, f1_score, balanced_accuracy_score
from sklearn.preprocessing import StandardScaler

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "output")
RES = os.path.join(OUT, "patient_split_v2_methods")
os.makedirs(RES, exist_ok=True)
SEEDS = [42, 43, 44, 45, 46]
N_BOOT = 2000

X = np.load(os.path.join(OUT, "nn_X.npy"), mmap_mode="r")
y = np.load(os.path.join(OUT, "nn_y.npy"))
pid = np.load(os.path.join(OUT, "nn_patient_ids.npy"), allow_pickle=True)
sp = pd.read_csv(os.path.join(OUT, "nn_train_test_split_patient.tsv"), sep="\t")
train_rows = sp["index"][sp["split"] == "train"].to_numpy()
test_rows = sp["index"][sp["split"] == "test"].to_numpy()
Xtr, ytr, ptr = np.asarray(X[train_rows]), y[train_rows], pid[train_rows]
Xte, yte, pte = np.asarray(X[test_rows]), y[test_rows], pid[test_rows]
assert not set(ptr) & set(pte), "patient overlap"

# grouped validation inside train (first fold, same as scripts 05/07/08)
tr_i, va_i = next(StratifiedGroupKFold(n_splits=10, shuffle=True, random_state=42)
                  .split(np.zeros(len(ytr)), ytr, ptr))
assert not set(ptr[tr_i]) & set(ptr[va_i])


def macro_f1(yt, yp):
    return f1_score(yt, yp, average="macro")


def fit_predict_pca(Xa, ya, Xb, C, seed, model):
    sc = StandardScaler().fit(Xa)
    pca = PCA(n_components=50, random_state=seed).fit(sc.transform(Xa))
    Za, Zb = pca.transform(sc.transform(Xa)), pca.transform(sc.transform(Xb))
    if model == "lr":
        clf = LogisticRegression(C=C, max_iter=2000)
    else:
        clf = RandomForestClassifier(n_estimators=300, random_state=seed, n_jobs=-1)
    clf.fit(Za, ya)
    return clf.predict(Zb)


def fit_predict_panel(Xa, ya, Xb, k, C, seed):
    sc = StandardScaler().fit(Xa)
    Xa_s, Xb_s = sc.transform(Xa), sc.transform(Xb)
    sel = LogisticRegression(penalty="l1", C=0.05, solver="saga", max_iter=3000, random_state=seed)
    sel.fit(Xa_s, ya)
    importance = np.abs(sel.coef_).sum(axis=0)
    genes = np.argsort(importance)[::-1][:k]
    clf = LogisticRegression(C=C, max_iter=2000)
    clf.fit(Xa_s[:, genes], ya)
    return clf.predict(Xb_s[:, genes])


def run_method(name, grid, fn):
    # select hyperparameters on grouped validation (seed 42), then refit + test once per seed
    best, best_f1 = None, -1.0
    val_table = []
    for hp in grid:
        pred = fn(Xtr[tr_i], ytr[tr_i], Xtr[va_i], hp, 42)
        f1v = macro_f1(ytr[va_i], pred)
        val_table.append({"hp": hp, "val_macro_f1": f1v})
        if f1v > best_f1:
            best, best_f1 = hp, f1v
    per_seed = []
    for s in SEEDS:
        pred = fn(Xtr, ytr, Xte, best, s)
        np.savez(os.path.join(RES, f"pred_seed{s}_{name}.npz"), y_true=yte, y_pred=pred, test_index=test_rows)
        per_seed.append({"seed": s, "accuracy": accuracy_score(yte, pred),
                         "macro_f1": macro_f1(yte, pred),
                         "balanced_accuracy": balanced_accuracy_score(yte, pred)})
    f1s = [r["macro_f1"] for r in per_seed]
    # patient-level bootstrap on the seed-42 prediction (fixed refit, resample test patients)
    pred42 = fn(Xtr, ytr, Xte, best, 42)
    test_patients = np.unique(pte)
    rng = np.random.default_rng(42)
    boots = []
    for _ in range(N_BOOT):
        pick = rng.choice(test_patients, size=len(test_patients), replace=True)
        idx = np.concatenate([np.where(pte == p)[0] for p in pick])
        boots.append(macro_f1(yte[idx], pred42[idx]))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    out = {"method": name, "selected_hp": best, "validation": val_table,
           "per_seed": per_seed,
           "macro_f1_mean": float(np.mean(f1s)), "macro_f1_min": float(np.min(f1s)),
           "macro_f1_max": float(np.max(f1s)),
           "bootstrap_macro_f1_95ci_seed42": [float(lo), float(hi)],
           "notes": "substitute method (see plan amendments); single test refit per seed; validation = 3 patients"}
    with open(os.path.join(RES, f"{name}.json"), "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"{name}: selected={best} test macro-F1 mean={np.mean(f1s):.3f} "
          f"[{np.min(f1s):.3f},{np.max(f1s):.3f}] CI42=({lo:.3f},{hi:.3f})", flush=True)
    return out


def lr_pca(Xa, ya, Xb, C, seed):
    return fit_predict_pca(Xa, ya, Xb, C, seed, "lr")


def rf_pca(Xa, ya, Xb, _hp, seed):
    return fit_predict_pca(Xa, ya, Xb, None, seed, "rf")


def panel(Xa, ya, Xb, hp, seed):
    k, C = hp
    return fit_predict_panel(Xa, ya, Xb, k, C, seed)


results = [
    run_method("M1p_pca_lr", [0.1, 1.0, 10.0], lr_pca),
    run_method("M1p_pca_rf", [None], rf_pca),
    run_method("M3p_panel_lr", [(50, 1.0), (200, 1.0)], panel),
]
with open(os.path.join(RES, "summary.json"), "w") as f:
    json.dump([{k: r[k] for k in ("method", "selected_hp", "macro_f1_mean", "macro_f1_min",
                                   "macro_f1_max", "bootstrap_macro_f1_95ci_seed42")}
               for r in results], f, indent=2, default=float)
print("DONE", RES)
