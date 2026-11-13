"""Method search on the TCGA development set (spent). Logs every try.

Protocol: paper/integration_protocol_methoddev_2026-11-12.md. Reads only TCGA features and TCGA RNA.
Score: repeated (10 x 5-fold) cross-validated Spearman rho; all tuning inside each training fold.

  python -I src/integration_methoddev_search.py
"""
from __future__ import annotations

import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import integration_bamd_gate1 as g1  # noqa: E402
from integration_bamd_features import feature_names  # noqa: E402

OUT = PROJECT_ROOT / "output" / "integration_methoddev"
VERHAAK = PROJECT_ROOT / "data" / "external" / "neftel2019" / "TCGA.Verhaak.GBM.txt"
SEED = 20261113


def verhaak_sets(path=VERHAAK) -> dict:
    t = pd.read_csv(path, sep="\t")
    return {c: tuple(t[c].dropna().astype(str)) for c in ("Neural", "Proneural", "Mesenchymal", "Classical")}


def targets(expr: pd.DataFrame) -> pd.DataFrame:
    S = g1.state_scores(expr, g1.load_modules())
    V = g1.state_scores(expr, verhaak_sets())
    T = pd.concat([S.add_prefix("neftel_"), V.add_prefix("verhaak_")], axis=1)
    T["neftel_MES_minus_NPC"] = T["neftel_MES"] - T["neftel_NPC"]
    T["neftel_AC_plus_MES"] = (T["neftel_AC"] + T["neftel_MES"]) / 2 - (T["neftel_OPC"] + T["neftel_NPC"]) / 2
    return T


def make_model(name):
    from sklearn.cross_decomposition import PLSRegression
    from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
    from sklearn.linear_model import ElasticNetCV, RidgeCV
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    alphas = np.logspace(-1, 4, 12)
    return {"ridge": lambda: make_pipeline(StandardScaler(), RidgeCV(alphas=alphas)),
            "enet": lambda: make_pipeline(StandardScaler(), ElasticNetCV(l1_ratio=[.2, .5, .9], cv=5, max_iter=5000)),
            "pls2": lambda: make_pipeline(StandardScaler(), PLSRegression(2)),
            "rf": lambda: RandomForestRegressor(300, min_samples_leaf=5, max_features=0.3, n_jobs=6, random_state=0),
            "et": lambda: ExtraTreesRegressor(300, min_samples_leaf=5, max_features=0.3, n_jobs=6, random_state=0),
            }[name]()


def cv_score(X, y, model, rng, repeats=10, k=5):
    rhos = []
    for _ in range(repeats):
        pred = np.empty(len(y))
        for test in g1.folds(len(y), k, rng):
            train = np.setdiff1d(np.arange(len(y)), test)
            Xtr, Xte = g1.impute(X[train], X[test])
            m = make_model(model).fit(Xtr, y[train])
            pred[test] = np.ravel(m.predict(Xte))
        rhos.append(spearmanr(pred, y)[0])
    return float(np.mean(rhos)), float(np.std(rhos))


def load():
    old = pd.read_csv(PROJECT_ROOT / "output" / "integration_bamd" / "features.csv")
    old = old[old["collection"] == "TCGA-GBM"]
    new = pd.read_csv(OUT / "features_TCGA-GBM.csv")
    with gzip.open(g1.EXPR, "rt") as fh:
        expr = pd.read_csv(fh, sep="\t", index_col=0)
    samples = g1.pick_samples(expr.columns)
    d = old.merge(new.drop(columns=["collection", "patient_id"]), on="brats_id", validate="one_to_one")
    d["patient_id"] = d["patient_id"].astype(str)
    d = d[d["patient_id"].isin(samples)].drop_duplicates("patient_id").sort_values("patient_id")
    T = targets(expr[[samples[p] for p in d["patient_id"]]])
    newcols = [c for c in new.columns if c not in ("collection", "patient_id", "brats_id")]
    return d, T, feature_names(), newcols


def main() -> int:
    d, T, oldcols, newcols = load()
    sets = {"old59": oldcols, "new": newcols, "both": oldcols + newcols,
            "location": [c for c in newcols if c.startswith("loc_")],
            "texture": [c for c in newcols if c.startswith("glcm_")],
            "intensity_norm": [c for c in newcols if c.startswith(("nmean_", "nsd_", "ring_", "hab_"))]}
    rng = np.random.Generator(np.random.PCG64DXSM(SEED))
    log = []
    for tname in T.columns:
        y = T[tname].to_numpy(float)
        for sname, cols in sets.items():
            X = d[cols].to_numpy(float)
            for mname in ("ridge", "enet", "pls2", "rf", "et"):
                mean, sd = cv_score(X, y, mname, rng)
                log.append({"target": tname, "features": sname, "n_features": len(cols), "model": mname,
                            "cv_rho_mean": mean, "cv_rho_sd": sd})
                print(f"{tname:24s} {sname:15s} {mname:6s} {mean:+.3f} ({sd:.3f})", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    res = {"n": len(d), "tries": len(log), "log": log, "seed": SEED,
           "note": "development only (TCGA spent); best of many tries is optimistic (selection bias)",
           "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (OUT / "log.json").write_text(json.dumps(res, indent=2))
    best = sorted(log, key=lambda r: -r["cv_rho_mean"])[:10]
    print(json.dumps(best, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
