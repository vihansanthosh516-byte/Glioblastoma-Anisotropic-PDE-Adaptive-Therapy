"""Track B v8: two-part collapse model on v5 features. Exploratory, try 8.

Amendment: paper/trackB_mu_v8_amendment_2026-11-14.md (locked before scoring).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression, QuantileRegressor
from sklearn.metrics import roc_auc_score

import track_b_mu_tcks as t
import track_b_mu_v3 as v3
import track_b_mu_v5 as v5

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "output" / "track_b_mu_v8_20261114"
TRY = 8


def run():
    d, _ = t.load()
    base_pred = pd.read_csv(v3.V2)
    assert (base_pred.patient_id.values == d.patient_id.values).all()
    comp = v5.features_table(d).drop(columns=["patient_id", "tp_in"]).reset_index(drop=True)
    collapse = (d.core_vox_out < 0.05 * d.core_vox_in).values
    pred, prob = np.full(len(d), np.nan), np.full(len(d), np.nan)
    for f in sorted(d.fold.unique()):
        trm, tem = (d.fold != f).values, (d.fold == f).values
        tr, te = d[trm], d[tem]
        Xtr = pd.concat([v3.design(tr, v3.latent(tr, tr)).reset_index(drop=True), comp[trm].reset_index(drop=True)], axis=1)
        Xte = pd.concat([v3.design(te, v3.latent(te, tr)).reset_index(drop=True), comp[tem].reset_index(drop=True)], axis=1)
        mu, sd = Xtr.mean(), Xtr.std().replace(0, 1.0)
        Ztr, Zte = (Xtr - mu) / sd, (Xte - mu) / sd
        c_tr = collapse[trm]
        clf = LogisticRegression(max_iter=1000).fit(Ztr, c_tr)
        p = clf.predict_proba(Zte)[:, 1]
        q = QuantileRegressor(quantile=0.5, alpha=0.0, solver="highs").fit(Ztr[~c_tr], tr.y_out.values[~c_tr])
        y_col = float(np.median(tr.y_out.values[c_tr]))
        pred[tem] = np.where(p > 0.5, y_col, q.predict(Zte))
        prob[tem] = p
    y = d.y_out.values
    errs = {k: np.abs(base_pred[k].values - y) for k in v3.BASE}
    errs["v8"] = np.abs(pred - y)
    pat = pd.DataFrame(errs).assign(pid=d.patient_id.values).groupby("pid").mean()
    mae = pat.mean().to_dict()
    best = min(v3.BASE, key=mae.get)
    level = 1 - 0.05 / TRY
    diff = (pat.v8 - pat[best]).values
    ci = v3.boot_ci(diff, np.random.default_rng(t.SEED), level)
    grow = (d.core_vox_out > 1.25 * d.core_vox_in).astype(int).values
    res = {
        "script": "src/track_b_mu_v8.py", "amendment": "paper/trackB_mu_v8_amendment_2026-11-14.md",
        "label": "exploratory, post hoc, try 8 on reused MU-Glioma-Post",
        "n_pairs": int(len(d)), "n_patients": int(len(pat)), "n_collapse": int(collapse.sum()),
        "n_predicted_collapse": int((prob > 0.5).sum()),
        "collapse_auc": float(roc_auc_score(collapse, prob)),
        "patient_mae_log_volume": mae, "best_baseline": best,
        "v8_minus_best": {"mean": float(diff.mean()), "ci_level": level, "ci": ci},
        "success": bool(ci[1] < 0),
        "secondary_auc_growth25_v8": float(roc_auc_score(grow, pred - d.y_in.values)),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "result.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    pd.DataFrame({"patient_id": d.patient_id, "fold": d.fold, "y_out": y, "v8": pred, "p_collapse": prob}).to_csv(OUT / "predictions.csv", index=False)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    run()
