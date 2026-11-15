"""Track B v5: input-scan tumor-composition latent model. Exploratory, try 5.

Amendment: paper/trackB_mu_v5_amendment_2026-11-14.md (locked before scoring).
Reads only the INPUT-timepoint tumor mask for features. Baselines read from the v2 run.
"""
from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from sklearn.linear_model import QuantileRegressor
from sklearn.metrics import roc_auc_score

import track_b_mu_tcks as t
import track_b_mu_v3 as v3

ROOT = Path(__file__).resolve().parent.parent
MU = ROOT / "data" / "tcia" / "MU-Glioma-Post"
OUT = ROOT / "output" / "track_b_mu_v5_20261114"
CACHE = OUT / "input_mask_features.csv"
TRY = 5


def mask_features(pid, tp):
    a = np.asarray(nib.load(MU / pid / f"Timepoint_{tp}" / f"{pid}_Timepoint_{tp}_tumorMask.nii.gz").dataobj)
    n1, n2, n3 = (int((a == k).sum()) for k in (1, 2, 3))
    core = n1 + n3
    box = 1
    if core:
        idx = np.argwhere((a == 1) | (a == 3))
        box = int(np.prod(idx.max(0) - idx.min(0) + 1))
    return {"log_edema": np.log(n2 + 1.0), "nec_frac": n1 / core if core else 0.0,
            "enh_frac": n3 / core if core else 0.0, "log_whole": np.log(n1 + n2 + n3 + 1.0),
            "compact": core / box if core else 0.0}


def features_table(d):
    if CACHE.exists():
        f = pd.read_csv(CACHE)
    else:
        rows = [{"patient_id": p, "tp_in": int(tp), **mask_features(p, int(tp))}
                for p, tp in d[["patient_id", "tp_in"]].drop_duplicates().itertuples(index=False)]
        f = pd.DataFrame(rows)
        OUT.mkdir(parents=True, exist_ok=True)
        f.to_csv(CACHE, index=False)
    return d[["patient_id", "tp_in"]].merge(f, on=["patient_id", "tp_in"], how="left")


def run():
    d, _ = t.load()
    base_pred = pd.read_csv(v3.V2)
    assert (base_pred.patient_id.values == d.patient_id.values).all()
    comp = features_table(d).drop(columns=["patient_id", "tp_in"]).reset_index(drop=True)
    pred = np.full(len(d), np.nan)
    for f in sorted(d.fold.unique()):
        trm, tem = (d.fold != f).values, (d.fold == f).values
        tr, te = d[trm], d[tem]
        Xtr = pd.concat([v3.design(tr, v3.latent(tr, tr)).reset_index(drop=True), comp[trm].reset_index(drop=True)], axis=1)
        Xte = pd.concat([v3.design(te, v3.latent(te, tr)).reset_index(drop=True), comp[tem].reset_index(drop=True)], axis=1)
        mu, sd = Xtr.mean(), Xtr.std().replace(0, 1.0)
        q = QuantileRegressor(quantile=0.5, alpha=0.0, solver="highs").fit((Xtr - mu) / sd, tr.y_out.values)
        pred[tem] = q.predict((Xte - mu) / sd)
    y = d.y_out.values
    errs = {k: np.abs(base_pred[k].values - y) for k in v3.BASE}
    errs["v5"] = np.abs(pred - y)
    pat = pd.DataFrame(errs).assign(pid=d.patient_id.values).groupby("pid").mean()
    mae = pat.mean().to_dict()
    best = min(v3.BASE, key=mae.get)
    level = 1 - 0.05 / TRY
    diff = (pat.v5 - pat[best]).values
    ci = v3.boot_ci(diff, np.random.default_rng(t.SEED), level)
    grow = (d.core_vox_out > 1.25 * d.core_vox_in).astype(int).values
    res = {
        "script": "src/track_b_mu_v5.py", "amendment": "paper/trackB_mu_v5_amendment_2026-11-14.md",
        "label": "exploratory, post hoc, try 5 on reused MU-Glioma-Post",
        "n_pairs": int(len(d)), "n_patients": int(len(pat)),
        "patient_mae_log_volume": mae, "best_baseline": best,
        "v5_minus_best": {"mean": float(diff.mean()), "ci_level": level, "ci": ci},
        "success": bool(ci[1] < 0),
        "secondary_auc_growth25_v5": float(roc_auc_score(grow, pred - d.y_in.values)),
    }
    (OUT / "result.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    pd.DataFrame({"patient_id": d.patient_id, "fold": d.fold, "y_out": y, "v5": pred}).to_csv(OUT / "predictions.csv", index=False)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    run()
