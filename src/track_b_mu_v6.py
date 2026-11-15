"""Track B v6: composition + edge shape latent model. Exploratory, try 6.

Amendment: paper/trackB_mu_v6_amendment_2026-11-14.md (locked before scoring).
Reads only INPUT-timepoint masks. Reuses v5 features and model; baselines from the v2 run.
"""
from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage
from sklearn.linear_model import QuantileRegressor
from sklearn.metrics import roc_auc_score

import track_b_mu_tcks as t
import track_b_mu_v3 as v3
import track_b_mu_v5 as v5

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "output" / "track_b_mu_v6_20261114"
CACHE = OUT / "input_shape_features.csv"
TRY = 6


def shape_features(pid, tp):
    a = np.asarray(nib.load(v5.MU / pid / f"Timepoint_{tp}" / f"{pid}_Timepoint_{tp}_tumorMask.nii.gz").dataobj)
    core = (a == 1) | (a == 3)
    v = int(core.sum())
    if v == 0:
        return {"surf_ratio": 0.0, "sphericity": 0.0, "edema_reach": 0.0}
    area = int((core & ~ndimage.binary_erosion(core)).sum())
    return {"surf_ratio": area / v,
            "sphericity": float(np.pi ** (1 / 3) * (6 * v) ** (2 / 3) / max(area, 1)),
            "edema_reach": float(np.log((int((a == 2).sum()) + v) / v))}


def shape_table(d):
    if CACHE.exists():
        f = pd.read_csv(CACHE)
    else:
        rows = [{"patient_id": p, "tp_in": int(tp), **shape_features(p, int(tp))}
                for p, tp in d[["patient_id", "tp_in"]].drop_duplicates().itertuples(index=False)]
        f = pd.DataFrame(rows)
        OUT.mkdir(parents=True, exist_ok=True)
        f.to_csv(CACHE, index=False)
    return d[["patient_id", "tp_in"]].merge(f, on=["patient_id", "tp_in"], how="left")


def run():
    d, _ = t.load()
    base_pred = pd.read_csv(v3.V2)
    assert (base_pred.patient_id.values == d.patient_id.values).all()
    comp = pd.concat([v5.features_table(d).drop(columns=["patient_id", "tp_in"]).reset_index(drop=True),
                      shape_table(d).drop(columns=["patient_id", "tp_in"]).reset_index(drop=True)], axis=1)
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
    errs["v6"] = np.abs(pred - y)
    pat = pd.DataFrame(errs).assign(pid=d.patient_id.values).groupby("pid").mean()
    mae = pat.mean().to_dict()
    best = min(v3.BASE, key=mae.get)
    level = 1 - 0.05 / TRY
    diff = (pat.v6 - pat[best]).values
    ci = v3.boot_ci(diff, np.random.default_rng(t.SEED), level)
    grow = (d.core_vox_out > 1.25 * d.core_vox_in).astype(int).values
    res = {
        "script": "src/track_b_mu_v6.py", "amendment": "paper/trackB_mu_v6_amendment_2026-11-14.md",
        "label": "exploratory, post hoc, try 6 on reused MU-Glioma-Post",
        "n_pairs": int(len(d)), "n_patients": int(len(pat)),
        "patient_mae_log_volume": mae, "best_baseline": best,
        "v6_minus_best": {"mean": float(diff.mean()), "ci_level": level, "ci": ci},
        "success": bool(ci[1] < 0),
        "secondary_auc_growth25_v6": float(roc_auc_score(grow, pred - d.y_in.values)),
    }
    (OUT / "result.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    pd.DataFrame({"patient_id": d.patient_id, "fold": d.fold, "y_out": y, "v6": pred}).to_csv(OUT / "predictions.csv", index=False)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    run()
