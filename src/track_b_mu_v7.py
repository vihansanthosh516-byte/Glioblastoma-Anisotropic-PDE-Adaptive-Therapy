"""Track B v6: MRI brightness + composition latent model. Exploratory, try 7.

Amendment: paper/trackB_mu_v7_amendment_2026-11-14.md (locked before scoring).
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
OUT = ROOT / "output" / "track_b_mu_v7_20261114"
CACHE = OUT / "input_intensity_features.csv"
TRY = 7


def shape_features(pid, tp):
    base = v5.MU / pid / f"Timepoint_{tp}" / f"{pid}_Timepoint_{tp}_"
    a = np.asarray(nib.load(str(base) + "tumorMask.nii.gz").dataobj)
    t1c = np.asarray(nib.load(str(base) + "brain_t1c.nii.gz").dataobj, dtype=np.float32)
    fl = np.asarray(nib.load(str(base) + "brain_t2f.nii.gz").dataobj, dtype=np.float32)
    norm = lambda im: im / max(float(np.median(im[im > 0])), 1e-6) if (im > 0).any() else im
    t1c, fl = norm(t1c), norm(fl)
    mean = lambda im, m: float(im[m].mean()) if m.any() else 0.0
    core, tumor = (a == 1) | (a == 3), a > 0
    shell = ndimage.binary_dilation(tumor, iterations=5) & ~tumor & (fl > 0)
    return {"t1c_enh": mean(t1c, a == 3), "t1c_core_sd": float(t1c[core].std()) if core.any() else 0.0,
            "flair_edema": mean(fl, a == 2), "flair_shell": mean(fl, shell)}


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
    errs["v7"] = np.abs(pred - y)
    pat = pd.DataFrame(errs).assign(pid=d.patient_id.values).groupby("pid").mean()
    mae = pat.mean().to_dict()
    best = min(v3.BASE, key=mae.get)
    level = 1 - 0.05 / TRY
    diff = (pat.v7 - pat[best]).values
    ci = v3.boot_ci(diff, np.random.default_rng(t.SEED), level)
    grow = (d.core_vox_out > 1.25 * d.core_vox_in).astype(int).values
    res = {
        "script": "src/track_b_mu_v7.py", "amendment": "paper/trackB_mu_v7_amendment_2026-11-14.md",
        "label": "exploratory, post hoc, try 7 on reused MU-Glioma-Post",
        "n_pairs": int(len(d)), "n_patients": int(len(pat)),
        "patient_mae_log_volume": mae, "best_baseline": best,
        "v7_minus_best": {"mean": float(diff.mean()), "ci_level": level, "ci": ci},
        "success": bool(ci[1] < 0),
        "secondary_auc_growth25_v7": float(roc_auc_score(grow, pred - d.y_in.values)),
    }
    (OUT / "result.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    pd.DataFrame({"patient_id": d.patient_id, "fold": d.fold, "y_out": y, "v7": pred}).to_csv(OUT / "predictions.csv", index=False)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    run()
