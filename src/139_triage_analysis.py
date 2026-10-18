#!/usr/bin/env python3
"""Script 139: triage analysis exactly as declared in Amendment 14 (A14.5). EXPLORATORY, leave-one-centre-out, spent sets.

L2 logistic regression, C = 1.0, class_weight balanced, features standardised on the training folds, no tuning.
Pooled out-of-fold predictions. Baseline B1 = log core volume only. Patient bootstrap 2000 draws, seed 20261018.
PRIMARY: pooled AUC for FAIL (standard coverage < 0.80); success = lower CI bound > 0.5.
S1 AUC(full) - AUC(B1); S2 AUC for FAIL90; S3 Spearman(risk, standard coverage); S4 per-centre AUC.
Output: output/triage/triage_results.json, triage_oof.csv. Run once; no re-runs with other settings (A14.6).
"""
from __future__ import annotations

import json
from pathlib import Path as _Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "triage"
FEATS = ["f1_log_core_ml", "f2_log1p_edema_ml", "f3_enh_share", "f4_lateral_mm", "f5_ap_pos", "f6_si_pos",
         "f7_dist_surface_mm", "f8_margin_truncation", "f9_max_extent_mm", "f10_n_components"]
BASE = ["f1_log_core_ml"]
SEED, NBOOT = 20261018, 2000


def loco(df, cols, target):
    oof = np.full(len(df), np.nan)
    for c in df["centre"].unique():
        te = (df["centre"] == c).to_numpy()
        tr = ~te
        ytr = df.loc[tr, target].to_numpy()
        if ytr.min() == ytr.max():
            continue
        sc = StandardScaler().fit(df.loc[tr, cols])
        m = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000).fit(sc.transform(df.loc[tr, cols]), ytr)
        oof[te] = m.predict_proba(sc.transform(df.loc[te, cols]))[:, 1]
    return oof


def auc(y, p):
    return float("nan") if len(np.unique(y)) < 2 else float(roc_auc_score(y, p))


def ci(vals):
    v = np.asarray(vals, float)
    v = v[~np.isnan(v)]
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))], int(len(v))


def main():
    df = pd.read_csv(OUT / "triage_table.csv")
    res = {"n": int(len(df)), "n_by_centre": df["centre"].value_counts().to_dict(),
           "prevalence_fail": float(df["fail"].mean()), "prevalence_fail90": float(df["fail90"].mean()),
           "prevalence_fail_by_centre": df.groupby("centre")["fail"].mean().round(3).to_dict(),
           "settings": {"C": 1.0, "class_weight": "balanced", "cv": "leave-one-centre-out", "nboot": NBOOT, "seed": SEED}}
    p_full, p_base = loco(df, FEATS, "fail"), loco(df, BASE, "fail")
    p90 = loco(df, FEATS, "fail90")
    y, y90 = df["fail"].to_numpy(), df["fail90"].to_numpy()
    rng = np.random.default_rng(SEED)
    n = len(df)
    bf, bd, b90 = [], [], []
    for _ in range(NBOOT):
        i = rng.integers(0, n, n)
        a, b = auc(y[i], p_full[i]), auc(y[i], p_base[i])
        bf.append(a)
        bd.append(a - b)
        b90.append(auc(y90[i], p90[i]))
    res["PRIMARY_auc_full"] = {"auc": auc(y, p_full), "ci95": ci(bf)[0], "success_lower_gt_0.5": bool(ci(bf)[0][0] > 0.5)}
    res["B1_auc_core_volume_only"] = auc(y, p_base)
    res["S1_auc_full_minus_B1"] = {"diff": auc(y, p_full) - auc(y, p_base), "ci95": ci(bd)[0], "success_lower_gt_0": bool(ci(bd)[0][0] > 0)}
    res["S2_auc_fail90"] = {"auc": auc(y90, p90), "ci95": ci(b90)[0]}
    rho, pv = spearmanr(p_full, df["standard_coverage"])
    res["S3_spearman_risk_vs_standard_coverage"] = {"rho": float(rho), "p": float(pv)}
    res["S4_auc_by_centre"] = {c: {"n": int((df["centre"] == c).sum()), "auc_full": auc(y[(df["centre"] == c).to_numpy()], p_full[(df["centre"] == c).to_numpy()]),
                                   "auc_B1": auc(y[(df["centre"] == c).to_numpy()], p_base[(df["centre"] == c).to_numpy()])} for c in df["centre"].unique()}
    df.assign(p_full=p_full, p_base=p_base, p_fail90=p90)[["pid", "centre", "standard_coverage", "fail", "fail90", "p_full", "p_base", "p_fail90"]].to_csv(OUT / "triage_oof.csv", index=False)
    (OUT / "triage_results.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
