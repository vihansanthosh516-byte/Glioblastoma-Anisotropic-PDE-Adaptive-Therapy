"""Frozen v8 external test on UCSD-PTGBM. Fit once on all MU primary pairs.

Protocol: paper/trackB_v8_ucsd_external_protocol_2026-11-14.md (locked before scoring).
"""
from __future__ import annotations

import glob
import json
import os
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import LogisticRegression, QuantileRegressor

import track_b_mu_tcks as t
import track_b_mu_v3 as v3
import track_b_mu_v5 as v5

ROOT = Path(__file__).resolve().parent.parent
UCSD = ROOT / "data" / "external" / "UCSD-PTGBM-v3"
UCLIN = ROOT / "research_notes" / "trackB_confirmatory_audit_ucsd" / "UCSD_PTGBM-clinical-information_v3_2026-12-Mar.xlsx"
OUT = ROOT / "output" / "track_b_v8_ucsd_external_20261114"
MIN_PATIENTS = 25


def comp_from_mask(path):
    a = np.asarray(nib.load(path).dataobj)
    n1, n2, n3 = (int((a == k).sum()) for k in (1, 2, 3))
    core = n1 + n3
    box = 1
    if core:
        idx = np.argwhere((a == 1) | (a == 3))
        box = int(np.prod(idx.max(0) - idx.min(0) + 1))
    return core, {"log_edema": np.log(n2 + 1.0), "nec_frac": n1 / core if core else 0.0,
                  "enh_frac": n3 / core if core else 0.0, "log_whole": np.log(n1 + n2 + n3 + 1.0),
                  "compact": core / box if core else 0.0}


def num(x):
    return pd.to_numeric(x, errors="coerce")


def ucsd_pairs():
    c = pd.read_excel(UCLIN, sheet_name="Cinical info")
    c.columns = [str(x).strip() for x in c.columns]
    col = lambda key: next(x for x in c.columns if key.lower() in x.lower())
    s = pd.DataFrame({
        "session": c["ID"].astype(str).str.strip(),
        "day": -num(c[col("Date of initial surgery, treatment or diagnosis")]),
        "rt_rel": num(c[col("Date of first radiation")]),
        "chemo_rel": num(c[col("Date of 1st chemo start")]),
        "pfs": num(c["Progression Free survival"]),
        "age": num(c["Patient's Age"]),
        "gbm": (c["Primary Diagnosis"].astype(str).str.strip() == "Glioblastoma").astype(float),
    })
    s["patient_id"] = s.session.str.rsplit("_", n=1).str[0]
    s["rt_s"], s["ad_s"] = s.day + s.rt_rel, s.day + s.chemo_rel
    segs = {}
    for p in sorted(glob.glob(str(UCSD / "**" / "*_BraTS_tumor_seg.nii.gz"), recursive=True)):
        segs.setdefault(os.path.basename(p).split("_BraTS")[0], p)
    s = s[s.session.isin(segs) & s.day.notna()].sort_values(["patient_id", "day"])
    feats = {}
    for sess in s.session:
        feats[sess] = comp_from_mask(segs[sess])
    rows = []
    for pid, g in s.groupby("patient_id"):
        g = g.reset_index(drop=True)
        for i in range(len(g) - 1):
            a, b = g.iloc[i], g.iloc[i + 1]
            dt = b.day - a.day
            core_in, comp = feats[a.session]
            core_out, _ = feats[b.session]
            r_hist = dt_prior = np.nan
            if i > 0:
                pv = g.iloc[i - 1]
                cp = feats[pv.session][0]
                if cp > 0 and a.day > pv.day:
                    dt_prior = a.day - pv.day
                    r_hist = (np.log(core_in + 1) - np.log(cp + 1)) / dt_prior
            rows.append({"patient_id": pid, "day_in": a.day, "day_out": b.day, "dt_days": dt,
                         "core_vox_in": core_in, "core_vox_out": core_out, "r_hist": r_hist, "dt_prior": dt_prior,
                         "rt_s": a.rt_s, "rt_e": a.rt_s + 42, "ad_s": a.ad_s, "ad_e": a.ad_s + 180,
                         "pfs": a.pfs, "age": a.age, "mgmt": np.nan, "gbm": a.gbm, **comp})
    d = pd.DataFrame(rows)
    d = d[(d.dt_days >= 14) & (d.dt_days <= 365) & (d.core_vox_in > 0)].reset_index(drop=True)
    d["y_in"], d["y_out"] = np.log(d.core_vox_in + 1.0), np.log(d.core_vox_out + 1.0)
    m = (d.day_in + d.day_out) / 2
    ph = np.full(len(d), "other", dtype=object)
    ph[((d.ad_s <= m) & (m <= d.ad_e)).values] = "adjuvant"
    ph[((d.rt_s <= m) & (m <= d.rt_e)).values] = "chemoradiation"
    ph[(d.pfs <= d.day_in).values] = "post_progression"
    d["phase"] = ph
    d["r_obs"] = (d.y_out - d.y_in) / d.dt_days
    return d


def fit_predict(tr, te, comp_tr, comp_te):
    out = {}
    r_pop = float(tr.r_obs.median())
    out["persistence"] = te.y_in.values
    out["two_point_exp"] = (te.y_in + te.r_hist.fillna(0.0) * te.dt_days).values
    out["pop_rate"] = (te.y_in + r_pop * te.dt_days).values
    t2, s2 = t.eb_params(tr, r_pop)
    out["eb_random_slope"] = np.asarray(t.shrink_pred(te, r_pop, t2, s2)[0])
    out["gompertz"] = t.gompertz_fit(tr)(te)
    g = HistGradientBoostingRegressor(random_state=0).fit(t.features(tr), tr.y_out)
    out["gbr"] = g.predict(t.features(te))
    Xtr = pd.concat([v3.design(tr, v3.latent(tr, tr)).reset_index(drop=True), comp_tr.reset_index(drop=True)], axis=1)
    Xte = pd.concat([v3.design(te, v3.latent(te, tr)).reset_index(drop=True), comp_te.reset_index(drop=True)], axis=1)
    mu, sd = Xtr.mean(), Xtr.std().replace(0, 1.0)
    Ztr, Zte = (Xtr - mu) / sd, (Xte - mu) / sd
    col = (tr.core_vox_out < 0.05 * tr.core_vox_in).values
    p = LogisticRegression(max_iter=1000).fit(Ztr, col).predict_proba(Zte)[:, 1]
    q = QuantileRegressor(quantile=0.5, alpha=0.0, solver="highs").fit(Ztr[~col], tr.y_out.values[~col])
    out["v8"] = np.where(p > 0.5, float(np.median(tr.y_out.values[col])), q.predict(Zte))
    return out, p


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT / "result.json").exists():
        raise SystemExit("locked test already run; refusing to overwrite")
    mu, _ = t.load()
    comp_mu = v5.features_table(mu).drop(columns=["patient_id", "tp_in"]).reset_index(drop=True)
    u = ucsd_pairs()
    n_pat = int(u.patient_id.nunique())
    res = {"script": "src/track_b_v8_ucsd_external.py", "protocol": "paper/trackB_v8_ucsd_external_protocol_2026-11-14.md",
           "label": "external test of frozen v8; UCSD previously exposed for a different target",
           "n_pairs": int(len(u)), "n_patients": n_pat, "phase_counts": u.phase.value_counts().to_dict()}
    if n_pat < MIN_PATIENTS:
        res["status"] = "no-go: fewer than 25 patients"
    else:
        comp_u = u[["log_edema", "nec_frac", "enh_frac", "log_whole", "compact"]]
        preds, p = fit_predict(mu, u, comp_mu, comp_u)
        y = u.y_out.values
        pat = pd.DataFrame({k: np.abs(v - y) for k, v in preds.items()}).assign(pid=u.patient_id.values).groupby("pid").mean()
        diff = (pat.v8 - pat.pop_rate).values
        ci = v3.boot_ci(diff, np.random.default_rng(t.SEED), 0.95)
        res.update({"status": "scored", "patient_mae_log_volume": pat.mean().to_dict(),
                    "comparator": "pop_rate", "v8_minus_pop_rate": {"mean": float(diff.mean()), "ci95": ci},
                    "success": bool(ci[1] < 0), "n_predicted_collapse": int((p > 0.5).sum()),
                    "n_collapse": int((u.core_vox_out < 0.05 * u.core_vox_in).sum())})
        pd.DataFrame(preds).assign(patient_id=u.patient_id, y_out=y, phase=u.phase).to_csv(OUT / "predictions.csv", index=False)
    (OUT / "result.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    run()
