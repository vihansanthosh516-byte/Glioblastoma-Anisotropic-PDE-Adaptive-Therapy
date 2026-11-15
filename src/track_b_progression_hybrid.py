"""Track B new question: 1-year progression from one post-treatment MRI. Dev on MU (CV), test once on UCSF-ALPTDG.

Protocol: paper/trackB_progression_hybrid_protocol_2026-11-14.md (locked before fitting).
"""
from __future__ import annotations

import gzip
import io
import json
import zipfile
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import track_b_mu_tcks as t
import track_b_mu_v5 as v5
import track_b_v8_ucsf_alptdg as ua

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "output" / "track_b_progression_hybrid_20261114"
COMP = ["log_edema", "nec_frac", "enh_frac", "log_whole", "compact"]
CLIN = ["age", "gbm", "log_days_dx", "on_tx"]
ALL = ["log_core", "log_enh"] + COMP + CLIN
TX_DAYS = 120.0  # assumed remaining on-treatment days in the projection


def mu_rows():
    d, allp = t.load()
    c = pd.read_excel(t.CLIN, sheet_name="MU Glioma Post")
    c.columns = [str(x).strip() for x in c.columns]
    death = pd.Series(pd.to_numeric(c[next(x for x in c.columns if "diagnosis to death" in x.lower())], errors="coerce").values,
                      index=c["Patient_ID"].astype(str).str.strip())
    last = allp.groupby("patient_id").day_out.max()
    comp = v5.features_table(d)
    s = pd.concat([d[["patient_id", "tp_in", "day_in", "core_vox_in", "fold", "prog", "age", "gbm",
                      "rt_s", "rt_e", "ad_s", "ad_e", "mgmt"]].reset_index(drop=True), comp[COMP].reset_index(drop=True)], axis=1)
    s = s.drop_duplicates(["patient_id", "tp_in"]).reset_index(drop=True)
    s["death"] = s.patient_id.map(death)
    s["last"] = s.patient_id.map(last)
    gap = s.prog - s.day_in
    follow = np.fmax(s.death.fillna(-1), s["last"].fillna(-1)) - s.day_in
    ev = (gap > 0) & (gap <= 365)
    non = (gap > 365) | (s.prog.isna() & (follow > 365))
    keep = (ev | non) & ~(s.prog <= s.day_in)
    s = s[keep].reset_index(drop=True)
    s["event"] = ev[keep].astype(int).values
    s["on_tx"] = (((s.rt_s <= s.day_in) & (s.day_in <= s.rt_e)) | ((s.ad_s <= s.day_in) & (s.day_in <= s.ad_e))).astype(float)
    s["log_core"] = np.log(s.core_vox_in + 1.0)
    s["log_enh"] = np.log(s.core_vox_in * s.enh_frac + 1.0)
    s["log_days_dx"] = np.log(s.day_in.clip(lower=0) + 1.0)
    s["y_in"] = s.log_core
    return s, d


def ucsf_rows():
    z = zipfile.ZipFile(ua.ZIP)
    x = pd.read_excel(io.BytesIO(z.read(ua.TABLE)), sheet_name=["Clinical Info", "UCSF-LPTDG", "UCSF-PDGM Mapping"])
    c = x["Clinical Info"]
    c.columns = [str(k).strip() for k in c.columns]
    col = lambda key: next(k for k in c.columns if key.lower() in k.lower())
    lp = x["UCSF-LPTDG"]
    age = lp[lp.Timepoint == 1].set_index("SubjectID")["Patient Age"]
    used = set(pd.read_csv(ua.USED).pid.astype(int)) | set(x["UCSF-PDGM Mapping"]["UCSF-LPTDG ID"].astype(int))
    rows = []
    for _, r in c.iterrows():
        sid = int(r["SubjectID"])
        prog = pd.to_numeric(r[col("1st scan to next progression")], errors="coerce")
        death = pd.to_numeric(r[col("death to 1st scan")], errors="coerce")
        death = abs(death) if pd.notna(death) else np.nan
        if pd.notna(prog) and prog <= 365:
            ev = 1
        elif (pd.notna(prog) and prog > 365) or (pd.isna(prog) and pd.notna(death) and death > 365):
            ev = 0
        else:
            continue
        a = np.asarray(nib.Nifti1Image.from_bytes(gzip.decompress(z.read(f"{sid}/{sid}_time1_seg.nii.gz"))).dataobj)
        n1, n2, n3 = (int((a == k).sum()) for k in (1, 2, 3))
        core = n1 + n3
        box = 1
        if core:
            idx = np.argwhere((a == 1) | (a == 3))
            box = int(np.prod(idx.max(0) - idx.min(0) + 1))
        days = pd.to_numeric(r[col("1st surgery/DX")], errors="coerce")
        rows.append({"patient_id": str(sid), "used_before": sid in used, "event": ev,
                     "log_core": np.log(core + 1.0), "log_enh": np.log(n3 + 1.0),
                     "log_edema": np.log(n2 + 1.0), "nec_frac": n1 / core if core else 0.0,
                     "enh_frac": n3 / core if core else 0.0, "log_whole": np.log(n1 + n2 + n3 + 1.0),
                     "compact": core / box if core else 0.0, "age": pd.to_numeric(age.get(sid), errors="coerce"),
                     "gbm": float("glioblastoma" in str(r[col("WHO 2021 Diagnosis")]).lower()),
                     "log_days_dx": np.log(max(days, 0) + 1.0) if pd.notna(days) else np.nan,
                     "on_tx": float(str(r[col("On tx at 1st scan")]).strip().lower() not in ("n", "nan"))})
    u = pd.DataFrame(rows)
    u["y_in"] = u.log_core
    return u


def mech_features(pairs, rows):
    """Gompertz + on-treatment rate projection fitted on MU pairs; applied to scan rows."""
    gom = t.gompertz_fit(pairs)
    on = pairs.phase.isin(["chemoradiation", "adjuvant"])
    r_tx, r_off = float(pairs.r_obs[on].median()), float(pairs.r_obs[~on].median())
    proj = pd.DataFrame({"y_in": rows.y_in.values, "dt_days": 365.0})
    delta = gom(proj) - rows.y_in.values + rows.on_tx.fillna(0).values * (r_tx - r_off) * TX_DAYS
    rate = delta / 365.0
    t40 = np.where(rate > 0, np.minimum(np.log(1.4) / np.maximum(rate, 1e-9), 730.0), 730.0)
    return np.c_[delta, t40, delta * rows.y_in.values, t40 * rows.y_in.values]


def models():
    lr = lambda: make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=2000))
    return {"B1_core": (["log_core"], lr), "B2_enh": (["log_enh"], lr), "B3_clinical": (CLIN, lr),
            "B4_clin_vol": (CLIN + ["log_core", "log_enh"], lr), "B5_composition": (ALL, lr),
            "B6_gbc": (ALL, lambda: HistGradientBoostingClassifier(random_state=0))}


def fit_all(tr_rows, tr_pairs, te_rows):
    med = tr_rows[ALL].median()
    Xtr, Xte = tr_rows[ALL].fillna(med), te_rows[ALL].fillna(med)
    out = {}
    for name, (cols, mk) in models().items():
        m = mk().fit(Xtr[cols], tr_rows.event)
        out[name] = m.predict_proba(Xte[cols])[:, 1]
    Htr = np.c_[Xtr.values, mech_features(tr_pairs, tr_rows.assign(on_tx=Xtr.on_tx))]
    Hte = np.c_[Xte.values, mech_features(tr_pairs, te_rows.assign(on_tx=Xte.on_tx))]
    h = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=2000)).fit(Htr, tr_rows.event)
    out["hybrid"] = h.predict_proba(Hte)[:, 1]
    return out


def boot_auc_diff(y, a, b, rng, n=5000):
    idx = np.arange(len(y))
    vals = []
    for _ in range(n):
        s = rng.choice(idx, len(idx))
        if y[s].min() == y[s].max():
            continue
        vals.append(roc_auc_score(y[s], a[s]) - roc_auc_score(y[s], b[s]))
    return [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))]


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT / "result.json").exists():
        raise SystemExit("locked run already done; refusing to overwrite")
    s, pairs = mu_rows()
    oof = {k: np.full(len(s), np.nan) for k in list(models()) + ["hybrid"]}
    for f in sorted(s.fold.unique()):
        trm, tem = (s.fold != f).values, (s.fold == f).values
        tp = pairs[pairs.fold != f]
        for k, v in fit_all(s[trm], tp, s[tem]).items():
            oof[k][tem] = v
    dev_auc = {k: float(roc_auc_score(s.event, v)) for k, v in oof.items()}
    comparator = max([k for k in dev_auc if k != "hybrid"], key=dev_auc.get)
    u = ucsf_rows()
    pt = fit_all(s, pairs, u)
    rng = np.random.default_rng(t.SEED)
    res = {"script": "src/track_b_progression_hybrid.py", "protocol": "paper/trackB_progression_hybrid_protocol_2026-11-14.md",
           "dev_mu": {"n_scans": int(len(s)), "n_patients": int(s.patient_id.nunique()), "n_events": int(s.event.sum()), "auc": dev_auc},
           "comparator_from_dev": comparator}
    for label, m in (("primary_untouched", ~u.used_before.values), ("all_eligible", np.ones(len(u), bool))):
        y = u.event.values[m]
        auc = {k: float(roc_auc_score(y, v[m])) for k, v in pt.items()}
        ci = boot_auc_diff(y, pt["hybrid"][m], pt[comparator][m], rng)
        res[label] = {"n_patients": int(m.sum()), "n_events": int(y.sum()), "auc": auc,
                      "hybrid_minus_comparator": {"mean": auc["hybrid"] - auc[comparator], "ci95": ci}}
    res["success"] = bool(res["primary_untouched"]["hybrid_minus_comparator"]["ci95"][0] > 0)
    (OUT / "result.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    run()
