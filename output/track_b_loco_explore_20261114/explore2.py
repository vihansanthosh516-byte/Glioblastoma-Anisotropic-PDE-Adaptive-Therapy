"""Exploration 2 on DEVELOPMENT cohorts only (MU, UCSD, UCSF). Post hoc, descriptive.

Adds days since first surgery and checks (1) time-since-surgery as a moderator,
(2) mask flicker after a zero outcome in MU, (3) leave-one-cohort-out ceiling of
simple median-change models. LUMIERE is not read here.
"""
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.linear_model import QuantileRegressor

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import track_b_mu_tcks as t
from explore1 import pool, pmae

num = lambda s: pd.to_numeric(s, errors="coerce")


def add_surgery_days(p):
    c = pd.read_excel(t.CLIN, sheet_name="MU Glioma Post"); c.columns = [str(k).strip() for k in c.columns]
    sx = dict(zip("MU_" + c["Patient_ID"].astype(str).str.strip(), num(c["Number of days from Diagnosis to First surgery or procedure"])))
    p["dsx"] = np.nan
    m = p.cohort == "MU"
    p.loc[m, "dsx"] = p.loc[m, "day_in"] - p.loc[m, "patient_id"].map(sx).fillna(0.0)
    m = p.cohort != "MU"          # UCSD and UCSF day_in is already days since first surgery / diagnosis
    p.loc[m, "dsx"] = p.loc[m, "day_in"]
    p.loc[(p.dsx < 0) | (p.dsx > 15000), "dsx"] = np.nan
    return p


def main():
    p = add_surgery_days(pool())
    p["ch"] = p.y_out - p.y_in
    p["zero"] = p.core_vox_out == 0
    p["edema_ratio"] = p.log_edema - p.y_in
    print("days since first surgery at input scan:", p.groupby("cohort").dsx.describe()[["count", "25%", "50%", "75%"]].round(0).to_dict("index"))
    p["sxb"] = pd.cut(p.dsx, [-1, 60, 180, 365, 730, 20000])
    print("\n1. by days since surgery: n | P(zero) | median change | mean |change| | median y_in")
    for c, g in p.groupby("cohort"):
        print(c, g.groupby("sxb", observed=True).agg(n=("ch", "size"), p0=("zero", "mean"), med=("ch", "median"), mad=("ch", lambda s: s.abs().mean()), yin=("y_in", "median")).round(2).to_dict("index"))
    print("\n   small lesions only (y_in < 7): median change by days since surgery")
    for c, g in p[p.y_in < 7].groupby("cohort"):
        print(c, g.groupby("sxb", observed=True).agg(n=("ch", "size"), p0=("zero", "mean"), med=("ch", "median")).round(2).to_dict("index"))

    # 2. flicker in MU: what is the core size at the scan after a zero outcome?
    pr = pd.read_csv(t.PAIRS).sort_values(["patient_id", "tp_in"])
    z = pr[(pr.core_vox_out == 0) & (pr.core_vox_in > 0)]
    nxt = pr.set_index(["patient_id", "tp_in"]).core_vox_out
    after = [nxt.get((r.patient_id, r.tp_out), np.nan) for r in z.itertuples()]
    print("\n2. MU zero outcomes (all pairs):", len(z), "| has a later scan:", int(np.sum(~np.isnan(after))), "| later scan core > 0:", int(np.nansum(np.array(after) > 0)), "| later core voxels:", [int(a) for a in after if not np.isnan(a)])
    emp = pr[pr.core_vox_in == 0]
    print("   MU pairs that START at zero:", len(emp), "| regrow (>0):", int((emp.core_vox_out > 0).sum()), "| median regrown voxels:", float(emp.core_vox_out[emp.core_vox_out > 0].median()) if (emp.core_vox_out > 0).any() else None)

    # 3. leave-one-cohort-out ceiling of simple median-change models
    p["z"] = p.y_in - p.groupby("cohort").y_in.transform("median")        # centred on own cohort inputs (no outcomes)
    p["ldt"] = np.log(p.dt_days / 60.0)
    p["lsx"] = np.log((p.dsx.fillna(p.dsx.median()) + 30) / 365.0)
    p["er"] = p.edema_ratio - p.groupby("cohort").edema_ratio.transform("median")
    p["small"] = np.clip(7 - p.y_in, 0, None)
    sets = {"const": [], "ldt": ["ldt"], "z": ["z"], "z+ldt": ["z", "ldt"], "small+ldt": ["small", "ldt"], "small+ldt+lsx": ["small", "ldt", "lsx"],
            "z+ldt+lsx": ["z", "ldt", "lsx"], "z+ldt+er": ["z", "ldt", "er"], "small+ldt+lsx+small:lsx": ["small", "ldt", "lsx", "sl"]}
    p["sl"] = p.small * p.lsx
    cohorts = ["MU", "UCSD", "UCSF"]
    print("\n3. LOCO patient MAE minus persistence (negative = better). Train on the other two cohorts.")
    rows = {}
    for name, cols in sets.items():
        out = {}
        for te_c in cohorts:
            tr, te = p[p.cohort != te_c], p[p.cohort == te_c]
            Xtr = tr[cols].values if cols else np.zeros((len(tr), 1)); Xte = te[cols].values if cols else np.zeros((len(te), 1))
            q = QuantileRegressor(quantile=0.5, alpha=0.0, solver="highs").fit(Xtr, tr.ch.values)
            pred = te.y_in.values + q.predict(Xte)
            out[te_c] = round(pmae(np.abs(pred - te.y_out.values), te.patient_id) - pmae(te.ch.abs(), te.patient_id), 3)
        rows[name] = out
    print(pd.DataFrame(rows).T.to_string())
    print("\n   same models, 5-fold patient CV inside each cohort (own-cohort ceiling)")
    rng = np.random.default_rng(0)
    rows = {}
    for name, cols in sets.items():
        out = {}
        for c in cohorts:
            g = p[p.cohort == c].reset_index(drop=True)
            pid = g.patient_id.unique(); fold = dict(zip(pid, rng.permutation(len(pid)) % 5))
            f = g.patient_id.map(fold).values
            pred = np.zeros(len(g))
            for k in range(5):
                Xtr = g.loc[f != k, cols].values if cols else np.zeros(((f != k).sum(), 1)); Xte = g.loc[f == k, cols].values if cols else np.zeros(((f == k).sum(), 1))
                q = QuantileRegressor(quantile=0.5, alpha=0.0, solver="highs").fit(Xtr, g.ch.values[f != k])
                pred[f == k] = g.y_in.values[f == k] + q.predict(Xte)
            out[c] = round(pmae(np.abs(pred - g.y_out.values), g.patient_id) - pmae(g.ch.abs(), g.patient_id), 3)
        rows[name] = out
    print(pd.DataFrame(rows).T.to_string())


if __name__ == "__main__":
    main()
