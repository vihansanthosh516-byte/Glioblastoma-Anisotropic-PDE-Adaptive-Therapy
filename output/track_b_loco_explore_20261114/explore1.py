"""Exploration on DEVELOPMENT cohorts only (MU, UCSD, UCSF; all spent). Post hoc, descriptive.

Builds one pooled pair table and prints where a transferable signal may be.
LUMIERE is not read here.
"""
import sys, zipfile
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import track_b_mu_tcks as t, track_b_mu_v5 as v5
import track_b_v8_ucsd_external as ext, track_b_v8_ucsf_alptdg as uf

HERE = Path(__file__).resolve().parent
POOL = HERE / "pooled_pairs.csv"
COMP = ["log_edema", "nec_frac", "enh_frac", "log_whole", "compact"]
KEEP = ["cohort", "patient_id", "day_in", "dt_days", "y_in", "y_out", "core_vox_in", "core_vox_out", "r_hist", "dt_prior", "phase", "age", "gbm"] + COMP


def pool():
    if POOL.exists():
        return pd.read_csv(POOL)
    mu, _ = t.load()
    mu = pd.concat([mu.reset_index(drop=True), v5.features_table(mu).drop(columns=["patient_id", "tp_in"]).reset_index(drop=True)], axis=1)
    mu["cohort"] = "MU"
    u = ext.ucsd_pairs(); u["cohort"] = "UCSD"
    f, _ = uf.build(zipfile.ZipFile(uf.ZIP)); f["cohort"] = "UCSF"
    p = pd.concat([d[KEEP] for d in (mu, u, f)], ignore_index=True)
    p["patient_id"] = p.cohort + "_" + p.patient_id.astype(str)
    p.to_csv(POOL, index=False)
    return p


def pmae(err, pid):
    return float(pd.Series(np.asarray(err)).groupby(np.asarray(pid)).mean().mean())


def main():
    p = pool()
    p["ch"] = p.y_out - p.y_in
    p["y_prev"] = p.y_in - p.r_hist * p.dt_prior
    p["edema_ratio"] = p.log_edema - p.y_in
    print("pairs/patients", p.groupby("cohort").agg(pairs=("y_in", "size"), pat=("patient_id", "nunique")).to_dict())
    bins = [0, 4, 5, 6, 7, 8, 9, 10, 20]
    p["bin"] = pd.cut(p.y_in, bins)
    print("\nA. by input-size bin: n | P(zero out) | median change | mean |change|")
    for c, g in p.groupby("cohort"):
        a = g.groupby("bin", observed=True).agg(n=("ch", "size"), p0=("core_vox_out", lambda s: (s == 0).mean()), med=("ch", "median"), mad=("ch", lambda s: s.abs().mean()))
        print(c); print(a.round(2).to_string())
    print("\nB. pairs with a prior scan: patient MAE of y_in + w*(y_prev - y_in)")
    for c, g in p[p.r_hist.notna()].groupby("cohort"):
        print(c, "n", len(g), {w: round(pmae(np.abs(g.y_in + w * (g.y_prev - g.y_in) - g.y_out), g.patient_id), 3) for w in (0, 0.1, 0.2, 0.3, 0.5)},
              "| max(y_in,y_prev)", round(pmae(np.abs(np.maximum(g.y_in, g.y_prev) - g.y_out), g.patient_id), 3),
              "| corr(last change, next change)", round(float(spearmanr(g.y_in - g.y_prev, g.ch)[0]), 3))
    print("\nC. Spearman of signed change with inputs (all pairs; then pairs with y_out>0)")
    for c, g in p.groupby("cohort"):
        nz = g[g.core_vox_out > 0]
        print(c, {k: (round(float(spearmanr(g[k], g.ch)[0]), 2), round(float(spearmanr(nz[k], nz.ch)[0]), 2)) for k in ["y_in", "dt_days", "edema_ratio", "nec_frac", "compact", "log_whole", "age"] if g[k].notna().sum() > 10})
    print("\nD. change by interval bin: n | median change | mean |change|")
    p["dtb"] = pd.cut(p.dt_days, [0, 45, 75, 105, 400])
    for c, g in p.groupby("cohort"):
        print(c, g.groupby("dtb", observed=True).agg(n=("ch", "size"), med=("ch", "median"), mad=("ch", lambda s: s.abs().mean())).round(2).to_dict("index"))
    print("\nE. change by phase label (as coded): n | median | mean |change|")
    for c, g in p.groupby("cohort"):
        print(c, g.groupby("phase").agg(n=("ch", "size"), med=("ch", "median"), mad=("ch", lambda s: s.abs().mean())).round(2).to_dict("index"))
    print("\nF. persistence patient MAE", {c: round(pmae(g.ch.abs(), g.patient_id), 3) for c, g in p.groupby("cohort")})
    print("   oracle: zero-outcome pairs forecast exactly", {c: round(pmae(np.where(g.core_vox_out == 0, 0, g.ch.abs()), g.patient_id), 3) for c, g in p.groupby("cohort")})
    print("   oracle: best constant shift per cohort (median change)", {c: round(pmae((g.ch - g.ch.median()).abs(), g.patient_id), 3) for c, g in p.groupby("cohort")})


if __name__ == "__main__":
    main()
