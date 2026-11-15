"""Read-only post hoc diagnostics (spent test sets; descriptive only). Writes nothing."""
import io, sys, zipfile
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(sys.argv[1]); sys.path.insert(0, str(ROOT / "src"))
import track_b_mu_tcks as t, track_b_mu_v3 as v3, track_b_mu_v5 as v5
import track_b_v8_ucsd_external as ext, track_b_v8_ucsf_alptdg as uf
from sklearn.linear_model import LogisticRegression, QuantileRegressor
pd.set_option("display.width", 250)

mu, _ = t.load()
comp_mu = v5.features_table(mu).drop(columns=["patient_id", "tp_in"]).reset_index(drop=True)
z = zipfile.ZipFile(uf.ZIP)
d, n_all = uf.build(z)
prim = (~d.used_before).values
te = d[prim].reset_index(drop=True)

def decomp(tr, comp_tr, te, comp_te, tag):
    Xtr = pd.concat([v3.design(tr, v3.latent(tr, tr)).reset_index(drop=True), comp_tr.reset_index(drop=True)], axis=1)
    Xte = pd.concat([v3.design(te, v3.latent(te, tr)).reset_index(drop=True), comp_te.reset_index(drop=True)], axis=1)
    m, sd = Xtr.mean(), Xtr.std().replace(0, 1.0)
    Ztr, Zte = (Xtr - m) / sd, (Xte - m) / sd
    col = (tr.core_vox_out < 0.05 * tr.core_vox_in).values
    q = QuantileRegressor(quantile=0.5, alpha=0.0, solver="highs").fit(Ztr[~col], tr.y_out.values[~col])
    tab = pd.DataFrame({"coef_std": q.coef_, "mu_mean": m.values, "mu_sd": sd.values, "test_mean": Xte.mean().values,
                        "z_shift": Zte.mean().values, "contrib": q.coef_ * Zte.mean().values}, index=Xtr.columns).round(3)
    print("==", tag, "intercept", round(float(q.intercept_), 3), "| MU mean y_in", round(tr.y_in.mean(), 3), "test mean y_in", round(te.y_in.mean(), 3))
    print(tab.to_string())
    pred = q.predict(Zte)
    print(" mean pred", round(pred.mean(), 3), "mean y_in", round(te.y_in.mean(), 3), "mean y_out", round(te.y_out.mean(), 3), "sum contrib", round(tab.contrib.sum(), 3))
    # effective slope on y_in in raw units (x_hat == y_in when no prior)
    print(" raw slope on y_in+x_hat", round(float(q.coef_[0] / sd.iloc[0] + q.coef_[1] / sd.iloc[1]), 3))

decomp(mu, comp_mu, te, te[uf.V5C], "v8 part-2 on UCSF untouched")
u = ext.ucsd_pairs()
decomp(mu, comp_mu, u, u[uf.V5C], "v8 part-2 on UCSD")

# patient-level SDs of paired differences (for power)
def sd_diff(path, a, b, filt=None):
    p = pd.read_csv(path)
    if filt is not None: p = p[filt(p)]
    return p
pm = pd.read_csv(ROOT / "output/track_b_mu_v8_20261114/predictions.csv")
pb = pd.read_csv(v3.V2)
e = pd.DataFrame({"v8": (pm.v8 - pm.y_out).abs(), "pop": (pb.pop_rate - pm.y_out).abs(), "pers": (pb.persistence - pm.y_out).abs(), "pid": pm.patient_id}).groupby("pid").mean()
print("MU patient diff v8-pop: mean", round((e.v8 - e["pop"]).mean(), 4), "sd", round((e.v8 - e["pop"]).std(), 3), "n", len(e), "| v8-pers mean", round((e.v8 - e.pers).mean(), 4), "sd", round((e.v8 - e.pers).std(), 3))
for name in ["track_b_v8_ucsd_external_20261114", "track_b_v8_ucsf_alptdg_20261114"]:
    p = pd.read_csv(ROOT / "output" / name / "predictions.csv")
    if "used_before" in p: p = p[~p.used_before]
    e = pd.DataFrame({k: (p[k] - p.y_out).abs() for k in ["v8", "pop_rate", "persistence"]}).assign(pid=p.patient_id.values).groupby("pid").mean()
    print(name, "n", len(e), "v8-pop mean", round((e.v8 - e.pop_rate).mean(), 4), "sd", round((e.v8 - e.pop_rate).std(), 3), "| v8-pers mean", round((e.v8 - e.persistence).mean(), 4), "sd", round((e.v8 - e.persistence).std(), 3))

# UCSF: surgery / treatment fields the scripts did not use
x = pd.read_excel(io.BytesIO(z.read(uf.TABLE)), sheet_name=None)
print(x["Clinical_Dictionary"].to_string())
c = x["Clinical Info"]; c.columns = [str(k).strip() for k in c.columns]
c["patient_id"] = c.SubjectID.astype(int).astype(str)
j = d.merge(c, on="patient_id", how="left")
sx2 = pd.to_numeric(j["Days from 1st scan to 2nd SX"], errors="coerce")
sx3 = pd.to_numeric(j["Days from 1st scan to 3rd SX"], errors="coerce")
between = ((sx2 > 0) & (sx2 <= j.dt_days)) | ((sx3 > 0) & (sx3 <= j.dt_days))
err = (j.y_out - j.y_in).abs()
zero = j.core_vox_out == 0
big = j.core_vox_out < 0.05 * j.core_vox_in
print("UCSF eligible", len(j), "day_in describe", j.day_in.describe().round(0).to_dict())
print("UCSF dt", j.dt_days.describe().round(0).to_dict())
print("UCSF surgery between scans: n", int(between.sum()), "| zero-out n", int(zero.sum()), "of which surgery between", int((zero & between).sum()),
      "| >95% drop n", int(big.sum()), "of which surgery between", int((big & between).sum()))
print("UCSF pers err: surgery-between mean", round(err[between].mean(), 3), "others", round(err[~between].mean(), 3), "share of err", round(err[between].sum() / err.sum(), 3))
for k in ["On tx at 1st scan (c = chemo, rt, n = none)", "On tx at 2nd scan (c = chemo, rt, n = none)", "1st Chemo type", "Number of surgeries"]:
    print(k, j[k].astype(str).str.strip().value_counts(dropna=False).head(8).to_dict())
print("UCSF 2nd SX sign sample", sx2.dropna().describe().round(0).to_dict())
print("UCSF chemo-start col", pd.to_numeric(j["Days from 1st chemo start to 1st scan"], errors="coerce").describe().round(0).to_dict())
print("UCSF RT col", pd.to_numeric(j["Days from 1st scan to 1st RT start (neg = RT first)"], errors="coerce").describe().round(0).to_dict())
print("UCSF next progression col", pd.to_numeric(j["Days from 1st scan to next progression"], errors="coerce").describe().round(0).to_dict())
print("UCSF y_in<log(500): n", int((j.core_vox_in < 500).sum()), "mean err", round(err[j.core_vox_in < 500].mean(), 3), "rest", round(err[j.core_vox_in >= 500].mean(), 3))

# UCSD: units and surgery/avastin fields
cu = pd.read_excel(ext.UCLIN, sheet_name="Cinical info"); cu.columns = [str(k).strip() for k in cu.columns]
for k in ["Progression Free survival", "Overall survival", "Days from Acquisition to Date of initial surgery, treatment or diagnosis",
          "Days from Acquisition to Date of last surgery prior to scan if different", "Days from Acquisition to Date of last Avastin treatment prior to scan",
          "Days from Acquistion to Date of last radiation prior to scan", "Number of surgeries", "Number of radiation courses"]:
    v = pd.to_numeric(cu[k], errors="coerce")
    print("UCSD", k[:60], "non-null", int(v.notna().sum()), v.describe().round(0).to_dict() if v.notna().any() else cu[k].astype(str).value_counts().head(5).to_dict())
print("UCSD Avastin col values", cu[[k for k in cu.columns if "Avastin" in k][0]].astype(str).value_counts().head(6).to_dict())

# MU clinical columns
cm = pd.read_excel(t.CLIN, sheet_name="MU Glioma Post"); print("MU clin cols", [str(k).strip() for k in cm.columns])
