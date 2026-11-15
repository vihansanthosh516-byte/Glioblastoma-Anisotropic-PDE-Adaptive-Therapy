"""Read-only post hoc diagnostics, part 3. Writes nothing."""
import io, sys, zipfile
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(sys.argv[1]); sys.path.insert(0, str(ROOT / "src"))
import track_b_mu_tcks as t, track_b_v8_ucsf_alptdg as uf
num = lambda s: pd.to_numeric(s, errors="coerce")

d, _ = t.load()
c = pd.read_excel(t.CLIN, sheet_name="MU Glioma Post"); c.columns = [str(k).strip() for k in c.columns]
print("MU clin cols tail", list(c.columns)[60:])
c["patient_id"] = c["Patient_ID"].astype(str).str.strip()
for k in ["Multiple surgeries", "Additional Therapy", "2nd_Additional Therapy", "Other Types of Therapy (LITT, more chemo, proton therapy)", "Treatment started after 2nd progression", "Immuno therapy", "Brachy therapy", "Name of Initial Chemo Therapy"]:
    if k in c: print("MU", k, c[k].astype(str).str.strip().str[:30].value_counts(dropna=False).head(7).to_dict())
starts = [k for k in c.columns if ("start" in k.lower() or "insertion" in k.lower() or "new treatment" in k.lower() or "further progression" in k.lower() or "first progression" in k.lower()) and "Radiation" not in k and "Initial Chemo" not in k]
print("event-date cols used:", starts)
j = d.merge(c[["patient_id"] + starts], on="patient_id", how="left")
E = j[starts].apply(num)
between = np.zeros(len(j), bool); before = np.zeros(len(j), bool)
for k in starts:
    between |= ((E[k] > j.day_in) & (E[k] <= j.day_out)).values
    before |= ((E[k] <= j.day_in) & (E[k] > j.day_in - 60)).values
col = (j.core_vox_out < 0.05 * j.core_vox_in).values
zero = (j.core_vox_out == 0).values
err = (j.y_out - j.y_in).abs().values
print("MU pairs", len(j), "collapse(>95% drop)", int(col.sum()), "zero-out", int(zero.sum()))
for name, m in [("event BETWEEN scans (future info)", between), ("event in 60 d BEFORE input scan", before)]:
    print(name, ": n", int(m.sum()), "| collapse rate in", round(col[m].mean(), 3), "out", round(col[~m].mean(), 3), "| collapse pairs covered", int((col & m).sum()), "of", int(col.sum()),
          "| mean pers err in", round(err[m].mean(), 3), "out", round(err[~m].mean(), 3))
print("MU collapse by phase", pd.crosstab(j.phase, col).to_dict())
print("MU collapse: median y_in collapse", round(float(np.median(j.y_in[col])), 2), "others", round(float(np.median(j.y_in[~col])), 2), "| median dt collapse", float(np.median(j.dt_days[col])), "others", float(np.median(j.dt_days[~col])))
# jumps up
up = (j.y_out - j.y_in).values > 3
print("MU jumps up >3 log:", int(up.sum()), "median y_in", round(float(np.median(j.y_in[up])), 2) if up.any() else None)
# label 4 (cavity) info is not in manifest; skip

# noise-vs-size on MU: abs change by input size tertile, short intervals
q = pd.qcut(j.core_vox_in, 3, labels=["small", "mid", "large"])
print("MU median |dy| by input size tertile", pd.Series(err).groupby(q.values, observed=True).median().round(3).to_dict(), "mean", pd.Series(err).groupby(q.values, observed=True).mean().round(3).to_dict())
print("MU signed mean change by tertile (regression to mean check)", pd.Series((j.y_out - j.y_in).values).groupby(q.values, observed=True).mean().round(3).to_dict())

# UCSF: opposite sign convention for 2nd/3rd surgery, corrupt dates
z = zipfile.ZipFile(uf.ZIP)
u, _ = uf.build(z)
x = pd.read_excel(io.BytesIO(z.read(uf.TABLE)), sheet_name="Clinical Info"); x.columns = [str(k).strip() for k in x.columns]
x["patient_id"] = x.SubjectID.astype(int).astype(str)
k = u.merge(x, on="patient_id", how="left")
s2, s3 = num(k["Days from 1st scan to 2nd SX"]), num(k["Days from 1st scan to 3rd SX"])
btw_neg = ((s2 < 0) & (s2 >= -k.dt_days)) | ((s3 < 0) & (s3 >= -k.dt_days))
zero = k.core_vox_out == 0
e = (k.y_out - k.y_in).abs()
print("UCSF alt sign (negative = after scan 1): surgery between n", int(btw_neg.sum()), "zero-out covered", int((zero & btw_neg).sum()), "of", int(zero.sum()), "mean err in", round(e[btw_neg].mean(), 3) if btw_neg.any() else None)
nsx = num(k["Number of surgeries"])
print("UCSF n_surg>=2:", int((nsx >= 2).sum()), "| 2nd SX date present", int(s2.notna().sum()), "| of those positive", int((s2 > 0).sum()), "negative", int((s2 < 0).sum()))
datecols = ["Days from 1st surgery/DX  to 1st scan", "Days from 1st scan to 1st RT start (neg = RT first)", "Days from 1st chemo start to 1st scan", "Days from 1st scan to 2nd SX", "Days from 1st scan to 3rd SX", "Days from 1st scan to next progression"]
allx = x[[q_ for q_ in datecols if q_ in x]].apply(num)
print("UCSF table rows with any |days| > 15000 (looks like raw date serial):", int((allx.abs() > 15000).any(axis=1).sum()), "of", len(x), "| in eligible:", int((k[[q_ for q_ in datecols if q_ in k]].apply(num).abs() > 15000).any(axis=1).sum()))
on1 = k["On tx at 1st scan (c = chemo, rt, n = none)"].astype(str).str.strip()
print("UCSF code phase vs dataset on-treatment flag at scan 1:", pd.crosstab(k.phase, on1).to_dict())
print("UCSF zero-out pairs: y_in", k.y_in[zero].round(2).tolist(), "dt", k.dt_days[zero].tolist(), "on tx", on1[zero].tolist())
prog = num(k["Days from 1st scan to next progression"])
pb = (prog > 0) & (prog <= k.dt_days)
print("UCSF progression called between scans: n", int(pb.sum()), "mean signed change", round((k.y_out - k.y_in)[pb].mean(), 3), "others", round((k.y_out - k.y_in)[~pb].mean(), 3))
