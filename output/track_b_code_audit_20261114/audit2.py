"""Code and data audit 2: UCSD, UCSF and LUMIERE pair builders, and the scoring code. Read-only.

All four cohorts are spent for this question; every number here is post hoc and descriptive.
"""
import io, sys, zipfile
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import track_b_mu_v3 as v3
import track_b_v8_ucsd_external as ext, track_b_v8_ucsf_alptdg as uf, track_b_vig_lumiere as vig

num = lambda s: pd.to_numeric(s, errors="coerce")
pd.set_option("display.width", 250)


def share(mask, err):
    return f"n={int(mask.sum())}, share of abs error={err[mask].sum() / err.sum():.3f}, mean |change|={err[mask].mean():.2f} vs {err[~mask].mean():.2f}"


def ucsd():
    print("== UCSD")
    c = pd.read_excel(ext.UCLIN, sheet_name="Cinical info"); c.columns = [str(x).strip() for x in c.columns]
    c["session"] = c["ID"].astype(str).str.strip(); c["pid"] = c.session.str.rsplit("_", n=1).str[0]; c["k"] = c.session.str.rsplit("_", n=1).str[1].astype(int)
    c["day"] = -num(c["Days from Acquisition to Date of initial surgery, treatment or diagnosis"])
    multi = c.groupby("pid").filter(lambda g: len(g) > 1).sort_values(["pid", "k"])
    bad = [p for p, g in multi.groupby("pid") if np.any(np.diff(g.day.values) <= 0) or g.day.isna().any()]
    print("sessions", len(c), "patients", c.pid.nunique(), "| patients with > 1 session", multi.pid.nunique(), "| consecutive pairs possible", len(multi) - multi.pid.nunique(),
          "| patients whose day does not increase with session number or is missing:", len(bad), bad[:6])
    u = ext.ucsd_pairs()
    print("eligible pairs", len(u), "patients", u.patient_id.nunique(), "| dt quantiles", u.dt_days.quantile([0, .5, 1]).to_dict())
    s = c.set_index("session")
    # which session is the outcome scan of each pair?
    idx = {(r.pid, r.day): r.session for r in c.itertuples()}
    u["sess_out"] = [idx.get((p, d)) for p, d in zip(u.patient_id, u.day_out)]
    sx = num(s["Days from Acquisition to Date of last surgery prior to scan if different"]); av = num(s["Days from Acquisition to Date of last Avastin treatment prior to scan"]); rt = num(s["Days from Acquistion to Date of last radiation prior to scan"])
    u["sx_between"] = [(-sx.get(k, np.nan)) < dt and (-sx.get(k, np.nan)) >= 0 for k, dt in zip(u.sess_out, u.dt_days)]
    u["av_between"] = [(-av.get(k, np.nan)) < dt and (-av.get(k, np.nan)) >= 0 for k, dt in zip(u.sess_out, u.dt_days)]
    u["rt_between"] = [(-rt.get(k, np.nan)) < dt and (-rt.get(k, np.nan)) >= 0 for k, dt in zip(u.sess_out, u.dt_days)]
    err = (u.y_out - u.y_in).abs()
    for k in ["sx_between", "av_between", "rt_between"]:
        print(f"   {k}: {share(u[k].values, err)} | median change {float((u.y_out - u.y_in)[u[k]].median()) if u[k].any() else None}")
    z = u.core_vox_out == 0
    print("   zero-outcome pairs:", int(z.sum()), "| with surgery between:", int((z & u.sx_between).sum()), "| with Avastin between:", int((z & u.av_between).sum()))


def ucsf():
    print("== UCSF")
    z = zipfile.ZipFile(uf.ZIP)
    x = pd.read_excel(io.BytesIO(z.read(uf.TABLE)), sheet_name="Clinical Info"); x.columns = [str(k).strip() for k in x.columns]
    dt = num(x["Days from 1st scan to 2nd scan"])
    lp = pd.read_excel(io.BytesIO(z.read(uf.TABLE)), sheet_name="UCSF-LPTDG")
    t1 = lp[lp.Timepoint == 1].set_index("SubjectID")
    core1 = num(t1["NCR_Volume_Label1"]).fillna(0) + num(t1["ET_Volume_Label3"]).fillna(0)
    x["core1_table"] = x.SubjectID.map(core1)
    print("patients", len(x), "| dt missing", int(dt.isna().sum()), "| dt < 14", int((dt < 14).sum()), "| dt > 365", int((dt > 365).sum()), "| table core volume at scan 1 == 0:", int((x.core1_table == 0).sum()))
    d, _ = uf.build(z)
    k = d.merge(x.assign(patient_id=x.SubjectID.astype(int).astype(str)), on="patient_id", how="left")
    ratio = (k.core_vox_in / k.core1_table.replace(0, np.nan))
    print("eligible", len(d), "| mask voxel count / table volume at scan 1: median", round(float(ratio.median()), 4), "| 5%-95%", ratio.quantile([.05, .95]).round(4).tolist(), "(1.0 = table is in mm3 and matches the mask)")


def lumiere():
    print("== LUMIERE (spent)")
    s = vig.lumiere_scans()
    r = pd.read_csv(vig.RATE_CSV); r.columns = ["Patient", "Time point", "lt3", "nm", "rating", "why"]
    r["rating"] = r.rating.astype(str).str.strip()
    dup = r.duplicated(["Patient", "Time point"], keep=False)
    print("rating rows", len(r), "| duplicated patient-time rows", int(dup.sum()), "| scans with mask and no rating row:", int((~pd.Series(list(zip(s["Patient"], s["Time point"]))).isin(set(zip(r["Patient"], r["Time point"])))).sum()))
    rate = r.drop_duplicates(["Patient", "Time point"]).set_index(["Patient", "Time point"]).rating
    s = s[~s.pre_op].sort_values(["Patient", "week", "sub"])
    rows = []
    for pid, g in s.groupby("Patient"):
        g = g.reset_index(drop=True)
        for i in range(len(g) - 1):
            a, b = g.iloc[i], g.iloc[i + 1]
            rows.append({"pid": pid, "dt": 7.0 * (b.week - a.week), "cin": a.core, "cout": b.core, "rate_in": rate.get((pid, a["Time point"]), "none"), "rate_out": rate.get((pid, b["Time point"]), "none")})
    d = pd.DataFrame(rows); d = d[(d.dt >= 14) & (d.dt <= 365) & (d.cin > 0)].reset_index(drop=True)
    d["ch"] = np.log(d.cout + 1) - np.log(d.cin + 1); err = d.ch.abs()
    print("eligible pairs", len(d), "| outcome-scan rating counts", d.rate_out.value_counts().to_dict())
    print("   input-scan rating counts", d.rate_in.value_counts().to_dict())
    po = d.rate_out.str.startswith("Post-Op").values
    print(f"   pairs that END in a post-operative scan (a surgery between the scans): {share(po, err)} | median change {float(d.ch[po].median()):.2f}")
    pi = d.rate_in.str.startswith("Post-Op").values
    print(f"   pairs that START at a post-operative scan: {share(pi, err)} | median change {float(d.ch[pi].median()):.2f}")
    print("   median change by outcome rating", d.groupby("rate_out").ch.median().round(2).to_dict())
    z = d.cout == 0
    print("   zero-outcome pairs", int(z.sum()), "| of those ending in a post-operative scan", int((z & po).sum()), "| ending in CR", int((z & (d.rate_out == 'CR').values).sum()))
    pm = lambda e: float(pd.Series(e).groupby(d.pid.values).mean().mean())
    print("   persistence patient MAE: all", round(pm(err.values), 3), "| excluding pairs that end in a post-operative scan", round(float(err[~po].groupby(d.pid[~po].values).mean().mean()), 3))


def scoring():
    print("== scoring code check (toy data with a known answer)")
    rng = np.random.default_rng(0)
    diff = rng.normal(-0.1, 1.0, 400)
    ci = v3.boot_ci(diff, np.random.default_rng(1), 0.95)
    se = diff.std(ddof=1) / np.sqrt(len(diff))
    print("   boot_ci", [round(x, 3) for x in ci], "| normal-theory CI", [round(diff.mean() - 1.96 * se, 3), round(diff.mean() + 1.96 * se, 3)])
    e = pd.DataFrame({"m": [1.0, 3.0, 5.0], "pid": ["a", "a", "b"]}).groupby("pid").mean()
    print("   patient-mean of errors [1,3 | 5] =", float(e.m.mean()), "(expected 3.5)")


if __name__ == "__main__":
    ucsd(); ucsf(); lumiere(); scoring()
