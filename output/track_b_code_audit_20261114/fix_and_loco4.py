"""Fixes from the audit + leave-one-cohort-out over FOUR cohorts. Exploratory, post hoc.

All four cohorts (MU, UCSD, UCSF, LUMIERE) are spent for this question. Nothing here is a test.

Fix 1 (UCSD): the locked builder drops a session when "days to initial surgery" is "Unknown".
       Here the scan interval is taken from ANY event that is fixed in time (initial surgery,
       death, last follow-up, first radiation, first chemo): dt = offset_a - offset_b.
Fix 2 (LUMIERE): an unrated scan labelled week-000-1 is treated as pre-operative.

Candidates (all fall back to persistence when a patient has no earlier scan, except where noted):
  rebound   : y_in - lam_fall * min(last change, 0)
  pk_gomp   : Gompertz pull toward the patient's own earlier peak size K_i = max(earlier sizes, y_in)
  vig       : the locked VIG (interval curve + rise/fall history terms)
  gomp_l1   : global Gompertz, absolute loss (a baseline; needs no history)
  blend     : mean of vig and gomp_l1
"""
import glob, os, sys, warnings, zipfile
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import minimize

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import track_b_mu_tcks as t
import track_b_v8_ucsd_external as ext, track_b_v8_ucsf_alptdg as uf, track_b_vig_lumiere as vig

HERE = Path(__file__).resolve().parent
num = lambda s: pd.to_numeric(s, errors="coerce")
FIXED = ["Days from Acquisition to Date of initial surgery, treatment or diagnosis", "Days from Acquisition to Date of death",
         "Days from Acquisition to Date of last follow-up", "Days from Acquistion to Date of first radiation", "Days from Acquisition to Date of 1st chemo start"]


def series(sizes, days, pid, cohort):
    """consecutive pairs with full earlier history (sizes are core voxels, days ascending)."""
    rows = []
    for i in range(len(sizes) - 1):
        dt = days[i + 1] - days[i]
        prev = [(days[j], sizes[j]) for j in range(i) if sizes[j] > 0 and days[j] < days[i]]
        r = {"cohort": cohort, "patient_id": f"{cohort}_{pid}", "dt_days": dt, "core_vox_in": sizes[i], "core_vox_out": sizes[i + 1],
             "n_prev": len(prev), "y_prev": np.log(prev[-1][1] + 1) if prev else np.nan, "dt_prior": days[i] - prev[-1][0] if prev else np.nan,
             "y_peak_prev": max(np.log(s + 1) for _, s in prev) if prev else np.nan}
        rows.append(r)
    return rows


def finish(d):
    d = d[(d.dt_days >= 14) & (d.dt_days <= 365) & (d.core_vox_in > 0)].reset_index(drop=True)
    d["y_in"], d["y_out"] = np.log(d.core_vox_in + 1.0), np.log(d.core_vox_out + 1.0)
    d["lastch"] = (d.y_in - d.y_prev).fillna(0.0)
    d["r_hist"] = (d.y_in - d.y_prev) / d.dt_prior
    d["r_obs"] = (d.y_out - d.y_in) / d.dt_days
    return d


def mu_pairs():
    c = pd.read_csv(HERE / "mu_mask_counts.csv")
    c = c[c.has_mask].copy(); c["core"] = c.n1 + c.n3
    cl = pd.read_excel(t.CLIN, sheet_name="MU Glioma Post"); cl.columns = [str(k).strip() for k in cl.columns]
    rows = []
    for _, r in cl.iterrows():
        pid = str(r["Patient_ID"]).strip(); g = c[c.patient_id == pid]
        pts = sorted((float(num(pd.Series([r[next(x for x in cl.columns if f"(Timepoint_{int(tp)})" in x)]]))[0]), float(core)) for tp, core in zip(g.tp, g.core))
        pts = [p for p in pts if not np.isnan(p[0])]
        rows += series([p[1] for p in pts], [p[0] for p in pts], pid, "MU")
    return finish(pd.DataFrame(rows))


def ucsd_pairs(fix=True):
    c = pd.read_excel(ext.UCLIN, sheet_name="Cinical info"); c.columns = [str(x).strip() for x in c.columns]
    c["session"] = c["ID"].astype(str).str.strip(); c["pid"] = c.session.str.rsplit("_", n=1).str[0]; c["k"] = c.session.str.rsplit("_", n=1).str[1].astype(int)
    segs = {}
    for p in sorted(glob.glob(str(ext.UCSD / "**" / "*_BraTS_tumor_seg.nii.gz"), recursive=True)):
        segs.setdefault(os.path.basename(p).split("_BraTS")[0], p)
    c = c[c.session.isin(segs)].sort_values(["pid", "k"])
    rows, disagree, n_rec = [], [], 0
    for pid, g in c.groupby("pid"):
        if len(g) < 2: continue
        off = g[FIXED].apply(num).values                       # sessions x anchors
        cols = [0] if not fix else list(range(len(FIXED)))
        days = [0.0]
        for i in range(1, len(g)):
            d = off[i - 1, cols] - off[i, cols]; d = d[~np.isnan(d)]
            if len(d) == 0: days = None; break
            if len(d) > 1: disagree.append(float(d.max() - d.min()))
            if np.isnan(off[i - 1, 0] - off[i, 0]): n_rec += 1
            days.append(days[-1] + float(np.median(d)))
        if days is None or np.any(np.diff(days) <= 0): continue
        sizes = [float(ext.comp_from_mask(segs[s])[0]) for s in g.session]
        rows += series(sizes, days, pid, "UCSD")
    d = finish(pd.DataFrame(rows))
    return d, {"intervals_recovered_from_other_anchors": n_rec, "anchor_disagreement_days_median": float(np.median(disagree)) if disagree else None,
               "anchor_disagreement_days_max": float(np.max(disagree)) if disagree else None, "n_checked": len(disagree)}


def ucsf_pairs():
    f, _ = uf.build(zipfile.ZipFile(uf.ZIP))
    f = f.assign(cohort="UCSF", patient_id="UCSF_" + f.patient_id.astype(str), n_prev=0, y_prev=np.nan, dt_prior=np.nan, y_peak_prev=np.nan)
    return finish(f[["cohort", "patient_id", "dt_days", "core_vox_in", "core_vox_out", "n_prev", "y_prev", "dt_prior", "y_peak_prev"]])


def lumiere_pairs():
    s = vig.lumiere_scans()
    r = pd.read_csv(vig.RATE_CSV); keys = set(zip(r.iloc[:, 0], r.iloc[:, 1]))
    unrated_pre = [(a, b) not in keys and b == "week-000-1" for a, b in zip(s["Patient"], s["Time point"])]
    s = s[~(s.pre_op | np.array(unrated_pre))].sort_values(["Patient", "week", "sub"])
    rows = []
    for pid, g in s.groupby("Patient"):
        g = g.drop_duplicates("week", keep="last")                 # same-week scans: keep the later one
        rows += series(g.core.astype(float).tolist(), (7.0 * g.week).tolist(), pid, "LUMIERE")
    return finish(pd.DataFrame(rows))


# ---- candidates -----------------------------------------------------------------------------
def gomp(th, y, dt):
    return th[1] - (th[1] - y) * np.exp(-th[0] * dt)


def fit_gomp_l1(tr):
    loss = lambda th: float(np.mean(np.abs(gomp(th, tr.y_in.values, tr.dt_days.values) - tr.y_out.values)))
    return min((minimize(loss, x0, method="Nelder-Mead", options={"xatol": 1e-7, "fatol": 1e-9, "maxiter": 10000}) for x0 in ([1e-3, tr.y_in.max() + 0.5], [1e-4, 9.0], [5e-3, 12.0])), key=lambda r: r.fun).x


def fit_rebound(tr):
    h = tr[tr.n_prev > 0]; res, lc = (h.y_out - h.y_in).values, np.minimum(h.lastch.values, 0.0)
    grid = np.linspace(0, 0.8, 41)
    return float(grid[int(np.argmin([np.mean(np.abs(res + g * lc)) for g in grid]))])


def fit_pk(tr):
    """patient-peak Gompertz: y_hat = K_i - (K_i - y_in) * exp(-a * dt), K_i = max(earlier peak, y_in); one rate a >= 0."""
    h = tr[tr.n_prev > 0]; K = np.maximum(h.y_peak_prev.values, h.y_in.values)
    grid = np.concatenate([[0.0], np.exp(np.linspace(np.log(1e-4), np.log(0.05), 60))])
    return float(grid[int(np.argmin([np.mean(np.abs(K - (K - h.y_in.values) * np.exp(-a * h.dt_days.values) - h.y_out.values)) for a in grid]))])


def predict_all(tr, te):
    out = {"persistence": te.y_in.values, "pop_rate": te.y_in.values + float(tr.r_obs.median()) * te.dt_days.values}
    out["gompertz"] = t.gompertz_fit(tr)(te)
    g1 = fit_gomp_l1(tr); out["gomp_l1"] = gomp(g1, te.y_in.values, te.dt_days.values)
    lam = fit_rebound(tr); out["rebound"] = te.y_in.values - lam * np.minimum(te.lastch.values, 0.0)
    a = fit_pk(tr); K = np.maximum(te.y_peak_prev.fillna(-np.inf).values, te.y_in.values)
    out["pk_gomp"] = K - (K - te.y_in.values) * np.exp(-a * te.dt_days.values)
    th, l2 = vig.fit_vig(tr)
    out["vig"] = te.y_in.values + vig.vig_curve(th, te.dt_days.values) - (l2[0] * np.maximum(te.lastch.values, 0) + l2[1] * np.minimum(te.lastch.values, 0))
    out["blend"] = 0.5 * (out["vig"] + out["gomp_l1"])
    return out, {"lam_rebound": lam, "pk_rate": a, "gomp_l1": [float(x) for x in g1]}


def ci(diff, n=20000, seed=20261114):
    rng = np.random.default_rng(seed); b = diff[rng.integers(0, len(diff), (n, len(diff)))].mean(1)
    return f"({np.percentile(b, 2.5):+.3f}, {np.percentile(b, 97.5):+.3f})"


def main():
    old, _ = ucsd_pairs(fix=False); new, info = ucsd_pairs(fix=True)
    print("Fix 1 UCSD: pairs", len(old), "->", len(new), "| patients", old.patient_id.nunique(), "->", new.patient_id.nunique(), "|", info)
    lum = lumiere_pairs(); mu = mu_pairs(); f = ucsf_pairs()
    print("Fix 2 LUMIERE: pairs", len(lum), "patients", lum.patient_id.nunique(), "(locked run: 400 pairs, 79 patients)")
    print("MU rebuilt from raw files: pairs", len(mu), "patients", mu.patient_id.nunique(), "(pair table: 334, 135)")
    p = pd.concat([mu, new, f, lum], ignore_index=True)
    p.to_csv(HERE / "pooled4_pairs.csv", index=False)
    print("pairs with an earlier scan:", p.groupby("cohort").apply(lambda g: f"{int((g.n_prev > 0).sum())}/{len(g)}").to_dict())
    cohorts = ["MU", "UCSD", "UCSF", "LUMIERE"]; pats = {}
    for te_c in cohorts:
        tr, te = p[p.cohort != te_c].reset_index(drop=True), p[p.cohort == te_c].reset_index(drop=True)
        preds, par = predict_all(tr, te)
        pats[te_c] = pd.DataFrame({k: np.abs(v - te.y_out.values) for k, v in preds.items()}).assign(pid=te.patient_id.values).groupby("pid").mean()
        print(f"\n== held out {te_c}: patients {len(pats[te_c])} | fitted on the other three: {par}")
        print("   patient MAE", pats[te_c].mean().round(3).to_dict())
    allp = pd.concat(pats.values())
    names = ["rebound", "pk_gomp", "vig", "blend", "gomp_l1", "gompertz", "pop_rate"]
    print("\nLOCO-4: method minus persistence, mean and 95% CI (negative = better)")
    for k in names:
        print(f"{k:9s}" + "".join(f" | {c} {float((pt[k] - pt.persistence).mean()):+.3f} {ci((pt[k] - pt.persistence).values)}" for c, pt in list(pats.items()) + [("POOLED", allp)]))
    print("\nLOCO-4: method minus the best baseline of that cohort (best of persistence, pop_rate, gompertz, gomp_l1)")
    for k in ["rebound", "pk_gomp", "vig", "blend"]:
        line = f"{k:9s}"
        for c, pt in pats.items():
            b = min(["persistence", "pop_rate", "gompertz", "gomp_l1"], key=lambda x: pt[x].mean())
            line += f" | {c} vs {b} {float((pt[k] - pt[b]).mean()):+.3f} {ci((pt[k] - pt[b]).values)}"
        print(line)


if __name__ == "__main__":
    main()
