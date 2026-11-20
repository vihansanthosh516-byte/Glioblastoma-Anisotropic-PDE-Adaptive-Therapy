"""Track B: VIG-2 development, leave-one-cohort-out over six cohorts. Pre-registration: paper/trackB_vig2_prereg_2026-11-19.md (committed first).
Development only. Not a test. Candidates H1, H2, H3 and references are exactly those of the pre-registration.
"""
import json, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import minimize

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parent.parent
AUDIT = ROOT / "output" / "track_b_code_audit_20261114"; OUT = ROOT / "output" / "track_b_vig2_dev_20261119"
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(AUDIT))
import track_b_vig_lumiere as vig
from fix_and_loco4 import fit_gomp_l1, gomp, fit_pk

SEED, NB = 20261119, 20000
COH = ["MU", "UCSD", "UCSF", "LUMIERE", "CFB", "WAW"]


def load():
    a = pd.read_csv(AUDIT / "pooled4_pairs.csv"); b = pd.read_csv(ROOT / "output" / "track_b_vig_cfb_waw_20261115" / "test_pairs.csv")
    d = pd.concat([a, b[a.columns]], ignore_index=True); assert set(d.cohort) == set(COH)
    return d


def lc(d): return (d.y_in - d.y_prev).fillna(0.0).values


def Kvec(d, th_K, peak):
    if not peak: return np.full(len(d), th_K)
    pk = d.y_peak_prev.values
    return np.where(np.isnan(pk), th_K, np.maximum(pk, d.y_in.values))


def hist_part(d, lam): return lam[0] * np.maximum(lc(d), 0) + lam[1] * np.minimum(lc(d), 0)


def change(d, th, kind, peak=False):
    y, dt = d.y_in.values, d.dt_days.values
    g = vig.vig_curve(th[:3], dt) if kind in ("H1", "H3") else 0.0
    K = Kvec(d, th[3], peak); a = th[4]
    return g + (1 - np.exp(-a * dt)) * (K - y)


def fit_h(tr, kind):
    peak = kind == "H3"; y, dt, ch = tr.y_in.values, tr.dt_days.values, (tr.y_out - tr.y_in).values
    def f(p):
        th = np.r_[p[:3], p[3], np.exp(p[4])]
        return float(np.mean(np.abs(ch - change(tr, th, kind, peak))))
    starts = [[np.log(0.01), 0.0, np.log(60), 9.4, np.log(1.6e-3)], [np.log(0.02), -0.1, np.log(30), 9.0, np.log(1e-3)], [np.log(0.005), 0.0, np.log(120), 10.0, np.log(3e-3)]]
    if kind == "H2": starts = [[0, 0, 0] + s[3:] for s in starts]
    best = None
    for s in starts:
        if kind == "H2":
            r = minimize(lambda q: f(np.r_[0, 0, 0, q]), s[3:], method="Nelder-Mead", options={"xatol": 1e-6, "fatol": 1e-9, "maxiter": 6000}); r.x = np.r_[0, 0, 0, r.x]
        else:
            r = minimize(f, s, method="Nelder-Mead", options={"xatol": 1e-6, "fatol": 1e-9, "maxiter": 8000})
        if best is None or r.fun < best.fun: best = r
    p = best.x; th = np.r_[p[:3], p[3], np.exp(p[4])]
    has = tr.r_hist.notna().values; res = (ch - change(tr, th, kind, peak))[has]; l = lc(tr)[has]
    X = np.c_[np.maximum(l, 0.0), np.minimum(l, 0.0)]; grid, lam = np.linspace(0.0, 0.6, 31), np.zeros(2)
    for _ in range(4):
        for j in range(2):
            c = []
            for gg in grid:
                b = lam.copy(); b[j] = gg; c.append(np.mean(np.abs(res + X @ b)))
            lam[j] = grid[int(np.argmin(c))]
    return th, lam


def predict_fold(tr, te):
    out = {"persistence": te.y_in.values, "pop_rate": te.y_in.values + float(tr.r_obs.median()) * te.dt_days.values}
    out["gomp_l1"] = gomp(fit_gomp_l1(tr), te.y_in.values, te.dt_days.values)
    a = fit_pk(tr); K = np.maximum(te.y_peak_prev.fillna(-np.inf).values, te.y_in.values); out["pk_gomp"] = K - (K - te.y_in.values) * np.exp(-a * te.dt_days.values)
    th, l2 = vig.fit_vig(tr); out["vig"] = te.y_in.values + vig.vig_curve(th, te.dt_days.values) - hist_part(te, l2)
    pars = {"vig": {"th": [float(x) for x in th], "lam": [float(x) for x in l2]}}
    for k in ("H1", "H2", "H3"):
        th, lam = fit_h(tr, k)
        out[k] = te.y_in.values + change(te, th, k, k == "H3") - hist_part(te, lam)
        pars[k] = {"th": [float(x) for x in th], "lam": [float(x) for x in lam]}
    return out, pars


def pat_err(te, out):
    e = pd.DataFrame({k: np.abs(v - te.y_out.values) for k, v in out.items()}); e["pid"] = te.patient_id.values
    return e.groupby("pid").mean()


def boot(diff, rng):
    x = np.asarray(diff); idx = rng.integers(0, len(x), (NB, len(x))); m = x[idx].mean(1)
    return float(x.mean()), [float(np.quantile(m, .025)), float(np.quantile(m, .975))]


def main():
    d = load(); rng = np.random.default_rng(SEED); res = {}; allpars = {}
    for c in COH:
        tr, te = d[d.cohort != c].reset_index(drop=True), d[d.cohort == c].reset_index(drop=True)
        out, pars = predict_fold(tr, te); e = pat_err(te, out); allpars[c] = pars
        res[c] = {"patients": int(len(e)), "mae": {k: float(e[k].mean()) for k in e}}
        res[c]["diff_vs_persistence"] = {k: dict(zip(("mean", "ci95"), boot(e[k] - e["persistence"], rng))) for k in e if k != "persistence"}
        print(c, len(e), {k: round(v, 3) for k, v in res[c]["mae"].items()}, flush=True)
    meths = [k for k in res[COH[0]]["mae"] if k != "persistence"]
    summ = {}
    for k in meths:
        dif = {c: res[c]["diff_vs_persistence"][k]["mean"] for c in COH}
        summ[k] = {"cohort_equal_mean_6": float(np.mean(list(dif.values()))), "n_better_than_persistence_of_6": int(sum(v < 0 for v in dif.values())),
                   "cohort_equal_mean_first4": float(np.mean([dif[c] for c in COH[:4]])), "per_cohort": dif}
    cand = [k for k in ("H1", "H2", "H3") if summ[k]["n_better_than_persistence_of_6"] >= 4 and summ[k]["cohort_equal_mean_6"] < summ["vig"]["cohort_equal_mean_6"] and summ[k]["cohort_equal_mean_6"] < summ["gomp_l1"]["cohort_equal_mean_6"]]
    pick = min(cand, key=lambda k: summ[k]["cohort_equal_mean_6"]) if cand else None
    json.dump({"folds": res, "summary": summ, "candidates_passing_rule": cand, "selected": pick, "fold_parameters": allpars, "seed": SEED, "n_boot": NB,
               "prereg": "paper/trackB_vig2_prereg_2026-11-19.md"}, open(OUT / "dev_result.json", "w"), indent=1)
    print(json.dumps({k: {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != "per_cohort"} for k, v in summ.items()}, indent=1)); print("passing rule:", cand, "selected:", pick)


if __name__ == "__main__":
    main()
