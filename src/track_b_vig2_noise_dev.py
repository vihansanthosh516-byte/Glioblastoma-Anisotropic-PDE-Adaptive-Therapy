"""Track B: noise-aware level-reversion, leave-one-cohort-out over six cohorts. Pre-registration: paper/trackB_vig2_noise_prereg_2026-11-19.md (committed first).
Development only. Candidates N0, N1, N2 and references exactly as pre-registered.
"""
import json, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import minimize

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
import track_b_vig2_dev as dv
from track_b_vig2_dev import vig, fit_gomp_l1, gomp, load, lc, hist_part, pat_err, boot, COH, SEED
OUT = ROOT / "output" / "track_b_vig3_noise_dev_20261119"
sig = lambda x: 1.0 / (1.0 + np.exp(-x))


def parts(d, p, kind, lo, hi):
    y, dt = d.y_in.values, d.dt_days.values; has = d.r_hist.notna().values.astype(float)
    K = lo + (hi - lo) * sig(p[3]); a0 = np.exp(p[4])
    if kind == "N0": a = a0
    else:
        a = a0 * sig((p[5] - y) / np.exp(p[6]))
        if kind == "N2": a = a * np.exp(-np.exp(p[7]) * has)
    return vig.vig_curve(p[:3], dt) + (1 - np.exp(-a * dt)) * (K - y)


def fit_n(tr, kind):
    lo, hi = float(tr.y_in.median()), float(tr.y_in.max()); ch = (tr.y_out - tr.y_in).values
    f = lambda p: float(np.mean(np.abs(ch - parts(tr, p, kind, lo, hi))))
    base = [[np.log(0.01), 0.0, np.log(60), 0.0, np.log(2e-3)], [np.log(0.02), -0.1, np.log(30), -1.0, np.log(1e-3)], [np.log(0.005), 0.0, np.log(120), 1.0, np.log(4e-3)]]
    ext = {"N0": [[], [], []], "N1": [[10.0, 0.0], [9.0, np.log(2.0)], [11.0, np.log(0.5)]], "N2": [[10.0, 0.0, np.log(0.5)], [9.0, np.log(2.0), np.log(0.3)], [11.0, np.log(0.5), np.log(1.0)]]}[kind]
    best = min((minimize(f, np.r_[b, e], method="Nelder-Mead", options={"xatol": 1e-6, "fatol": 1e-9, "maxiter": 12000}) for b, e in zip(base, ext)), key=lambda r: r.fun)
    p = best.x; has = tr.r_hist.notna().values; res = (ch - parts(tr, p, kind, lo, hi))[has]; l = lc(tr)[has]
    X = np.c_[np.maximum(l, 0.0), np.minimum(l, 0.0)]; grid, lam = np.linspace(0.0, 0.6, 31), np.zeros(2)
    for _ in range(4):
        for j in range(2):
            c = []
            for g in grid:
                b = lam.copy(); b[j] = g; c.append(np.mean(np.abs(res + X @ b)))
            lam[j] = grid[int(np.argmin(c))]
    return p, lam, lo, hi


def predict_fold(tr, te):
    out = {"persistence": te.y_in.values}
    out["gomp_l1"] = gomp(fit_gomp_l1(tr), te.y_in.values, te.dt_days.values)
    th, l2 = vig.fit_vig(tr); out["vig"] = te.y_in.values + vig.vig_curve(th, te.dt_days.values) - hist_part(te, l2)
    th, lam = dv.fit_h(tr, "H1"); out["H1"] = te.y_in.values + dv.change(te, th, "H1") - hist_part(te, lam)
    pars = {}
    for k in ("N0", "N1", "N2"):
        p, lam, lo, hi = fit_n(tr, k); out[k] = te.y_in.values + parts(te, p, k, lo, hi) - hist_part(te, lam)
        pars[k] = {"p": [float(x) for x in p], "lam": [float(x) for x in lam], "K_range": [lo, hi], "K": float(lo + (hi - lo) * sig(p[3]))}
    return out, pars


def main():
    d = load(); rng = np.random.default_rng(SEED); res, allp = {}, {}
    for c in COH:
        tr, te = d[d.cohort != c].reset_index(drop=True), d[d.cohort == c].reset_index(drop=True)
        out, pars = predict_fold(tr, te); e = pat_err(te, out); allp[c] = pars
        res[c] = {"patients": int(len(e)), "mae": {k: float(e[k].mean()) for k in e}, "diff_vs_persistence": {k: dict(zip(("mean", "ci95"), boot(e[k] - e["persistence"], rng))) for k in e if k != "persistence"}}
        print(c, len(e), {k: round(v, 3) for k, v in res[c]["mae"].items()}, flush=True)
    meths = [k for k in res[COH[0]]["mae"] if k != "persistence"]; summ = {}
    for k in meths:
        dif = {c: res[c]["diff_vs_persistence"][k]["mean"] for c in COH}
        summ[k] = {"mean6": float(np.mean(list(dif.values()))), "mean_first4": float(np.mean([dif[c] for c in COH[:4]])), "n_better_of_6": int(sum(v < 0 for v in dif.values())),
                   "worst_cohort": float(max(dif.values())), "per_cohort": dif}
    ok = [k for k in ("N1", "N2") if summ[k]["mean6"] < summ["vig"]["mean6"] and summ[k]["mean6"] < summ["H1"]["mean6"] and summ[k]["mean_first4"] <= summ["vig"]["mean_first4"]
          and summ[k]["worst_cohort"] <= 0.03 and summ[k]["n_better_of_6"] >= 5]
    pick = min(ok, key=lambda k: summ[k]["mean6"]) if ok else None
    json.dump({"folds": res, "summary": summ, "passing_rule": ok, "selected": pick, "fold_parameters": allp, "seed": SEED, "prereg": "paper/trackB_vig2_noise_prereg_2026-11-19.md"}, open(OUT / "noise_dev_result.json", "w"), indent=1)
    print(json.dumps({k: {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != "per_cohort"} for k, v in summ.items()}, indent=1)); print("per cohort:", {k: {c: round(x, 3) for c, x in summ[k]["per_cohort"].items()} for k in ("vig", "H1", "N0", "N1", "N2")}); print("passing rule:", ok, "selected:", pick)


if __name__ == "__main__":
    main()
