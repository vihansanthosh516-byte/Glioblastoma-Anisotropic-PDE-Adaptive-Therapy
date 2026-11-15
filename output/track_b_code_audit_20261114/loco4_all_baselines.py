"""LOCO over four cohorts: VIG against ALL six standard baselines (+ absolute-loss variants). Exploratory, post hoc.

Reads pooled4_pairs.csv written by fix_and_loco4.py. All four cohorts are spent; nothing here is a test.
Gradient boosting here uses the inputs that exist in all four cohorts (y_in, dt, earlier rate, has-earlier-scan).
"""
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(Path(__file__).resolve().parent))
import track_b_mu_tcks as t, track_b_vig_lumiere as vig
from fix_and_loco4 import fit_gomp_l1, gomp, ci

HERE = Path(__file__).resolve().parent
BASE = ["persistence", "two_point_exp", "pop_rate", "eb_random_slope", "gompertz", "gbr", "gompertz_l1", "gbr_l1"]


def feats(d):
    return pd.DataFrame({"y_in": d.y_in, "dt": d.dt_days, "r_hist": d.r_hist, "has_prior": d.r_hist.notna().astype(float)})


def main():
    p = pd.read_csv(HERE / "pooled4_pairs.csv")
    pats = {}
    for te_c in ["MU", "UCSD", "UCSF", "LUMIERE"]:
        tr, te = p[p.cohort != te_c].reset_index(drop=True), p[p.cohort == te_c].reset_index(drop=True)
        r_pop = float(tr.r_obs.median()); o = {}
        o["persistence"] = te.y_in.values
        o["two_point_exp"] = (te.y_in + te.r_hist.fillna(0.0) * te.dt_days).values
        o["pop_rate"] = (te.y_in + r_pop * te.dt_days).values
        t2, s2 = t.eb_params(tr, r_pop); o["eb_random_slope"] = np.asarray(t.shrink_pred(te, r_pop, t2, s2)[0])
        o["gompertz"] = t.gompertz_fit(tr)(te)
        o["gbr"] = HistGradientBoostingRegressor(random_state=0).fit(feats(tr), tr.y_out).predict(feats(te))
        o["gompertz_l1"] = gomp(fit_gomp_l1(tr), te.y_in.values, te.dt_days.values)
        o["gbr_l1"] = HistGradientBoostingRegressor(loss="absolute_error", random_state=0).fit(feats(tr), tr.y_out).predict(feats(te))
        th, lam = vig.fit_vig(tr)
        o["vig"] = te.y_in.values + vig.vig_curve(th, te.dt_days.values) - (lam[0] * np.maximum(te.lastch.values, 0) + lam[1] * np.minimum(te.lastch.values, 0))
        pats[te_c] = pd.DataFrame({k: np.abs(v - te.y_out.values) for k, v in o.items()}).assign(pid=te.patient_id.values).groupby("pid").mean()
        print(f"held out {te_c} (n={len(pats[te_c])}) patient MAE:", pats[te_c].mean().round(3).to_dict())
    allp = pd.concat(pats.values())
    print(f"\npooled patients: {len(allp)}. Pooled patient MAE:", allp.mean().round(3).to_dict())
    print("\nVIG minus each comparator, mean and 95% CI (negative = VIG better)")
    for b in BASE:
        print(f"{b:16s}" + "".join(f" | {c} {float((pt.vig - pt[b]).mean()):+.3f} {ci((pt.vig - pt[b]).values)}" for c, pt in list(pats.items()) + [("POOLED", allp)]))
    worst = {b: max(float((pt[b] - pt.persistence).mean()) for pt in pats.values()) for b in BASE + ["vig"]}
    print("\nworst cohort result of each method against persistence (positive = worse than persistence somewhere):", {k: round(v, 3) for k, v in worst.items()})
    # equal weight per cohort instead of per patient
    eq = {k: float(np.mean([pt[k].mean() for pt in pats.values()])) for k in BASE + ["vig"]}
    print("mean of the four cohort MAEs (each cohort equal weight):", {k: round(v, 3) for k, v in eq.items()})


if __name__ == "__main__":
    main()
