"""Track B image-feature development, leave-one-cohort-out over six cohorts. Pre-registration: paper/trackB_imgfeat_prereg_2026-11-19.md (committed first).
Development only. References GBR-std, GBR-img; candidates VF1, VF2; VIG refitted each fold."""
import json, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import QuantileRegressor

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from track_b_vig2_dev import vig, hist_part, pat_err, boot, COH, SEED
OUT = ROOT / "output" / "track_b_imgfeat_dev_20261119"
IMG = ["log_edema", "nec_frac", "enh_frac", "log_cavity"]
STD = lambda d: pd.DataFrame({"y_in": d.y_in, "dt": d.dt_days, "r_hist": d.r_hist, "has_prior": d.r_hist.notna().astype(float)})
FULL = lambda d: pd.concat([STD(d), d[IMG].reset_index(drop=True)], axis=1)


def vig_pred(th, lam, d): return d.y_in.values + vig.vig_curve(th, d.dt_days.values) - hist_part(d, lam)


def predict_fold(tr, te):
    out = {"persistence": te.y_in.values}; ch = (tr.y_out - tr.y_in).values
    th, lam = vig.fit_vig(tr); out["vig"] = vig_pred(th, lam, te); res = tr.y_out.values - vig_pred(th, lam, tr)
    for name, F in (("GBR_std", STD), ("GBR_img", FULL)):
        out[name] = te.y_in.values + HistGradientBoostingRegressor(loss="absolute_error", random_state=0).fit(F(tr), ch).predict(F(te))
    Z = lambda d: d[["y_in"] + IMG].reset_index(drop=True)
    gb = HistGradientBoostingRegressor(loss="absolute_error", max_depth=2, max_iter=100, learning_rate=0.05, min_samples_leaf=40, random_state=0).fit(Z(tr), res)
    out["VF1"] = vig_pred(th, lam, te) + gb.predict(Z(te))
    med = Z(tr).median(); mu, sd = Z(tr).fillna(med).mean(), Z(tr).fillna(med).std().replace(0, 1)
    P = lambda d: ((Z(d).fillna(med) - mu) / sd).values
    qr = QuantileRegressor(quantile=0.5, alpha=1e-3, solver="highs").fit(P(tr), res)
    out["VF2"] = vig_pred(th, lam, te) + qr.predict(P(te))
    return out, {"vig": {"th": [float(x) for x in th], "lam": [float(x) for x in lam]}, "VF2_coef": dict(zip(["y_in"] + IMG, [float(x) for x in qr.coef_]))}


def main():
    d = pd.read_csv(OUT / "pairs_with_features.csv"); rng = np.random.default_rng(SEED); res, allp = {}, {}
    for c in COH:
        tr, te = d[d.cohort != c].reset_index(drop=True), d[d.cohort == c].reset_index(drop=True)
        out, pars = predict_fold(tr, te); e = pat_err(te, out); allp[c] = pars
        res[c] = {"patients": int(len(e)), "mae": {k: float(e[k].mean()) for k in e}, "diff_vs_persistence": {k: dict(zip(("mean", "ci95"), boot(e[k] - e["persistence"], rng))) for k in e if k != "persistence"}}
        print(c, len(e), {k: round(v, 3) for k, v in res[c]["mae"].items()}, flush=True)
    summ = {}
    for k in [k for k in res[COH[0]]["mae"] if k != "persistence"]:
        dif = {c: res[c]["diff_vs_persistence"][k]["mean"] for c in COH}
        summ[k] = {"mean6": float(np.mean(list(dif.values()))), "mean_first4": float(np.mean([dif[c] for c in COH[:4]])), "n_better_of_6": int(sum(v < 0 for v in dif.values())), "worst_cohort": float(max(dif.values())), "per_cohort": dif}
    ok = [k for k in ("VF1", "VF2") if summ[k]["mean6"] <= summ["vig"]["mean6"] - 0.010 and summ[k]["mean6"] < summ["GBR_img"]["mean6"] and summ[k]["mean_first4"] <= summ["vig"]["mean_first4"]
          and summ[k]["worst_cohort"] <= 0.03 and summ[k]["n_better_of_6"] >= 5]
    pick = min(ok, key=lambda k: summ[k]["mean6"]) if ok else None
    json.dump({"folds": res, "summary": summ, "passing_rule": ok, "selected": pick, "fold_parameters": allp, "seed": SEED, "prereg": "paper/trackB_imgfeat_prereg_2026-11-19.md"}, open(OUT / "imgfeat_result.json", "w"), indent=1)
    print(json.dumps({k: {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != "per_cohort"} for k, v in summ.items()}, indent=1))
    print("per cohort:", {k: {c: round(x, 3) for c, x in summ[k]["per_cohort"].items()} for k in summ}); print("passing rule:", ok, "selected:", pick)


if __name__ == "__main__":
    main()
