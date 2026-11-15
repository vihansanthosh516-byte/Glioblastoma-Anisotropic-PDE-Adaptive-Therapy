"""Track B: VIG hybrid on CFB-GBM + WAW. Protocol: paper/trackB_vig_cfb_waw_protocol_2026-11-15.md (locked first).

Modes follow the lock chain of the protocol:
  --freeze       step 2: fit VIG and the eight comparators on the 976 spent pairs; write frozen_params.json
  --build-pairs  step 5: build the sealed test_pairs.csv from tool-made CFB volumes and released WAW masks (prints no volume)
  --inventory    step 6: counts of eligible pairs and patients only
  --score        step 7: ONE scoring run; refuses to overwrite result.json
  --selftest     toy data only; writes nothing

VIG code (fit_vig, vig_curve) is imported unchanged from track_b_vig_lumiere.py.
"""
import argparse, gzip, hashlib, json, re, sys, warnings, zipfile
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parent.parent
AUDIT = ROOT / "output" / "track_b_code_audit_20261114"
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(AUDIT))
import track_b_mu_tcks as t, track_b_vig_lumiere as vig
from fix_and_loco4 import series, finish, fit_gomp_l1, gomp

OUT = ROOT / "output" / "track_b_vig_cfb_waw_20261115"
TRAIN = AUDIT / "pooled4_pairs.csv"
CFB_TABLE = ROOT / "research_notes" / "GBM patient analysis readiness 3" / "source_tables" / "CFB-GBM_mri_availability_v02_20260129.tsv"
CFB_VOLUMES = OUT / "cfb_tool_volumes.csv"            # written by the labelling step: id_patient, temporality, core_mm3
WAW = ROOT / "data" / "external" / "waw"
WAW_MASKS = WAW / "3_2_postop_masks.zip"
BASE = ["persistence", "two_point_exp", "pop_rate", "eb_random_slope", "gompertz", "gbr", "gompertz_l1", "gbr_l1"]
SEED, N_BOOT, LEVEL, MIN_PATIENTS = 20261115, 50000, 0.975, 100
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def feats(d):
    return pd.DataFrame({"y_in": d.y_in, "dt": d.dt_days, "r_hist": d.r_hist, "has_prior": d.r_hist.notna().astype(float)})


def fit_all(tr):
    r_pop = float(tr.r_obs.median()); t2, s2 = t.eb_params(tr, r_pop)
    gom, g1 = t.gompertz_fit(tr), fit_gomp_l1(tr)
    gb = HistGradientBoostingRegressor(random_state=0).fit(feats(tr), tr.y_out)
    gb1 = HistGradientBoostingRegressor(loss="absolute_error", random_state=0).fit(feats(tr), tr.y_out)
    th, lam = vig.fit_vig(tr)

    def predict(te):
        y, dt, lc = te.y_in.values, te.dt_days.values, te.lastch.values
        nh = y + vig.vig_curve(th, dt)
        return {"persistence": y, "two_point_exp": (te.y_in + te.r_hist.fillna(0.0) * te.dt_days).values, "pop_rate": y + r_pop * dt,
                "eb_random_slope": np.asarray(t.shrink_pred(te, r_pop, t2, s2)[0], float), "gompertz": np.asarray(gom(te), float),
                "gbr": gb.predict(feats(te)), "gompertz_l1": gomp(g1, y, dt), "gbr_l1": gb1.predict(feats(te)),
                "vig_nohist": nh, "vig": nh - (lam[0] * np.maximum(lc, 0.0) + lam[1] * np.minimum(lc, 0.0))}

    par = {"vig_rho_per_day": float(np.exp(th[0])), "vig_drift": float(th[1]), "vig_tau_days": float(np.exp(th[2])), "vig_lam_rise": float(lam[0]),
           "vig_lam_fall": float(lam[1]), "pop_rate_per_day": r_pop, "gompertz_l1": [float(x) for x in g1]}
    return predict, par


def probe():
    """fixed synthetic pairs; their predictions are the fingerprint of the frozen fit."""
    rows = series([4000., 9000., 2500., 2600.], [0., 30., 95., 250.], "a", "PROBE") + series([150., 20., 900.], [10., 70., 130.], "b", "PROBE") \
        + series([30000., 31000.], [0., 21.], "c", "PROBE") + series([600., 0., 50., 4000.], [0., 60., 150., 500.], "d", "PROBE")
    return finish(pd.DataFrame(rows))


def train():
    tr = pd.read_csv(TRAIN)
    return tr, fit_all(tr)


def freeze():
    f = OUT / "frozen_params.json"
    assert not f.exists(), "frozen_params.json exists; the freeze is done once"
    tr, (predict, par) = train(); pr = probe()
    import sklearn, scipy
    out = {"protocol": "paper/trackB_vig_cfb_waw_protocol_2026-11-15.md", "train_file": str(TRAIN.relative_to(ROOT)).replace("\\", "/"), "train_sha256": sha(TRAIN),
           "train_pairs": int(len(tr)), "train_patients": int(tr.patient_id.nunique()), "train_pairs_by_cohort": tr.cohort.value_counts().to_dict(),
           "parameters": par, "probe_pairs": int(len(pr)), "probe_predictions": {k: [float(x) for x in v] for k, v in predict(pr).items()},
           "versions": {"python": sys.version.split()[0], "numpy": np.__version__, "pandas": pd.__version__, "sklearn": sklearn.__version__, "scipy": scipy.__version__}}
    OUT.mkdir(parents=True, exist_ok=True); f.write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ("train_sha256", "train_pairs", "train_patients", "train_pairs_by_cohort", "parameters", "probe_pairs")}, indent=1))


def frozen_predict():
    fz = json.loads((OUT / "frozen_params.json").read_text())
    assert sha(TRAIN) == fz["train_sha256"], "training file changed after the freeze"
    _, (predict, _) = train(); got = predict(probe())
    for k, v in fz["probe_predictions"].items():
        assert np.allclose(got[k], v, atol=1e-6), f"refit of {k} does not reproduce the frozen fit"
    return predict


# ---- sealed test pairs (step 5) --------------------------------------------------------------
def cfb_pairs():
    m = pd.read_csv(CFB_TABLE, sep="\t"); v = pd.read_csv(CFB_VOLUMES)
    assert {"id_patient", "temporality", "core_mm3"} <= set(v.columns) and not v.duplicated(["id_patient", "temporality"]).any()
    day = {(r.id_patient, r.temporality): 7.0 * r._3 for r in m[["id_patient", "temporality", "time_diff_t0 (weeks)"]].itertuples()}
    vol = {(r.id_patient, r.temporality): float(r.core_mm3) for r in v.itertuples()}
    rows = []
    for pid in sorted({k[0] for k in vol}):
        run = []                                               # adjacent visits only: a missing visit breaks the chain
        for tp in ("t0", "t1", "t2", None):
            if tp is not None and (pid, tp) in vol and (pid, tp) in day:
                run.append(tp); continue
            if len(run) > 1:
                r = series([vol[(pid, k)] for k in run], [day[(pid, k)] for k in run], pid, "CFB")
                for a, x in zip(run[:-1], r): x["step"] = a + run[run.index(a) + 1]
                rows += r
            run = []
    d = pd.DataFrame(rows); d["grade"] = "4"
    return d


def waw_pairs():
    import nibabel as nib
    meta = pd.read_csv(WAW / "2_2_postop_scans_coreg_meta.csv", encoding="utf-8-sig")
    key = meta.file_name.str.extract(r"^(\d+)_(\d+)_"); days = {(a, b): float(d) for a, b, d in zip(key[0], key[1], meta.days_postop) if pd.notna(d)}
    known = set(zip(key[0], key[1]))                         # an exam with no day in the table cannot form a pair: it is dropped before pairing
    cl = pd.read_csv(WAW / "1_0_clinical_molecular_data.csv", sep=";", encoding="utf-8-sig"); grade = dict(zip(cl.id.astype(str), cl.who_grade.astype(str).str.strip()))
    vol = {}
    with zipfile.ZipFile(WAW_MASKS) as z:
        for n in sorted(x for x in z.namelist() if x.endswith(".nii.gz")):
            k = re.match(r"^(\d+)_(\d+)_", Path(n).name)
            assert k and k.groups() in known, f"mask file does not map to a post-operative exam in the table: {n}"
            assert k.groups() not in vol, f"two masks for one exam: {n}"
            if k.groups() not in days: continue
            img = nib.Nifti1Image.from_bytes(gzip.decompress(z.read(n))); a = np.asarray(img.dataobj)
            assert set(np.unique(a).astype(int)) <= {0, 1, 2, 3, 4}, f"unexpected label in {n}"
            vol[k.groups()] = float(((a == 1) | (a == 3)).sum()) * float(np.prod(img.header.get_zooms()[:3]))
    rows = []
    for pid in sorted({k[0] for k in vol}):
        ex = sorted((days[k], int(k[1]), vol[k]) for k in vol if k[0] == pid)
        r = series([e[2] for e in ex], [e[0] for e in ex], pid, "WAW")
        for x in r: x["step"], x["grade"] = "post", grade.get(pid, "")
        rows += r
    return pd.DataFrame(rows)


def build_pairs():
    f = OUT / "test_pairs.csv"
    assert (OUT / "frozen_params.json").exists(), "freeze first (step 2)"
    for g in ("gate_g1.json", "gate_g2.json"):
        assert json.loads((OUT / g).read_text()).get("pass") is True, f"{g} did not pass; no test pairs are built"
    assert not f.exists(), "test_pairs.csv exists"
    d = pd.concat([cfb_pairs(), waw_pairs()], ignore_index=True); meta = d[["step", "grade"]]
    d = finish(d.drop(columns=["step", "grade"]).assign(step=meta.step, grade=meta.grade))
    d.to_csv(f, index=False); print("sealed test_pairs.csv written; sha256", sha(f))


def inventory():
    d = pd.read_csv(OUT / "test_pairs.csv")
    out = {"test_pairs_sha256": sha(OUT / "test_pairs.csv"), "pairs": int(len(d)), "patients": int(d.patient_id.nunique()),
           "by_cohort": {c: {"pairs": int(len(g)), "patients": int(g.patient_id.nunique()), "pairs_with_earlier_scan": int((g.n_prev > 0).sum())} for c, g in d.groupby("cohort")},
           "by_step": d.step.value_counts().to_dict(), "interval_days_quantiles_10_50_90": d.dt_days.quantile([.1, .5, .9]).round(1).tolist(),
           "minimum_patients": MIN_PATIENTS, "enough_patients": bool(d.patient_id.nunique() >= MIN_PATIENTS)}
    (OUT / "inventory.json").write_text(json.dumps(out, indent=1)); print(json.dumps(out, indent=1))


# ---- scoring (step 7) ------------------------------------------------------------------------
def strat_boot(diff, strata, n=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed); tot = np.zeros(n)
    for s in sorted(set(strata)):
        x = diff[strata == s]; tot += x[rng.integers(0, len(x), (n, len(x)))].sum(1)
    return tot / len(diff)


def compare(pt, a, b):
    """a minus b on the patient table pt (index = patient, column cohort = stratum)."""
    diff = (pt[a] - pt[b]).values; bt = strat_boot(diff, pt.cohort.values)
    lo, hi = (1 - LEVEL) / 2 * 100, (1 + LEVEL) / 2 * 100
    flips = np.random.default_rng(SEED + 1).choice([-1.0, 1.0], (20000, len(diff)))
    return {"mean": float(diff.mean()), "ci_primary": [float(np.percentile(bt, lo)), float(np.percentile(bt, hi))], "ci_level": LEVEL,
            "ci95": [float(np.percentile(bt, 2.5)), float(np.percentile(bt, 97.5))], "p_signflip": float((np.abs((flips * diff).mean(1)) >= abs(diff.mean())).mean()),
            "n_patients": int(len(diff))}


def patient_table(d, pred):
    e = pd.DataFrame({k: np.abs(v - d.y_out.values) for k, v in pred.items()})
    e["pid"], e["cohort"] = d.patient_id.values, d.cohort.values
    return e.groupby("pid").agg({**{k: "mean" for k in pred}, "cohort": "first"})


def evaluate(d, pred):
    pt = patient_table(d, pred)
    res = {"patients": int(len(pt)), "pairs": int(len(d)), "patient_mae": {k: float(pt[k].mean()) for k in pred},
           "vig_minus": {b: compare(pt, "vig", b) for b in BASE}}
    res["passed"] = [b for b in BASE if res["vig_minus"][b]["ci_primary"][1] < 0]
    res["primary_success"] = len(res["passed"]) == len(BASE)
    sub = {"cohort_" + c: d.cohort == c for c in sorted(d.cohort.unique())}
    sub.update({"step_" + s: d.step == s for s in sorted(d.step.astype(str).unique())})
    sub.update({"waw_grade4": (d.cohort == "WAW") & (d.grade.astype(str) == "4"), "interval_lt60": d.dt_days < 60,
                "interval_60_120": (d.dt_days >= 60) & (d.dt_days <= 120), "interval_gt120": d.dt_days > 120})
    res["secondary_no_claim"] = {"vig_nohist_minus_persistence": compare(pt, "vig_nohist", "persistence"), "vig_minus_vig_nohist": compare(pt, "vig", "vig_nohist")}
    for k, m in sub.items():
        if m.sum() < 5: continue
        p2 = patient_table(d[m.values], {a: v[m.values] for a, v in pred.items()})
        res["secondary_no_claim"][k] = {"patients": int(len(p2)), "patient_mae": {a: round(float(p2[a].mean()), 4) for a in pred}, "vig_minus_persistence": compare(p2, "vig", "persistence")}
    return res


def score():
    f = OUT / "result.json"
    assert not f.exists(), "result.json exists; the test is run once"
    inv = json.loads((OUT / "inventory.json").read_text())
    assert inv["test_pairs_sha256"] == sha(OUT / "test_pairs.csv") and inv["enough_patients"], "inventory missing, stale, or below the minimum patient count"
    predict = frozen_predict(); d = pd.read_csv(OUT / "test_pairs.csv"); pred = predict(d)
    res = evaluate(d, pred); res.update({"protocol": "paper/trackB_vig_cfb_waw_protocol_2026-11-15.md", "seed": SEED, "n_boot": N_BOOT, "test_pairs_sha256": inv["test_pairs_sha256"]})
    d.assign(**{"pred_" + k: v for k, v in pred.items()}).to_csv(OUT / "predictions.csv", index=False)
    f.write_text(json.dumps(res, indent=1)); print(json.dumps({k: res[k] for k in ("patients", "pairs", "patient_mae", "passed", "primary_success")}, indent=1))


def selftest():
    rng = np.random.default_rng(0); rows = []
    for c, n in (("CFB", 60), ("WAW", 20)):
        for i in range(n):
            k = int(rng.integers(2, 5)); days = np.cumsum(rng.integers(20, 200, k)).astype(float); y = 7 + np.cumsum(rng.normal(0, 1, k))
            r = series(list(np.round(np.exp(y))), list(days), i, c)
            for x in r: x["step"], x["grade"] = "toy", "4"
            rows += r
    d = pd.DataFrame(rows); meta = d[["step", "grade"]]; d = finish(d.drop(columns=["step", "grade"]).assign(step=meta.step, grade=meta.grade))
    _, (predict, par) = train(); res = evaluate(d, predict(d))
    x = rng.normal(-0.1, 1.0, 300); s = np.array(["a"] * 200 + ["b"] * 100); bt = strat_boot(x, s, n=20000); se = x.std(ddof=1) / np.sqrt(len(x))
    print("toy pairs", len(d), "patients", res["patients"], "| toy patient MAE", {k: round(v, 3) for k, v in res["patient_mae"].items()})
    print("toy passed:", res["passed"], "| primary_success", res["primary_success"], "(random-walk toy data: no method should pass all eight)")
    print("stratified bootstrap 95% CI", [round(float(np.percentile(bt, q)), 3) for q in (2.5, 97.5)], "| normal theory", [round(x.mean() - 1.96 * se, 3), round(x.mean() + 1.96 * se, 3)])
    print("secondary keys:", sorted(res["secondary_no_claim"]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); g = ap.add_mutually_exclusive_group(required=True)
    for k in ("freeze", "build-pairs", "inventory", "score", "selftest"): g.add_argument("--" + k, action="store_true")
    a = ap.parse_args()
    {"freeze": freeze, "build_pairs": build_pairs, "inventory": inventory, "score": score, "selftest": selftest}[next(k for k, v in vars(a).items() if v)]()
