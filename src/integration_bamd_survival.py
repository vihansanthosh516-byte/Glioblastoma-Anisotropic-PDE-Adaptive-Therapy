"""BAMD step 5: UPENN survival analysis (primary H1, secondary H2/H3, all benchmarks).

Protocol: paper/integration_protocol_bamd_2026-11-12.md + note 1. Run once, after Gate 1
and after frozen_map.json is committed. Refuses to run if Gate 1 passed no state, unless
--residual-only (then only H2 and benchmarks are reported, and H1 is "not tested").

  python -I src/integration_bamd_survival.py [--boot 1000] [--deepsurv]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
from integration_bamd_features import feature_names  # noqa: E402

OUT = PROJECT_ROOT / "output" / "integration_bamd"
CLIN = PROJECT_ROOT / "data" / "external" / "upenn_meta" / "clin.csv"
SEED = 20261115
VOLS = ["logvol_ET", "logvol_NCR", "logvol_ED"]


# ---------- survival utilities ----------
def cox_fit(X, t, e, ridge=1e-4, iters=50):
    """Breslow Cox by Newton-Raphson with a small ridge. Returns coefficients."""
    n, p = X.shape
    order = np.argsort(-t, kind="stable")
    X, t, e = X[order], t[order], e[order]
    beta = np.zeros(p)
    for _ in range(iters):
        eta = np.clip(X @ beta, -30, 30)
        w = np.exp(eta)
        cw = np.cumsum(w)
        cwx = np.cumsum(w[:, None] * X, axis=0)
        cwxx = np.cumsum(w[:, None, None] * X[:, :, None] * X[:, None, :], axis=0)
        # tie handling (Breslow): use the risk set at the last index of each tied time block
        last = np.r_[np.nonzero(np.diff(t))[0], n - 1]
        idx = np.searchsorted(last, np.arange(n))
        risk_end = last[idx]
        S0, S1, S2 = cw[risk_end], cwx[risk_end], cwxx[risk_end]
        ev = e.astype(bool)
        mean = S1[ev] / S0[ev, None]
        grad = (X[ev] - mean).sum(0) - ridge * beta
        H = (S2[ev] / S0[ev, None, None] - mean[:, :, None] * mean[:, None, :]).sum(0) + ridge * np.eye(p)
        step = np.linalg.solve(H, grad)
        beta = beta + step
        if np.abs(step).max() < 1e-8:
            break
    return beta


def lasso_cox(X, t, e, lam, iters=500, lr=None):
    """Proximal gradient Cox with L1 penalty (Breslow, no tie correction beyond ordering)."""
    n, p = X.shape
    order = np.argsort(-t, kind="stable")
    X, t, e = X[order], t[order], e[order].astype(bool)
    beta = np.zeros(p)
    lr = lr or 1.0 / max(1.0, e.sum())
    for _ in range(iters):
        w = np.exp(np.clip(X @ beta, -30, 30))
        cw = np.cumsum(w)
        cwx = np.cumsum(w[:, None] * X, axis=0)
        grad = -(X[e] - cwx[e] / cw[e, None]).sum(0)
        z = beta - lr * grad
        beta = np.sign(z) * np.maximum(np.abs(z) - lr * lam, 0.0)
    return beta


def cindex(risk, t, e):
    """Harrell's C: higher risk should mean shorter time."""
    t, e, r = np.asarray(t), np.asarray(e, bool), np.asarray(risk)
    comparable = (t[:, None] < t[None, :]) & e[:, None]
    conc = (r[:, None] > r[None, :]) & comparable
    ties = (r[:, None] == r[None, :]) & comparable
    den = comparable.sum()
    return float((conc.sum() + 0.5 * ties.sum()) / den) if den else np.nan


def standardize(tr, te):
    mu, sd = tr.mean(0), tr.std(0)
    sd = np.where(sd > 0, sd, 1.0)
    return (tr - mu) / sd, (te - mu) / sd


# ---------- data ----------
def load_upenn():
    feats = pd.read_csv(OUT / "features.csv")
    up = feats[feats["collection"] == "UPENN-GBM"].copy()
    clin = pd.read_csv(CLIN).rename(columns={"ID": "patient_id"})
    df = up.merge(clin, on="patient_id", how="inner", validate="one_to_one")
    df = df[df["IDH1"] != "Mutated"].copy()
    dead = df["Survival_Status"] == "Deceased"
    df["event"] = dead.astype(int)
    t_dead = pd.to_numeric(df["Survival_from_surgery_days_UPDATED"], errors="coerce")
    t_cens = pd.to_numeric(df["Survival_Censor"], errors="coerce")
    df["time"] = np.where(dead, t_dead, t_cens)
    df = df[np.isfinite(df["time"]) & (df["time"] > 0)].copy()
    df["age"] = pd.to_numeric(df["Age_at_scan_years"], errors="coerce")
    df["male"] = (df["Gender"] == "M").astype(float)
    df["gtr_yes"] = (df["GTR_over90percent"] == "Y").astype(float)
    df["gtr_unknown"] = (~df["GTR_over90percent"].isin(["Y", "N"])).astype(float)
    df = df[np.isfinite(df["age"])].reset_index(drop=True)
    return df


def bamd_parts(df, frozen):
    names = frozen["feature_names"]
    X = df[names].to_numpy(float)
    med = np.array(frozen["impute_median"])
    X = np.where(np.isnan(X), med, X)
    Xs = (X - np.array(frozen["mu"])) / np.array(frozen["sd"])
    P = Xs @ np.array(frozen["W"])
    R = Xs - P @ np.array(frozen["B"])
    pcs = (R - np.array(frozen["R_mean"])) @ np.array(frozen["pcs"]).T
    return P, pcs


# ---------- models ----------
CLINICAL = ["age", "male", "gtr_yes", "gtr_unknown"]


def design(df, frozen, which):
    base = df[CLINICAL].to_numpy(float)
    vols = df[VOLS].to_numpy(float)
    P, pcs = bamd_parts(df, frozen) if frozen else (np.zeros((len(df), 0)), np.zeros((len(df), 0)))
    allmri = df[feature_names()].to_numpy(float)
    allmri = np.where(np.isnan(allmri), np.nanmedian(allmri, axis=0), allmri)
    return {
        "clinical": base,
        "B0": np.column_stack([base, vols]),
        "H1_B0+states": np.column_stack([base, vols, P]),
        "H2_B0+residual": np.column_stack([base, vols, pcs]),
        "BAMD_B0+states+residual": np.column_stack([base, vols, P, pcs]),
        "all_mri": np.column_stack([base, allmri]),
    }[which]


def fit_predict(model, Xtr, ttr, etr, Xte, rng):
    Xtr, Xte = standardize(Xtr, Xte)
    if model in ("lasso_all", "ridge_all"):
        lams = [0.01, 0.03, 0.1, 0.3, 1.0] if model == "lasso_all" else [0.1, 1, 10, 100, 1000]
        best, best_c = lams[0], -1
        idx = rng.permutation(len(Xtr))
        for lam in lams:
            cs = []
            for k in range(5):
                te = idx[k::5]
                tr = np.setdiff1d(idx, te)
                b = lasso_cox(Xtr[tr], ttr[tr], etr[tr], lam * len(tr)) if model == "lasso_all" \
                    else cox_fit(Xtr[tr], ttr[tr], etr[tr], ridge=lam)
                cs.append(cindex(Xtr[te] @ b, ttr[te], etr[te]))
            if np.nanmean(cs) > best_c:
                best, best_c = lam, np.nanmean(cs)
        b = lasso_cox(Xtr, ttr, etr, best * len(Xtr)) if model == "lasso_all" else cox_fit(Xtr, ttr, etr, ridge=best)
        return Xte @ b
    if model == "deepsurv_all":
        return deepsurv(Xtr, ttr, etr, Xte, rng)
    return Xte @ cox_fit(Xtr, ttr, etr)


def deepsurv(Xtr, ttr, etr, Xte, rng):
    import torch
    torch.manual_seed(int(rng.integers(0, 2**31 - 1)))
    order = np.argsort(-ttr, kind="stable")
    X = torch.tensor(Xtr[order], dtype=torch.float32)
    e = torch.tensor(etr[order], dtype=torch.float32)
    net = torch.nn.Sequential(torch.nn.Linear(X.shape[1], 32), torch.nn.ReLU(), torch.nn.Dropout(0.2),
                              torch.nn.Linear(32, 1))
    opt = torch.optim.Adam(net.parameters(), lr=1e-3, weight_decay=1e-3)
    for _ in range(200):
        opt.zero_grad()
        r = net(X).squeeze(1)
        loss = -((r - torch.logcumsumexp(r, 0)) * e).sum() / e.sum().clamp(min=1)
        loss.backward()
        opt.step()
    net.eval()
    with torch.no_grad():
        return net(torch.tensor(Xte, dtype=torch.float32)).squeeze(1).numpy()


MODELS = {  # name -> (design key, fitter)
    "clinical_cox": ("clinical", "cox"),
    "B0_clinical+volumes": ("B0", "cox"),
    "lasso_cox_all_mri": ("all_mri", "lasso_all"),
    "ridge_cox_all_mri": ("all_mri", "ridge_all"),
    "deepsurv_all_mri": ("all_mri", "deepsurv_all"),
    "H1_B0+states": ("H1_B0+states", "cox"),
    "H2_B0+residual": ("H2_B0+residual", "cox"),
    "BAMD_B0+states+residual": ("BAMD_B0+states+residual", "cox"),
}


def cv_scores(df, frozen, models, rng, repeats, k=5, rows=None):
    rows = np.arange(len(df)) if rows is None else rows
    t, e = df["time"].to_numpy(float)[rows], df["event"].to_numpy(int)[rows]
    Xs = {m: design(df, frozen, MODELS[m][0])[rows] for m in models}
    out = {m: [] for m in models}
    for _ in range(repeats):
        perm = rng.permutation(len(rows))
        for f in range(k):
            te = perm[f::k]
            tr = np.setdiff1d(perm, te)
            for m in models:
                risk = fit_predict(MODELS[m][1], Xs[m][tr], t[tr], e[tr], Xs[m][te], rng)
                out[m].append(cindex(risk, t[te], e[te]))
    return {m: np.array(v) for m, v in out.items()}


def lr_test(df, frozen, small, big):
    from scipy.stats import chi2
    t, e = df["time"].to_numpy(float), df["event"].to_numpy(int)

    def loglik(X):
        X = (X - X.mean(0)) / np.where(X.std(0) > 0, X.std(0), 1)
        b = cox_fit(X, t, e)
        order = np.argsort(-t, kind="stable")
        eta = (X @ b)[order]
        ee = e[order].astype(bool)
        return float((eta - np.log(np.cumsum(np.exp(eta))))[ee].sum())
    Xs, Xb = design(df, frozen, small), design(df, frozen, big)
    stat = 2 * (loglik(Xb) - loglik(Xs))
    dfree = Xb.shape[1] - Xs.shape[1]
    return {"chi2": stat, "df": dfree, "p": float(chi2.sf(max(stat, 0), dfree))}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--repeats", type=int, default=10)
    ap.add_argument("--deepsurv", action="store_true")
    ap.add_argument("--residual-only", action="store_true")
    args = ap.parse_args(argv)
    gate = json.loads((OUT / "gate1.json").read_text())
    frozen = None
    if gate["passed"]:
        body = (OUT / "frozen_map.json").read_bytes()
        if hashlib.sha256(body).hexdigest() != gate["frozen_map_sha256"]:
            raise SystemExit("frozen map hash mismatch")
        frozen = json.loads(body)
    elif not args.residual_only:
        raise SystemExit("Gate 1 passed no state; H1 is not tested (protocol). Use --residual-only only if declared.")
    models = [m for m in MODELS if (args.deepsurv or "deepsurv" not in m)]
    if frozen is None:
        models = [m for m in models if m.startswith(("clinical", "B0", "lasso", "ridge", "deepsurv"))]
    rng = np.random.Generator(np.random.PCG64DXSM(SEED))
    df = load_upenn()
    t0 = time.time()
    point = cv_scores(df, frozen, models, rng, args.repeats)
    point_mean = {m: float(np.nanmean(v)) for m, v in point.items()}
    boot = {m: [] for m in models}
    n = len(df)
    for _ in range(args.boot):
        rows = rng.integers(0, n, n)
        s = cv_scores(df, frozen, models, rng, 1, rows=rows)
        for m in models:
            boot[m].append(float(np.nanmean(s[m])))
    boot = {m: np.array(v) for m, v in boot.items()}

    def delta(a, b):
        d = boot[a] - boot[b]
        return {"point": point_mean[a] - point_mean[b], "ci95": [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]}

    res = {"n": n, "events": int(df["event"].sum()), "models": models, "c_index_cv": point_mean,
           "c_index_boot_ci95": {m: [float(np.percentile(boot[m], 2.5)), float(np.percentile(boot[m], 97.5))] for m in models},
           "repeats": args.repeats, "boot": args.boot, "runtime_s": round(time.time() - t0, 1),
           "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           "gate1_passed": gate["passed"]}
    if frozen:
        res["H1"] = {"delta_vs_B0": delta("H1_B0+states", "B0_clinical+volumes"),
                     "lr_test": lr_test(df, frozen, "B0", "H1_B0+states")}
        res["H1"]["passes"] = res["H1"]["delta_vs_B0"]["ci95"][0] > 0 and res["H1"]["lr_test"]["p"] < 0.05
        res["H2"] = {"delta_vs_B0": delta("H2_B0+residual", "B0_clinical+volumes"),
                     "lr_test": lr_test(df, frozen, "B0", "H2_B0+residual")}
        res["H3"] = {"delta_vs_B0": delta("BAMD_B0+states+residual", "B0_clinical+volumes"),
                     "lr_test": lr_test(df, frozen, "B0", "BAMD_B0+states+residual")}
        others = [m for m in models if m not in ("BAMD_B0+states+residual", "H1_B0+states", "H2_B0+residual")]
        res["BAMD_vs_baselines"] = {m: delta("BAMD_B0+states+residual", m) for m in others}
        res["BAMD_better_than_all"] = all(v["ci95"][0] > 0 for v in res["BAMD_vs_baselines"].values())
    res["not_run"] = ["random_survival_forest", "gradient_boosted_survival"] + ([] if args.deepsurv else ["deepsurv_all_mri"])
    (OUT / "survival_result.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(json.dumps({k: res.get(k) for k in ("n", "events", "c_index_cv", "H1", "BAMD_better_than_all")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
