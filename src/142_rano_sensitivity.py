#!/usr/bin/env python3
"""Script 142: Amendment 18 sensitivity analyses. EXPLORATORY / SIMULATION. Declared in analysis_plan_v1.md
(commit 17c95dc) before this file was written. Replaces nothing: A15 (script 140) and A16 (script 141) stand as run.

Part A: finer constant-noise grid, s in {0.05, 0.10, 0.15, 0.20} (+ 0.25 as a reproduction check of A16), N=100, D=8.
Part B: AR(1) correlated scan noise, sigma in {0.10, 0.25, 0.50}, rho in {0, 0.5, 0.8}.
Part C: independent exponential death (median 78 weeks, assumed) as a PFS event, s in {0, 0.25, 0.5}.
Part D: real-data sensitivity of A15 (exclusion flow, floor 0/250/500/1000, nadir reference).
Operational choices fixed here before the first run (they only make A18 executable):
  * Part A/B/C reuse the A16 model, seed [20261018, N, D] and draw order, so s=0.25 reproduces A16 exactly.
  * AR(1): e_0 = u_0; e_j = rho e_{j-1} + sqrt(1 - rho^2) u_j (marginal SD sigma), per patient.
  * Death: T_d ~ Exp(median 78 weeks), separate stream seed [20261018, N, D, 7]; death at T_d <= last scan time is an event
    placed at scan index ceil(T_d / D) (first scan at or after death); PFS event = earlier of progression call and death.
  * Part D nadir: per tool, progression ratio = V_i / min(earlier volumes of that tool >= 500 mm3 in the post-operative
    sequence); response ratio = V_i / V_first (first post-operative scan, must be >= 500). Class PD if nadir ratio >= 1.40,
    else PR if baseline ratio <= 0.35, else SD. Pair eligible if both tools have V_i >= 500, a nadir and a baseline.
    Floor 0 means "all four volumes > 0".
  * Materiality: E1 |mean| >= 4.345 weeks with CI excluding 0; E2 <= -0.05 with CI excluding 0.
Run once. No re-run with other settings.
Output: output/rano_sensitivity/ (cells csv, summary json).
"""
from __future__ import annotations

import importlib.util
import json
import time
from pathlib import Path as _Path

import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
SRC = PROJECT_ROOT / "src"
OUT = PROJECT_ROOT / "output" / "rano_sensitivity"


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, SRC / fname)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


T = _load("trial141", "141_trial_sim.py")
P = _load("pairs140", "140_rano_tool_dependence.py")

SEED = 20261018
R = 2000
N, D = 100, 8
HORIZON = 104
MEDIAN_DEATH_WK = 78.0
E1_MATERIAL, E2_MATERIAL = 4.345, -0.05


# ------------------------------------------------------------------ simulation core
def arm_calls(lnv, eps, rule, death):
    """Return idx, ev for one arm. eps = additive ln noise array (R,N,J+1) ; death = (idx_d, ev_d) or None."""
    idx, ev = T.event_index(lnv + eps, rule)
    if death is not None:
        idx_d, ev_d = death
        J = lnv.shape[2] - 1
        big = J + 10
        ip = np.where(ev, idx, big)
        idd = np.where(ev_d, idx_d, big)
        mn = np.minimum(ip, idd)
        ev = mn < big
        idx = np.where(ev, mn, J)
    return idx, ev


def make_world(times, J):
    rng = np.random.default_rng([SEED, N, D])
    Zc = rng.standard_normal((R, N)); Ze = rng.standard_normal((R, N))
    L0c = T.LNV0_MU + T.LNV0_SD * rng.standard_normal((R, N))
    L0e = T.LNV0_MU + T.LNV0_SD * rng.standard_normal((R, N))
    Uc = rng.standard_normal((R, N, J + 1)); Ue = rng.standard_normal((R, N, J + 1))
    gc = T.G_MED * np.exp(T.G_SD * Zc)
    ge = T.G_MED * np.exp(T.G_SD * Ze)
    rd = np.random.default_rng([SEED, N, D, 7])
    lam = np.log(2.0) / MEDIAN_DEATH_WK
    Tdc = rd.exponential(1.0 / lam, (R, N)); Tde = rd.exponential(1.0 / lam, (R, N))
    tj = float(times[-1])

    def dpack(Td):
        return (np.minimum(np.ceil(Td / D).astype(int), J), Td <= tj)

    return dict(
        lnv_c=L0c[:, :, None] + gc[:, :, None] * times[None, None, :],
        lnv_e={1.0: L0e[:, :, None] + ge[:, :, None] * times[None, None, :],
               0.7: L0e[:, :, None] + (0.7 * ge)[:, :, None] * times[None, None, :]},
        Uc=Uc, Ue=Ue, death_c=dpack(Tdc), death_e=dpack(Tde))


def ar1(U, sigma, rho):
    e = np.empty_like(U)
    e[:, :, 0] = U[:, :, 0]
    k = np.sqrt(1.0 - rho ** 2)
    for j in range(1, U.shape[2]):
        e[:, :, j] = rho * e[:, :, j - 1] + k * U[:, :, j]
    return sigma * e


def evaluate(W, J, times, sigma_scale, rho, rule, use_death):
    """Per-trial outputs for one noise scenario: control median, chi2/lhr at theta 0.7 and 1.0."""
    if sigma_scale == 0:
        ec = ee = 0.0
    elif rho is None:  # A16 convention: constant change-SD s, per-scan SD s / sqrt(2), independent
        ec = sigma_scale / np.sqrt(2.0) * W["Uc"]
        ee = sigma_scale / np.sqrt(2.0) * W["Ue"]
    else:
        ec = ar1(W["Uc"], sigma_scale, rho)
        ee = ar1(W["Ue"], sigma_scale, rho)
    dc = W["death_c"] if use_death else None
    de = W["death_e"] if use_death else None
    ic, vc = arm_calls(W["lnv_c"], ec, rule, dc)
    med = T.km_median_grid(ic, vc, times)
    out = dict(med=med)
    for th in (0.7, 1.0):
        ie, ve = arm_calls(W["lnv_e"][th], ee, rule, de)
        out[th] = T.logrank(ic, vc, ie, ve, J)
    return out


def summarise(cur, base):
    ok = ~np.isnan(cur["med"]) & ~np.isnan(base["med"])
    e1 = T.mc_ci((cur["med"] - base["med"])[ok]) if ok.sum() > 2 else (np.nan,) * 3
    chi2, lhr, valid = cur[0.7]
    chi2b, lhrb, validb = base[0.7]
    d = (chi2 > T.CHI2_CRIT).astype(float) - (chi2b > T.CHI2_CRIT).astype(float)
    e2 = T.mc_ci(d)
    okh = valid & validb
    e4 = T.mc_ci(np.where(okh, lhr - lhrb, np.nan))
    typ1 = float((cur[1.0][0] > T.CHI2_CRIT).mean())
    return dict(E1=e1[0], E1_lo=e1[1], E1_hi=e1[2], E1_n=int(ok.sum()),
                E2=e2[0], E2_lo=e2[1], E2_hi=e2[2], E4=e4[0], E4_lo=e4[1], E4_hi=e4[2],
                type1=typ1, ctrl_median=float(np.nanmean(cur["med"])))


def mat_e1(r):
    return bool(abs(r["E1"]) >= E1_MATERIAL and (r["E1_lo"] > 0 or r["E1_hi"] < 0))


def mat_e2(r):
    return bool(r["E2"] <= E2_MATERIAL and r["E2_hi"] < 0)


def simulate_parts():
    times = np.arange(0, HORIZON + 1, D, dtype=float)
    J = len(times) - 1
    W = make_world(times, J)
    rows = []
    for rule in T.CONFIRMS:
        base = evaluate(W, J, times, 0, None, rule, False)
        base_d = evaluate(W, J, times, 0, None, rule, True)
        # Part A
        for s in (0.05, 0.10, 0.15, 0.20, 0.25):
            r = summarise(evaluate(W, J, times, s, None, rule, False), base)
            rows.append(dict(part="A", rule=rule, s=s, sigma=np.nan, rho=np.nan, death=False, **r))
        # Part B
        for sigma in (0.10, 0.25, 0.50):
            for rho in (0.0, 0.5, 0.8):
                r = summarise(evaluate(W, J, times, sigma, rho, rule, False), base)
                rows.append(dict(part="B", rule=rule, s=np.nan, sigma=sigma, rho=rho, death=False, **r))
        # Part C
        for s in (0.0, 0.25, 0.5):
            cur = evaluate(W, J, times, s, None, rule, True)
            r = summarise(cur, base_d)
            rows.append(dict(part="C", rule=rule, s=s, sigma=np.nan, rho=np.nan, death=True, **r))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ Part D
def build_sequences():
    hd, dbt = P.load("hdglioauto"), P.load("deepbratumia")
    m = hd.merge(dbt, on=["Patient", "tp"], how="outer", suffixes=("_hd", "_dbt"))
    rat = pd.read_csv(P.LUM / "LUMIERE-ExpertRating-v202211.csv")
    rat.columns = ["Patient", "tp", "lt3", "nonmeas", "rating", "rationale"]
    rat["rating"] = rat["rating"].astype(str).str.strip()
    pre = rat[rat["rating"] == "Pre-Op"].set_index("Patient")["tp"].map(P.tp_key).to_dict()
    m["k"] = m["tp"].map(P.tp_key)
    m = m[m["Patient"].isin(pre)]
    m = m[[k > pre[p] for p, k in zip(m["Patient"], m["k"])]]
    rating_of = {(p, t): r for p, t, r in zip(rat["Patient"], rat["tp"], rat["rating"])}
    seqs = {}
    for pid, g in m.sort_values(["Patient", "k"]).groupby("Patient"):
        g = g.reset_index(drop=True)
        seqs[pid] = dict(tp=list(g["tp"]), wk=[k[0] for k in g["k"]], hd=g["vol_hd"].to_numpy(float),
                         dbt=g["vol_dbt"].to_numpy(float), rating=[rating_of.get((pid, t), "") for t in g["tp"]])
    return seqs


def pair_table(seqs):
    rows = []
    for pid, s in seqs.items():
        n = len(s["tp"])
        for i in range(1, n):
            v = [s["hd"][i - 1], s["dbt"][i - 1], s["hd"][i], s["dbt"][i]]
            rec = dict(Patient=pid, i=i, hd0=v[0], dbt0=v[1], hd1=v[2], dbt1=v[3],
                       weeks_since_first=s["wk"][i] - s["wk"][0], expert_later=s["rating"][i])
            # nadir / baseline classes
            cls = {}
            for tool in ("hd", "dbt"):
                vols = s[tool]
                cur = vols[i]
                earlier = vols[:i]
                earlier = earlier[~np.isnan(earlier) & (earlier >= P.FLOOR)]
                base = vols[0]
                if (not np.isnan(cur)) and cur >= P.FLOOR and len(earlier) > 0 and (not np.isnan(base)) and base >= P.FLOOR:
                    rn, rb = cur / earlier.min(), cur / base
                    cls[tool] = "PD" if rn >= P.PD_LINE else ("PR" if rb <= P.PR_LINE else "SD")
                else:
                    cls[tool] = None
            rec["nadir_cls_hd"], rec["nadir_cls_dbt"] = cls["hd"], cls["dbt"]
            rows.append(rec)
    return pd.DataFrame(rows)


def flip_stats(df, rng):
    df = df.copy()
    df["r_hd"], df["r_dbt"] = df["hd1"] / df["hd0"], df["dbt1"] / df["dbt0"]
    df["flip"] = (P.cls(df["r_hd"]) != P.cls(df["r_dbt"])).astype(int)
    ci = P.boot_ci(df, lambda d: d["flip"].mean(), rng)
    return dict(n_pairs=int(len(df)), n_patients=int(df["Patient"].nunique()), flip_rate=float(df["flip"].mean()), ci95=ci)


def part_d():
    seqs = build_sequences()
    pt = pair_table(seqs)
    rng = np.random.default_rng(SEED)
    vols = pt[["hd0", "dbt0", "hd1", "dbt1"]]
    res = {"n_consecutive_pairs": int(len(pt))}
    # (i) exclusion flow
    miss = vols.isna(); zero = (vols == 0); low = (vols > 0) & (vols < P.FLOOR)
    any_miss, any_zero, any_low = miss.any(axis=1), zero.any(axis=1), low.any(axis=1)
    primary = np.where(any_miss, "missing", np.where(any_zero, "zero", np.where(any_low, "below_floor", "eligible")))
    res["exclusion_flow"] = dict(
        primary_reason_mutually_exclusive=pd.Series(primary).value_counts().to_dict(),
        pairs_with_any_missing=int(any_miss.sum()), pairs_with_any_zero=int(any_zero.sum()), pairs_with_any_below_floor=int(any_low.sum()),
        by_volume_column={c: dict(missing=int(miss[c].sum()), zero=int(zero[c].sum()), below_floor=int(low[c].sum())) for c in vols.columns})
    elig = primary == "eligible"
    rated = pt["expert_later"].isin(["PD", "SD", "PR", "CR"])

    def grp(mask):
        sub = pt[mask]
        rr = sub[sub["expert_later"].isin(["PD", "SD", "PR", "CR"])]
        return dict(n_pairs=int(len(sub)), median_weeks_since_first_postop=float(sub["weeks_since_first"].median()),
                    n_rated=int(len(rr)), share_expert_PD_among_rated=float((rr["expert_later"] == "PD").mean()) if len(rr) else None)

    res["eligible_vs_excluded"] = dict(eligible=grp(elig), excluded=grp(~elig))
    # (ii) floors
    fl = {}
    for f in (0, 250, 500, 1000):
        mask = (vols > 0).all(axis=1) if f == 0 else (vols >= f).all(axis=1)
        fl[str(f)] = flip_stats(pt[mask], rng)
    res["floor_sensitivity"] = fl
    # (iii) nadir reference
    ok = pt["nadir_cls_hd"].notna() & pt["nadir_cls_dbt"].notna()
    nd = pt[ok].copy()
    nd["flip"] = (nd["nadir_cls_hd"] != nd["nadir_cls_dbt"]).astype(int)
    res["nadir_reference"] = dict(n_pairs=int(len(nd)), n_patients=int(nd["Patient"].nunique()), flip_rate=float(nd["flip"].mean()),
                                  ci95=P.boot_ci(nd, lambda d: d["flip"].mean(), rng),
                                  class_counts_hd=nd["nadir_cls_hd"].value_counts().to_dict(),
                                  class_counts_dbt=nd["nadir_cls_dbt"].value_counts().to_dict())
    return res, pt


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    cells = simulate_parts()
    cells.to_csv(OUT / "rano_sensitivity_cells.csv", index=False)
    dres, pt = part_d()
    pt.to_csv(OUT / "rano_sensitivity_pairs.csv", index=False)

    # reproduction check against A16 at s = 0.25, single-scan, N=100, D=8, theta 0.7
    a16 = pd.read_csv(PROJECT_ROOT / "output" / "trial_sim" / "trial_sim_cells.csv")
    a16 = a16[(a16.N == N) & (a16.D == D) & (a16.theta == 0.7)]
    chk = {}
    for rule in T.CONFIRMS:
        mine = cells[(cells.part == "A") & (cells.rule == rule) & (np.isclose(cells.s, 0.25))].iloc[0]
        ref = a16[(a16.s == "0.25") & (a16.confirm == rule)].iloc[0]
        chk[rule] = dict(E1_here=float(mine.E1), E1_A16=float(ref.E1_mean_wk), E2_here=float(mine.E2), E2_A16=float(ref.E2_power_change),
                         match=bool(abs(mine.E1 - ref.E1_mean_wk) < 1e-6 and abs(mine.E2 - ref.E2_power_change) < 1e-9))

    # Part A reading: smallest s (A grid plus A16 grid) with material effects
    a = cells[cells.part == "A"]
    first = {}
    for rule in T.CONFIRMS:
        rows = [dict(s=float(r.s), E1=r.E1, E1_lo=r.E1_lo, E1_hi=r.E1_hi, E2=r.E2, E2_hi=r.E2_hi) for r in a[a.rule == rule].itertuples()]
        for s in ("0.5", "0.75", "1.0"):
            q = a16[(a16.s == s) & (a16.confirm == rule)].iloc[0]
            rows.append(dict(s=float(s), E1=q.E1_mean_wk, E1_lo=q.E1_lo, E1_hi=q.E1_hi, E2=q.E2_power_change, E2_hi=q.E2_hi))
        rows.sort(key=lambda r: r["s"])
        s1 = next((r["s"] for r in rows if mat_e1(r)), None)
        s2 = next((r["s"] for r in rows if mat_e2(r)), None)
        first[rule] = dict(smallest_s_material_E1=s1, smallest_s_material_E2_power_loss=s2)
    pa1 = first["none"]["smallest_s_material_E1"]
    pa2 = first["none"]["smallest_s_material_E2_power_loss"]
    # Part B reading
    b = cells[(cells.part == "B") & (cells.rule == "none")]
    pb = {}
    for sg in (0.10, 0.25, 0.50):
        e0 = b[(np.isclose(b.sigma, sg)) & (b.rho == 0.0)].iloc[0].E1
        e8 = b[(np.isclose(b.sigma, sg)) & (b.rho == 0.8)].iloc[0].E1
        pb[str(sg)] = dict(abs_E1_rho0=abs(float(e0)), abs_E1_rho08=abs(float(e8)), smaller_at_rho08=bool(abs(e8) < abs(e0)))
    summary = dict(
        label="sensitivity analysis; simulation parts use assumed parameters; Part D is post hoc on LUMIERE; not evidence about any real trial",
        amendment="A18 in analysis_plan_v1.md", seed=SEED, R=R, N=N, D=D, runtime_s=round(time.time() - t0, 1),
        reproduction_check_vs_A16=chk, partA_first_material=first,
        PA1_supported=bool(pa1 is not None and pa1 <= 0.15), PA2_supported=bool(pa2 is not None and pa2 <= 0.15),
        partB_single_scan=pb, PB1_supported=bool(all(v["smaller_at_rho08"] for v in pb.values())),
        partD=dres,
    )
    lows = [dres["floor_sensitivity"][k]["ci95"][0] for k in ("250", "500", "1000")] + [dres["nadir_reference"]["ci95"][0]]
    summary["PD1_lower_CI_bounds_floor250_500_1000_nadir"] = lows
    summary["PD1_supported"] = bool(all(x >= 0.20 for x in lows))
    (OUT / "rano_sensitivity_summary.json").write_text(json.dumps(summary, indent=2, default=float))
    print(json.dumps({k: summary[k] for k in summary if k != "partD"}, indent=2, default=float))
    print(json.dumps(dres, indent=2, default=float))


if __name__ == "__main__":
    main()
