"""Script 141: trial simulation, Amendment 16 (analysis_plan_v1.md).

SIMULATION ONLY. All parameters are assumed or taken from the A15 noise bracket.
Implements A16.2 to A16.5 exactly. Run once; no re-runs with other settings (A16.6).
Writes output/trial_sim/trial_sim_cells.csv and trial_sim_summary.json.
"""
import json
import time
from pathlib import Path as _Path

import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "trial_sim"

SEED = 20261018
R = 2000
LN14 = float(np.log(1.40))
G_MED = LN14 / 34.8          # per week, control median g
G_SD = 0.6
LNV0_MU, LNV0_SD = float(np.log(4000.0)), 0.8
HORIZON = 104
S_GRID = [0.0, 0.25, 0.5, 0.75, 1.0, "vol"]
CONFIRMS = ["none", "two-scan"]
THETAS = [1.0, 0.7]
NS = [50, 100]
DS = [8, 12]
VOL_EDGES = (2503.39, 5118.34)       # mm3, A15 size bands
VOL_SD = (1.286, 1.133, 0.718)       # A15 upper-bound change SD per band
CHI2_CRIT = 3.8415


def change_sd(s, true_lnv):
    """SD of the CHANGE error of one tool, per scan, same shape as true_lnv."""
    if s == "vol":
        v = np.exp(true_lnv)
        return np.where(v < VOL_EDGES[0], VOL_SD[0], np.where(v < VOL_EDGES[1], VOL_SD[1], VOL_SD[2]))
    return np.full(true_lnv.shape, float(s))


def event_index(lnvhat, confirm):
    """First progression scan index (1..J) or J if censored. Returns (idx, event)."""
    nt, J1 = lnvhat.shape[1], lnvhat.shape[2]
    J = J1 - 1
    ref = np.minimum.accumulate(lnvhat, axis=2)          # min over scans 0..j
    crit = np.zeros(lnvhat.shape, dtype=bool)
    # criterion at scan j against the nadir of earlier scans (0..j-1)
    crit[:, :, 1:] = (lnvhat[:, :, 1:] - ref[:, :, :-1]) >= LN14
    if confirm == "two-scan":
        conf = np.zeros(lnvhat.shape, dtype=bool)
        # confirm at j+1 vs the SAME ref_j (ref[:, :, j-1])
        conf[:, :, 1:J] = crit[:, :, 1:J] & ((lnvhat[:, :, 2:] - ref[:, :, :-1][:, :, :J - 1]) >= LN14)
        hit = conf
    else:
        hit = crit
    anyhit = hit.any(axis=2)
    first = np.argmax(hit, axis=2)
    idx = np.where(anyhit, first, J)
    return idx, anyhit


def km_median_grid(idx, ev, times):
    """Kaplan-Meier median on the scan grid per trial; NaN if not reached."""
    Rr = idx.shape[0]
    J = len(times) - 1
    S = np.ones(Rr)
    med = np.full(Rr, np.nan)
    for j in range(1, J + 1):
        n = (idx >= j).sum(1)
        d = ((idx == j) & ev).sum(1)
        with np.errstate(divide="ignore", invalid="ignore"):
            S = S * np.where(n > 0, 1.0 - d / np.maximum(n, 1), 1.0)
        newly = np.isnan(med) & (S <= 0.5)
        med[newly] = times[j]
    return med


def km_median_continuous(T, ev):
    """Kaplan-Meier median for continuous times, per trial (rows)."""
    order = np.argsort(T, axis=1)
    Ts = np.take_along_axis(T, order, axis=1)
    es = np.take_along_axis(ev, order, axis=1)
    n = T.shape[1]
    at_risk = n - np.arange(n)
    S = np.cumprod(1.0 - es / at_risk, axis=1)
    below = S <= 0.5
    has = below.any(1)
    first = np.argmax(below, axis=1)
    med = np.where(has, Ts[np.arange(T.shape[0]), first], np.nan)
    return med


def logrank(idx_c, ev_c, idx_e, ev_e, J):
    """Log-rank on the scan grid. Group 1 = experimental. Returns chi2, peto log HR, valid."""
    Rr = idx_c.shape[0]
    O = np.zeros(Rr)
    E = np.zeros(Rr)
    V = np.zeros(Rr)
    for j in range(1, J + 1):
        n1 = (idx_e >= j).sum(1).astype(float)
        n0 = (idx_c >= j).sum(1).astype(float)
        n = n0 + n1
        d1 = ((idx_e == j) & ev_e).sum(1).astype(float)
        d0 = ((idx_c == j) & ev_c).sum(1).astype(float)
        d = d0 + d1
        ok = (n > 1) & (d > 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            e1 = np.where(ok, d * n1 / n, 0.0)
            v = np.where(ok, d * (n1 / n) * (n0 / n) * (n - d) / (n - 1), 0.0)
        O += d1
        E += e1
        V += v
    valid = V > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = np.where(valid, (O - E) ** 2 / V, 0.0)
        lhr = np.where(valid, (O - E) / V, np.nan)
    return chi2, lhr, valid


def mc_ci(x):
    x = x[~np.isnan(x)]
    m = float(x.mean())
    se = float(x.std(ddof=1) / np.sqrt(len(x)))
    return m, m - 1.96 * se, m + 1.96 * se


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    e0_rows = []
    for D in DS:
        times = np.arange(0, HORIZON + 1, D, dtype=float)
        J = len(times) - 1
        for N in NS:
            rng = np.random.default_rng([SEED, N, D])
            # common random numbers for this (N, D)
            Zc = rng.standard_normal((R, N)); Ze = rng.standard_normal((R, N))
            L0c = LNV0_MU + LNV0_SD * rng.standard_normal((R, N))
            L0e = LNV0_MU + LNV0_SD * rng.standard_normal((R, N))
            Uc = rng.standard_normal((R, N, J + 1)); Ue = rng.standard_normal((R, N, J + 1))
            gc = G_MED * np.exp(G_SD * Zc)
            ge_base = G_MED * np.exp(G_SD * Ze)
            lnv_c = L0c[:, :, None] + gc[:, :, None] * times[None, None, :]
            Ttrue_c = LN14 / gc
            res = {}  # (s, confirm) -> dict
            for conf in CONFIRMS:
                for s in S_GRID:
                    sd_c = change_sd(s, lnv_c) / np.sqrt(2.0)
                    idx_c, ev_c = event_index(lnv_c + sd_c * Uc, conf)
                    med_c = km_median_grid(idx_c, ev_c, times)
                    t_obs = times[idx_c]
                    early = (ev_c & (t_obs < (Ttrue_c - D))).mean()
                    lr = {}
                    for th in THETAS:
                        ge = ge_base * th
                        lnv_e = L0e[:, :, None] + ge[:, :, None] * times[None, None, :]
                        sd_e = change_sd(s, lnv_e) / np.sqrt(2.0)
                        idx_e, ev_e = event_index(lnv_e + sd_e * Ue, conf)
                        lr[th] = logrank(idx_c, ev_c, idx_e, ev_e, J)
                    res[(s, conf)] = dict(med=med_c, early=float(early), lr=lr, idx=idx_c, ev=ev_c)
            # discretisation check E0 (s = 0) against continuous KM median, same censoring
            tj = float(times[-1])
            Tc = np.minimum(Ttrue_c, tj); evc = Ttrue_c <= tj
            med_cont = km_median_continuous(Tc, evc.astype(float))
            for conf in CONFIRMS:
                m0 = res[(0.0, conf)]["med"]
                ok = ~np.isnan(m0) & ~np.isnan(med_cont)
                e0_rows.append(dict(N=N, D=D, confirm=conf, E0_mean_weeks=float((m0[ok] - med_cont[ok]).mean()),
                                    n_used=int(ok.sum()), expected_order_D_over_2=D / 2))
            # assemble cells
            for conf in CONFIRMS:
                base = res[(0.0, conf)]
                for s in S_GRID:
                    cur = res[(s, conf)]
                    ok = ~np.isnan(cur["med"]) & ~np.isnan(base["med"])
                    e1 = (cur["med"] - base["med"])[ok]
                    e1m = mc_ci(e1) if ok.sum() > 2 else (np.nan, np.nan, np.nan)
                    for th in THETAS:
                        chi2, lhr, valid = cur["lr"][th]
                        chi2b, lhrb, validb = base["lr"][th]
                        rej = chi2 > CHI2_CRIT
                        rejb = chi2b > CHI2_CRIT
                        row = dict(N=N, D=D, s=str(s), confirm=conf, theta=th,
                                   reject_rate=float(rej.mean()), reject_rate_s0=float(rejb.mean()),
                                   E1_mean_wk=e1m[0], E1_lo=e1m[1], E1_hi=e1m[2],
                                   E1_n_used=int(ok.sum()), E1_n_excluded=int(R - ok.sum()),
                                   ctrl_median_mean_wk=float(np.nanmean(cur["med"])),
                                   E5_early_share=cur["early"])
                        if th == 0.7:
                            d = rej.astype(float) - rejb.astype(float)
                            row["E2_power_change"], row["E2_lo"], row["E2_hi"] = mc_ci(d)
                            okh = valid & validb
                            dh = np.where(okh, lhr - lhrb, np.nan)
                            row["E4_atten"], row["E4_lo"], row["E4_hi"] = mc_ci(dh)
                            row["E4_n_used"] = int(okh.sum())
                        rows.append(row)
    cells = pd.DataFrame(rows)
    cells.to_csv(OUT / "trial_sim_cells.csv", index=False)
    pd.DataFrame(e0_rows).to_csv(OUT / "trial_sim_E0.csv", index=False)

    # ---- declared readings (A16.5) ----
    def pick(**kw):
        d = cells
        for k, v in kw.items():
            d = d[d[k] == v]
        return d

    sel = ["0.5", "0.75", "1.0"]
    h1 = {}
    for s in sel:
        r = pick(N=100, D=8, s=s, confirm="none", theta=1.0).iloc[0]
        h1[s] = dict(E1=r.E1_mean_wk, lo=r.E1_lo, hi=r.E1_hi, supported=bool(r.E1_mean_wk < 0 and r.E1_hi < 0))
    h2 = {}
    for s in sel:
        r = pick(N=100, D=8, s=s, confirm="none", theta=0.7).iloc[0]
        h2[s] = dict(E4=r.E4_atten, lo=r.E4_lo, hi=r.E4_hi, supported=bool(r.E4_atten > 0 and r.E4_lo > 0))
    a = pick(N=100, D=8, s="0.75", confirm="none", theta=1.0).iloc[0]
    b = pick(N=100, D=8, s="0.75", confirm="two-scan", theta=1.0).iloc[0]
    m1 = pick(N=100, D=8, s="0.0", confirm="none", theta=1.0).iloc[0]
    m2 = pick(N=100, D=8, s="0.0", confirm="two-scan", theta=1.0).iloc[0]
    h3 = dict(abs_E1_single=abs(a.E1_mean_wk), abs_E1_two=abs(b.E1_mean_wk),
              ctrl_median_s0_single=m1.ctrl_median_mean_wk, ctrl_median_s0_two=m2.ctrl_median_mean_wk,
              supported=bool(abs(b.E1_mean_wk) < abs(a.E1_mean_wk) and m2.ctrl_median_mean_wk > m1.ctrl_median_mean_wk))
    t1 = cells[cells.theta == 1.0]
    out_of_band = t1[(t1.reject_rate < 0.04) | (t1.reject_rate > 0.06)]
    h4 = dict(n_cells=int(len(t1)), n_outside=int(len(out_of_band)),
              outside=out_of_band[["N", "D", "s", "confirm", "reject_rate"]].to_dict("records"),
              supported=bool(len(out_of_band) == 0))
    material = cells[(cells.theta == 0.7)]
    summary = dict(
        label="simulation, assumed parameters; not evidence about any real trial",
        amendment="A16 in analysis_plan_v1.md", seed=SEED, R=R, n_cells=int(len(cells)),
        runtime_s=round(time.time() - t0, 1),
        H_S1=h1, H_S2=h2, H_S3=h3, H_S4=h4,
        E0_discretisation=e0_rows,
        E1_material_cells=int(((cells.E1_mean_wk.abs() >= 4.345) & ((cells.E1_lo > 0) | (cells.E1_hi < 0))).sum() / 2),
        E2_material_loss_cells=int(((material.E2_power_change <= -0.05) & (material.E2_hi < 0)).sum()),
    )
    (OUT / "trial_sim_summary.json").write_text(json.dumps(summary, indent=2, default=float))
    print(json.dumps({k: summary[k] for k in ("n_cells", "runtime_s", "H_S1", "H_S2", "H_S3")}, indent=2, default=float))
    print("H_S4:", h4["supported"], "outside:", h4["n_outside"], "of", h4["n_cells"])


if __name__ == "__main__":
    main()
