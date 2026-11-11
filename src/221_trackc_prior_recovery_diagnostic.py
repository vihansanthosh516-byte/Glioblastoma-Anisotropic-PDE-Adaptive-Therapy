"""Outcome-blind particle-filter parameter-recovery diagnostic (report-only).

Frozen in paper/trackc_amendment_pf_prior_support_2026-11-10.md, revision 3, section 4.
Compares PF_PRIOR_NARROW_V2 with PF_PRIOR_INTEGRATED_V3 on synthetic patients under a
fixed treatment and scan schedule. No controller runs and no utility is computed.
"""
import hashlib
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

import trackc_integrated_sampler as smp
import trackc_sim as tc
from trackc_filter_priors import PF_PRIOR_INTEGRATED_V3, PF_PRIOR_NARROW_V2

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = PROJECT_ROOT / "output" / "trackc_prior_recovery_diagnostic_20261110"

WORLDS_PER_FAMILY = 50
PATIENTS_PER_WORLD = 4
WORLD_SEED0, COHORT_SEED0, FILTER_SEED0, NOISE_KEY0, BOOT_KEY = 20263000, 20263100, 20263200, 2630000, 2630100
SIGMA_CYCLE = (0.05, 0.10, 0.20)
TREAT_CYCLES = (0, 2, 4, 6, 8, 10)
SCAN_CYCLES = (0, 3, 6, 9, 12)
LAST_CYCLE = 12
N_PARTICLES = 500
N_BOOT = 2000
PRIORS = {"narrow_v2": PF_PRIOR_NARROW_V2, "v3": PF_PRIOR_INTEGRATED_V3}


def weighted_quantile(values, weights, q):
    order = np.argsort(values, 1)
    v = np.take_along_axis(values, order, 1)
    cw = np.cumsum(np.take_along_axis(weights, order, 1), 1)
    return v[np.arange(v.shape[0]), np.argmax(cw >= q, 1)]


def simulate_world(index):
    family = "single" if index < WORLDS_PER_FAMILY else "mixture"
    world = smp.draw_world_parameters(family, WORLD_SEED0 + index)
    scen = smp.draw_patient_cohort(world, PATIENTS_PER_WORLD, COHORT_SEED0 + index)
    n = PATIENTS_PER_WORLD
    scen["sigma"] = np.array([SIGMA_CYCLE[k % 3] for k in range(n)])
    truth = {k: np.asarray(scen[k], float) for k in tc.PARAMS}
    p = dict(truth, a_sr=np.ones(n), a_rs=np.ones(n))
    noise = np.random.default_rng(NOISE_KEY0 + index).standard_normal((len(SCAN_CYCLES), n))
    S, R = (1 - truth["f0"]) * truth["B0"], truth["f0"] * truth["B0"]
    obs = {}
    for j in range(LAST_CYCLE + 1):
        if j in SCAN_CYCLES:
            obs[j] = np.log(S + R) + scen["sigma"] * noise[SCAN_CYCLES.index(j)]
        if j < LAST_CYCLE:
            S, R = tc.run_cycle(S, R, p, np.full(n, j in TREAT_CYCLES), j * tc.CYCLE_DAYS)

    result = {"family": family, "world_index": index, "truth": truth, "priors": {}}
    for label, ranges in PRIORS.items():
        start = time.perf_counter()
        pf = tc.ParticleFilter(scen["sigma"], N_PARTICLES, np.random.default_rng(FILTER_SEED0 + index), ranges=ranges)
        ess = []
        for j in range(LAST_CYCLE + 1):
            if j in SCAN_CYCLES:
                B = np.maximum(pf.S + pf.R, 1e-300)
                lw = pf.logw - 0.5 * ((obs[j][:, None] - np.log(B)) / pf.sigma) ** 2
                w = np.exp(lw - lw.max(1, keepdims=True))
                w /= w.sum(1, keepdims=True)
                ess.append(1.0 / (w ** 2).sum(1))
                pf.update(j, obs[j])
            if j < LAST_CYCLE:
                pf.propagate(j, np.full(n, j in TREAT_CYCLES))
        w = pf.w
        lo = {k: weighted_quantile(pf.theta[k], w, 0.05) for k in tc.PARAMS}
        hi = {k: weighted_quantile(pf.theta[k], w, 0.95) for k in tc.PARAMS}
        ess = np.array(ess)
        result["priors"][label] = {
            "lo": lo, "hi": hi, "ess_min": ess.min(0), "ess_median": np.median(ess, 0),
            "n_resample_events": int(pf.n_resample), "seconds": time.perf_counter() - start,
        }
    return result


def per_patient_table(results):
    rows = []
    for res in results:
        for i in range(PATIENTS_PER_WORLD):
            t = {k: float(res["truth"][k][i]) for k in tc.PARAMS}
            row = {"family": res["family"], "world_index": res["world_index"], "patient": i, "truth": t,
                   "truth_outside_narrow": bool(t["r"] < PF_PRIOR_NARROW_V2["r"][0] or t["r"] > PF_PRIOR_NARROW_V2["r"][1]
                                                or t["B0"] < PF_PRIOR_NARROW_V2["B0"][0] or t["B0"] > PF_PRIOR_NARROW_V2["B0"][1])}
            for label, pr in res["priors"].items():
                kind = PRIORS[label]
                cov = {k: bool(pr["lo"][k][i] <= t[k] <= pr["hi"][k][i]) for k in tc.PARAMS}
                width = {k: float(np.log(pr["hi"][k][i]) - np.log(pr["lo"][k][i])) if kind[k][2] == "log"
                         else float(pr["hi"][k][i] - pr["lo"][k][i]) for k in tc.PARAMS}
                row[label] = {"covered": cov, "joint_covered": all(cov.values()), "width": width,
                              "ess_min": float(pr["ess_min"][i]), "ess_median": float(pr["ess_median"][i])}
            rows.append(row)
    return rows


def summarize(rows):
    out = {}
    for label in PRIORS:
        out[label] = {
            "marginal_coverage": {k: float(np.mean([r[label]["covered"][k] for r in rows])) for k in tc.PARAMS},
            "joint_coverage": float(np.mean([r[label]["joint_covered"] for r in rows])),
            "median_width": {k: float(np.median([r[label]["width"][k] for r in rows])) for k in tc.PARAMS},
            "ess_min_median": float(np.median([r[label]["ess_min"] for r in rows])),
            "ess_median_median": float(np.median([r[label]["ess_median"] for r in rows])),
        }
    out["n_patients"] = len(rows)
    return out


def paired_stats(rows):
    """Per-patient paired differences V3 minus narrow."""
    stats = {}
    for k in tc.PARAMS:
        stats[f"coverage_{k}"] = [r["v3"]["covered"][k] - r["narrow_v2"]["covered"][k] for r in rows]
        stats[f"width_{k}"] = [r["v3"]["width"][k] - r["narrow_v2"]["width"][k] for r in rows]
    stats["joint_coverage"] = [r["v3"]["joint_covered"] - r["narrow_v2"]["joint_covered"] for r in rows]
    return {key: np.asarray(v, float) for key, v in stats.items()}


def cluster_bootstrap(rows, rng):
    worlds = sorted({r["world_index"] for r in rows})
    by_world = {w: [i for i, r in enumerate(rows) if r["world_index"] == w] for w in worlds}
    stats = paired_stats(rows)
    idx_lists = [by_world[w] for w in worlds]
    draws = {key: np.empty(N_BOOT) for key in stats}
    for b in range(N_BOOT):
        pick = rng.integers(0, len(worlds), len(worlds))
        idx = np.concatenate([idx_lists[p] for p in pick])
        for key, v in stats.items():
            draws[key][b] = v[idx].mean()
    return {key: {"mean": float(v.mean()), "ci95": [float(np.percentile(draws[key], 2.5)),
                                                    float(np.percentile(draws[key], 97.5))]}
            for key, v in stats.items()}


def git_head():
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, capture_output=True,
                              text=True, check=True).stdout.strip()
    except Exception:
        return None


def main():
    t0 = time.perf_counter()
    results = [simulate_world(i) for i in range(2 * WORLDS_PER_FAMILY)]
    rows = per_patient_table(results)
    boot_rng = np.random.default_rng(BOOT_KEY)
    report = {"families": {}}
    for family in ("single", "mixture"):
        fam = [r for r in rows if r["family"] == family]
        report["families"][family] = {
            "all": summarize(fam),
            "truth_outside_narrow": summarize([r for r in fam if r["truth_outside_narrow"]]),
            "truth_inside_narrow": summarize([r for r in fam if not r["truth_outside_narrow"]]),
            "paired_v3_minus_narrow_world_cluster_bootstrap": cluster_bootstrap(fam, boot_rng),
        }
    seconds = {label: float(sum(res["priors"][label]["seconds"] for res in results)) for label in PRIORS}
    resamples = {label: int(sum(res["priors"][label]["n_resample_events"] for res in results)) for label in PRIORS}
    doc = {
        "status": "completed_report_only",
        "protocol": "paper/trackc_amendment_pf_prior_support_2026-11-10.md revision 3 section 4",
        "no_policy_or_utility": True,
        "reporting_rule": "report-only; no pass threshold; results cannot change bounds or settings",
        "git_head_at_run": git_head(),
        "python": sys.version.split()[0], "numpy": np.__version__, "platform": platform.platform(),
        "settings": {"worlds_per_family": WORLDS_PER_FAMILY, "patients_per_world": PATIENTS_PER_WORLD,
                     "n_particles": N_PARTICLES, "treat_cycles": TREAT_CYCLES, "scan_cycles": SCAN_CYCLES,
                     "sigma_cycle": SIGMA_CYCLE, "n_boot": N_BOOT,
                     "width_difference_statistic": "mean of paired per-patient width differences (log width for log parameters)"},
        "seeds": {"world": [WORLD_SEED0, WORLD_SEED0 + 99], "cohort": [COHORT_SEED0, COHORT_SEED0 + 99],
                  "filter_policy_stream": [FILTER_SEED0, FILTER_SEED0 + 99],
                  "noise_keys": [NOISE_KEY0, NOISE_KEY0 + 99], "bootstrap_key": BOOT_KEY},
        "priors": {k: {p: list(v) for p, v in pr.items()} for k, pr in PRIORS.items()},
        "filter_seconds_total": seconds, "resample_events_total": resamples,
        "wall_seconds": time.perf_counter() - t0,
        "report": report,
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "diagnostic.json").write_text(json.dumps(doc, indent=2), encoding="utf-8")
    with open(OUT_DIR / "per_patient.jsonl", "w", encoding="utf-8") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
    digest = hashlib.sha256((OUT_DIR / "diagnostic.json").read_bytes()).hexdigest()
    print(f"wrote {OUT_DIR / 'diagnostic.json'} sha256={digest}")


if __name__ == "__main__":
    main()
