"""Track B compute feasibility pilot: call-count formulas plus a budget-bounded SYNTHETIC timing run.

Label for every result: "projection, not a bound; this PC; synthetic". No patient data is read. Each timed item
runs in a child process with a hard timeout; the whole run has a global deadline (default 25 min). A timed-out
item is recorded as censored (lower bound only) and never replaced silently.

PDE timing uses the LUMIERE full grid (91x109x91, 2 mm) but REDUCED posterior settings so that a call finishes
inside the budget (a full-settings call was recorded at 640-19,500 s in
output/track_b_timing_pilot_20261031/lumiere_history_gap_stress.json). The reduced call is scaled to the full
settings by the ratio of solver fields (count of rho values passed to TensorFK.run), an approximation.

Usage: python -m src.fk_track_b_feasibility_pilot --out <json path>
"""
from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
GRID_SHAPE = (91, 109, 91)
SPACING_MM = 2.0
LABEL = "projection, not a bound; this PC; synthetic"
FULL_SETTINGS = {"n_coarse": 21, "n_refine": 11, "n_draws": 200, "n_anchors": 256, "mass_fraction": 0.9,
                 "rho_bounds": (1e-5, 1e-1), "D_bounds": (1e-5, 10.0)}
REDUCED_SETTINGS = {**FULL_SETTINGS, "n_coarse": 5, "n_refine": 3, "n_draws": 20, "n_anchors": 64}
REPRESENTATIVE_GAP = {"history_weeks": 6, "forecast_weeks": 13}
EXTREME_GAP = {"history_weeks": 23, "forecast_weeks": 13}
SEED = 20261102


# ------------------------------------------------------------------ call-count formulas (pure, no timing)

def call_counts(n: int, *, n_central: int | None = None, outer_folds: int = 5, inner_folds: int = 5,
                prior_scenarios: int = 3, refits: int = 1, epochs: int = 300, n_seeds: int = 3,
                n_bootstrap: int = 10_000, n_rates: int = 20, film_reuse: bool = True) -> dict:
    """Counts for one full protocol (P prior-scenario analyses on n patients + R central-gap refits on n_central).

    Per analysis on m eligible patients with K outer folds (sum over folds of training size = (K-1) m):
      full forecasts            = m                      (one per held-out patient)
      inner calibration-only    = (K-1) m                (each training patient forecast once, cross-fitted;
                                                          independent of the inner fold count I)
      PDE calls                 = K m                    (upper bound: inner records need 3 usable masks)
      FiLM optimiser steps      = (K-1) m * epochs * seeds   per FiLM training set
      FiLM trainings            = K * seeds models per training set
      fixed-growth forecasts    = (K-1) m * n_rates      (rate selection) + m (test forecasts)
      volume-nuisance fits      = K (outer) + K * I (inner, cheap)
      bootstrap resamples       = B
    Prior scenarios share one FiLM cache (fingerprint ignores the provider); the central-gap refit has its own.
    """
    m_c = n if n_central is None else n_central
    K, P, R = outer_folds, prior_scenarios, refits

    def per_analysis(m):
        return {"full_forecasts": m, "inner_calibration_only": (K - 1) * m, "pde_calls": K * m,
                "film_steps_per_training_set": (K - 1) * m * epochs * n_seeds,
                "film_models_per_training_set": K * n_seeds,
                "fixed_growth_forecasts": (K - 1) * m * n_rates + m, "nuisance_fits": K + K * inner_folds,
                "bootstrap_resamples": n_bootstrap}

    a, c = per_analysis(n), per_analysis(m_c)
    film_sets = (1 if film_reuse else P) + R
    film_steps = (a["film_steps_per_training_set"] * (1 if film_reuse else P)
                  + c["film_steps_per_training_set"] * R)
    return {"n": n, "n_central": m_c, "outer_folds": K, "inner_folds": inner_folds, "prior_scenarios": P,
            "refits": R, "per_analysis_primary_size": a, "per_analysis_central_gap_size": c,
            "total": {"analyses": P + R,
                      "pde_calls": P * a["pde_calls"] + R * c["pde_calls"],
                      "full_forecasts": P * a["full_forecasts"] + R * c["full_forecasts"],
                      "inner_calibration_only": P * a["inner_calibration_only"] + R * c["inner_calibration_only"],
                      "film_training_sets": film_sets, "film_steps": film_steps,
                      "film_predictions": (P * n + R * m_c),
                      "fixed_growth_forecasts": P * a["fixed_growth_forecasts"] + R * c["fixed_growth_forecasts"],
                      "bootstrap_resamples": (P + R) * n_bootstrap,
                      "nuisance_fits": (P + R) * a["nuisance_fits"]}}


# ------------------------------------------------------------------ synthetic fixtures

def _domain():
    g = np.ogrid[tuple(slice(0, s) for s in GRID_SHAPE)]
    return sum(((g[i] - (GRID_SHAPE[i] - 1) / 2) / (0.45 * GRID_SHAPE[i])) ** 2 for i in range(3)) <= 1.0


def _ball(center, radius_mm):
    g = np.ogrid[tuple(slice(0, s) for s in GRID_SHAPE)]
    return sum(((g[i] - center[i]) * SPACING_MM) ** 2 for i in range(3)) <= radius_mm ** 2


def _masks(rng, domain):
    c = tuple(s // 2 + int(rng.integers(-8, 9)) for s in GRID_SHAPE)
    r = float(rng.uniform(9, 13))
    return c, _ball(c, r) & domain, _ball(c, r + 3) & domain, _ball(c, r + 6) & domain


# ------------------------------------------------------------------ child items (each prints one JSON line)

def item_pde(args):
    from src import fk_mask_inverse as fi
    from src.fk_hybrid_forecast import hybrid_posterior_predictive_forecast
    counter = {"fields": 0, "run_calls": 0}
    orig = fi.ria.TensorFK.run

    def counted(self, u0, rhos, *a, **k):
        counter["fields"] += int(np.size(rhos))
        counter["run_calls"] += 1
        return orig(self, u0, rhos, *a, **k)
    fi.ria.TensorFK.run = counted
    rng = np.random.default_rng(SEED)
    domain = _domain()
    _, m1, m2, _ = _masks(rng, domain)
    d_h = 7.0 * args.history_weeks
    t0 = time.perf_counter()
    try:
        hybrid_posterior_predictive_forecast(
            masks=[m1, m2], days=np.array([0.0, d_h]), target_day=d_h + 7.0 * args.forecast_weeks,
            prior={"distribution": "log_uniform"}, h=SPACING_MM,
            boundary_covariance={"sigma_b_mm": 0.8, "tau_b_mm": 0.3, "ell_mm": 6.0},
            population_log_volume_slope_per_day=0.002, shrinkage_time_constant_days=90.0,
            volume_variance_fit={"intercept": 0.02, "slope_per_day": 0.0001}, domain=domain, seed=SEED,
            calibration_only=args.calibration_only, **REDUCED_SETTINGS)
        status = "ok"
    except ValueError as exc:
        status = f"model failure: {exc}"
    return {"seconds": time.perf_counter() - t0, "status": status, "solver_fields": counter["fields"],
            "solver_run_calls": counter["run_calls"], "settings": {k: list(v) if isinstance(v, tuple) else v
                                                                   for k, v in REDUCED_SETTINGS.items()}}


def item_film(_args):
    import torch
    from dataclasses import replace
    from src import fk_time_aware_baselines as ta
    from src.fk_learned_baselines import UNET_RECIPE, unet_inputs
    rng = np.random.default_rng(SEED)
    domain = _domain()
    pairs = []
    for _ in range(3):
        _, m1, m2, m3 = _masks(rng, domain)
        pairs.append((unet_inputs(m1, m2, domain), m3))
    recipe = replace(UNET_RECIPE, epochs=1, n_seeds=1, seed_offsets=(0,))

    def train(k):
        t = time.perf_counter()
        out = ta.train_film_models(pairs[:k], [(13, 6)] * k, fold=0, recipe=recipe)
        return time.perf_counter() - t, out
    t_warm, _ = train(1)
    t1, _ = train(1)
    t3, trained = train(3)
    step = max((t3 - t1) / 2.0, 0.0)
    pred = None
    if trained.get("status") == "ok":
        t = time.perf_counter()
        ta.predict_film(trained, pairs[0][0], (13, 6), domain=domain)
        pred = time.perf_counter() - t
    return {"step_seconds": step, "train_1_pair_seconds": t1, "train_3_pair_seconds": t3,
            "warmup_train_1_pair_seconds": t_warm, "predict_one_model_seconds": pred,
            "film_status": trained.get("status"), "padded_shape": list(trained.get("padded_shape", [])),
            "torch_threads": torch.get_num_threads(), "torch_version": torch.__version__,
            "method": "(t(3 pairs) - t(1 pair)) / 2 with epochs=1, 1 seed; fixed overhead cancels"}


def item_fixed_growth(args):
    from src import fk_time_aware_baselines as ta
    rng = np.random.default_rng(SEED)
    domain = _domain()
    cases = []
    for i in range(args.n_train):
        _, m1, m2, m3 = _masks(rng, domain)
        cases.append(ta.TimedGrowthCase(f"s{i}", m2, domain, (SPACING_MM,) * 3, 13, m3))
    t = time.perf_counter()
    ta.select_rate(cases)
    s = time.perf_counter() - t
    return {"n_train": args.n_train, "n_rates": len(ta.RATE_GRID_MM_PER_DAY), "select_rate_seconds": s,
            "seconds_per_training_patient": s / args.n_train}


def item_bootstrap(args):
    from src.fk_lumiere_stats import N_BOOTSTRAP, max_t_pass_test
    rng = np.random.default_rng(SEED)
    comps = ("copy_last", "time_aware_fixed_growth", "pde_only", "film_unet")
    ids = [f"s{i:03d}" for i in range(args.n_train)]
    dice = {k: {p: float(rng.uniform(0.2, 0.6)) for p in ids} for k in comps + ("candidate",)}
    t = time.perf_counter()
    out = max_t_pass_test(dice, comps, n_bootstrap=N_BOOTSTRAP)
    return {"n_patients": len(ids), "n_bootstrap": N_BOOTSTRAP, "seconds": time.perf_counter() - t,
            "status": out.get("status")}


ITEMS = {"pde": item_pde, "film": item_film, "fixed_growth": item_fixed_growth, "bootstrap": item_bootstrap}


def child_main(argv):
    p = argparse.ArgumentParser()
    p.add_argument("--item", required=True)
    p.add_argument("--history-weeks", type=int, default=6)
    p.add_argument("--forecast-weeks", type=int, default=13)
    p.add_argument("--calibration-only", action="store_true")
    p.add_argument("--n-train", type=int, default=22)
    a = p.parse_args(argv)
    print("RESULT_JSON " + json.dumps(ITEMS[a.item](a)))


# ------------------------------------------------------------------ orchestrator and projection

def _run_item(name, extra, cap_seconds, deadline):
    remaining = deadline - time.monotonic()
    timeout = max(1.0, min(cap_seconds, remaining))
    cmd = [sys.executable, "-m", "src.fk_track_b_feasibility_pilot", "--child"] + ["--item", name] + extra
    t = time.perf_counter()
    try:
        proc = subprocess.run(cmd, cwd=PROJECT_ROOT, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"status": "censored: exceeded budget", "timeout_seconds": timeout,
                "note": "lower bound only: true time > timeout"}
    line = [x for x in proc.stdout.splitlines() if x.startswith("RESULT_JSON ")]
    if proc.returncode != 0 or not line:
        return {"status": f"error rc={proc.returncode}", "stderr_tail": proc.stderr[-600:]}
    out = json.loads(line[0][len("RESULT_JSON "):])
    out["wall_including_import_seconds"] = time.perf_counter() - t
    out["status_item"] = "measured"
    return out


def full_field_upper_bound(settings=FULL_SETTINGS) -> int:
    """Coarse grid + refine grid + at most n_draws unique draws (one solver field each); an upper-ish scale."""
    return int(settings["n_coarse"] ** 2 + settings["n_refine"] ** 2 + settings["n_draws"])


def project(measured: dict, ns=(20, 30, 40, 61), *, central_fraction=0.8, bootstrap_seconds_at=None) -> list:
    rows = []
    film = measured.get("film", {})
    fg = measured.get("fixed_growth", {})
    bs = measured.get("bootstrap", {})
    ok = lambda d: isinstance(d, dict) and d.get("status_item") == "measured"
    pde_s = {}
    for label in ("representative", "extreme"):
        e = measured.get(f"pde_{label}_full")
        c = measured.get(f"pde_{label}_cal")
        if ok(e) and ok(c) and e.get("solver_fields"):
            scale = full_field_upper_bound() / e["solver_fields"]
            pde_s[label] = {"full": e["seconds"] * scale, "cal": c["seconds"] * (
                full_field_upper_bound() / max(c["solver_fields"], 1)), "scale": scale}
    for n in ns:
        nc = int(round(central_fraction * n))
        cc = call_counts(n, n_central=nc)
        t = cc["total"]
        row = {"N": n, "N_central_gap": nc, "counts": t}
        film_h = (t["film_steps"] * film["step_seconds"] / 3600.0) if ok(film) else None
        # seconds_per_training_patient already covers all rates for one patient -> forecasts/n_rates patients
        fg_h = (t["fixed_growth_forecasts"] / 20.0 * fg["seconds_per_training_patient"] / 3600.0) if ok(fg) else None
        bs_h = (t["bootstrap_resamples"] / bs["n_bootstrap"] * bs["seconds"] * (n / bs["n_patients"]) / 3600.0
                if ok(bs) else None)
        row["hours"] = {"film_training": film_h, "fixed_growth": fg_h, "bootstrap": bs_h}
        for label, v in pde_s.items():
            pde_h = (t["full_forecasts"] * v["full"] + t["inner_calibration_only"] * v["cal"]) / 3600.0
            row["hours"][f"pde_{label}"] = pde_h
            row["hours"][f"total_{label}_gap"] = pde_h + sum(x for x in (film_h, fg_h, bs_h) if x is not None)
        rows.append(row)
    return rows


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True)
    p.add_argument("--budget-minutes", type=float, default=25.0)
    a = p.parse_args(argv)
    out_path = Path(a.out)
    if out_path.exists():
        raise FileExistsError(f"refusing to overwrite {out_path}")
    start = time.monotonic()
    deadline = start + a.budget_minutes * 60.0
    measured = {}
    plan = [
        ("film", "film", [], 360),
        ("fixed_growth", "fixed_growth", ["--n-train", "22"], 120),
        ("bootstrap", "bootstrap", ["--n-train", "40"], 120),
        ("pde_representative_full", "pde", ["--history-weeks", "6", "--forecast-weeks", "13"], 360),
        ("pde_representative_cal", "pde", ["--history-weeks", "6", "--forecast-weeks", "13", "--calibration-only"], 360),
        ("pde_extreme_full", "pde", ["--history-weeks", "23", "--forecast-weeks", "13"], 420),
        ("pde_extreme_cal", "pde", ["--history-weeks", "23", "--forecast-weeks", "13", "--calibration-only"], 420),
    ]
    for key, item, extra, cap in plan:
        measured[key] = _run_item(item, extra, cap, deadline)
    record = {
        "schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(), "label": LABEL,
        "evidence_scope": "synthetic full-grid cases only; no patient data read; not patient-result evidence",
        "grid_shape": list(GRID_SHAPE), "spacing_mm": SPACING_MM,
        "budget_minutes": a.budget_minutes, "elapsed_seconds": time.monotonic() - start,
        "machine": {"platform": platform.platform(), "python": platform.python_version(),
                    "processor": platform.processor()},
        "pde_settings_full": {k: list(v) if isinstance(v, tuple) else v for k, v in FULL_SETTINGS.items()},
        "pde_settings_timed_reduced": {k: list(v) if isinstance(v, tuple) else v for k, v in REDUCED_SETTINGS.items()},
        "gaps": {"representative": REPRESENTATIVE_GAP, "extreme": EXTREME_GAP,
                 "note": "155-week history gap was recorded at 16,098-19,549 s per full-settings call in "
                         "output/track_b_timing_pilot_20261031/lumiere_history_gap_stress.json; not re-timed here"},
        "scale_to_full_settings": {"method": "reduced-call seconds x (n_coarse^2 + n_refine^2 + n_draws) / counted "
                                             "solver fields; approximation (cost per field depends on D and horizon)",
                                   "full_field_upper_bound": full_field_upper_bound()},
        "measured": measured,
        "formulas_example_N28": call_counts(28, n_central=22),
        "projection_assumptions": [
            "central-gap subset = 0.8 N (LUMIERE metadata subset was 22 of 28)",
            "K=5 outer folds, 5 inner folds, 3 prior scenarios + 1 central-gap refit, 300 epochs x 3 seeds, B=10,000",
            "FiLM reused across the 3 prior scenarios (one training set), central-gap refit trains its own",
            "all PDE calls timed at one gap (representative or extreme), serial, no caching, this PC",
            "count model, FiLM prediction, PDE-only cutoff selection and file I/O not timed (small)",
        ],
        "projection": project(measured),
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(record, indent=1, default=str))
    print(f"wrote {out_path}")


if __name__ == "__main__":
    if "--child" in sys.argv:
        sys.argv.remove("--child")
        child_main(sys.argv[1:])
    else:
        main()
