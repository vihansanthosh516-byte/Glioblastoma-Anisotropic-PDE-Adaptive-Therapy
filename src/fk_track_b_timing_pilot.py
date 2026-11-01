"""Synthetic-only runtime pilot for the Track B MU runner arms. No patient data is read.

Each arm runs in its own subprocess so peak memory is per arm. The U-Net arm times real training steps with the
frozen recipe (architecture, optimizer, augmentation, loss) on full-size synthetic volumes; it does not run all
300 epochs x 3 seeds, which would take many hours. Synthetic timing is not patient-result evidence.

Usage:  python -m src.fk_track_b_timing_pilot            (all arms, writes the record)
        python -m src.fk_track_b_timing_pilot --arm unet (one arm, prints JSON)
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output" / "track_b_timing_pilot_20261031"

GRID_SHAPE = (120, 120, 78)  # SRI24 240x240x155 at 1 mm, resampled to 2 mm
SPACING_MM = 2.0
N_PATIENTS = 18              # frozen MU development cohort size
FOLD_TRAIN_SIZES = (14, 14, 14, 15, 15)  # 18 patients in folds of 4,4,4,3,3
UNET_TIMED_STEPS = 12
SEED = 20261031


def _ball(center, radius_mm, shape=GRID_SHAPE):
    g = np.ogrid[tuple(slice(0, s) for s in shape)]
    return sum(((g[i] - center[i]) * SPACING_MM) ** 2 for i in range(3)) <= radius_mm ** 2


def synthetic_cohort(n=N_PATIENTS):
    rng = np.random.default_rng(SEED)
    domain = sum(((np.ogrid[tuple(slice(0, s) for s in GRID_SHAPE)][i] - (GRID_SHAPE[i] - 1) / 2)
                  / (0.45 * GRID_SHAPE[i])) ** 2 for i in range(3)) <= 1.0
    patients = []
    for i in range(n):
        c = (60 + int(rng.integers(-12, 13)), 60 + int(rng.integers(-12, 13)), 39 + int(rng.integers(-6, 7)))
        r1 = float(rng.uniform(8, 14))
        r2 = r1 + float(rng.uniform(0, 4))
        r3 = r2 + float(rng.uniform(-2, 6))
        gap = 76.0 + float(rng.integers(0, 29))
        patients.append({"id": f"syn{i:02d}", "days": [0.0, 30.0, 30.0 + gap],
                         "masks": [_ball(c, r1) & domain, _ball(c, r2) & domain, _ball(c, max(r3, 4)) & domain]})
    return patients, domain


def _peak_gb():
    import psutil
    info = psutil.Process().memory_info()
    return float(getattr(info, "peak_wset", info.rss)) / 1e9


def arm_unet(crop_voxels=None, timed_steps=UNET_TIMED_STEPS):
    import torch
    from src import fk_learned_baselines as lb
    patients, domain = synthetic_cohort()
    r = lb.UNET_RECIPE
    pairs = []
    for p in patients[:2]:
        x, y = lb.unet_inputs(p["masks"][0], p["masks"][1], domain), p["masks"][2]
        if crop_voxels:
            box = lb.history_crop_box(p["masks"][1], GRID_SHAPE, crop_voxels)
            x, y = x[(slice(None),) + box], y[box]
        pairs.append((x, y))
    padded = lb._padded_shape(pairs[0][1].shape, len(r.widths))
    xs = [torch.from_numpy(np.stack([lb._pad(c, padded) for c in x]))[None] for x, _ in pairs]
    ys = [torch.from_numpy(lb._pad(y.astype(np.float32), padded))[None, None] for _, y in pairs]
    torch.use_deterministic_algorithms(True)
    torch.manual_seed(r.base_seed)
    gen = np.random.default_rng(r.base_seed)
    model = lb.UNet3D(r.widths, r.leaky_slope)
    opt = torch.optim.Adam(model.parameters(), lr=r.learning_rate, betas=r.betas, weight_decay=r.weight_decay)
    model.train()
    times = []
    for step in range(timed_steps + 2):
        x, y = xs[step % 2], ys[step % 2]
        for axis in (2, 3, 4):
            if gen.random() < 0.5:
                x, y = torch.flip(x, (axis,)), torch.flip(y, (axis,))
        k = int(gen.integers(4))
        x, y = torch.rot90(x, k, (2, 3)), torch.rot90(y, k, (2, 3))
        t = time.perf_counter()
        opt.zero_grad()
        loss = lb._loss(model(x), y, r.dice_smoothing)
        loss.backward()
        opt.step()
        times.append(time.perf_counter() - t)
    steady = times[2:]  # first two steps are warm-up
    model.eval()
    t = time.perf_counter()
    with torch.no_grad():
        for _ in range(r.n_seeds):
            torch.sigmoid(model(xs[0]))
    predict_s = time.perf_counter() - t
    return {"crop_voxels": crop_voxels, "padded_shape": list(padded), "n_parameters": int(sum(p.numel() for p in model.parameters())),
            "timed_steps": len(steady), "step_seconds_mean": float(np.mean(steady)),
            "step_seconds_min": float(np.min(steady)), "step_seconds_max": float(np.max(steady)),
            "warmup_step_seconds": [float(v) for v in times[:2]],
            "predict_three_seeds_seconds_one_patient": predict_s,
            "torch_threads": int(torch.get_num_threads()), "torch_version": str(torch.__version__)}


def arm_fixed_growth():
    from src.fk_fixed_growth import FixedGrowthCase, fixed_growth_forecast, select_radius
    patients, domain = synthetic_cohort()
    cases = [FixedGrowthCase(p["id"], p["masks"][1], domain, (SPACING_MM,) * 3, p["masks"][2])
             for p in patients[:max(FOLD_TRAIN_SIZES)]]
    t = time.perf_counter()
    sel = select_radius(cases)
    select_s = time.perf_counter() - t
    t = time.perf_counter()
    fixed_growth_forecast(patients[-1]["masks"][1], domain, spacing_mm=(SPACING_MM,) * 3, radius_mm=sel["radius_mm"])
    return {"select_seconds_15_training": select_s, "forecast_seconds_one_patient": time.perf_counter() - t,
            "chosen_radius_mm_synthetic": sel["radius_mm"]}


def arm_count():
    from src import fk_learned_baselines as lb
    patients, domain = synthetic_cohort()
    t = time.perf_counter()
    counts = [lb.new_label_pieces(p["masks"][1], p["masks"][2], spacing_mm=(SPACING_MM,) * 3, domain=domain)
              for p in patients]
    pieces_s = time.perf_counter() - t
    lv = [np.log(p["masks"][1].sum() * SPACING_MM ** 3) for p in patients]
    gaps = [p["days"][2] - p["days"][1] for p in patients]
    t = time.perf_counter()
    model = lb.fit_count_model(counts, lv, gaps)
    return {"new_piece_seconds_18_patients": pieces_s, "fit_seconds": time.perf_counter() - t,
            "fit_status_synthetic": model["status"]}


def arm_pde():
    from src.fk_hybrid_forecast import hybrid_posterior_predictive_forecast
    patients, domain = synthetic_cohort()
    p = patients[0]
    settings = {"rho_bounds": (1e-5, 1e-1), "D_bounds": (1e-5, 10.0), "n_coarse": 21, "mass_fraction": 0.9,
                "n_refine": 11, "n_draws": 200, "n_anchors": 256}
    out = {}
    for calibration_only in (True, False):
        t = time.perf_counter()
        try:
            hybrid_posterior_predictive_forecast(
                masks=p["masks"][:2], days=p["days"][:2], target_day=p["days"][1] + 90.0,
                prior={"distribution": "log_uniform"}, h=SPACING_MM,
                boundary_covariance={"sigma_b_mm": 0.8, "tau_b_mm": 0.3, "ell_mm": 6.0},
                population_log_volume_slope_per_day=0.002, shrinkage_time_constant_days=90.0,
                volume_variance_fit={"intercept": 0.02, "slope_per_day": 0.0001}, domain=domain,
                seed=SEED, calibration_only=calibration_only, **settings)
            status = "ok"
        except ValueError as exc:
            status = f"model failure: {exc}"
        key = "calibration_only" if calibration_only else "full_forecast"
        out[f"{key}_seconds_one_patient"] = time.perf_counter() - t
        out[f"{key}_status"] = status
    out["settings"] = {k: list(v) if isinstance(v, tuple) else v for k, v in settings.items()}
    return out


def arm_unet_full_rerun():
    return arm_unet(None, timed_steps=40)


def arm_unet_crop64():
    return arm_unet(64, timed_steps=40)


CROP_ARMS = {"unet_full_rerun": arm_unet_full_rerun, "unet_crop64": arm_unet_crop64}
ARMS = {"unet": arm_unet, "fixed_growth": arm_fixed_growth, "count": arm_count, "pde_hybrid": arm_pde}


def project(arms):
    n_pairs_total = sum(FOLD_TRAIN_SIZES)
    from src.fk_learned_baselines import UNET_RECIPE as r
    unet_steps = n_pairs_total * r.epochs * r.n_seeds
    unet_h = unet_steps * arms["unet"]["step_seconds_mean"] / 3600
    unet_h_hi = unet_steps * arms["unet"]["step_seconds_max"] / 3600
    pde = arms["pde_hybrid"]
    # Assumes one full forecast per patient and one calibration-only forecast per outer-training patient per fold.
    pde_h = (N_PATIENTS * pde["full_forecast_seconds_one_patient"]
             + n_pairs_total * pde["calibration_only_seconds_one_patient"]) / 3600
    one_run_h = unet_h + pde_h
    return {
        "assumptions": [
            "U-Net cost = (sum of outer-training sizes) x 300 epochs x 3 seeds x measured mean step time",
            "step time is constant across epochs, seeds and folds; no other load on the machine",
            "PDE cost = 18 full forecasts + one calibration-only forecast per outer-training patient per fold; "
            "inner cross-fitting may call more forecasts than this, so the PDE figure is a lower bound",
            "fixed-growth and count arms are negligible next to the U-Net",
            "the protocol-complete wrapper runs four nested analyses (three prior supports + central-gap refit); "
            "each retrains the U-Net",
        ],
        "unet_training_steps_one_analysis": unet_steps,
        "unet_hours_one_analysis_mean_step": unet_h,
        "unet_hours_one_analysis_max_step": unet_h_hi,
        "pde_hours_one_analysis_lower_bound": pde_h,
        "total_hours_one_analysis": one_run_h,
        "total_hours_protocol_complete_four_analyses": 4 * one_run_h,
        "one_outer_fold_unet_hours": FOLD_TRAIN_SIZES[-1] * r.epochs * r.n_seeds
                                     * arms["unet"]["step_seconds_mean"] / 3600,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", choices=sorted({**ARMS, **CROP_ARMS}))
    parser.add_argument("--crop-pilot", action="store_true",
                        help="time the full grid and the proposed (not adopted) 64^3 crop; write a separate record")
    args = parser.parse_args()
    if args.arm:
        t = time.perf_counter()
        result = {**ARMS, **CROP_ARMS}[args.arm]()
        result["arm_wall_seconds"] = time.perf_counter() - t
        result["peak_memory_gb"] = _peak_gb()
        print(json.dumps(result))
        return
    if args.crop_pilot:
        crop_main()
        return
    arms = {}
    for name in ARMS:
        proc = subprocess.run([sys.executable, "-m", "src.fk_track_b_timing_pilot", "--arm", name],
                              cwd=PROJECT_ROOT, capture_output=True, text=True, check=True)
        arms[name] = json.loads(proc.stdout.strip().splitlines()[-1])
    record = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "evidence_scope": "synthetic-only runtime pilot; no patient data read; not patient-result evidence",
        "grid_shape": list(GRID_SHAPE), "spacing_mm": SPACING_MM, "n_synthetic_patients": N_PATIENTS,
        "fold_train_sizes": list(FOLD_TRAIN_SIZES),
        "machine": {"platform": platform.platform(), "processor": platform.processor(),
                    "cpu_count": os.cpu_count(), "python": platform.python_version()},
        "arms": arms,
        "projection": project(arms),
    }
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_DIR / "timing_record.json"
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    path.write_text(json.dumps(record, indent=2))
    print(json.dumps(record["projection"], indent=2))


def _run_arm(name):
    proc = subprocess.run([sys.executable, "-m", "src.fk_track_b_timing_pilot", "--arm", name],
                          cwd=PROJECT_ROOT, capture_output=True, text=True, check=True)
    return json.loads(proc.stdout.strip().splitlines()[-1])


def crop_main():
    """Second pilot: full grid rerun and 64^3 crop, same session, same machine. Proposal only, not adopted."""
    from src.fk_learned_baselines import UNET_RECIPE as r
    arms = {name: _run_arm(name) for name in CROP_ARMS}
    first = json.loads((OUTPUT_DIR / "timing_record.json").read_text())
    pde_h = first["projection"]["pde_hours_one_analysis_lower_bound"]
    steps = sum(FOLD_TRAIN_SIZES) * r.epochs * r.n_seeds
    gap_steps = steps  # central-gap subset is smaller; using the full count keeps this an upper figure
    out = {}
    for name, arm in arms.items():
        unet_h = steps * arm["step_seconds_mean"] / 3600
        out[name] = {
            "unet_hours_one_analysis_mean_step": unet_h,
            "unet_hours_one_analysis_max_step": steps * arm["step_seconds_max"] / 3600,
            "protocol_complete_hours_no_reuse": 4 * unet_h + 4 * pde_h,
            "protocol_complete_hours_with_reuse": unet_h + gap_steps * arm["step_seconds_mean"] / 3600 + 4 * pde_h,
        }
    record = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "evidence_scope": "synthetic-only runtime pilot; no patient data read; not patient-result evidence",
        "status": "crop and U-Net reuse are PROPOSED, NOT ADOPTED; frozen recipe unchanged",
        "grid_shape": list(GRID_SHAPE), "spacing_mm": SPACING_MM, "n_synthetic_patients": N_PATIENTS,
        "fold_train_sizes": list(FOLD_TRAIN_SIZES),
        "machine": {"platform": platform.platform(), "processor": platform.processor(),
                    "cpu_count": os.cpu_count(), "python": platform.python_version()},
        "arms": arms,
        "pde_hours_one_analysis_lower_bound_from_first_pilot": pde_h,
        "projection": {
            "assumptions": [
                "U-Net steps per analysis = 72 pairs x 300 epochs x 3 seeds, constant step time",
                "PDE hours are the first pilot's lower bound, times four analyses; PDE is not reused",
                "with reuse: one U-Net training for the three prior-support runs plus one for the central-gap "
                "refit, which is counted at full size (upper figure)",
                "prediction time is small next to training and is not added",
            ],
            **out,
        },
    }
    path = OUTPUT_DIR / "crop_timing_record.json"
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    path.write_text(json.dumps(record, indent=2))
    print(json.dumps(record["projection"], indent=2))


if __name__ == "__main__":
    main()
