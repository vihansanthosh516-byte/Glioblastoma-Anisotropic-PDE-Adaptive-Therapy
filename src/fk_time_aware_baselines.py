"""PROPOSED, NOT ADOPTED time-aware Track B baselines (LUMIERE memo r5). Synthetic prototype only.

Not wired into any runner. The frozen required baselines and the scoring protocol are unchanged.
- FiLMUNet3D: the frozen U-Net backbone plus FiLM after every InstanceNorm layer. A small network maps
  z = (log(forecast weeks / 13), log(history weeks / 13)) to a separate scale and shift for each norm channel.
  The last FiLM layer starts at zero, so the untrained model equals the frozen U-Net with the same backbone.
- Time-aware fixed growth: radius = k x nominal days (7 x week-label gap), k chosen per outer fold from
  outer-training patients only, ties to the smaller k; the grid keeps every radius under the frozen 20 mm cap.

Usage:  python -m src.fk_time_aware_baselines --time-film   (synthetic step timing at the LUMIERE grid)
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path as _Path

import numpy as np
import torch

from src.fk_fixed_growth import TIE_TOLERANCE, fixed_growth_forecast
from src.fk_learned_baselines import UNET_RECIPE, UNet3D
from src.fk_mu_scoring import dice_with_empty_rule

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output" / "track_b_timing_pilot_20261031"

INTERVAL_SCALE_WEEKS = 13.0
FILM_HIDDEN = 16
RADIUS_CAP_MM = 20.0
DAYS_PER_WEEK = 7.0
RATE_GRID_MM_PER_DAY = tuple(round(0.01 * i, 2) for i in range(20))  # 0.00 ... 0.19
MAX_GAP_WEEKS = 15


def interval_features(forecast_weeks: float, history_weeks: float) -> torch.Tensor:
    """z = (log(forecast/13), log(history/13)) from published week-label differences; fixed scaling."""
    f, h = float(forecast_weeks), float(history_weeks)
    if not (f > 0 and h > 0):
        raise ValueError("week-label intervals must be positive")
    return torch.tensor([np.log(f / INTERVAL_SCALE_WEEKS), np.log(h / INTERVAL_SCALE_WEEKS)], dtype=torch.float32)


class FiLMUNet3D(torch.nn.Module):
    def __init__(self, widths=UNET_RECIPE.widths, slope=UNET_RECIPE.leaky_slope, in_channels=3):
        super().__init__()
        self.backbone = UNet3D(widths, slope, in_channels)
        blocks = list(self.backbone.down) + list(self.backbone.up_blocks)
        self.norm_channels = [m.num_features for b in blocks for m in b if isinstance(m, torch.nn.InstanceNorm3d)]
        self.film = torch.nn.Sequential(torch.nn.Linear(2, FILM_HIDDEN), torch.nn.LeakyReLU(slope),
                                        torch.nn.Linear(FILM_HIDDEN, 2 * sum(self.norm_channels)))
        torch.nn.init.zeros_(self.film[-1].weight)
        torch.nn.init.zeros_(self.film[-1].bias)

    def _modulations(self, z):
        out = self.film(z.reshape(-1, 2))
        gammas, betas = out.split(sum(self.norm_channels), dim=1)
        return gammas.split(self.norm_channels, dim=1), betas.split(self.norm_channels, dim=1)

    def _block(self, block, x, gammas, betas, k):
        for m in block:
            x = m(x)
            if isinstance(m, torch.nn.InstanceNorm3d):
                g, b = gammas[k][:, :, None, None, None], betas[k][:, :, None, None, None]
                x = x * (1 + g) + b
                k += 1
        return x, k

    def forward(self, x, z):
        gammas, betas = self._modulations(z)
        net, k, skips = self.backbone, 0, []
        for i, block in enumerate(net.down):
            x, k = self._block(block, x, gammas, betas, k)
            if i < len(net.down) - 1:
                skips.append(x)
                x = net.pool(x)
        for up, block, skip in zip(net.up, net.up_blocks, reversed(skips)):
            x, k = self._block(block, torch.cat([up(x), skip], dim=1), gammas, betas, k)
        return net.head(x)


# ---------------------------------------------------------------- time-aware fixed growth

@dataclass(frozen=True)
class TimedGrowthCase:
    patient_id: str
    last_mask: np.ndarray
    support: np.ndarray
    spacing_mm: tuple[float, float, float]
    gap_weeks: int
    target_mask: np.ndarray | None = None


def _check_grid(grid, max_gap_weeks):
    grid = sorted(float(k) for k in grid)
    if not grid or grid[0] != 0.0 or len(set(grid)) != len(grid):
        raise ValueError("rate grid must be unique, nonnegative and include 0")
    if grid[-1] * DAYS_PER_WEEK * max_gap_weeks > RADIUS_CAP_MM + 1e-9:
        raise ValueError(f"rate grid exceeds the {RADIUS_CAP_MM} mm cap at {max_gap_weeks} weeks")
    return grid


def timed_radius_mm(rate_mm_per_day: float, gap_weeks: int) -> float:
    if not 0 < gap_weeks <= MAX_GAP_WEEKS:
        raise ValueError(f"gap must be 1-{MAX_GAP_WEEKS} week labels")
    return float(rate_mm_per_day) * DAYS_PER_WEEK * int(gap_weeks)


def select_rate(training_cases: Sequence[TimedGrowthCase], *, rate_grid=RATE_GRID_MM_PER_DAY) -> dict:
    """Pick k with the best patient-macro Dice on training cases; ties within TIE_TOLERANCE go to the smaller k."""
    grid = _check_grid(rate_grid, MAX_GAP_WEEKS)
    if not training_cases:
        raise ValueError("rate selection needs at least one training case")
    ids = [c.patient_id for c in training_cases]
    if len(set(ids)) != len(ids):
        raise ValueError("each training patient must appear once")
    mean_dice = []
    for k in grid:
        scores = []
        for c in training_cases:
            if c.target_mask is None:
                raise ValueError(f"training case {c.patient_id} has no target mask")
            f = fixed_growth_forecast(c.last_mask, c.support, spacing_mm=c.spacing_mm,
                                      radius_mm=timed_radius_mm(k, c.gap_weeks))
            scores.append(dice_with_empty_rule(f, np.asarray(c.target_mask, bool)))
        mean_dice.append(float(np.mean(scores)))
    best = max(mean_dice)
    chosen = next(k for k, d in zip(grid, mean_dice) if d >= best - TIE_TOLERANCE)
    return {"rate_mm_per_day": chosen, "grid": grid, "training_mean_dice": mean_dice,
            "n_training": len(training_cases), "boundary_hit": chosen == grid[-1], "equals_copy_last": chosen == 0.0}


def run_outer_timed_growth(cases: Sequence[TimedGrowthCase], fold_by_patient: Mapping[str, int]) -> dict:
    """Per outer fold: select k on that fold's training patients only, forecast its test patients."""
    folds = sorted({fold_by_patient[c.patient_id] for c in cases})
    out = {}
    for f in folds:
        train = [c for c in cases if fold_by_patient[c.patient_id] != f]
        test = [c for c in cases if fold_by_patient[c.patient_id] == f]
        sel = select_rate(train)
        out[f] = {"selection": sel, "training_ids": [c.patient_id for c in train],
                  "forecasts": {c.patient_id: fixed_growth_forecast(
                      c.last_mask, c.support, spacing_mm=c.spacing_mm,
                      radius_mm=timed_radius_mm(sel["rate_mm_per_day"], c.gap_weeks)) for c in test}}
    return out


# ---------------------------------------------------------------- synthetic step timing

def _time_steps(model, xs, ys, zs, timed_steps, r):
    from src import fk_learned_baselines as lb
    gen = np.random.default_rng(r.base_seed)
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
        logits = model(x) if zs is None else model(x, zs[step % 2])
        loss = lb._loss(logits, y, r.dice_smoothing)
        loss.backward()
        opt.step()
        times.append(time.perf_counter() - t)
    steady = times[2:]
    return {"n_parameters": int(sum(p.numel() for p in model.parameters())), "timed_steps": len(steady),
            "step_seconds_mean": float(np.mean(steady)), "step_seconds_min": float(np.min(steady)),
            "step_seconds_max": float(np.max(steady)), "warmup_step_seconds": [float(v) for v in times[:2]]}


def time_film_main(timed_steps: int = 40):
    """Same method as the earlier LUMIERE pilot (2 warm-up + timed steps, same augmentation and seeds)."""
    os.environ["TRACK_B_PILOT_COHORT"] = "lumiere"
    from src import fk_learned_baselines as lb
    from src import fk_track_b_timing_pilot as tp
    path = OUTPUT_DIR / "film_unet_timing_record.json"
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    r = lb.UNET_RECIPE
    patients, domain = tp.synthetic_cohort()
    pairs = [(lb.unet_inputs(p["masks"][0], p["masks"][1], domain), p["masks"][2]) for p in patients[:2]]
    padded = lb._padded_shape(pairs[0][1].shape, len(r.widths))
    xs = [torch.from_numpy(np.stack([lb._pad(c, padded) for c in x]))[None] for x, _ in pairs]
    ys = [torch.from_numpy(lb._pad(y.astype(np.float32), padded))[None, None] for _, y in pairs]
    zs = [interval_features(13, 13)[None], interval_features(15, 155)[None]]
    prior = torch.are_deterministic_algorithms_enabled()
    torch.use_deterministic_algorithms(True)
    try:
        torch.manual_seed(r.base_seed)
        frozen = _time_steps(lb.UNet3D(r.widths, r.leaky_slope), xs, ys, None, timed_steps, r)
        torch.manual_seed(r.base_seed)
        film = _time_steps(FiLMUNet3D(r.widths, r.leaky_slope), xs, ys, zs, timed_steps, r)
    finally:
        torch.use_deterministic_algorithms(prior)
    sizes = tp.FOLD_TRAIN_SIZES
    steps = sum(sizes) * r.epochs * r.n_seeds
    record = {
        "schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
        "evidence_scope": "synthetic-only step timing of the PROPOSED, NOT ADOPTED FiLM U-Net at the LUMIERE "
                          "2 mm grid; no patient data read; not patient-result evidence",
        "grid_shape": list(tp.GRID_SHAPE), "padded_shape": list(padded), "spacing_mm": tp.SPACING_MM,
        "fold_train_sizes": list(sizes),
        "machine": {"platform": platform.platform(), "processor": platform.processor(),
                    "cpu_count": os.cpu_count(), "python": platform.python_version()},
        "environment": lb.unet_runtime_environment(),
        "frozen_unet_same_session": frozen, "film_unet": film,
        "film_over_frozen_step_ratio": film["step_seconds_mean"] / frozen["step_seconds_mean"],
        "projection_inputs": {"unet_steps_per_training_set": steps,
                              "film_hours_per_training_set": steps * film["step_seconds_mean"] / 3600,
                              "frozen_hours_per_training_set_same_session": steps * frozen["step_seconds_mean"] / 3600},
    }
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2))
    print(json.dumps({k: record[k] for k in ("film_over_frozen_step_ratio", "projection_inputs")}, indent=2))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--time-film", action="store_true")
    if ap.parse_args().time_film:
        time_film_main()
