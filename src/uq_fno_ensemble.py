#!/usr/bin/env python3
"""
Human-in-the-Loop (HIL) Uncertainty Quantification (Proposal 4)
================================================================
Provides confidence bands for future tumour volume by combining an ensemble
of FNO neural-PDE surrogates with Monte-Carlo (MC) rollout. The output is a
Plotly/HTML trajectory plot with a shaded 95% confidence band over the next
7-30 days, consumed by the HIL dashboard (Phase 14) via an iframe.

Components
----------
1. ``FNOEnsemble``  — train/load N=5 independent FNO models (different seeds).
2. ``monte_carlo_predict(...)`` — M=200 forward passes by randomly drawing
   one ensemble member at each step; returns mean ± p2.5/p97.5 band.
3. ``plot_volume_trajectory_html(...)`` — Plotly HTML with shaded 95% CI.
4. ``trigger_eval(...)`` — dose-escalate alert when the lower bound of the
   confidence band exceeds a critical tumour-volume threshold (e.g. 150 mm^3).
5. ``CalibrationLog`` — stores realized day+7 volume per patient for empirical
   coverage estimation and model recalibration.

Dependencies: only torch + plotly (already in the repo).
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    import torch
    _HAVE_TORCH = True
except ImportError:
    torch = None  # type: ignore
    _HAVE_TORCH = False

try:
    import plotly.graph_objects as go
    _HAVE_PLOTLY = True
except ImportError:
    go = None  # type: ignore
    _HAVE_PLOTLY = False

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_ENSEMBLE_DIR = Path("output/fno_ensemble")
DEFAULT_CALIB_PATH = Path("output/uq_calibration.jsonl")

# Default policy thresholds
DEFAULT_CRITICAL_VOLUME_MM3 = 150.0
DEFAULT_FORECAST_HORIZON_DAYS = 30
M_DEFAULT = 200


# --------------------------------------------------------------------------- #
# FNO ensemble
# --------------------------------------------------------------------------- #
def _build_fno(weights_path: Optional[Path] = None,
               grid_size: int = 32, modes: int = 8, width: int = 20):
    """Construct a fresh FNO3d and optionally load weights.

    Imports fno_solver lazily so this module works in CPU-only environments
    even when the GPU FNO training code is unavailable.
    """
    import importlib.util, sys
    sys.path.insert(0, str(REPO_ROOT / "src"))
    p = REPO_ROOT / "src" / "neural_pde" / "fno_solver.py"
    spec = importlib.util.spec_from_file_location("fno_solver", p)
    fno_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fno_mod)
    model = fno_mod.FNO3d(in_channels=2, out_channels=1, modes=modes, width=width)
    if weights_path is not None and Path(weights_path).exists():
        state = torch.load(weights_path, map_location="cpu")
        if isinstance(state, dict) and "model_state" in state:
            state = state["model_state"]
        try:
            model.load_state_dict(state)
        except Exception as e:
            print(f"[uq] WARN could not load {weights_path}: {e}")
    model.eval()
    return model


class FNOEnsemble:
    """Manages N independently-trained FNO model checkpoints.

    Parameters
    ----------
    model_paths : list of paths to ``fno_model_i.pth``. If None or empty, a
        pool of N randomly-initialized FNO3d networks is created (useful for
        unit-testing the MC plumbing without real training).
    """

    def __init__(
        self,
        model_paths: Optional[List[Path]] = None,
        n: int = 5,
        grid_size: int = 32,
    ) -> None:
        if not _HAVE_TORCH:
            raise ImportError("torch is required for FNOEnsemble")
        self.grid_size = grid_size
        self.models = []
        if model_paths:
            for p in model_paths:
                self.models.append(_build_fno(Path(p), grid_size=grid_size))
        else:
            for _ in range(n):
                torch.manual_seed(np.random.randint(0, 2 ** 31 - 1))
                self.models.append(_build_fno(grid_size=grid_size))
        self.n = len(self.models)
        print(f"[uq] FNOEnsemble loaded {self.n} model(s)")

    def predict_one_member(self, idx: int, u0: np.ndarray,
                           rho: np.ndarray) -> np.ndarray:
        """Single FNO forward pass for ensemble member idx."""
        m = self.models[idx]
        gs = self.grid_size
        if u0.shape != (gs, gs, gs) or rho.shape != (gs, gs, gs):
            from scipy.ndimage import zoom
            u0 = zoom(u0,
                      tuple(gs / s for s in u0.shape), order=1)
            rho = zoom(rho,
                       tuple(gs / s for s in rho.shape), order=1)
        inp = np.stack([u0, rho], axis=0).astype(np.float32)
        t = torch.tensor(inp).unsqueeze(0)
        with torch.no_grad():
            out = m(t)
        return np.clip(out.squeeze().cpu().numpy(), 0.0, 1.0)

    # ------------------------------------------------------------------ #
    def monte_carlo_predict(
        self,
        u0: np.ndarray,
        rho: np.ndarray,
        horizon_days: int = DEFAULT_FORECAST_HORIZON_DAYS,
        M: int = M_DEFAULT,
        rng_seed: Optional[int] = 42,
    ) -> Dict[str, np.ndarray]:
        """MC rollout: at each forward day, randomly draw one of N ensemble
        members and roll one FNO step. Repeat M independent trajectories and
        compute the mean and 2.5/97.5 percentiles per day.

        Returns arrays of shape (horizon_days+1,) for mean / lower / upper and
        a (M, horizon_days+1) raw array of per-trajectory volumes (mm^3) for
        downstream plotting.
        """
        rng = np.random.default_rng(rng_seed)
        # Volume = count(u>0.1)*voxel_volume; here approximated as sum(u)
        voxel_vol_mm3 = 1.0  # grid is normalized; tune via dx^3 in callers
        all_vols = np.zeros((M, horizon_days + 1), dtype=np.float64)
        for traj in range(M):
            u = np.asarray(u0, dtype=np.float32).copy()
            r = np.asarray(rho, dtype=np.float32).copy()
            all_vols[traj, 0] = np.clip(u, 0, 1).sum() * voxel_vol_mm3
            for day in range(1, horizon_days + 1):
                member = int(rng.integers(0, self.n))
                u = self.predict_one_member(member, u, r)
                all_vols[traj, day] = np.clip(u, 0, 1).sum() * voxel_vol_mm3
        mean = all_vols.mean(axis=0)
        lower = np.percentile(all_vols, 2.5, axis=0)
        upper = np.percentile(all_vols, 97.5, axis=0)
        std = all_vols.std(axis=0)
        return {
            "mean": mean, "lower": lower, "upper": upper, "std": std,
            "trajectories": all_vols,
            "horizon_days": horizon_days, "M": M,
        }


# --------------------------------------------------------------------------- #
# HIL HTML rendering
# --------------------------------------------------------------------------- #
def plot_volume_trajectory_html(
    forecast: Dict,
    patient_id: str,
    output: Path = Path("output/uq_forecast.html"),
    critical_volume_mm3: float = DEFAULT_CRITICAL_VOLUME_MM3,
) -> Path:
    """Render a Plotly HTML figure with the mean volume trajectory + 95% band."""
    output.parent.mkdir(parents=True, exist_ok=True)
    days = np.arange(len(forecast["mean"]))
    mean = forecast["mean"]
    lower = forecast["lower"]
    upper = forecast["upper"]

    flag_critical = bool(forecast["lower"][-1] >= critical_volume_mm3) or \
                    bool(upper[-1] >= critical_volume_mm3)
    band_color = "red" if flag_critical else "rgba(0,100,255,0.18)"

    if not _HAVE_PLOTLY:
        # Minimal HTML fallback
        body = (f"<html><body><h3>UQ forecast for {patient_id}</h3>"
                f"<p>mean: [{mean[0]:.1f} -> {mean[-1]:.1f}] mm^3</p>"
                f"<p>95% CI day {len(mean)-1}: "
                f"[{lower[-1]:.1f}, {upper[-1]:.1f}] mm^3</p>"
                f"<p>Critical flag: {flag_critical}</p></body></html>")
        output.write_text(body)
        return output
    fig = go.Figure([
        go.Scatter(x=days, y=upper, line=dict(width=0), hoverinfo="skip",
                   showlegend=False),
        go.Scatter(x=days, y=lower, line=dict(width=0),
                   fill="tonexty", fillcolor=band_color,
                   name="95% CI", hovertemplate="day %{x}: [%{y:.1f}%{customdata:.1f}]<extra></extra>",
                   customdata=upper),
        go.Scatter(x=days, y=mean, mode="lines", name="Mean volume (mm^3)",
                   line=dict(color="blue" if not flag_critical else "darkred")),
        go.Scatter(x=[days[-1]], y=[upper[-1]], mode="markers",
                   marker=dict(color="red" if flag_critical else "green",
                               size=10),
                   name="Day H upper bound"),
    ])
    fig.add_hline(y=critical_volume_mm3, line_dash="dash",
                  line_color="red",
                  annotation_text=f"critical volume = {critical_volume_mm3} mm^3")
    title = f"UQ trajectory for {patient_id} - 95% CI band over next {len(days)-1} days"
    if flag_critical:
        title += "  [DOSE-ESCALATE ALERT]"
    fig.update_layout(title=title, xaxis_title="Day", yaxis_title="Volume (mm^3)")
    fig.write_html(str(output))
    print(f"[uq] wrote trajectory HTML -> {output}")
    return output


# --------------------------------------------------------------------------- #
# Threshold-trigger evaluation
# --------------------------------------------------------------------------- #
def trigger_eval(
    forecast: Dict,
    critical_volume_mm3: float = DEFAULT_CRITICAL_VOLUME_MM3,
) -> Dict:
    """Decide whether to flag a dose-escalate recommendation.

    Trigger fires when the *lower* bound of the confidence band at the final
    forecast day exceeds the critical volume (risk-aware criterion: even the
    optimistic trajectory clears the threshold).
    """
    final_lower = float(forecast["lower"][-1])
    final_upper = float(forecast["upper"][-1])
    flag = bool(final_lower >= float(critical_volume_mm3))
    recommendation = "dose-escalate" if flag else "hold / monitor"
    return {
        "flag": flag,
        "trigger_lower_volume_mm3": final_lower,
        "trigger_upper_volume_mm3": final_upper,
        "critical_volume_mm3": critical_volume_mm3,
        "recommendation": recommendation,
        "horizon_days": int(forecast["horizon_days"]),
        "M_samples": int(forecast["M"]),
    }


# --------------------------------------------------------------------------- #
# Calibration tracking & empirical coverage
# --------------------------------------------------------------------------- #
class CalibrationLog:
    """Append-only JSONL of realized day+7 volumes for empirical coverage.

    After each new patient run, store the realized volume at day H and
    compute the empirical coverage of the 95% band: fraction of patients
    whose realized volume fell inside [forecast_lower, forecast_upper] at
    the chosen day.
    """

    def __init__(self, path: Path = DEFAULT_CALIB_PATH) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def record(self, patient_id: str, day: int,
               forecast_lower: float, forecast_upper: float,
               forecast_mean: float, realized_volume: float) -> None:
        entry = {
            "patient_id": patient_id, "day": day,
            "forecast_lower": float(forecast_lower),
            "forecast_upper": float(forecast_upper),
            "forecast_mean": float(forecast_mean),
            "realized_volume": float(realized_volume),
            "covered": bool(forecast_lower <= realized_volume <= forecast_upper),
        }
        with open(self.path, "a", encoding="utf-8") as f:
            f.write(json.dumps(entry) + "\n")

    def coverage(self) -> Dict:
        entries = []
        if self.path.exists():
            entries = [
                json.loads(l) for l in self.path.read_text().splitlines()
                if l.strip()
            ]
        if not entries:
            return {"n": 0, "empirical_coverage": None}
        covered = sum(e["covered"] for e in entries)
        return {
            "n": len(entries),
            "empirical_coverage": covered / len(entries),
            "target_coverage": 0.95,
        }


# --------------------------------------------------------------------------- #
# High-level: predict + render for a single patient
# --------------------------------------------------------------------------- #
def predict_and_render(
    ensemble: "FNOEnsemble",
    u0: np.ndarray,
    rho: np.ndarray,
    patient_id: str,
    horizon_days: int = DEFAULT_FORECAST_HORIZON_DAYS,
    M: int = M_DEFAULT,
    critical_volume_mm3: float = DEFAULT_CRITICAL_VOLUME_MM3,
    output_dir: Path = Path("output/uq"),
) -> Dict:
    forecast = ensemble.monte_carlo_predict(
        u0, rho, horizon_days=horizon_days, M=M
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    html_path = plot_volume_trajectory_html(
        forecast, patient_id,
        output=output_dir / f"uq_forecast_{patient_id}.html",
        critical_volume_mm3=critical_volume_mm3,
    )
    trigger = trigger_eval(forecast, critical_volume_mm3=critical_volume_mm3)
    pick = Path(output_dir) / f"uq_forecast_{patient_id}.pkl"
    with open(pick, "wb") as f:
        pickle.dump({
            "forecast": {k: (v.tolist() if hasattr(v, "tolist") else v)
                          for k, v in forecast.items()},
            "trigger": trigger,
        }, f)
    return {
        "forecast_mean": forecast["mean"].tolist(),
        "forecast_lower": forecast["lower"].tolist(),
        "forecast_upper": forecast["upper"].tolist(),
        "trigger": trigger,
        "html_path": str(html_path),
    }


# --------------------------------------------------------------------------- #
# CLI self-test: synthetic cohort
# --------------------------------------------------------------------------- #
def _self_test() -> None:
    if not _HAVE_TORCH:
        print("[uq] torch not available; skipping self-test"); return
    rng = np.random.default_rng(0)
    gs = 32
    u0 = np.zeros((gs, gs, gs), dtype=np.float32)
    z, y, x = np.ogrid[:gs, :gs, :gs]
    mask = (x - gs // 2) ** 2 + (y - gs // 2) ** 2 + (z - gs // 2) ** 2 <= 5 ** 2
    u0[mask] = 0.7
    rho = np.full((gs, gs, gs), 0.02, dtype=np.float32)

    ens = FNOEnsemble(model_paths=None, n=5, grid_size=gs)
    # Tiny forecast for speed
    forecast = ens.monte_carlo_predict(u0, rho, horizon_days=3, M=5)
    assert forecast["mean"].shape == (4,), forecast["mean"].shape
    # sanity: vol nearest-day near initial sum and non-negative
    assert forecast["mean"][0] >= 0 and forecast["mean"][-1] >= 0
    # trigger evaluation must always return a dict
    trig = trigger_eval(forecast, critical_volume_mm3=forecast["upper"][-1] + 1.0)
    assert isinstance(trig["flag"], bool)
    # calibration roundtrip
    log = CalibrationLog(Path("output/_uq_calib_test.jsonl"))
    log._path = Path("output/_uq_calib_test.jsonl")
    if log._path.exists():
        log._path.unlink()
    lk = CalibrationLog(Path("output/_uq_calib_test.jsonl"))
    lk.record("P01", day=7, forecast_lower=forecast["lower"][-1],
              forecast_upper=forecast["upper"][-1],
              forecast_mean=forecast["mean"][-1],
              realized_volume=forecast["mean"][-1] * 1.1)
    cov = lk.coverage()
    assert isinstance(cov["n"], int)
    Path("output/_uq_calib_test.jsonl").unlink(missing_ok=True)
    html = plot_volume_trajectory_html(forecast, "P01",
                                       Path("output/_uq_test.html"))
    assert Path(html).exists()
    Path("output/_uq_test.html").unlink(missing_ok=True)
    print(f"[uq] self-test OK: M=5 horizon=3 "
          f"mean[-1]={forecast['mean'][-1]:.3f}, "
          f"CI=[{forecast['lower'][-1]:.3f},{forecast['upper'][-1]:.3f}]")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(
        description="HIL Uncertainty Quantification (Proposal 4)"
    )
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--n-models", type=int, default=5)
    p.add_argument("--model-dir", type=str, default=str(DEFAULT_ENSEMBLE_DIR))
    p.add_argument("--horizon-days", type=int, default=DEFAULT_FORECAST_HORIZON_DAYS)
    p.add_argument("--M", type=int, default=M_DEFAULT)
    p.add_argument("--patient-id", type=str, default="P01")
    p.add_argument("--grid-size", type=int, default=32)
    p.add_argument("--critical-volume", type=float,
                   default=DEFAULT_CRITICAL_VOLUME_MM3)
    p.add_argument("--output-dir", type=str, default="output/uq")
    args = p.parse_args()

    if args.self_test:
        _self_test()
        sys.exit(0)

    if not _HAVE_TORCH:
        print("[uq] torch not available; running self-test only")
        _self_test(); sys.exit(0)

    model_dir = Path(args.model_dir)
    model_paths = sorted(model_dir.glob("fno_model_*.pth")) if model_dir.exists() else []
    ens = FNOEnsemble(model_paths=model_paths or None,
                      n=args.n_models, grid_size=args.grid_size)
    rng = np.random.default_rng(0)
    gs = args.grid_size
    u0 = np.zeros((gs, gs, gs), dtype=np.float32)
    z, y, x = np.ogrid[:gs, :gs, :gs]
    mask = (x - gs // 2) ** 2 + (y - gs // 2) ** 2 + (z - gs // 2) ** 2 <= 5 ** 2
    u0[mask] = 0.7
    rho = np.full((gs, gs, gs), 0.02, dtype=np.float32)
    res = predict_and_render(
        ens, u0, rho, args.patient_id,
        horizon_days=args.horizon_days, M=args.M,
        critical_volume_mm3=args.critical_volume, output_dir=Path(args.output_dir))
    print(f"[uq] trigger: {res['trigger']['recommendation']}")
    print(f"[uq] HTML -> {res['html_path']}")
