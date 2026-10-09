#!/usr/bin/env python3
"""Script 110: 3D tumour-forecast viewer, research prototype of a clinical tool (GRAND_PLAN section 11; Amendment 7 A7.5).

Not for patient care. Descriptive only: no endpoint is computed here.
For one MU forecast pair it builds a single self-contained HTML page:
  brain outline (see-through) | tumour core today (input scan) | model forecast (visible-tumour surface u >= U_VISIBLE)
  at several days across the interval (slider) | true core at the next scan | panel with dates, volumes, Dice.
Model: the fold's selected cell for the registered primary arm (aniso r=1, `output/pde_manifest_monotone/selected_cells.json`,
chosen on other folds), run with the exact positivity-preserving solver (ledger 2ae, 2ak), same treatment schedule as
script 100. The final-day Dice is checked against the script 100 monotone cache for that pair (must match to 1e-6 for
later pairs; first pairs use the grid cache with the same cell and may differ by < 0.002 because the grid run
shares one time step across a batch of rho values).
Default patient: the pair with the MEDIAN thresholded PDE Dice among primary pairs (no cherry-picking); any pair can be
chosen with --patient/--tp-in/--tp-out.
Surfaces: marching cubes (scikit-image) on the 2 mm masks, light smoothing of the brain only. Page: Plotly mesh3d,
plotly.js embedded so the file opens offline.
Output: output/viewer3d/<patient>_<tp_in>_<tp_out>.html and .json (numbers shown on the page)
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path as _Path

import numpy as np

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "viewer3d"
PM = PROJECT_ROOT / "output" / "pde_manifest_monotone"
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import run_improved_aniso as ria  # noqa: E402
from solver_monotone import TensorFKMonotone  # noqa: E402

spec = importlib.util.spec_from_file_location("m100", PROJECT_ROOT / "src" / "100_pde_manifest.py")
m100 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m100)

ARM = "aniso|r=1"
N_FRAMES = 6
H = 2.0


def median_pair():
    import pandas as pd
    df = pd.read_csv(PM / "oof_pairs.csv")
    df = df.dropna(subset=[ARM]).sort_values(ARM).reset_index(drop=True)
    r = df.iloc[len(df) // 2]
    return r["patient_id"], int(r["tp_in"]), int(r["tp_out"])


def cached_dice(pair, cell):
    key = f"{pair['patient_id']}_{pair['tp_in']}_{pair['tp_out']}.json"
    if pair["has_prior_scan"]:
        c = json.loads((PM / "cache_sel" / key).read_text())
        return c["dice_sel"][ARM]
    c = json.loads((PM / "cache_h2" / key).read_text())
    return c["dice"]["grid"][cell]


def forecast_frames(pair, cell):
    S, T, brain = m100.load_inputs(pair["patient_id"], pair["tp_in"], pair["tp_out"], H)
    A = ria._ATLAS
    kv = dict(x.split("=") for x in cell.split("|")[1:])
    arm = cell.split("|")[0]
    margin = ria.MARGIN_VOX
    lo = np.maximum(np.argwhere(S).min(0) - margin, 0)
    hi = np.minimum(np.argwhere(S).max(0) + margin + 1, S.shape)
    box = tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))
    dom = (brain[box] & A["tissue"][box]) | S[box]
    T0 = ria.arm_tensor(A, arm, float(kv["r"]), box)   # as script 100 tensor_for()
    T0[S[box] & ~A["tissue"][box]] = np.array([1, 1, 1, 0, 0, 0], np.float32)
    solver = TensorFKMonotone(T0 * float(kv["d"]), dom, h=H)
    days = np.linspace(0, pair["dt_days"], N_FRAMES + 1)[1:]
    frames = []
    for t in days:   # each frame integrated from the input scan, so the last frame is exactly the script 100 run
        uu, _ = solver.run(S[box].astype(np.float32), [float(kv["rho"])], float(t),
                           kill_schedule=ria.build_schedule(pair.get("treatment")), t_start=pair["day_in"])
        u = np.zeros(S.shape, np.float32)
        u[box] = uu[0]
        frames.append((float(t), u))
    return S, T, brain, frames


def mesh(mask, smooth=0.0, step=1):
    from scipy import ndimage
    from skimage import measure
    v = mask.astype(np.float32)
    if smooth:
        v = ndimage.gaussian_filter(v, smooth)
    if v.max() < 0.5:
        return None
    verts, faces, _, _ = measure.marching_cubes(np.pad(v, 1), 0.5, step_size=step)
    verts = (verts - 1) * H
    return verts, faces


def mesh3d(m, name, color, opacity, visible=True):
    import plotly.graph_objects as go
    if m is None:
        return go.Mesh3d(x=[], y=[], z=[], name=name, visible=visible)
    v, f = m
    return go.Mesh3d(x=v[:, 0], y=v[:, 1], z=v[:, 2], i=f[:, 0], j=f[:, 1], k=f[:, 2], name=name, color=color,
                     opacity=opacity, showlegend=True, visible=visible, flatshading=False,
                     lighting=dict(ambient=0.5, diffuse=0.8, specular=0.2), hoverinfo="name")


def build_page(pair, cell, S, T, brain, frames, info, path):
    import plotly.graph_objects as go
    brain_m = mesh(brain, smooth=1.2, step=2)
    traces = [mesh3d(brain_m, "Brain", "#9aa3ad", 0.08), mesh3d(mesh(S), "Tumour today (input scan)", "#3b82f6", 0.55),
              mesh3d(mesh(T), f"Real tumour at next scan (day {pair['dt_days']:.0f})", "#22c55e", 0.35)]
    fc = [mesh3d(mesh(u >= ria.U_VISIBLE), f"Forecast day {t:.0f}", "#ef4444", 0.45, visible=(k == len(frames) - 1))
          for k, (t, u) in enumerate(frames)]
    fig = go.Figure(traces + fc)
    n0 = len(traces)
    steps = []
    for k, (t, _) in enumerate(frames):
        vis = [True] * n0 + [j == k for j in range(len(frames))]
        steps.append(dict(method="restyle", args=[{"visible": vis}], label=f"day {t:.0f}"))
    fig.update_layout(
        title=dict(text=(f"GBM forecast viewer: {pair['patient_id']} scan {pair['tp_in']} -> {pair['tp_out']} "
                         f"| RESEARCH PROTOTYPE - NOT FOR PATIENT CARE"), font=dict(size=15)),
        sliders=[dict(active=len(frames) - 1, steps=steps, currentvalue=dict(prefix="Forecast: "), pad=dict(t=30))],
        scene=dict(aspectmode="data", xaxis_visible=False, yaxis_visible=False, zaxis_visible=False,
                   bgcolor="#0b1220"),
        paper_bgcolor="#0b1220", font=dict(color="#e5e7eb"), legend=dict(itemclick="toggle", x=0.01, y=0.95),
        margin=dict(l=0, r=0, t=60, b=0), height=780)
    panel = (f"<div style='font:14px system-ui;color:#e5e7eb;background:#0b1220;padding:12px 16px;line-height:1.5'>"
             f"<b>Patient</b> {pair['patient_id']} &nbsp; <b>Interval</b> {pair['dt_days']:.0f} days &nbsp; "
             f"<b>Model</b> exact anisotropic Fisher-Kolmogorov, atlas tensor, population cell {cell} (chosen on other folds)<br>"
             f"<b>Core volume</b> today {info['v_in_ml']:.1f} mL &nbsp; forecast {info['v_forecast_ml']:.1f} mL &nbsp; "
             f"real {info['v_true_ml']:.1f} mL<br>"
             f"<b>Dice</b> forecast vs real {info['dice_forecast']:.3f} &nbsp; 'no change' vs real {info['dice_persistence']:.3f}"
             f" &nbsp; (cohort result: no PDE model beats 'no change' with confidence; ledger 2ak)<br>"
             f"<span style='color:#fca5a5'>Retrospective research prototype. Not validated for clinical use. "
             f"Masks only; 2 mm grid.</span></div>")
    html = fig.to_html(include_plotlyjs=True, full_html=True)
    html = html.replace("<body>", "<body style='margin:0;background:#0b1220'>" + panel, 1)
    path.write_text(html, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient")
    ap.add_argument("--tp-in", type=int)
    ap.add_argument("--tp-out", type=int)
    a = ap.parse_args()
    pid, tin, tout = (a.patient, a.tp_in, a.tp_out) if a.patient else median_pair()
    ria._forecast_init()
    pairs = m100.get_pairs(2.0, stage="forecast") + m100.get_pairs(2.0, stage="selected")
    pair = next(p for p in pairs if p["patient_id"] == pid and int(p["tp_in"]) == tin and int(p["tp_out"]) == tout)
    sel = json.loads((PM / "selected_cells.json").read_text())
    cell = sel[str(int(pair["fold"]))][ARM]
    ria.TensorFK = TensorFKMonotone
    S, T, brain, frames = forecast_frames(pair, cell)
    u_last = frames[-1][1]
    vv = H ** 3 / 1000.0
    info = {"patient_id": pid, "tp_in": tin, "tp_out": tout, "dt_days": float(pair["dt_days"]), "cell": cell,
            "solver": "TensorFKMonotone (selling)", "v_in_ml": float(S.sum() * vv), "v_true_ml": float(T.sum() * vv),
            "v_forecast_ml": float((u_last >= ria.U_VISIBLE).sum() * vv),
            "dice_forecast": float(ria.dice(u_last > ria.U_VISIBLE, T)), "dice_persistence": float(ria.dice(S, T)),
            "frame_days": [t for t, _ in frames]}
    info["dice_cache_script100"] = float(cached_dice(pair, cell))
    info["dice_matches_cache"] = bool(abs(info["dice_forecast"] - info["dice_cache_script100"]) < 1e-6)
    OUT.mkdir(parents=True, exist_ok=True)
    stem = f"{pid}_{tin}_{tout}"
    build_page(pair, cell, S, T, brain, frames, info, OUT / f"{stem}.html")
    (OUT / f"{stem}.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
