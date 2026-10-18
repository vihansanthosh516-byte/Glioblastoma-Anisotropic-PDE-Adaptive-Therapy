#!/usr/bin/env python3
"""Script 137: 4D viewer (3D space + time), research prototype for the science fair. Not for patient care.

Two modes in one self-contained HTML page (plotly.js embedded, opens offline):
  1. Recurrence-risk map (no time axis). Per patient: pre-op tumour core, the standard 15 mm plan, the model plan at the SAME
     volume (frozen H-3 model C1, which was never trained on UCSF), and the real recurrence (ALPTDG time1, labels 1-3).
     Coverage numbers are recomputed here and checked against output/predictgbm/h3_final/ucsf_test_pairs.csv.
  2. Growth simulation (time slider). Fisher-Kolmogorov growth from the pre-op core on the 2 mm grid with ASSUMED round
     parameters (rho, D_wm, D_gm are not fitted to the patient). Scenarios: no treatment, and an ASSUMED constant drug
     exposure whose kill rate is a fraction of the growth rate. This is a what-if simulation, not a forecast.
Patients: UCSF locked-test patients at the 75th, 90th and 95th percentile of all-recurrence gain (model minus standard):
  upper-part examples, labelled as such on the page together with the cohort mean gain.
Output: output/viewer4d/viewer4d.html and viewer4d_info.json.
"""
from __future__ import annotations

import gzip
import json
import os
import sys
import zipfile
from pathlib import Path as _Path

import nibabel as nib
import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "viewer4d"
H3OUT = PROJECT_ROOT / "output" / "predictgbm" / "h3_final"
H3 = _Path(os.environ.get("UCSF_H3_DIR", r"C:\Users\vihan\data\ucsf_h3"))
ZIP = _Path(os.environ.get("UCSF_ALPTDG_ZIP", r"C:\Users\vihan\Downloads\UCSF_POSTOP_GLIOMA_DATASET_FINAL_v1.0.zip"))
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import predictgbm_io as pg  # noqa: E402

QUANTILES = (0.75, 0.9, 0.95)
RHO, D_WM, D_GM = 0.02, 0.10, 0.02          # per day, mm^2/day, mm^2/day: ASSUMED round values, not fitted
KILL_FRACTIONS = (0.0, 0.2, 0.4)            # ASSUMED drug kill rate as a fraction of RHO
DAYS = (0, 30, 60, 90, 120, 150, 180)
U_VISIBLE = 0.5
H = 2.0


def up1(a):
    u = np.repeat(np.repeat(np.repeat(a, 2, 0), 2, 1), 2, 2)
    return np.concatenate([u, u[:, :, -1:]], 2)


def label(z, sid, tp):
    a = np.asarray(nib.Nifti1Image.from_bytes(gzip.decompress(z.read(f"{sid}/{sid}_{tp}_seg.nii.gz"))).dataobj)
    return np.ascontiguousarray(np.rint(a).astype(np.uint8)[::-1, ::-1])


def down2(a):
    a = a.astype(np.float32)
    pad = [(0, s % 2) for s in a.shape]
    a = np.pad(a, pad)
    s = a.shape
    return a.reshape(s[0] // 2, 2, s[1] // 2, 2, s[2] // 2, 2).mean((1, 3, 5))


def simulate(core, brain, wm, gm, kill_frac):
    """Explicit Fisher-Kolmogorov on the 2 mm grid, dt = 1 day, no flux out of the brain mask. Returns {day: u}."""
    D = ((D_WM * wm + D_GM * gm) * brain).astype(np.float32)
    from scipy import ndimage
    reg = ndimage.binary_dilation(core > 0.5, iterations=2) & (brain > 0)
    D[reg] = np.maximum(D[reg], D_GM)   # tissue maps are empty inside the tumour: give it grey-matter diffusivity
    u = (core * brain).astype(np.float32)
    kill = kill_frac * RHO
    snaps = {0: u.copy()}
    for day in range(1, max(DAYS) + 1):
        du = np.zeros_like(u)
        for ax in range(3):
            a, b = [slice(None)] * 3, [slice(None)] * 3
            a[ax], b[ax] = slice(0, -1), slice(1, None)
            a, b = tuple(a), tuple(b)
            df = 2 * D[a] * D[b] / (D[a] + D[b] + 1e-12)
            flux = df * (u[b] - u[a]) / H ** 2
            du[a] += flux
            du[b] -= flux
        u = np.clip(u + du + RHO * u * (1 - u) * brain - kill * u, 0, 1)
        if day in DAYS:
            snaps[day] = u.copy()
    return snaps


def mesh(mask, spacing, smooth=0.0, step=1):
    from scipy import ndimage
    from skimage import measure
    v = mask.astype(np.float32)
    if smooth:
        v = ndimage.gaussian_filter(v, smooth)
    if v.max() < 0.5:
        return None
    verts, faces, _, _ = measure.marching_cubes(np.pad(v, 1), 0.5, step_size=step)
    return ((verts - 1) * spacing).astype(np.float32), faces.astype(np.int32)


def tr(m, name, color, opacity, group):
    import plotly.graph_objects as go
    if m is None:
        v, f = np.zeros((3, 3), np.float32), np.array([[0, 1, 2]], np.int32)
    else:
        v, f = m
    return go.Mesh3d(x=v[:, 0], y=v[:, 1], z=v[:, 2], i=f[:, 0], j=f[:, 1], k=f[:, 2], name=name, color=color,
                     opacity=opacity, showlegend=True, visible=False, flatshading=False, legendgroup=group,
                     lighting=dict(ambient=0.5, diffuse=0.8, specular=0.2), hoverinfo="name")


def pick_patients(df):
    d = (df["time1_final_all"] - df["time1_standard_all"]).to_numpy()
    order = np.argsort(d, kind="stable")
    return [str(df["pid"].iloc[order[int(round(q * (len(d) - 1)))]]) for q in QUANTILES]


def main():
    import plotly.graph_objects as go
    pairs = pd.read_csv(H3OUT / "ucsf_test_pairs.csv")
    pairs["pid"] = pairs["pid"].astype(str)
    pids = pick_patients(pairs)
    cohort_gain = float((pairs['time1_final_all'] - pairs['time1_standard_all']).mean())
    L = np.load(H3OUT / "ucsf_logits_2mm.npz")
    z = zipfile.ZipFile(ZIP)
    traces, meta, info = [], [], []
    for pi, sid in enumerate(pids):
        x = np.load(H3 / "inputs" / f"{sid}.npz")
        seg, brain = np.rint(x["seg"]).astype(np.uint8), x["brain"] > 0
        core = pg.core(seg)
        std = pg.standard_plan(seg, brain)
        k = int(std.sum())
        score = 1.0 / (1.0 + np.exp(-up1(L[f"ucsf_{sid}"].astype(np.float32))))
        plan = pg.model_plan(score, seg, brain, k=k)
        lab = label(z, sid, "time1")
        rec = np.isin(lab, (1, 2, 3)) & brain
        cs, cm = pg.coverage(rec, std), pg.coverage(rec, plan)
        row = pairs[pairs["pid"] == sid].iloc[0]
        assert abs(cs - row["time1_standard_all"]) < 1e-6 and abs(cm - row["time1_final_all"]) < 1e-6, sid
        ml = 1e-3
        info.append({"pid": sid, "gain_percentile": float(100 * QUANTILES[pi]), "cohort_mean_gain": cohort_gain, "plan_volume_ml": k * ml, "core_ml": float(core.sum() * ml),
                     "recurrence_ml": float(rec.sum() * ml), "coverage_standard_all": cs, "coverage_model_all": cm,
                     "gain": cm - cs, "matches_script133_csv": True})
        statics = [("Brain", mesh(brain, 1.0, smooth=1.5, step=3), "#9aa3ad", 0.07, "brain"),
                   ("Tumour core (pre-op)", mesh(core, 1.0), "#3b82f6", 0.7, "core"),
                   ("Standard 15 mm plan", mesh(std, 1.0, step=2), "#cbd5e1", 0.10, "std"),
                   ("Model plan (same volume)", mesh(plan, 1.0, step=2), "#f59e0b", 0.35, "model"),
                   ("Real recurrence (later scan)", mesh(rec, 1.0), "#22c55e", 0.35, "rec")]
        for name, m, col, op, grp in statics:
            traces.append(tr(m, name, col, op, f"{pi}{grp}"))
            meta.append({"p": pi, "kind": "static", "layer": grp})
        d2 = dict(core=down2(core) >= 0.25, brain=down2(brain) >= 0.5, wm=down2(x["wm"]), gm=down2(x["gm"]))
        for kf in KILL_FRACTIONS:
            snaps = simulate(d2["core"], d2["brain"], d2["wm"], d2["gm"], kf)
            scen = "none" if kf == 0 else f"drug{int(kf * 100)}"
            sname = "No treatment" if kf == 0 else f"Assumed drug, kill = {int(kf * 100)}% of growth rate"
            vols = []
            for d in DAYS:
                vis = snaps[d] >= U_VISIBLE
                vols.append(float(vis.sum() * H ** 3 * ml))
                traces.append(tr(mesh(vis, H), f"Simulated tumour, day {d}", "#ef4444", 0.55, f"{pi}sim{scen}{d}"))
                meta.append({"p": pi, "kind": "sim", "scen": scen, "day": d})
            info[-1].setdefault("sim_volume_ml", {})[sname] = dict(zip(map(str, DAYS), vols))
        print(sid, round(cs, 3), "->", round(cm, 3), flush=True)

    fig = go.Figure(traces)
    fig.update_layout(scene=dict(aspectmode="data", xaxis_visible=False, yaxis_visible=False, zaxis_visible=False,
                                 bgcolor="#0b1220"),
                      paper_bgcolor="#0b1220", font=dict(color="#e5e7eb"), legend=dict(x=0.01, y=0.98),
                      margin=dict(l=0, r=0, t=0, b=0), height=700)
    script = f"""
var META={json.dumps(meta)}, INFO={json.dumps(info)}, DAYS={json.dumps(list(DAYS))};
var RISK=['brain','core','std','model','rec'];
function gd(){{return document.getElementsByClassName('plotly-graph-div')[0];}}
function upd(){{
 var p=+document.getElementById('pat').value, mode=document.getElementById('mode').value,
     scen=document.getElementById('scen').value, di=+document.getElementById('day').value, i=INFO[p];
 var vis=META.map(function(m){{ if(m.p!==p) return false;
   if(m.kind==='static') return mode==='risk' ? true : (m.layer==='brain'||m.layer==='core');
   return mode==='sim' && m.scen===scen && m.day===DAYS[di]; }});
 Plotly.restyle(gd(),{{visible:vis}});
 document.getElementById('simctl').style.display = mode==='sim' ? 'inline' : 'none';
 var t;
 if(mode==='risk') t='<b>Recurrence-risk map from the pre-op MRI.</b> Orange = model plan, grey = standard 15 mm plan, both '
   +i.plan_volume_ml.toFixed(0)+' mL. Green = where the tumour really came back. Share of the recurrence inside the plan: '
   +'standard '+(100*i.coverage_standard_all).toFixed(1)+'%, model '+(100*i.coverage_model_all).toFixed(1)+'%. '
   +'Frozen model, never trained on this centre. All recurrence labels counted. '
   +'This is an upper-part example (gain percentile '+i.gain_percentile.toFixed(0)+'); mean gain over all 59 patients is +'+(100*i.cohort_mean_gain).toFixed(1)+' points.';
 else {{ var s=document.getElementById('scen').selectedOptions[0].text, v=i.sim_volume_ml[s];
   t='<b>SIMULATION with assumed parameters, not a forecast.</b> Growth 0.02 per day, diffusion 0.10 (white) and 0.02 (grey) mm2 per day, '
   +'round values, not fitted to this patient. Drug scenarios are what-ifs with an assumed kill rate. Day '+DAYS[di]+': '
   +'simulated tumour '+v[String(DAYS[di])].toFixed(1)+' mL (start '+v['0'].toFixed(1)+' mL).'; }}
 document.getElementById('note').innerHTML=t; document.getElementById('dl').textContent='day '+DAYS[di];
}}
window.addEventListener('load',function(){{ ['pat','mode','scen','day'].forEach(function(id){{document.getElementById(id).addEventListener('input',upd);}}); setTimeout(upd,300); }});
"""
    controls = ("<div style='font:14px system-ui;color:#e5e7eb;background:#0b1220;padding:10px 16px;line-height:1.7'>"
                "<b style='font-size:16px'>GBM 4D viewer</b> &nbsp;<span style='color:#fca5a5'>Research prototype. Not for patient care.</span><br>"
                "Patient <select id='pat'>" + "".join(f"<option value='{i}'>Patient {i + 1} (UCSF {p})</option>" for i, p in enumerate(pids)) + "</select> "
                "&nbsp;View <select id='mode'><option value='risk'>Recurrence-risk map</option><option value='sim'>Growth simulation over time</option></select> "
                "<span id='simctl' style='display:none'>&nbsp;Scenario <select id='scen'>"
                "<option value='none'>No treatment</option>"
                "<option value='drug20'>Assumed drug, kill = 20% of growth rate</option>"
                "<option value='drug40'>Assumed drug, kill = 40% of growth rate</option></select> "
                "&nbsp;Time <input id='day' type='range' min='0' max='" + str(len(DAYS) - 1) + "' value='0' style='vertical-align:middle;width:260px'> "
                "<b id='dl'>day 0</b></span>"
                "<div id='note' style='margin-top:4px;max-width:1100px'></div></div>")
    html = fig.to_html(include_plotlyjs=True, full_html=True)
    html = html.replace("<body>", "<body style='margin:0;background:#0b1220'>" + controls, 1)
    html = html.replace("</body>", "<script>" + script + "</script></body>", 1)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "viewer4d.html").write_text(html, encoding="utf-8")
    (OUT / "viewer4d_info.json").write_text(json.dumps({"patients": info, "assumed": {
        "rho_per_day": RHO, "D_wm_mm2_day": D_WM, "D_gm_mm2_day": D_GM, "kill_fractions_of_rho": KILL_FRACTIONS,
        "note": "assumed round values, not fitted; simulation only"}, "selection": f"quantiles {QUANTILES} of all-recurrence gain among 59 UCSF patients"}, indent=1))
    print("wrote", OUT / "viewer4d.html", round((OUT / "viewer4d.html").stat().st_size / 1e6, 1), "MB")


if __name__ == "__main__":
    main()
