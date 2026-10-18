#!/usr/bin/env python3
"""Script 137: 4D viewer (3D space + time), research prototype for the science fair. Not for patient care.

Two modes in one self-contained HTML page (plotly.js embedded, opens offline):
  1. Recurrence-risk map (no time axis). Per patient: pre-op tumour core, the standard 15 mm plan, the model plan at the SAME
     volume (frozen H-3 model C1, which was never trained on UCSF), and the real recurrence (ALPTDG time1, labels 1-3).
     Coverage numbers are recomputed here and checked against output/predictgbm/h3_final/ucsf_test_pairs.csv.
  2. Growth simulation (time slider). Fisher-Kolmogorov growth from the pre-op core on the 2 mm grid with ASSUMED round
     parameters (rho, D_wm, D_gm are not fitted to the patient). Scenarios: no treatment, and an ASSUMED constant drug
     exposure whose kill rate is a fraction of the growth rate. This is a what-if simulation, not a forecast.
Patients: UCSF locked-test patients at the 60th to 95th percentile of all-recurrence gain (model minus standard):
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

QUANTILES = (0.6, 0.7, 0.8, 0.85, 0.9, 0.95)
RHO, D_WM, D_GM = 0.02, 0.10, 0.02          # per day, mm^2/day, mm^2/day: ASSUMED round values, not fitted
KILL_FRACTIONS = (0.0, 0.2, 0.4)            # ASSUMED drug kill rate as a fraction of RHO
DAYS = (0, 45, 90, 135, 180)
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
    fd = np.uint16 if len(verts) < 65535 else np.int32   # smaller page
    return np.round((verts - 1) * spacing, 1).astype(np.float32), faces.astype(fd)


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
        statics = [("Brain", mesh(brain, 1.0, smooth=2.0, step=6), "#8b9bb4", 0.06, "brain"),
                   ("Tumour core (pre-op)", mesh(core, 1.0, step=2), "#60a5fa", 0.75, "core"),
                   ("Standard 15 mm plan", mesh(std, 1.0, smooth=0.8, step=3), "#94a3b8", 0.10, "std"),
                   ("Model plan (same volume)", mesh(plan, 1.0, smooth=0.8, step=3), "#fbbf24", 0.38, "model"),
                   ("Real recurrence (later scan)", mesh(rec, 1.0, smooth=0.6, step=2), "#34d399", 0.35, "rec")]
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
                traces.append(tr(mesh(vis, H) if d else None, f"Simulated tumour, day {d}", "#f87171", 0.55, f"{pi}sim{scen}{d}"))
                meta.append({"p": pi, "kind": "sim", "scen": scen, "day": d})
            info[-1].setdefault("sim_volume_ml", {})[sname] = dict(zip(map(str, DAYS), vols))
        print(sid, round(cs, 3), "->", round(cm, 3), flush=True)

    fig = go.Figure(traces)
    fig.update_layout(scene=dict(aspectmode="data", xaxis_visible=False, yaxis_visible=False, zaxis_visible=False,
                                 bgcolor="#0a0f1c", camera=dict(eye=dict(x=1.5, y=-1.5, z=0.7))),
                      paper_bgcolor="#0a0f1c", font=dict(color="#e5e7eb", family="system-ui, sans-serif"),
                      legend=dict(x=0.01, y=0.99, bgcolor="rgba(10,15,28,0.6)", font=dict(size=13)),
                      margin=dict(l=0, r=0, t=0, b=0), height=680)
    script = f"""
var META={json.dumps(meta)}, INFO={json.dumps(info)}, DAYS={json.dumps(list(DAYS))}, P=0, MODE='risk';
function gd(){{return document.getElementsByClassName('plotly-graph-div')[0];}}
function $(id){{return document.getElementById(id);}}
function setP(i){{P=i; document.querySelectorAll('.pt').forEach(function(b,k){{b.classList.toggle('on',k===i);}}); upd();}}
function setM(m){{MODE=m; $('bm_risk').classList.toggle('on',m==='risk'); $('bm_sim').classList.toggle('on',m==='sim'); upd();}}
function upd(){{
 var scen=$('scen').value, di=+$('day').value, i=INFO[P];
 var vis=META.map(function(m){{ if(m.p!==P) return false;
   if(m.kind==='static') return MODE==='risk' ? true : (m.layer==='brain'||m.layer==='core');
   return MODE==='sim' && m.scen===scen && m.day===DAYS[di]; }});
 Plotly.restyle(gd(),{{visible:vis}});
 $('simctl').style.display = MODE==='sim' ? 'flex' : 'none';
 $('cards_risk').style.display = MODE==='risk' ? 'flex' : 'none';
 $('cards_sim').style.display = MODE==='sim' ? 'flex' : 'none';
 var cs=100*i.coverage_standard_all, cm=100*i.coverage_model_all;
 $('v_std').textContent=cs.toFixed(0)+'%'; $('v_mod').textContent=cm.toFixed(0)+'%'; $('v_gain').textContent='+'+(cm-cs).toFixed(0)+' pts';
 $('b_std').style.width=cs+'%'; $('b_mod').style.width=cm+'%';
 $('v_vol').textContent=i.plan_volume_ml.toFixed(0)+' mL';
 var s=$('scen').selectedOptions[0].text, v=i.sim_volume_ml[s], d=String(DAYS[di]);
 $('s_start').textContent=v['0'].toFixed(1)+' mL'; $('s_now').textContent=v[d].toFixed(1)+' mL';
 var ch=100*(v[d]/v['0']-1); $('s_ch').textContent=(ch>=0?'+':'')+ch.toFixed(0)+'%'; $('dl').textContent='day '+DAYS[di];
 $('note').innerHTML = MODE==='risk'
  ? 'Both plans cover the same volume ('+i.plan_volume_ml.toFixed(0)+' mL). Green is where the tumour really came back on a later scan. The model never saw this centre in training. Examples are from the upper part of the cohort (gain percentile '+i.gain_percentile.toFixed(0)+'); the mean gain over all 59 patients is +'+(100*i.cohort_mean_gain).toFixed(1)+' points.'
  : '<b>Simulation with assumed parameters. Not a forecast.</b> Growth 0.02 per day; diffusion 0.10 (white matter) and 0.02 (grey matter) mm\\u00b2 per day; round values, not fitted to this patient. Drug scenarios are what-ifs with an assumed kill rate.';
}}
window.addEventListener('load',function(){{ $('scen').addEventListener('input',upd); $('day').addEventListener('input',upd); setTimeout(upd,300); }});
"""
    pills = "".join(f"<button class='pt{' on' if i == 0 else ''}' onclick='setP({i})'>Patient {i + 1}<small>UCSF {p}</small></button>" for i, p in enumerate(pids))
    controls = f"""<style>
:root{{--bg:#0a0f1c;--panel:#111a2e;--line:#22304d;--text:#e5e7eb;--mut:#9aa8c2;--acc:#fbbf24}}
body{{margin:0;background:var(--bg);color:var(--text);font:14px/1.5 system-ui,sans-serif}}
.top{{padding:14px 20px 6px;display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;gap:8px}}
.top h1{{font-size:20px;margin:0;letter-spacing:.2px}}.pill{{font-size:12px;color:#fca5a5;border:1px solid #7f1d1d;border-radius:99px;padding:2px 10px}}
.row{{display:flex;gap:8px;flex-wrap:wrap;align-items:center;padding:4px 20px}}
button{{background:var(--panel);color:var(--text);border:1px solid var(--line);border-radius:8px;padding:6px 12px;cursor:pointer;font:inherit}}
button small{{display:block;font-size:11px;color:var(--mut)}}button.on{{border-color:var(--acc);box-shadow:0 0 0 1px var(--acc) inset}}
select,input{{background:var(--panel);color:var(--text);border:1px solid var(--line);border-radius:8px;padding:5px 8px;font:inherit}}
.cards{{display:flex;gap:10px;padding:8px 20px;flex-wrap:wrap}}.card{{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:8px 14px;min-width:150px}}
.card .k{{font-size:12px;color:var(--mut)}}.card .v{{font-size:24px;font-weight:600}}.bar{{height:6px;border-radius:3px;background:#1e293b;margin-top:4px}}.bar i{{display:block;height:6px;border-radius:3px}}
#note{{padding:2px 20px 8px;color:var(--mut);max-width:1100px}}
@media(max-width:700px){{.card{{min-width:120px}}.top h1{{font-size:17px}}}}
</style>
<div class='top'><h1>GBM 4D viewer</h1><span class='pill'>Research prototype. Not for patient care.</span></div>
<div class='row'>{pills}</div>
<div class='row'><button id='bm_risk' class='on' onclick="setM('risk')">Recurrence-risk map</button><button id='bm_sim' onclick="setM('sim')">Growth simulation over time</button>
<span id='simctl' style='display:none;gap:8px;align-items:center'><select id='scen'>
<option value='none'>No treatment</option><option value='drug20'>Assumed drug, kill = 20% of growth rate</option><option value='drug40'>Assumed drug, kill = 40% of growth rate</option></select>
<input id='day' type='range' min='0' max='{len(DAYS) - 1}' value='0' style='width:240px'><b id='dl'>day 0</b></span></div>
<div class='cards' id='cards_risk'>
<div class='card'><div class='k'>Recurrence inside standard 15 mm plan</div><div class='v' id='v_std'></div><div class='bar'><i id='b_std' style='background:#94a3b8'></i></div></div>
<div class='card'><div class='k'>Recurrence inside model plan, same volume</div><div class='v' id='v_mod' style='color:var(--acc)'></div><div class='bar'><i id='b_mod' style='background:#fbbf24'></i></div></div>
<div class='card'><div class='k'>Difference</div><div class='v' id='v_gain'></div></div><div class='card'><div class='k'>Plan volume</div><div class='v' id='v_vol'></div></div></div>
<div class='cards' id='cards_sim' style='display:none'><div class='card'><div class='k'>Simulated tumour at start</div><div class='v' id='s_start'></div></div>
<div class='card'><div class='k'>Simulated tumour now</div><div class='v' id='s_now'></div></div><div class='card'><div class='k'>Change</div><div class='v' id='s_ch'></div></div></div>
<div id='note'></div>"""
    html = fig.to_html(include_plotlyjs=True, full_html=True, config=dict(displaylogo=False, displayModeBar=False))
    html = html.replace("<body>", "<body>" + controls, 1)
    html = html.replace("</body>", "<script>" + script + "</script></body>", 1)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "viewer4d.html").write_text(html, encoding="utf-8")
    (OUT / "viewer4d_info.json").write_text(json.dumps({"patients": info, "assumed": {
        "rho_per_day": RHO, "D_wm_mm2_day": D_WM, "D_gm_mm2_day": D_GM, "kill_fractions_of_rho": KILL_FRACTIONS,
        "note": "assumed round values, not fitted; simulation only"}, "selection": f"quantiles {QUANTILES} of all-recurrence gain among 59 UCSF patients"}, indent=1))
    print("wrote", OUT / "viewer4d.html", round((OUT / "viewer4d.html").stat().st_size / 1e6, 1), "MB")


if __name__ == "__main__":
    main()
