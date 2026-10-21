#!/usr/bin/env python3
"""Script 144: patient timeline viewer (Amendment 20, analysis_plan_v1.md, commit 390fac2).

Builds one self-contained HTML file from the LUMIERE tool volumes and src/rano_calc.py. The viewer does no statistics of
its own: every probability is computed by rano_calc.rano_calc. Display only; not a diagnosis; not for clinical use.
Output: output/viewer_rano/viewer_rano.html, viewer_rano_info.json
"""
from __future__ import annotations

import importlib.util
import json
import math
import random
from pathlib import Path as _Path

import numpy as np

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "viewer_rano"
SRC = PROJECT_ROOT / "src"


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, SRC / fname)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


S = _load("sens142", "142_rano_sensitivity.py")
RC = _load("rano_calc", "rano_calc.py")
TOOLS = {"hd": "HD-GLIO-AUTO", "dbt": "DeepBraTumIA"}
MIN_SCANS = 3
LABELS = ["model-based proxy sensitivity scenarios from one centre (LUMIERE, 2 tools), not formal noise bounds", "not a diagnosis", "not for clinical use"]


def build_data():
    seqs = S.build_sequences()
    data, pts = {}, []
    for pid, s in seqs.items():
        w0 = s["wk"][0]
        for tool in TOOLS:
            vols = s[tool]
            elig = [(not np.isnan(v)) and v >= RC.FLOOR_MM3 for v in vols]
            if sum(elig) < MIN_SCANS:
                continue
            scans = []
            for i, v in enumerate(vols):
                rec = {"wk": int(s["wk"][i] - w0), "tp": s["tp"][i], "rating": s["rating"][i],
                       "vol": None if np.isnan(v) else float(v), "eligible": bool(elig[i])}
                if elig[i]:
                    lo, hi = RC.bracket(RC.size_band(v, v))
                    f = 1.96 / math.sqrt(2.0)
                    rec["in_lo"], rec["in_hi"] = v * math.exp(-f * lo), v * math.exp(f * lo)
                    rec["out_lo"], rec["out_hi"] = v * math.exp(-f * hi), v * math.exp(f * hi)
                earlier = [x for x, e in zip(vols[:i], elig[:i]) if e]
                rec["nadir"] = float(min(earlier)) if (elig[i] and earlier) else None
                if rec["nadir"] is not None:
                    c = RC.rano_calc(rec["nadir"], float(v))
                    rec["det"] = c["deterministic_class"]
                    rec["band"] = c["size_band"]
                    rec["s_lo"], rec["s_hi"] = c["s_bracket"]
                    rec["p_pd_lo"], rec["p_pd_hi"] = c["backward_at_low_s"]["P_PD"], c["backward_at_high_s"]["P_PD"]
                    rec["p_pr_lo"], rec["p_pr_hi"] = c["backward_at_low_s"]["P_PR"], c["backward_at_high_s"]["P_PR"]
                scans.append(rec)
            data.setdefault(pid, {})[tool] = scans
            pts.append((pid, tool))
    return data, pts


def checks(data, pts, html_path):
    rnd = random.Random(20261018)
    cand = [(p, t, i) for p, t in pts for i, r in enumerate(data[p][t]) if r.get("nadir") is not None]
    pick = rnd.sample(cand, min(30, len(cand)))
    maxdiff = 0.0
    for p, t, i in pick:
        r = data[p][t][i]
        c = RC.rano_calc(r["nadir"], r["vol"])
        maxdiff = max(maxdiff, abs(c["backward_at_low_s"]["P_PD"] - r["p_pd_lo"]), abs(c["backward_at_high_s"]["P_PD"] - r["p_pd_hi"]))
    text = html_path.read_text(encoding="utf-8")
    size_mb = html_path.stat().st_size / 1e6
    each_ok = all(sum(r["eligible"] for r in data[p][t]) >= MIN_SCANS for p, t in pts)
    return dict(check1_prob_integrity_max_abs_diff=maxdiff, check1_pass=bool(maxdiff < 1e-9), n_points_checked=len(pick),
                check2_size_mb=round(size_mb, 2), check2_pass=bool(0 < size_mb < 8),
                check3_pass=bool(each_ok), check4_labels_present={l: (l in text) for l in LABELS},
                check4_pass=bool(all(l in text for l in LABELS)))


TEMPLATE = r"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Probabilistic RANO timeline</title>
<style>
:root{--bg:#fff;--fg:#1d2733;--mut:#5b6b7b;--card:#f4f6f8;--bd:#d5dbe1}
@media (prefers-color-scheme:dark){:root{--bg:#12171d;--fg:#e6ebf0;--mut:#9aa8b6;--card:#1b222b;--bd:#2c3744}}
body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.45 system-ui,Segoe UI,Arial,sans-serif}
header{padding:14px 18px;border-bottom:1px solid var(--bd)}h1{font-size:19px;margin:0 0 4px}
.warn{color:var(--mut);font-size:13px}.bar{display:flex;gap:12px;flex-wrap:wrap;padding:12px 18px;align-items:center}
select{background:var(--card);color:var(--fg);border:1px solid var(--bd);border-radius:6px;padding:6px 8px;font-size:14px}
.cards{display:flex;gap:10px;flex-wrap:wrap;padding:0 18px}.card{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:8px 12px;font-size:13px}
#a,#b{margin:0 8px}footer{padding:14px 18px;color:var(--mut);font-size:12.5px;border-top:1px solid var(--bd)}
</style></head><body>
<header><h1>Probabilistic RANO: volume timeline with uncertainty bands</h1>
<div class="warn">__WARN__</div></header>
<div class="bar"><label>Patient <select id="pat"></select></label><label>Tool <select id="tool"></select></label><span id="cnt" class="warn"></span></div>
<div class="cards" id="cards"></div>
<div id="a" style="height:430px"></div><div id="b" style="height:250px"></div>
<footer>Bands show low and high ends of two proxy sensitivity scenarios from LUMIERE: cross-tool disagreement and AI-anchored expert correction. Neither is a formal bound on scan-rescan noise; the expert-correction proxy is near zero for larger tumours, so the inner band may collapse. Probabilities: chance that the true change from the nadir (smallest earlier volume of the same tool, at least 500 mm3) is +40% or more, flat prior, normal error. The 500 mm3 floor is our own operational rule, not a RANO 2.0 rule. Expert ratings are LUMIERE 2D RANO ratings. __WARN__</footer>
<script>__PLOTLY__</script>
<script>
const DATA=__DATA__, TOOLS=__TOOLS__;
const pat=document.getElementById('pat'), tool=document.getElementById('tool');
const pids=Object.keys(DATA).sort();
pids.forEach(p=>{const o=document.createElement('option');o.value=p;o.textContent=p;pat.appendChild(o)});
function toolsOf(p){return Object.keys(DATA[p])}
function fillTools(){tool.innerHTML='';const ts=toolsOf(pat.value);ts.forEach(t=>{const o=document.createElement('option');o.value=t;o.textContent=TOOLS[t];tool.appendChild(o)});if(ts.length>1){const o=document.createElement('option');o.value='both';o.textContent='Both tools';tool.appendChild(o)}}
const COL={hd:'#1f77b4',dbt:'#d95f02'};
function render(){
 const p=pat.value, sel=tool.value, ts=sel==='both'?toolsOf(p):[sel];
 const tr=[], bars=[]; let cards='';
 ts.forEach(t=>{const S=DATA[p][t], c=COL[t]; const el=S.filter(r=>r.eligible);
  const x=el.map(r=>r.wk);
  const rgba=(a)=>c+a;
  tr.push({x:x.concat(x.slice().reverse()),y:el.map(r=>r.out_hi).concat(el.map(r=>r.out_lo).reverse()),fill:'toself',fillcolor:c+'26',line:{width:0},hoverinfo:'skip',name:TOOLS[t]+' outer band',showlegend:false,type:'scatter',mode:'lines'});
  tr.push({x:x.concat(x.slice().reverse()),y:el.map(r=>r.in_hi).concat(el.map(r=>r.in_lo).reverse()),fill:'toself',fillcolor:c+'55',line:{width:0},hoverinfo:'skip',name:TOOLS[t]+' inner band',showlegend:false,type:'scatter',mode:'lines'});
  tr.push({x:S.map(r=>r.wk),y:S.map(r=>r.vol),mode:'lines+markers',name:TOOLS[t]+' volume',line:{color:c,width:2},marker:{color:S.map(r=>r.eligible?c:'#999'),size:S.map(r=>r.eligible?8:6),symbol:S.map(r=>r.eligible?'circle':'x')},
   text:S.map(r=>'week '+r.wk+' ('+r.tp+')<br>'+(r.vol==null?'missing':Math.round(r.vol)+' mm3')+(r.eligible?'':' (below 500 mm3 / missing: not measurable)')+(r.rating?'<br>expert: '+r.rating:'')),hoverinfo:'text'});
  const nx=S.filter(r=>r.nadir!=null);
  tr.push({x:S.filter(r=>r.eligible).map(r=>r.wk),y:S.filter(r=>r.eligible).map(r=>{const e=S.filter(q=>q.eligible&&q.wk<=r.wk).map(q=>q.vol);return Math.min(...e)}),mode:'lines',line:{color:c,dash:'dot',width:1.5,shape:'hv'},name:TOOLS[t]+' nadir',hoverinfo:'skip'});
  tr.push({x:S.filter(r=>r.eligible).map(r=>r.wk),y:S.filter(r=>r.eligible).map(r=>{const e=S.filter(q=>q.eligible&&q.wk<=r.wk).map(q=>q.vol);return 1.4*Math.min(...e)}),mode:'lines',line:{color:c,dash:'dash',width:1,shape:'hv'},name:TOOLS[t]+' +40% line',hoverinfo:'skip'});
  S.filter(r=>r.rating&&r.vol!=null).forEach(r=>{tr.push({x:[r.wk],y:[r.vol],mode:'text',text:[r.rating],textposition:'top center',showlegend:false,hoverinfo:'skip',textfont:{size:11,color:'#888'}})});
  bars.push({x:nx.map(r=>r.wk),y:nx.map(r=>r.p_pd_hi),type:'bar',name:TOOLS[t]+' P(PD) higher proxy SD',marker:{color:c+'66'},text:nx.map(r=>'det. class '+r.det+'<br>P(PD) lower/higher proxy '+r.p_pd_lo.toFixed(2)+' / '+r.p_pd_hi.toFixed(2)+'<br>P(PR) lower/higher proxy '+r.p_pr_lo.toFixed(2)+' / '+r.p_pr_hi.toFixed(2)+'<br>change-SD proxy span '+r.s_lo.toFixed(2)+' to '+r.s_hi.toFixed(2)),hoverinfo:'text'});
  bars.push({x:nx.map(r=>r.wk),y:nx.map(r=>r.p_pd_lo),type:'bar',name:TOOLS[t]+' P(PD) lower proxy SD',marker:{color:c},hoverinfo:'skip'});
  const pd=nx.filter(r=>r.det==='PD').length, unsure=nx.filter(r=>r.p_pd_hi>0.05&&r.p_pd_hi<0.95).length;
  cards+='<div class="card"><b>'+TOOLS[t]+'</b>: '+el.length+' measurable scans, '+nx.length+' with a nadir; deterministic PD '+pd+'; calls uncertain at the high end (P between 0.05 and 0.95): '+unsure+'</div>';
 });
 document.getElementById('cards').innerHTML=cards;
 const dark=matchMedia('(prefers-color-scheme:dark)').matches, fg=dark?'#e6ebf0':'#1d2733', gr=dark?'#2c3744':'#e3e8ed';
 const lay=(yt,extra)=>Object.assign({paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(0,0,0,0)',font:{color:fg,size:12},margin:{l:60,r:20,t:20,b:44},xaxis:{title:'weeks from first post-operative scan',gridcolor:gr},yaxis:{title:yt,gridcolor:gr},legend:{orientation:'h',y:-0.22}},extra);
 Plotly.react('a',tr,lay('enhancing volume (mm3, log)',{yaxis:{title:'enhancing volume (mm3, log)',type:'log',gridcolor:gr}}),{responsive:true,displaylogo:false});
 Plotly.react('b',bars,lay('P(true change from nadir >= +40%)',{barmode:'overlay',yaxis:{title:'P(true change >= +40%)',range:[0,1],gridcolor:gr}}),{responsive:true,displaylogo:false});
}
pat.onchange=()=>{fillTools();render()};tool.onchange=render;
document.getElementById('cnt').textContent=pids.length+' patients ('+__HDN__+' with HD-GLIO-AUTO, '+__DBN__+' with DeepBraTumIA)';
fillTools();render();
</script></body></html>"""


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data, pts = build_data()
    from plotly.offline import get_plotlyjs
    warn = "Display of model-based quantities only. " + "; ".join(LABELS) + "."
    hdn = sum(1 for p, t in pts if t == "hd"); dbn = sum(1 for p, t in pts if t == "dbt")
    html = (TEMPLATE.replace("__WARN__", warn).replace("__PLOTLY__", get_plotlyjs())
            .replace("__DATA__", json.dumps(data, separators=(",", ":"))).replace("__TOOLS__", json.dumps(TOOLS))
            .replace("__HDN__", str(hdn)).replace("__DBN__", str(dbn)))
    path = OUT / "viewer_rano.html"
    path.write_text(html, encoding="utf-8")
    chk = checks(data, pts, path)
    info = dict(label=warn, amendment="A20", n_patients=len(data), n_patient_tool_series=len(pts),
                n_hd=hdn, n_dbt=dbn, min_scans=MIN_SCANS, floor_mm3=RC.FLOOR_MM3, **chk)
    (OUT / "viewer_rano_info.json").write_text(json.dumps(info, indent=2, default=float))
    print(json.dumps(info, indent=2, default=float))


if __name__ == "__main__":
    main()
