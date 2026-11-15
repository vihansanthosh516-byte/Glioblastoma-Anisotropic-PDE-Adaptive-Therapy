"""Fixed data, folds and score for the Track B hybrid search. Candidates must not edit this file.

Protocol: paper/trackB_hybrid_search_protocol_2026-11-14.md
Usage: py -3.14 -I src/track_b_hybrid_harness.py <name>   (imports src/track_b_hybrid_<name>.py)
A candidate module defines fit_predict(tr, te, comp_tr, comp_te) -> np.ndarray of predicted y_out for te.
tr/te: rows of the MU pair table (columns include y_in, y_out (tr only is fair), dt_days, r_hist, dt_prior,
phase, r_obs, core_vox_in, core_vox_out, age, mgmt, gbm, day_in, rt_s, rt_e, ad_s, ad_e, prog, patient_id).
comp_*: input-scan-only mask features (composition, edge shape, MRI brightness).
"""
from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track_b_mu_tcks as t  # noqa: E402
import track_b_mu_v3 as v3  # noqa: E402
import track_b_mu_v5 as v5  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
N_TRIES = 14
V6F = ROOT / "output" / "track_b_mu_v6_20261114" / "input_shape_features.csv"
V7F = ROOT / "output" / "track_b_mu_v7_20261114" / "input_intensity_features.csv"
V8P = ROOT / "output" / "track_b_mu_v8_20261114" / "predictions.csv"


def load_dev():
    d, _ = t.load()
    key = d[["patient_id", "tp_in"]]
    comp = v5.features_table(d).drop(columns=["patient_id", "tp_in"]).reset_index(drop=True)
    for f in (V6F, V7F):
        extra = key.merge(pd.read_csv(f), on=["patient_id", "tp_in"], how="left").drop(columns=["patient_id", "tp_in"])
        comp = pd.concat([comp, extra.reset_index(drop=True)], axis=1)
    return d, comp


def evaluate(name):
    mod = importlib.import_module(f"track_b_hybrid_{name}")
    d, comp = load_dev()
    pred = np.full(len(d), np.nan)
    for f in sorted(d.fold.unique()):
        trm, tem = (d.fold != f).values, (d.fold == f).values
        te = d[tem].drop(columns=["y_out", "core_vox_out", "r_obs"])
        p = np.asarray(mod.fit_predict(d[trm].reset_index(drop=True), te.reset_index(drop=True),
                                       comp[trm].reset_index(drop=True), comp[tem].reset_index(drop=True)), float)
        assert p.shape == (tem.sum(),) and np.isfinite(p).all(), "bad prediction shape or non-finite values"
        pred[tem] = p
    y = d.y_out.values
    base = pd.read_csv(v3.V2)
    v8 = pd.read_csv(V8P)
    errs = {k: np.abs(base[k].values - y) for k in v3.BASE}
    errs["v8"] = np.abs(v8["v8"].values - y)
    errs["candidate"] = np.abs(pred - y)
    pat = pd.DataFrame(errs).assign(pid=d.patient_id.values).groupby("pid").mean()
    mae = pat.mean().to_dict()
    best = min(v3.BASE, key=mae.get)
    level = 1 - 0.05 / N_TRIES
    rng = np.random.default_rng(t.SEED)
    ci_b = v3.boot_ci((pat.candidate - pat[best]).values, rng, level)
    ci_v8 = v3.boot_ci((pat.candidate - pat.v8).values, rng, level)
    res = {"candidate": name, "module": f"src/track_b_hybrid_{name}.py",
           "protocol": "paper/trackB_hybrid_search_protocol_2026-11-14.md",
           "label": "exploratory dev result on reused MU-Glioma-Post",
           "n_pairs": int(len(d)), "n_patients": int(len(pat)), "patient_mae_log_volume": mae,
           "best_baseline": best, "ci_level": level,
           "candidate_minus_best": {"mean": float((pat.candidate - pat[best]).mean()), "ci": ci_b},
           "candidate_minus_v8": {"mean": float((pat.candidate - pat.v8).mean()), "ci": ci_v8},
           "dev_pass": bool(ci_b[1] < 0)}
    out = ROOT / "output" / f"track_b_hybrid_{name}"
    out.mkdir(parents=True, exist_ok=True)
    (out / "result.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    pd.DataFrame({"patient_id": d.patient_id, "fold": d.fold, "y_out": y, "pred": pred}).to_csv(out / "predictions.csv", index=False)
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    evaluate(sys.argv[1])
