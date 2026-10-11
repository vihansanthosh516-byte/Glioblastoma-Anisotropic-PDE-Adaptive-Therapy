#!/usr/bin/env python3
"""Script 123: H-2 locked TEST evaluation, run ONCE after the freeze (configs/h2_frozen.yaml + git tag h2-frozen).

Reads the frozen config, checks that the logits file hash matches the one recorded in it, then scores the final model
on the LUMIERE + RHUH test patients exactly as the benchmark: logits (2 mm) -> nearest upsample to 1 mm -> sigmoid ->
predictgbm_io.model_plan (top-k inside the brain at the standard-plan volume). The test recurrence is read through the
loader lock (refuses before the freeze).
Primary (frozen): mean paired difference in enhancing-recurrence coverage, final model minus the 15 mm standard plan,
all test patients with enhancing recurrence; patient bootstrap 95% CI (patient_stats) and Wilcoxon signed-rank p.
Secondary (frozen): all-recurrence coverage; per dataset (LUMIERE, RHUH); X2 volume: coverage at 90% of the standard
volume minus the standard plan at 100%, non-inferiority margin -0.01 (lower CI bound > -0.01); context only: released
published plans (U-Net, GliODIL, LOTI, PINN) on the same patients.
Output: output/predictgbm/final/test_results.json, test_pairs.csv
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path as _Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import wilcoxon

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "predictgbm" / "final"
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import predictgbm_io as pg  # noqa: E402
from src import patient_stats as ps  # noqa: E402


def up1(a):
    u = np.repeat(np.repeat(np.repeat(a, 2, 0), 2, 1), 2, 2)
    return np.concatenate([u, u[:, :, -1:]], 2)


def main():
    if not pg.is_frozen():
        raise SystemExit("not frozen: configs/h2_frozen.yaml and tag h2-frozen are required")
    cfg = yaml.safe_load(pg.FROZEN.read_text())
    f = OUT / "test_logits_2mm.npz"
    sha = hashlib.sha256(f.read_bytes()).hexdigest()
    if sha != cfg["final_model"]["test_logits_sha256"]:
        raise SystemExit(f"logits hash {sha} does not match the frozen config")
    z = np.load(f)
    rows = []
    for pid in sorted(z.files):
        c = pg.load_case(pid, with_recurrence=True)
        std = pg.standard_plan(c["seg"], c["brain"])
        k = int(std.sum())
        score = 1.0 / (1.0 + np.exp(-up1(z[pid].astype(np.float32))))
        re, ra = pg.rec_enhancing(c["rec"]), pg.rec_all(c["rec"])
        plan = pg.model_plan(score, c["seg"], c["brain"], k=k)
        plan90 = pg.topk_plan(score, int(round(0.9 * k)), c["brain"])
        r = {"pid": pid, "dataset": "rhuh" if pid.upper().startswith("RHUH") else "lumiere", "n_rec_enh": int(re.sum()),
             "standard_enh": pg.coverage(re, std), "final_enh": pg.coverage(re, plan), "final90_enh": pg.coverage(re, plan90),
             "standard_all": pg.coverage(ra, std), "final_all": pg.coverage(ra, plan)}
        for m in ("unet", "gliodil", "sbtc", "pinngbm"):
            p = pg.DATA / pid / f"{m}_plan.nii.gz"
            if p.exists():
                r[f"{m}_released_enh"] = pg.coverage(re, pg._load(pid, f"{m}_plan.nii.gz") > 0)
        rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "test_pairs.csv", index=False)
    use = df[df["n_rec_enh"] > 0]

    def contrast(d, a, b):
        x = (d[a] - d[b]).dropna().to_numpy()
        s = ps.summarize_delta(x)
        s["wilcoxon_p_two_sided"] = float(wilcoxon(x).pvalue) if np.any(x != 0) else 1.0
        return s

    res = {"script": "123_predictgbm_test_score", "frozen_config_sha256": hashlib.sha256(pg.FROZEN.read_bytes()).hexdigest(),
           "logits_sha256": sha, "n_test": int(len(df)), "n_with_enhancing_recurrence": int(len(use)),
           "means": {c: float(use[c].mean()) for c in use.columns if c.endswith("_enh") and c != "n_rec_enh"},
           "PRIMARY_final_minus_standard_enh": contrast(use, "final_enh", "standard_enh"),
           "secondary": {"all_recurrence": contrast(df[df["standard_all"].notna()], "final_all", "standard_all"),
                         "by_dataset": {g: contrast(d, "final_enh", "standard_enh") for g, d in use.groupby("dataset")},
                         "volume90_vs_standard100": contrast(use, "final90_enh", "standard_enh"),
                         "context_released": {m: contrast(use, f"{m}_released_enh", "standard_enh")
                                              for m in ("unet", "gliodil", "sbtc", "pinngbm") if f"{m}_released_enh" in use}}}
    lo = res["secondary"]["volume90_vs_standard100"]["mean_ci95"][0]
    res["secondary"]["volume90_noninferior_margin_minus_0.01"] = bool(lo > -0.01)
    (OUT / "test_results.json").write_text(json.dumps(res, indent=1, default=float))
    p = res["PRIMARY_final_minus_standard_enh"]
    print(json.dumps(res["means"], indent=1))
    print("PRIMARY", round(p["mean"], 4), p["mean_ci95"], "wilcoxon p", p["wilcoxon_p_two_sided"])


if __name__ == "__main__":
    main()
