#!/usr/bin/env python3
"""Script 104: measured segmentation noise floor (Amendment 6 A6.5; GRAND_PLAN #12, #18; weak points W16, W24).

Pre-registered rule (A6.5, committed 4d8da2e before this script was written):
  For each LUMIERE scan with both automated segmentations, disagreement = |ln(V_HD-GLIO-AUTO / V_DeepBraTumIA)|
  for the contrast-enhancing label (CT1 row). Floor = median disagreement / 2 per volume tertile.
  H-1 SNR tables are reported with the assumed floor (script 99) and this measured floor side by side.
  Neither floor is tuned.
Volume for the tertiles = geometric mean of the two segmenters' volumes (mm3). Scans where either volume is 0 or
missing are excluded from the floor and counted.
The floor is applied to MU primary pairs by input volume (v_in voxels x 8 mm3 at 2 mm), using the LUMIERE tertile
edges; SNR = |ln((V_out + 1) / (V_in + 1))| / floor, as in script 99.

Exploratory addition (not pre-registered): change noise. For consecutive scans of a patient, the two segmenters'
ln-volume changes are compared, |dln_HD - dln_DBT|. This is closer to what matters for a forecast (same pipeline at
both scans, so a constant segmenter bias cancels).

Input : data/external/lumiere/LUMIERE-pyradiomics-{hdglioauto,deepbratumia}-features.csv, output/forecastability_pairs.csv,
        output/forecastability.json
Output: output/seg_noise_floor.json, output/seg_noise_floor.manifest.json
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path as _Path

import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output"
LUM = PROJECT_ROOT / "data" / "external" / "lumiere"
sys.path.insert(0, str(PROJECT_ROOT))

from src import patient_stats as ps  # noqa: E402
from src.run_manifest import write_run_manifest  # noqa: E402

VOXEL_MM3_MU = 8.0   # MU pairs are scored at 2 mm isotropic
LABEL = "Contrast-enhancing"


def tp_key(tp):
    m = re.match(r"week-(\d+)(?:-(\d+))?$", tp)
    return int(m.group(1)), int(m.group(2) or 0)


def load(name):
    p = pd.read_csv(LUM / f"LUMIERE-pyradiomics-{name}-features.csv",
                    usecols=["Patient", "Time point", "Label name", "Sequence", "original_shape_VoxelVolume"])
    p = p[(p["Label name"] == LABEL) & (p["Sequence"] == "CT1")]
    return p.rename(columns={"Time point": "tp", "original_shape_VoxelVolume": "vol"})[["Patient", "tp", "vol"]]


def main():
    hd, dbt = load("hdglioauto"), load("deepbratumia")
    m = hd.merge(dbt, on=["Patient", "tp"], how="outer", suffixes=("_hd", "_dbt"))
    n_scans = len(m)
    both = m[(m["vol_hd"] > 0) & (m["vol_dbt"] > 0)].copy()
    both["ln_ratio"] = np.log(both["vol_hd"] / both["vol_dbt"])
    both["dis"] = both["ln_ratio"].abs()
    both["vgeo"] = np.sqrt(both["vol_hd"] * both["vol_dbt"])
    edges = list(both["vgeo"].quantile([1 / 3, 2 / 3]))
    both["tert"] = np.digitize(both["vgeo"], edges)

    res = {"script": "104_measured_noise_floor", "rule": "Amendment 6 A6.5", "label": LABEL,
           "n_scans_any": int(n_scans), "n_scans_both_positive": int(len(both)),
           "n_patients_both_positive": int(both["Patient"].nunique()),
           "n_excluded_zero_or_missing": {"hd_zero_or_missing": int(((m["vol_hd"].fillna(0)) <= 0).sum()),
                                          "dbt_zero_or_missing": int(((m["vol_dbt"].fillna(0)) <= 0).sum())},
           "tertile_edges_mm3": edges,
           "bias_median_ln_hd_over_dbt": float(both["ln_ratio"].median()),
           "disagreement_median_ln_all": float(both["dis"].median()),
           "floor_by_tertile": {}}
    for t, g in both.groupby("tert"):
        lo, hi = ps.cluster_bootstrap(g, lambda d: d["dis"].median() / 2, patient="Patient")
        res["floor_by_tertile"][str(int(t))] = {
            "n_scans": int(len(g)), "n_patients": int(g["Patient"].nunique()),
            "vgeo_median_mm3": float(g["vgeo"].median()),
            "floor_ln": float(g["dis"].median() / 2), "floor_ln_ci95_patient_bootstrap": [float(lo), float(hi)]}

    # exploratory: change noise between consecutive scans (same segmenter at both scans)
    both["k"] = both["tp"].map(tp_key)
    rows = []
    for pid, g in both.sort_values(["Patient", "k"]).groupby("Patient"):
        g = g.reset_index(drop=True)
        for i in range(1, len(g)):
            d_hd = np.log(g.loc[i, "vol_hd"] / g.loc[i - 1, "vol_hd"])
            d_dbt = np.log(g.loc[i, "vol_dbt"] / g.loc[i - 1, "vol_dbt"])
            rows.append({"Patient": pid, "vgeo_in": g.loc[i - 1, "vgeo"], "dis_change": abs(d_hd - d_dbt),
                         "abs_change_hd": abs(d_hd)})
    ch = pd.DataFrame(rows)
    ch["tert"] = np.digitize(ch["vgeo_in"], edges)
    res["exploratory_change_noise"] = {
        "n_pairs": int(len(ch)), "n_patients": int(ch["Patient"].nunique()),
        "median_abs_change_disagreement_ln": float(ch["dis_change"].median()),
        "median_abs_change_hd_ln": float(ch["abs_change_hd"].median()),
        "by_tertile": {str(int(t)): {"n_pairs": int(len(g)), "median_abs_change_disagreement_ln": float(g["dis_change"].median())}
                       for t, g in ch.groupby("tert")}}

    # exploratory, added after the first run: the pair set above includes pre-op -> post-op steps (surgery), so
    # repeat on RANO-rated follow-up scans only (ratings PD/SD/PR/CR, as script 91). Both versions are kept.
    rat = pd.read_csv(LUM / "LUMIERE-ExpertRating-v202211.csv")
    rat.columns = ["Patient", "tp", "lt3", "nonmeas", "rating", "rationale"]
    rated = set(map(tuple, rat[rat["rating"].isin({"PD", "SD", "PR", "CR"})][["Patient", "tp"]].to_numpy()))
    fu = both[[(p, t) in rated for p, t in zip(both["Patient"], both["tp"])]]
    rows = []
    for pid, g in fu.sort_values(["Patient", "k"]).groupby("Patient"):
        g = g.reset_index(drop=True)
        for i in range(1, len(g)):
            d_hd = np.log(g.loc[i, "vol_hd"] / g.loc[i - 1, "vol_hd"])
            d_dbt = np.log(g.loc[i, "vol_dbt"] / g.loc[i - 1, "vol_dbt"])
            rows.append({"Patient": pid, "dis_change": abs(d_hd - d_dbt), "abs_change_hd": abs(d_hd)})
    cf = pd.DataFrame(rows)
    res["exploratory_change_noise_rated_followups"] = {
        "n_scans": int(len(fu)), "n_pairs": int(len(cf)), "n_patients": int(cf["Patient"].nunique()),
        "median_abs_change_disagreement_ln": float(cf["dis_change"].median()),
        "median_abs_change_hd_ln": float(cf["abs_change_hd"].median())}

    # apply both floors to the MU primary pairs (script 99 table)
    fp = pd.read_csv(OUT / "forecastability_pairs.csv")
    fp["v_in_mm3"] = fp["v_in"] * VOXEL_MM3_MU
    fp["tert"] = np.digitize(fp["v_in_mm3"], edges)
    floor_map = {int(k): v["floor_ln"] for k, v in res["floor_by_tertile"].items()}
    fp["floor_measured"] = fp["tert"].map(floor_map)
    fp["snr_measured"] = fp["abs_ln_change"] / fp["floor_measured"]
    bins = [0, 1, 2, 4, 1e9]
    labels = ["<1", "1-2", "2-4", ">4"]
    side = {}
    for name, col in (("assumed_floor_script99", "snr"), ("measured_floor_script104", "snr_measured")):
        b = pd.cut(fp[col], bins, labels=labels)
        side[name] = {"share_pairs_snr_lt_1": float((fp[col] < 1).mean()),
                      "share_pairs_snr_lt_2": float((fp[col] < 2).mean()),
                      "floor_median_ln": float(fp["noise_floor_ln" if col == "snr" else "floor_measured"].median()),
                      "persistence_dice_by_snr_bin": {str(k): {"n_pairs": int(len(g)), "mean": float(g["dice_persistence"].mean())}
                                                      for k, g in fp.groupby(b, observed=True)}}
    res["mu_primary_pairs"] = {"n_pairs": int(len(fp)), "n_patients": int(fp["patient_id"].nunique()),
                               "mu_pairs_per_lumiere_tertile": {str(int(k)): int(v) for k, v in fp["tert"].value_counts().sort_index().items()},
                               "side_by_side": side}

    # sensitivity, added 2026-10-08 (GRAND_PLAN L7, L8; post hoc, the A6.5 floor above stays primary)
    # L7: continuous floor = median regression of |ln ratio| on ln volume, / 2 (no tertile steps)
    import statsmodels.api as sm
    X = sm.add_constant(np.log(both["vgeo"].to_numpy()))
    qr = sm.QuantReg(both["dis"].to_numpy(), X).fit(q=0.5)
    vmin, vmax = both["vgeo"].min(), both["vgeo"].max()
    v_mu = fp["v_in_mm3"].clip(vmin, vmax)   # no extrapolation outside the LUMIERE range
    fp["floor_cont"] = np.maximum(qr.params[0] + qr.params[1] * np.log(v_mu), 1e-6) / 2
    # L8: keep scans where one tool finds nothing, with a +8 mm3 (one 2 mm voxel) offset, tertiles as above
    # an empty label is a blank volume in the CSV row (tool ran, found nothing) -> 0; a missing row (tool not run) is excluded
    m8 = hd.assign(vol=hd["vol"].fillna(0)).merge(dbt.assign(vol=dbt["vol"].fillna(0)), on=["Patient", "tp"],
                                                   how="inner", suffixes=("_hd", "_dbt"))
    m8 = m8[(m8["vol_hd"] > 0) | (m8["vol_dbt"] > 0)].copy()
    m8["dis"] = np.abs(np.log((m8["vol_hd"] + VOXEL_MM3_MU) / (m8["vol_dbt"] + VOXEL_MM3_MU)))
    m8["tert"] = np.digitize(np.sqrt((m8["vol_hd"] + VOXEL_MM3_MU) * (m8["vol_dbt"] + VOXEL_MM3_MU)), edges)
    floor8 = {int(t): float(g["dis"].median() / 2) for t, g in m8.groupby("tert")}
    fp["floor_incl_zero"] = fp["tert"].map(floor8)
    sens = {}
    for name, col in (("continuous_floor_L7", "floor_cont"), ("including_empty_scans_L8", "floor_incl_zero")):
        snr = fp["abs_ln_change"] / fp[col]
        sens[name] = {"floor_median_ln": float(fp[col].median()), "share_pairs_snr_lt_1": float((snr < 1).mean()),
                      "share_pairs_snr_lt_2": float((snr < 2).mean())}
    sens["continuous_floor_L7"].update({"quantreg_intercept": float(qr.params[0]), "quantreg_slope_per_ln_mm3": float(qr.params[1]),
                                        "lumiere_volume_range_mm3": [float(vmin), float(vmax)]})
    sens["including_empty_scans_L8"].update({"n_scans": int(len(m8)), "floor_by_tertile": {str(k): v for k, v in floor8.items()}})
    res["sensitivity_post_hoc"] = sens
    res["weakest_points"] = [
        "Inter-method disagreement (two automated tools) is not test-retest noise; a constant bias between tools adds to it.",
        "LUMIERE contrast-enhancing label vs MU core (labels 1+3, includes necrosis): different targets and pipelines.",
        "Floor is per volume tertile, not per patient; MU volumes outside the LUMIERE range get the end tertile.",
        "Scans where one tool finds nothing are excluded, so the floor ignores the hardest small-lesion cases."]

    out = OUT / "seg_noise_floor.json"
    out.write_text(json.dumps(res, indent=1, default=float))
    write_run_manifest("seg_noise_floor_104", OUT / "seg_noise_floor.manifest.json", script="src/104_measured_noise_floor.py",
                       seed=ps.DEFAULT_SEED, config={"label": LABEL, "voxel_mm3_mu": VOXEL_MM3_MU, "n_boot": ps.DEFAULT_N_BOOT},
                       inputs=[LUM / "LUMIERE-pyradiomics-hdglioauto-features.csv",
                               LUM / "LUMIERE-pyradiomics-deepbratumia-features.csv", OUT / "forecastability_pairs.csv"],
                       dataset="LUMIERE (floor); MU-Glioma-Post (applied)", patient_split="none (descriptive)",
                       primary_endpoint="median |ln(V_HD/V_DBT)| / 2 per volume tertile")
    print(json.dumps({k: res[k] for k in ("n_scans_both_positive", "tertile_edges_mm3", "bias_median_ln_hd_over_dbt",
                                          "floor_by_tertile", "exploratory_change_noise")}, indent=1, default=float))
    print(json.dumps(side, indent=1, default=float))


if __name__ == "__main__":
    main()
