"""Track C seed registry (SIMULATION bookkeeping only; runs no simulation).

Every scenario seed, observation-noise key, and policy-RNG seed that a Track C script has used is
listed here, with its status. The point: new work must not reuse a seed that fed a sealed test or a
tuning draw.

Status values
    spent      scored once as a sealed test, or the tuning draw that chose a sealed arm. Never reuse.
    dev_used   used in development runs. Do not reuse as a test.
    reserved   set aside for later work, not yet used. Spend only for the stated purpose.

Ranges are inclusive (low, high). A script that uses "base + k" is given a block of 100 scenario
seeds and 1,000 noise keys, which is more than it needs, so nobody reuses a nearby integer.

Reserved on 2026-11-03 by paper/trackc_restart_plan_2026-11-03.md:
    stage2_dev     scenario 20261800-20262099, noise keys 1,800,000-1,899,999
    stage3_sealed  scenario 20262300-20262499, noise keys 2,300,000-2,399,999
                   (spend only after the Stage 3 gate and the user's OK)
The stage2_dev block shrinks as seeds are spent: a dev_used entry is added after the run's manifest is written,
and the reserved stage2_dev ranges are computed from the dev_used entries (_reserved_stage2_dev), so a reserved
block never touches a used seed.

Run:  python src/trackc_seed_registry.py            check and print
      python src/trackc_seed_registry.py --write    also write output/trackc_seed_registry.json
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REGISTRY_PATH = ROOT / "output" / "trackc_seed_registry.json"
DATE = "2026-11-06"


def _entry(script, label, status, scen, noise, evidence=(), purpose=None, note="", policy=()):
    return {"script": script, "label": label, "status": status, "purpose": purpose,
            "scenario_seeds": [list(r) for r in scen], "noise_keys": [list(r) for r in noise],
            "policy_seeds": [list(r) for r in policy],
            "evidence": [list(e) for e in evidence], "note": note}


STAGE2_DEV_SCENARIO_BLOCK = (20261800, 20262099)
STAGE2_DEV_NOISE_BLOCK = (1800000, 1899999)


def _remainder(block, used):
    """Inclusive ranges inside `block` that no range in `used` touches, sorted."""
    lo, hi = block
    out, cur = [], lo
    for a, b in sorted((max(a, lo), min(b, hi)) for a, b in used if a <= hi and b >= lo):
        if a > cur:
            out.append((cur, a - 1))
        cur = max(cur, b + 1)
    if cur <= hi:
        out.append((cur, hi))
    return out


def _reserved_stage2_dev(entries):
    """The reserved stage2_dev entry: the block minus every range that a used entry already holds."""
    used = [e for e in entries if e["status"] != "reserved"]
    scen = _remainder(STAGE2_DEV_SCENARIO_BLOCK, [tuple(r) for e in used for r in e["scenario_seeds"]])
    noise = _remainder(STAGE2_DEV_NOISE_BLOCK, [tuple(r) for e in used for r in e["noise_keys"]])
    return _entry("restart", "Stage 2 dev runs (restart plan)", "reserved", scen, noise, purpose="stage2_dev",
                  note="Reserved 2026-11-03; shrinks as dev seeds are spent (computed from the dev_used entries). "
                       "Dev runs only. Never a test.")


ENTRIES = [
    _entry("150", "pilot, first run", "dev_used", [(20261023, 20261025)],
           [(1023001, 1023001), (1024001, 1024002), (1025001, 1025001)],
           [("output/trackc_pilot/pilot_results.json", "20261023")]),
    _entry("150", "pilot, rerun (pilot2)", "dev_used", [(20261031, 20261033)],
           [(1031001, 1031001), (1032001, 1032002), (1033001, 1033001)],
           [("output/trackc_pilot2/pilot_results.json", "20261031")],
           note="Seed 20261031 is the shared 'script-150 training draw' reused by 153, 154 and 156."),
    _entry("153", "headroom (oracle)", "dev_used", [(20261031, 20261031)], [(1031001, 1031001)],
           [("output/trackc_headroom/headroom_results.json", "20261031")],
           note="Same draw as the 150 training draw, by design."),
    _entry("154", "containment, train and eval", "dev_used", [(20261031, 20261031), (20261042, 20261042)],
           [(1031001, 1031001), (1042001, 1042002)],
           [("output/trackc_containment/containment_results.json", "20261042")]),
    _entry("155", "wrong models (dev)", "dev_used", [(20261060, 20261069)], [(1060000, 1060999)],
           note="Seed rule from the script docstring: 20261060 + k, noise 1_060_00k*."),
    _entry("156", "guarded containment", "dev_used",
           [(20261031, 20261031), (20261071, 20261071), (20261080, 20261089)],
           [(1031001, 1031001), (1071001, 1071001), (1080000, 1080999)],
           [("output/trackc_guarded/guarded_results.json", "20261080 + k")]),
    _entry("158", "prior trade-off", "dev_used", [(20261090, 20261099)], [(1090000, 1090999)],
           note="Seed rule from the script docstring: 20261090 + k, noise 1_090_0k1 and 1_090_0k2."),
    _entry("159", "numerics check", "dev_used", [(20261101, 20261103)], [(1101001, 1101001)],
           [("output/trackc_numerics_check/numerics_results.json", "20261101")]),
    _entry("166", "filter calibration, development", "dev_used", [(20261111, 20261111)], [],
           [("output/trackc_filter_calibration/development_results.json", "20261111")],
           note="Noise comes from rng(seed) inside the script; no separate noise key."),
    _entry("166", "filter calibration, sealed", "spent", [(20261112, 20261112)], [],
           [("output/trackc_filter_calibration/sealed_results.json", "20261112")]),
    _entry("167", "tuning draw shared by 167, 169, 171, 173", "spent", [(20261121, 20261121)],
           [(1121001, 1121002)], [("output/trackc_sealed_containment/frozen_protocol.json", "20261121")]),
    _entry("167", "sealed draw", "spent", [(20261122, 20261122)], [(1122001, 1122004)],
           [("output/trackc_sealed_containment/frozen_protocol.json", "20261122")]),
    _entry("168", "sealed wrong models", "spent", [(20261200, 20261299)], [(1200000, 1200999)],
           [("output/trackc_sealed_wrong_models/sealed_results.json", "20261200 + k")]),
    _entry("169", "sealed new axes", "spent", [(20261300, 20261399)], [(1300000, 1300999)],
           [("output/trackc_sealed_new_axes/sealed_results.json", "scenario 20261300 + k")]),
    _entry("170", "sealed missed scans", "spent", [(20261400, 20261499)], [(1400000, 1400999)],
           [("output/trackc_sealed_missed_scans/sealed_results.json", "scenario 20261400 + k")]),
    _entry("171", "sealed tail hybrid", "spent", [(20261500, 20261599)], [(1500000, 1500999)],
           [("output/trackc_sealed_tail_hybrid/sealed_results.json", "scenario 20261500 + k")]),
    _entry("173", "sealed tail trigger", "spent", [(20261700, 20261799)], [(1700000, 1700999)],
           [("output/trackc_sealed_tail_trigger/sealed_results.json", "scenario 20261700 + k")]),
    _entry("176", "dev plasticity", "dev_used", [(20262100, 20262109)],
           [(1300001, 1300002), (1300011, 1300012), (1300021, 1300022)],
           [("output/trackc_dev_plasticity/dev_results.json", "scenario 20262100 + k")],
           note="Script 176 says its seeds were 'unused before'. That is true for the scenario seeds. "
                "Its noise keys 1_300_000 + 10k + rep equal the keys that sealed script 169 used for "
                "k = 0-2, on different scenarios (20261300 + k)."),
    _entry("177", "dev plasticity sensitivity", "dev_used", [(20262200, 20262200)], [(1620000, 1620001)],
           [("output/trackc_dev_plasticity_sensitivity/dev_sensitivity_results.json", "20262200")]),
    _entry("190", "Stage 2A tail estimator", "dev_used", [(20261800, 20261800)], [(1800001, 1800016)],
           [("output/trackc_stage2_A_tail_estimator/results.json", "20261800")],
           note="Dev run started 2026-11-03. Noise key = [1_800_000 + replicate, chunk]."),
    _entry("192", "Stage 2C calibration of the standard arm", "dev_used", [(20261801, 20261802)], [(1810000, 1810000)],
           [("output/trackc_stage2_C_calibration/results.json", "20261801")],
           note="Dev run 2026-11-04. Worlds: Latin hypercube seed 20261801. Patients: seed [20261802, world index]. "
                "Noise key [1_810_000, chunk]; the open-loop arm ignores readings."),
    _entry("193", "Stage 2D negative controls on the kept worlds", "dev_used", [(20261803, 20261803)], [(1820000, 1820000)],
           [("output/trackc_stage2_D_negative_controls/results.json", "20261803")],
           note="Dev run 2026-11-04. Fresh patients for the worlds kept by script 192: seed [20261803, world index]; "
                "noise key [1_820_000, world index]; the open-loop arms ignore readings."),
    _entry("194", "Stage 2E regime map on the controlled worlds", "dev_used", [(20261804, 20261805)], [(1830000, 1830000)],
           [("output/trackc_stage2_E_regime_map/results.json", "20261804")],
           note="Dev run 2026-11-04. Evaluation patients seed [20261804, world index], selection patients seed "
                "[20261805, world index], noise key [1_830_000, world index], the same key in every arm."),
    _entry("198", "Stage 2I integration-step check on the controlled worlds", "dev_used", [(20261806, 20261806)],
           [(1840000, 1840000)], [("output/trackc_stage2_I_numerics/results.json", "20261806")],
           note="Dev run 2026-11-04. Patients seed [20261806, world index], noise key [1_840_000, world index]; "
                "the open-loop arms ignore readings."),
    _entry("199", "Stage 2J response-depth check on the controlled worlds", "dev_used", [(20261807, 20261807)], [],
           [("output/trackc_stage2_J_response_depth/results.json", "20261807")],
           note="Dev run 2026-11-04. Patients seed [20261807, world index]; open loop, no noise key."),
    _entry("192", "Stage 2C calibration, fresh-seed replicate (config single_rep)", "dev_used", [(20261810, 20261811)],
           [(1850000, 1850000)], [("output/trackc_stage2_C_calibration_rep/results.json", "20261810")],
           note="Dev run 2026-11-04. Same pipeline as 192 with fresh seeds: worlds LHS seed 20261810, patients seed "
                "[20261811, world index], noise key [1_850_000, chunk]; adds the 3- and 4-course variants. "
                "The matching negative controls (193) use 20261812 and noise key 1_860_000, not yet spent."),
    _entry("192", "Stage 2C calibration of the mixture family (config mixture)", "dev_used", [(20261820, 20261821)],
           [(1890000, 1890000)], [("output/trackc_stage2_C_calibration_mixture/results.json", "20261820")],
           note="Dev run 2026-11-04. Mixture family (responders keep the earlier kill range, non-responders weak): worlds "
                "LHS seed 20261820, patients seed [20261821, world index], noise key [1_890_000, chunk]. Later runs of "
                "this config use 20261822 (193) with noise 1_890_100, 20261823-20261824 (194) with noise 1_890_200, "
                "20261825 (198) with noise 1_890_300 and 20261826 (199); not yet spent."),
    _entry("193", "Stage 2D negative controls of the single_rep replicate", "dev_used", [(20261812, 20261812)],
           [(1860000, 1860000)], [("output/trackc_stage2_D_negative_controls_rep/results.json", "20261812")],
           note="Dev run 2026-11-04. Fresh patients for the worlds kept by the single_rep calibration: seed "
                "[20261812, world index], noise key [1_860_000, world index]; the open-loop arms ignore readings."),
    _entry("193", "Stage 2D negative controls of the mixture family", "dev_used", [(20261822, 20261822)],
           [(1890100, 1890100)], [("output/trackc_stage2_D_negative_controls_mixture/results.json", "20261822")],
           note="Dev run 2026-11-04. Fresh patients for the worlds kept by the mixture calibration: seed "
                "[20261822, world index], noise key [1_890_100, world index]; the open-loop arms ignore readings."),
    _entry("192", "Stage 2C calibration with a mixed course cap (config single_mix)", "dev_used", [(20261830, 20261831)],
           [(1891000, 1891000)], [("output/trackc_stage2_C_calibration_mix36/results.json", "20261830")],
           note="Dev run 2026-11-04. Same pipeline as 192 with std_courses=mix36 (cap 3 for 63% of the patients, cap 6 for 37%): "
                "worlds LHS seed 20261830, patients seed [20261831, world index], noise key [1_891_000, chunk]. Later runs of "
                "this config use 20261832 (193) with noise 1_891_100, 20261833-20261834 (194) with noise 1_891_200, "
                "20261835 (198) with noise 1_891_300 and 20261836 (199); not yet spent."),
    _entry("192", "Stage 2C calibration of the mixture family with a mixed course cap (config mixture_mix)", "dev_used",
           [(20261840, 20261841)], [(1892000, 1892000)],
           [("output/trackc_stage2_C_calibration_mixture_mix36/results.json", "20261840")],
           note="Dev run 2026-11-04. Mixture family with std_courses=mix36: worlds LHS seed 20261840, patients seed "
                "[20261841, world index], noise key [1_892_000, chunk]. Later runs of this config use 20261842 (193) with "
                "noise 1_892_100, 20261843-20261844 (194) with noise 1_892_200, 20261845 (198) with noise 1_892_300, "
                "20261846 (199) and 20261847 (201) with noise 1_892_400; not yet spent."),
    _entry("201", "Stage 2K stratified MGMT anchor of the mixture family", "dev_used", [(20261827, 20261827)],
           [(1890400, 1890400)], [("output/trackc_stage2_K_stratified_anchor_mixture/results.json", "20261827")],
           note="Dev run 2026-11-04. Fresh patients seed [20261827, world index], noise key [1_890_400, world index]; the "
                "open-loop arms ignore readings."),
    _entry("198", "Stage 2I integration-step check of the mixture family", "dev_used", [(20261825, 20261825)],
           [(1890300, 1890300)], [("output/trackc_stage2_I_numerics_mixture/results.json", "20261825")],
           note="Dev run 2026-11-04. Patients seed [20261825, world index], noise key [1_890_300, world index]; the "
                "open-loop arms ignore readings."),
    _entry("199", "Stage 2J response-depth check of the mixture family", "dev_used", [(20261826, 20261826)], [],
           [("output/trackc_stage2_J_response_depth_mixture/results.json", "20261826")],
           note="Dev run 2026-11-04. Patients seed [20261826, world index]; open loop, no noise key."),
    _entry("198", "Stage 2I integration-step check of the single_rep replicate", "dev_used", [(20261815, 20261815)],
           [(1880000, 1880000)], [("output/trackc_stage2_I_numerics_rep/results.json", "20261815")],
           note="Dev run 2026-11-04. Patients seed [20261815, world index], noise key [1_880_000, world index]; the "
                "open-loop arms ignore readings."),
    _entry("199", "Stage 2J response-depth check of the single_rep replicate", "dev_used", [(20261816, 20261816)], [],
           [("output/trackc_stage2_J_response_depth_rep/results.json", "20261816")],
           note="Dev run 2026-11-04. Patients seed [20261816, world index]; open loop, no noise key."),
    _entry("194", "Stage 2E regime map of the mixture family", "dev_used", [(20261823, 20261824)], [(1890200, 1890200)],
           [("output/trackc_stage2_E_regime_map_mixture/results.json", "20261823")],
           note="Dev run 2026-11-04. Evaluation patients seed [20261823, world index], selection patients seed "
                "[20261824, world index], noise key [1_890_200, world index], the same key in every arm. A first launch "
                "of this run was stopped a few minutes after its manifest was written, before any result existed, so "
                "that the response-guided arms could be added; the run was relaunched with the same seeds."),
    _entry("194", "Stage 2E regime map of the single_rep replicate", "dev_used", [(20261813, 20261814)], [(1870000, 1870000)],
           [("output/trackc_stage2_E_regime_map_rep/results.json", "20261813")],
           note="Dev run 2026-11-04. Evaluation patients seed [20261813, world index], selection patients seed "
                "[20261814, world index], noise key [1_870_000, world index], the same key in every arm. A first launch "
                "of this run was stopped a few minutes after its manifest was written, before any result existed, so "
                "that the response-guided arms could be added; the run was relaunched with the same seeds."),
    _entry("193", "Stage 2D negative controls of the single_mix configuration", "dev_used", [(20261832, 20261832)],
           [(1891100, 1891100)], [("output/trackc_stage2_D_negative_controls_mix36/results.json", "20261832")],
           note="Dev run 2026-11-04. Fresh patients for the worlds kept by the single_mix calibration (mixed course cap, "
                "63% of the patients allowed 3 courses and 37% allowed 6): seed [20261832, world index], noise key "
                "[1_891_100, world index]; the open-loop arms ignore readings."),
    _entry("193", "Stage 2D negative controls of the mixture_mix configuration", "dev_used", [(20261842, 20261842)],
           [(1892100, 1892100)], [("output/trackc_stage2_D_negative_controls_mixture_mix36/results.json", "20261842")],
           note="Dev run 2026-11-04. Fresh patients for the worlds kept by the mixture_mix calibration (mixed course cap): "
                "seed [20261842, world index], noise key [1_892_100, world index]; the open-loop arms ignore readings."),
    _entry("194", "Stage 2E regime map of the single_mix configuration", "dev_used", [(20261833, 20261834)], [(1891200, 1891200)],
           [("output/trackc_stage2_E_regime_map_mix36/results.json", "20261833")],
           note="Dev run 2026-11-04. Mixed course cap (63% of the patients allowed 3 courses, 37% allowed 6). Evaluation patients "
                "seed [20261833, world index], selection patients seed [20261834, world index], noise key [1_891_200, world index], "
                "the same key in every arm."),
    _entry("194", "Stage 2E regime map of the mixture_mix configuration", "dev_used", [(20261843, 20261844)], [(1892200, 1892200)],
           [("output/trackc_stage2_E_regime_map_mixture_mix36/results.json", "20261843")],
           note="Dev run 2026-11-04. Mixed course cap. Evaluation patients seed [20261843, world index], selection patients seed "
                "[20261844, world index], noise key [1_892_200, world index], the same key in every arm."),
    _entry("198", "Stage 2I integration-step check of the single_mix configuration", "dev_used", [(20261835, 20261835)],
           [(1891300, 1891300)], [("output/trackc_stage2_I_numerics_mix36/results.json", "20261835")],
           note="Dev run 2026-11-04. Mixed course cap. Patients seed [20261835, world index], noise key [1_891_300, world index]; "
                "the open-loop arms ignore readings."),
    _entry("198", "Stage 2I integration-step check of the mixture_mix configuration", "dev_used", [(20261845, 20261845)],
           [(1892300, 1892300)], [("output/trackc_stage2_I_numerics_mixture_mix36/results.json", "20261845")],
           note="Dev run 2026-11-04. Mixed course cap. Patients seed [20261845, world index], noise key [1_892_300, world index]; "
                "the open-loop arms ignore readings."),
    _entry("199", "Stage 2J response-depth check of the single_mix configuration", "dev_used", [(20261836, 20261836)], [],
           [("output/trackc_stage2_J_response_depth_mix36/results.json", "20261836")],
           note="Dev run 2026-11-04. Mixed course cap. Patients seed [20261836, world index]; open loop, no noise key."),
    _entry("199", "Stage 2J response-depth check of the mixture_mix configuration", "dev_used", [(20261846, 20261846)], [],
           [("output/trackc_stage2_J_response_depth_mixture_mix36/results.json", "20261846")],
           note="Dev run 2026-11-04. Mixed course cap. Patients seed [20261846, world index]; open loop, no noise key."),
    _entry("201", "Stage 2K stratified MGMT anchor of the mixture_mix configuration", "dev_used", [(20261847, 20261847)],
           [(1892400, 1892400)], [("output/trackc_stage2_K_stratified_anchor_mixture_mix36/results.json", "20261847")],
           note="Dev run 2026-11-04. Mixed course cap. Fresh patients seed [20261847, world index], noise key [1_892_400, world "
                "index]; the open-loop arms ignore readings."),
    _entry("203", "Stage 2H value of knowing a patient-level quantity, mixture family", "dev_used", [(20261854, 20261855)], [],
           [("output/trackc_stage2_H_value_of_information_mixture/results.json", "20261854")],
           note="Dev run 2026-11-04. Selection patients seed [20261854, world index], evaluation patients seed [20261855, world "
                "index]; open loop, no noise key. The other configurations have their own entries."),
    _entry("203", "Stage 2H value of knowing a patient-level quantity, single_mix configuration", "dev_used", [(20261856, 20261857)], [],
           [("output/trackc_stage2_H_value_of_information_mix36/results.json", "20261856")],
           note="Dev run 2026-11-04. Mixed course cap in the worlds of the configuration. Selection patients seed [20261856, world "
                "index], evaluation patients seed [20261857, world index]; open loop, no noise key."),
    _entry("203", "Stage 2H value of knowing a patient-level quantity, single configuration", "dev_used", [(20261850, 20261851)], [],
           [("output/trackc_stage2_H_value_of_information/results.json", "20261850")],
           note="Dev run 2026-11-04. Selection patients seed [20261850, world index], evaluation patients seed [20261851, world "
                "index]; open loop, no noise key."),
    _entry("203", "Stage 2H value of knowing a patient-level quantity, single_rep configuration", "dev_used", [(20261852, 20261853)], [],
           [("output/trackc_stage2_H_value_of_information_rep/results.json", "20261852")],
           note="Dev run 2026-11-04. Selection patients seed [20261852, world index], evaluation patients seed [20261853, world "
                "index]; open loop, no noise key."),
    _entry("203", "Stage 2H value of knowing a patient-level quantity, mixture_mix configuration", "dev_used", [(20261858, 20261859)], [],
           [("output/trackc_stage2_H_value_of_information_mixture_mix36/results.json", "20261858")],
           note="Dev run 2026-11-04. Mixed course cap in the worlds of the configuration. Selection patients seed [20261858, world "
                "index], evaluation patients seed [20261859, world index]; open loop, no noise key."),
    _entry("204", "Stage 2M what the RANO-like headroom consists of, single configuration", "dev_used", [(20261860, 20261861)], [],
           [("output/trackc_stage2_M_endpoint_mechanism/results.json", "20261860")],
           note="Dev run 2026-11-04. Selection patients seed [20261860, world index], evaluation patients seed [20261861, world "
                "index]; open loop, no noise key."),
    _entry("204", "Stage 2M what the RANO-like headroom consists of, single_rep configuration", "dev_used", [(20261862, 20261863)], [],
           [("output/trackc_stage2_M_endpoint_mechanism_rep/results.json", "20261862")],
           note="Dev run 2026-11-04. Selection patients seed [20261862, world index], evaluation patients seed [20261863, world "
                "index]; open loop, no noise key."),
    _entry("204", "Stage 2M what the RANO-like headroom consists of, mixture configuration", "dev_used", [(20261864, 20261865)], [],
           [("output/trackc_stage2_M_endpoint_mechanism_mixture/results.json", "20261864")],
           note="Dev run 2026-11-04. Selection patients seed [20261864, world index], evaluation patients seed [20261865, world "
                "index]; open loop, no noise key."),
    _entry("205", "Stage 2G finer search for the informed arm and the fixed frontier, mixture configuration", "dev_used",
           [(20261874, 20261875)], [],
           [("output/trackc_stage2_G_informed_search_mixture/results.json", "20261874")],
           note="Dev run 2026-11-04. Selection patients seed [20261874, world index], evaluation patients seed [20261875, world "
                "index]; open loop, no noise key."),
    _entry("205", "Stage 2G finer search for the informed arm and the fixed frontier, single configuration", "dev_used",
           [(20261870, 20261871)], [],
           [("output/trackc_stage2_G_informed_search/results.json", "20261870")],
           note="Dev run 2026-11-04. Selection patients seed [20261870, world index], evaluation patients seed [20261871, world "
                "index]; open loop, no noise key."),
    _entry("205", "Stage 2G finer search for the informed arm and the fixed frontier, single_rep configuration", "dev_used",
           [(20261872, 20261873)], [],
           [("output/trackc_stage2_G_informed_search_rep/results.json", "20261872")],
           note="Dev run 2026-11-04. Selection patients seed [20261872, world index], evaluation patients seed [20261873, world "
                "index]; open loop, no noise key."),
    _entry("205", "Stage 2G finer search for the informed arm and the fixed frontier, single_mix configuration", "dev_used",
           [(20261876, 20261877)], [],
           [("output/trackc_stage2_G_informed_search_mix36/results.json", "20261876")],
           note="Dev run 2026-11-04. Mixed course cap in the worlds of the configuration. Selection patients seed [20261876, world "
                "index], evaluation patients seed [20261877, world index]; open loop, no noise key."),
    _entry("205", "Stage 2G finer search for the informed arm and the fixed frontier, mixture_mix configuration", "dev_used",
           [(20261878, 20261879)], [],
           [("output/trackc_stage2_G_informed_search_mixture_mix36/results.json", "20261878")],
           note="Dev run 2026-11-04. Mixed course cap in the worlds of the configuration. Selection patients seed [20261878, world "
                "index], evaluation patients seed [20261879, world index]; open loop, no noise key."),
    _entry("213", "finite-sample signal recovery (development)", "dev_used",
           [(20261880, 20261904)], [],
           [("output/trackc_finite_sample_signal_recovery/diagnostic.json", '"seed_root": 20261880')],
           note="Development-only synthetic parameter recovery, run 2026-11-06; 25 roots keyed by f0 range, truth anchor, and "
                "sample size. No treatment outcomes or patient data."),
    _entry("215", "LOO paired-signal predictive check (development; gate not passed)", "dev_used",
           [(20261905, 20261929)], [],
           [("output/trackc_loo_predictive_check/diagnostic.json", '"seed_root": 20261905')],
           note="Frozen synthetic signal-prediction diagnostic, run 2026-11-06. Gate failed at n=18 and n=28; do not reuse "
                "these roots to tune the method. No treatment outcomes or patient data."),
    _entry("216", "ambiguity MPC initial timing attempt (discarded; seed collision)", "dev_used",
           [(20261109, 20261109), (20261111, 20261111)], [],
           [("output/trackc_robust_signal_ambiguity_20261106/timing.json", '"discarded_timing_attempt"')],
           note="Timing only; all timing values discarded and no policy outcomes saved or examined. Scenario seed 20261111 "
                "overlaps script 166 development. Controller seed 20261112 numerically matches script 166 sealed scenario "
                "seed, but belongs to a separate RNG stream. Recorded for audit; never use these seeds as evaluation data.",
           policy=[(20261110, 20261110), (20261112, 20261112)]),
    _entry("216", "ambiguity MPC timing attempts (fresh; superseded timings)", "dev_used",
           [(20262600, 20262600), (20262602, 20262602)], [],
           [("output/trackc_robust_signal_ambiguity_20261106/timing.json", '"additional_fresh_timing_attempts"')],
           note="Timing only, no policy outcomes. These single timing observations were superseded by the final timing "
                "run and are not treated as performance evidence.",
           policy=[(20262601, 20262601), (20262603, 20262603)]),
    _entry("216", "ambiguity MPC final timing run", "dev_used",
           [(20262604, 20262604)], [],
           [("output/trackc_robust_signal_ambiguity_20261106/timing.json", '"scenario_seed": 20262604')],
           note="Final synthetic timing only; no outcomes. Controller RNG seed is recorded separately as a policy seed.",
           policy=[(20262605, 20262605)]),
]
ENTRIES.append(_reserved_stage2_dev(ENTRIES))
ENTRIES.append(_entry("restart", "Stage 3 sealed run (restart plan)", "reserved", [(20262300, 20262499)],
                      [(2300000, 2399999)], purpose="stage3_sealed",
                      note="Reserved 2026-11-03. Spend only after the Stage 3 gate and the user's OK."))

# Numeric seed overlaps that exist on purpose or are recorded as a flaw. The test compares the
# computed same-stream and cross-stream overlaps with this list, so a new accidental collision fails.
KNOWN_OVERLAPS = [
    {"kind": "scenario", "range": [20261031, 20261031], "scripts": ["150", "153", "154", "156"],
     "why": "shared script-150 training draw, by design (documented in each script)"},
    {"kind": "noise", "range": [1031001, 1031001], "scripts": ["150", "153", "154", "156"],
     "why": "shared script-150 training draw, by design (documented in each script)"},
    {"kind": "noise", "range": [1300001, 1300022], "scripts": ["169", "176"],
     "why": "176 noise keys equal 169 sealed noise keys for k = 0-2, on different scenarios; "
            "176 calls its seeds 'unused before'. Recorded, no result changed."},
    {"kind": "scenario", "range": [20261111, 20261111], "scripts": ["166", "216"],
     "why": "Script 216 used this only in a discarded timing attempt; no policy outcome was saved or examined."},
    {"kind": "scenario-policy", "range": [20261112, 20261112], "scripts": ["166", "216"],
     "why": "Script 216's controller seed numerically matched script 166's sealed scenario seed; the timing attempt "
           "was discarded and no policy outcome was saved or examined."},
]


def _overlap(a, b):
    lo, hi = max(a[0], b[0]), min(a[1], b[1])
    return (lo, hi) if lo <= hi else None


def validate_entries(entries=ENTRIES):
    """Return a list of problems with the table itself (empty list = fine)."""
    problems = []
    for e in entries:
        for field in ("scenario_seeds", "noise_keys", "policy_seeds"):
            for lo, hi in e[field]:
                if not (isinstance(lo, int) and isinstance(hi, int) and lo <= hi):
                    problems.append(f"{e['script']} {e['label']}: bad range {lo, hi} in {field}")
        if e["status"] not in ("spent", "dev_used", "reserved"):
            problems.append(f"{e['script']} {e['label']}: bad status {e['status']}")
        if (e["status"] == "reserved") != (e["purpose"] is not None):
            problems.append(f"{e['script']} {e['label']}: purpose is set only on reserved entries")
    return problems


def find_overlaps(entries=ENTRIES):
    """Numeric seed overlaps between used entries; includes cross-stream collisions.
    Reserved entries are checked by check_free. Returns stable {kind, range, scripts} rows."""
    used = [e for e in entries if e["status"] != "reserved"]
    cells = {}
    streams = (("scenario", "scenario_seeds"), ("noise", "noise_keys"), ("policy", "policy_seeds"))
    for i, a in enumerate(used):
        for b in used[i + 1:]:
            for ia, (kind_a, field_a) in enumerate(streams):
                for kind_b, field_b in streams[ia:]:
                    kind = kind_a if kind_a == kind_b else f"{kind_a}-{kind_b}"
                    for ra in a[field_a]:
                        for rb in b[field_b]:
                            ov = _overlap(ra, rb)
                            if ov:
                                cells.setdefault((kind, ov), set()).update((a["script"], b["script"]))
    out = [{"kind": k, "range": list(r), "scripts": sorted(s)} for (k, r), s in cells.items()]
    return sorted(out, key=lambda o: (o["kind"], o["range"]))


def unexpected_overlaps(entries=ENTRIES, known=KNOWN_OVERLAPS):
    """Computed overlaps that are not covered by a known entry."""
    bad = []
    for o in find_overlaps(entries):
        covered = any(k["kind"] == o["kind"] and k["range"][0] <= o["range"][0]
                      and o["range"][1] <= k["range"][1] and set(o["scripts"]) <= set(k["scripts"])
                      for k in known)
        if not covered:
            bad.append(o)
    return bad


def check_free(scenario=None, noise=None, purpose=None, entries=ENTRIES, policy=None):
    """Return a list of problems with using these ranges. Empty list = free to use.

    scenario / noise / policy: inclusive (low, high) or None. Policy covers controller RNG seeds.
    purpose: 'stage2_dev' or 'stage3_sealed' to spend a reserved scenario/noise block. The request
    must then lie inside that reserved block. No policy-RNG block is reserved; check policy seeds
    without a purpose. Without a purpose, any overlap in the requested stream is a problem.
    """
    problems = []
    requests = [("scenario", "scenario_seeds", scenario), ("noise", "noise_keys", noise),
                ("policy", "policy_seeds", policy)]
    streams = (("scenario", "scenario_seeds"), ("noise", "noise_keys"), ("policy", "policy_seeds"))
    for kind, field, req in requests:
        if req is None:
            continue
        req = (int(req[0]), int(req[1]))
        if req[0] > req[1]:
            problems.append(f"{kind}: low {req[0]} is above high {req[1]}")
            continue
        inside_block = False
        for e in entries:
            for existing_kind, existing_field in streams:
                for rng in e[existing_field]:
                    ov = _overlap(req, rng)
                    if not ov:
                        continue
                    if (existing_field == field and e["status"] == "reserved" and e["purpose"] == purpose):
                        inside_block = inside_block or (rng[0] <= req[0] and req[1] <= rng[1])
                    else:
                        problems.append(f"{kind} {ov[0]}-{ov[1]} numerically overlaps {existing_kind} seed "
                                        f"in {e['script']} {e['label']} ({e['status']}, {rng[0]}-{rng[1]})")
        if purpose and not inside_block:
            problems.append(f"{kind} {req[0]}-{req[1]} is not inside the reserved block for {purpose}")
    return problems


def verify_evidence(root=ROOT, entries=ENTRIES):
    """Check that each saved output file contains the text the table cites.
    Returns (failures, skipped). A missing file is skipped (for example in a scratch copy)."""
    failures, skipped = [], []
    for e in entries:
        for rel, needle in e["evidence"]:
            path = Path(root) / rel
            if not path.exists():
                skipped.append(rel)
                continue
            if needle not in path.read_text(encoding="utf-8", errors="replace"):
                failures.append(f"{rel}: text {needle!r} not found")
    return failures, skipped


def registry_document(entries=ENTRIES):
    body = {"date": DATE, "label": "SIMULATION bookkeeping only; no result is stored here",
            "entries": entries, "known_overlaps": KNOWN_OVERLAPS, "computed_overlaps": find_overlaps(entries)}
    digest = hashlib.sha256(json.dumps(body, sort_keys=True).encode("utf-8")).hexdigest()
    return dict(body, sha256=digest)


def write_registry(path=REGISTRY_PATH):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    doc = registry_document()
    path.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    return path, doc["sha256"]


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    problems = validate_entries() + [f"unexpected overlap {o}" for o in unexpected_overlaps()]
    failures, skipped = verify_evidence()
    problems += failures
    for e in ENTRIES:
        lo = ", ".join(f"{a}-{b}" for a, b in e["scenario_seeds"])
        policy = ", ".join(f"{a}-{b}" for a, b in e["policy_seeds"]) or "-"
        print(f"{e['status']:9s} {e['script']:8s} {e['label']:48s} scenario {lo} policy {policy}")
    for o in find_overlaps():
        print(f"overlap {o['kind']} {o['range']} scripts {o['scripts']}")
    print(f"evidence files skipped (missing): {len(skipped)}")
    if "--write" in argv and not problems:
        path, digest = write_registry()
        print(f"wrote {path} sha256 {digest[:12]}")
    print("PROBLEMS:" if problems else "registry ok")
    for p in problems:
        print("  ", p)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
