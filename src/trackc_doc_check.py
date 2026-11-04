"""Read-only check of the Track C numbers marked [V] in section 6 of
paper/trackc_final_status_2026-10-31.md (addendum of 2026-11-03).

SIMULATION RECORDS ONLY. This script writes nothing. It reads saved JSON and CSV files, the
simulator constants and a few source lines, then prints PASS or FAIL per check. Exit code 1 if
any check fails. If a source file changes, a failing line check means the addendum needs an update.

Run:  python src/trackc_doc_check.py        (needs numpy and scipy, because it imports trackc_sim)
"""
from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "output"
sys.path.insert(0, str(ROOT / "src"))

FAILS: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(("PASS  " if ok else "FAIL  ") + name + (f"   [{detail}]" if detail else ""))
    if not ok:
        FAILS.append(name)


def near(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


def cvar(values, q: float = 0.10) -> float:
    k = math.ceil(q * len(values))
    return sum(sorted(values)[:k]) / k


def load_json(path: Path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def git_subject(commit: str) -> str:
    try:
        res = subprocess.run(["git", "show", "-s", "--format=%s", commit], cwd=ROOT,
                             capture_output=True, text=True, timeout=30)
        return res.stdout.strip()
    except Exception as exc:  # git missing or repo unavailable
        return f"<git error: {exc}>"


def main() -> int:
    import trackc_sim as tc

    # --- items 12 and 14: constants and solver mode
    consts = (tc.CYCLE_DAYS, tc.COURSE_DAYS, tc.N_CYCLES, tc.HORIZON, tc.BUDGET, tc.PROG_FACTOR)
    check("sim constants: 28 d cycle, 5 d course, 52 cycles, 1456 d, 12 courses, 1.2 progression",
          consts == (28, 5, 52, 1456, 12, 1.2), str(consts))
    check("default EVENT_MODE is refine", tc.EVENT_MODE == "refine", str(tc.EVENT_MODE))

    d167 = OUT / "trackc_sealed_containment"
    p167 = load_json(d167 / "frozen_protocol.json")
    check("167 protocol records event_mode refine", p167["protocol"].get("event_mode") == "refine")

    # --- item 2: one tuning draw shared by 167, 169, 171, 173
    tune = {"seed": 20261121, "n": 300, "noise": [1121001, 1121002]}
    check("167 tuning draw is seed 20261121, n 300", p167["protocol"]["tune"] == tune)
    dirs = {"169": "trackc_sealed_new_axes", "171": "trackc_sealed_tail_hybrid",
            "173": "trackc_sealed_tail_trigger"}
    sealed = {k: load_json(OUT / d / "sealed_results.json") for k, d in dirs.items()}
    for k, res in sealed.items():
        check(f"{k} protocol records the same tuning draw", res["protocol"]["tune"] == tune)
    f171 = load_json(OUT / dirs["171"] / "frozen_hybrid.json")
    f173 = load_json(OUT / dirs["173"] / "frozen_trigger.json")
    a = p167["tuning_mean_ttp"]["['risk', 0.05, 2, 0.0, True]"]
    b = f171["tuning_mean_ttp"]["['riskf', 0.05, 2, 0]"]
    check("167 wide arm and 171 K=0 arm have identical tuning mean TTP", near(a, b, 1e-9), f"{a} vs {b}")
    a = f171["tuning_mean_ttp"]["['riskf', 0.05, 2, 4]"]
    b = f173["tuning_mean_ttp"]["nominal ['riskf', 0.05, 2, 4]"]
    check("171 and 173 tuning mean TTP of the chosen arm are identical", near(a, b, 1e-9), f"{a} vs {b}")

    # --- item 5 (part): repeats used by 171 and 173
    check("171 and 173 used 2 noise repeats",
          sealed["171"]["protocol"]["n_rep"] == 2 and sealed["173"]["protocol"]["n_rep"] == 2)

    # --- item 8: comparator grid edge
    grid = [spec[1] for spec in p167["protocol"]["families"]["raw_single"]]
    tm = p167["tuning_mean_ttp"]
    vals = [tm[f"['raw', {t}, {t}]"] for t in grid]
    best = grid[max(range(len(vals)), key=lambda i: vals[i])]
    check("raw single: grid minimum is 0.6 and 0.6 is the tuning winner", min(grid) == 0.6 and best == 0.6,
          f"grid {min(grid)}..{max(grid)}, winner {best}")
    check("raw single tuning mean TTP: 683.1 at 0.6, 682.1 at 0.7",
          near(tm["['raw', 0.6, 0.6]"], 683.11, 0.01) and near(tm["['raw', 0.7, 0.7]"], 682.07, 0.01))

    # --- items 5, 6, 7, 14: the 167 per-scenario file
    path = d167 / "sealed_ttp_by_scenario_replicate.csv"
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    n_s = max(int(r["scenario"]) for r in rows) + 1
    n_r = max(int(r["replicate"]) for r in rows) + 1
    check("167 file has 400 scenarios x 4 repeats", (n_s, n_r, len(rows)) == (400, 4, 1600), f"{n_s}x{n_r}")
    arms = ["winner_containment_wide", "winner_containment_narrow", "winner_raw_single",
            "winner_periodic", "reference_zhang"]
    data = {a: [[0.0] * n_r for _ in range(n_s)] for a in arms}
    for r in rows:
        for a in arms:
            data[a][int(r["scenario"])][int(r["replicate"])] = float(r[a])

    def smean(a: str, n: int):
        return [sum(data[a][s][:n]) / n for s in range(n_s)]

    expect = {"winner_containment_wide": (171.1, 174.5, 176.0),
              "winner_raw_single": (121.1, 125.5, 132.6),
              "winner_periodic": (182.2, 182.2, 182.2)}
    for a, exp in expect.items():
        got = tuple(cvar(smean(a, n)) for n in (1, 2, 4))
        check(f"CVaR10 of scenario means at 1, 2, 4 repeats: {a}",
              all(near(g, e, 0.05) for g, e in zip(got, exp)), ", ".join(f"{g:.1f}" for g in got))
    check("front-loaded repeats are identical (open loop)",
          all(near(data["winner_periodic"][s][0], data["winner_periodic"][s][k], 1e-9)
              for s in range(n_s) for k in range(1, n_r)))
    w, p = smean("winner_containment_wide", 4), smean("winner_periodic", 4)
    diff = [w[s] - p[s] for s in range(n_s)]
    share = sum(1 for x in diff if x < 0) / n_s
    worst = sum(sorted(diff)[:math.ceil(0.1 * n_s)]) / math.ceil(0.1 * n_s)
    check("wide minus front-loaded: lower in 21.0% of scenarios", near(share, 0.210, 0.0006), f"{share:.4f}")
    check("wide minus front-loaded: worst-10% mean paired difference -42.1 d", near(worst, -42.1, 0.05), f"{worst:.2f}")
    check("wide minus front-loaded: mean difference +34.2 d", near(sum(diff) / n_s, 34.2, 0.05),
          f"{sum(diff) / n_s:.2f}")
    horizon = tc.HORIZON - 0.01
    for a, exp in (("winner_containment_narrow", 0.059), ("winner_containment_wide", 0.061),
                   ("winner_raw_single", 0.034), ("winner_periodic", 0.022)):
        frac = sum(1 for s in range(n_s) for k in range(n_r) if data[a][s][k] >= horizon) / (n_s * n_r)
        check(f"share with no progression by the horizon: {a}", near(frac, exp, 0.0006), f"{frac:.4f}")
    means = {a: sum(sum(row) for row in data[a]) / (n_s * n_r) for a in arms}
    check("arm means run from 596.0 (Zhang) to 701.9 (narrow containment)",
          near(means["reference_zhang"], 596.02, 0.01) and near(means["winner_containment_narrow"], 701.89, 0.01),
          f"{means['reference_zhang']:.2f}, {means['winner_containment_narrow']:.2f}")

    # --- item 10: per-scenario outcome files exist only for script 167
    exts = {".csv", ".npz", ".npy", ".parquet", ".pkl", ".feather"}
    folders = {"167": "trackc_sealed_containment", "168": "trackc_sealed_wrong_models",
               "169": "trackc_sealed_new_axes", "170": "trackc_sealed_missed_scans",
               "171": "trackc_sealed_tail_hybrid", "173": "trackc_sealed_tail_trigger"}
    for k, d in folders.items():
        found = sorted(q.name for q in (OUT / d).iterdir() if q.suffix.lower() in exts)
        if k == "167":
            check("167 saved a per-scenario file", len(found) >= 1, ", ".join(found))
        else:
            check(f"{k} saved no per-scenario file", not found, ", ".join(found))

    # --- item 9: oracle label
    head = "\n".join((ROOT / "src" / "trackc_oracle.py").read_text(encoding="utf-8").splitlines()[:30])
    check("oracle docstring says it is not a proven upper bound and a lower bound on the optimum",
          "NOT a proven upper bound" in head and "lower bound on the true optimum" in head)

    # --- item 11: information given to the filter (line numbers as of 2026-11-03)
    lines = (ROOT / "src" / "trackc_sim.py").read_text(encoding="utf-8").splitlines()
    check("sim line 14 says sigma is known to the filter", "sigma known to the filter" in lines[13])
    check("ParticleFilter is built from scen['sigma'] at lines 340 and 375",
          'ParticleFilter(scen["sigma"]' in lines[339] and 'ParticleFilter(scen["sigma"]' in lines[374])

    # --- item 13: commit pins are freeze/result commits, JSON pins are run-time HEADs
    pins = {"8762ff0": "frozen ablation", "f2d16bf": "sealed result", "313017f": "sealed result"}
    for commit, words in pins.items():
        subject = git_subject(commit)
        check(f"commit {commit} is a freeze/result commit", words in subject, subject[:80])
    run_heads = {"169 frozen_ablation.json": (OUT / dirs["169"] / "frozen_ablation.json", "d55dcb2"),
                 "171 sealed_results.json": (OUT / dirs["171"] / "sealed_results.json", "60234ab"),
                 "173 sealed_results.json": (OUT / dirs["173"] / "sealed_results.json", "2f7c373")}
    for label, (fpath, short) in run_heads.items():
        got = load_json(fpath).get("git_commit_at_run", "")
        check(f"{label} git_commit_at_run starts with {short}", got.startswith(short), got[:10])

    # --- item 15: the two superseded documents carry a banner and say what the addendum says
    draft = (ROOT / "paper" / "research_plan_trackC_draft.md").read_text(encoding="utf-8")
    check("research_plan_trackC_draft.md reports 169-170 and carries the superseded banner",
          "results of scripts 169-170" in draft and "SUPERSEDED for planning on 2026-11-03" in draft)
    report = ROOT / "Claude outputs" / "trackC_results_report_claude.md"
    if report.exists():  # untracked file, so it may be absent in a fresh clone
        text = report.read_text(encoding="utf-8")
        check("trackC_results_report_claude.md says 'First positive result' and carries the banner",
              "## First positive result" in text and "SUPERSEDED 2026-11-03" in text)

    # --- item 16: seed registry agrees with the saved outputs and reserved blocks are free
    import trackc_seed_registry as reg
    check("registry table is well formed", reg.validate_entries() == [])
    check("registry has no overlap other than the known ones", reg.unexpected_overlaps() == [])
    check("169 and 176 share noise keys 1300001-1300022 (recorded)",
          any(o["kind"] == "noise" and o["scripts"] == ["169", "176"] for o in reg.find_overlaps()))
    failures, skipped = reg.verify_evidence()
    check("registry evidence text is found in the saved output files", not failures,
          "; ".join(failures) or f"{len(skipped)} files skipped")
    check("stage 2 dev block is free", reg.check_free((20261800, 20262099), (1800000, 1899999), "stage2_dev") == [])
    check("stage 3 sealed block is free", reg.check_free((20262300, 20262499), (2300000, 2399999),
                                                         "stage3_sealed") == [])

    print(f"\n{len(FAILS)} failed" if FAILS else "\nall checks passed")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
