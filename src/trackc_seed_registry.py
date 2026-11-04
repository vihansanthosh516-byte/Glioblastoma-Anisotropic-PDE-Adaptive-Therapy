"""Track C seed registry (SIMULATION bookkeeping only; runs no simulation).

Every scenario seed and observation-noise key that a Track C script has used is listed here, with
its status. The point: new work must not reuse a seed that fed a sealed test or a tuning draw.

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
DATE = "2026-11-03"


def _entry(script, label, status, scen, noise, evidence=(), purpose=None, note=""):
    return {"script": script, "label": label, "status": status, "purpose": purpose,
            "scenario_seeds": [list(r) for r in scen], "noise_keys": [list(r) for r in noise],
            "evidence": [list(e) for e in evidence], "note": note}


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
    _entry("restart", "Stage 2 dev runs (restart plan)", "reserved", [(20261800, 20262099)],
           [(1800000, 1899999)], purpose="stage2_dev",
           note="Reserved 2026-11-03. Dev runs only. Never a test."),
    _entry("restart", "Stage 3 sealed run (restart plan)", "reserved", [(20262300, 20262499)],
           [(2300000, 2399999)], purpose="stage3_sealed",
           note="Reserved 2026-11-03. Spend only after the Stage 3 gate and the user's OK."),
]

# Overlaps that exist on purpose or are recorded as a flaw. The test compares the computed
# overlaps with this list, so a new accidental overlap fails the test.
KNOWN_OVERLAPS = [
    {"kind": "scenario", "range": [20261031, 20261031], "scripts": ["150", "153", "154", "156"],
     "why": "shared script-150 training draw, by design (documented in each script)"},
    {"kind": "noise", "range": [1031001, 1031001], "scripts": ["150", "153", "154", "156"],
     "why": "shared script-150 training draw, by design (documented in each script)"},
    {"kind": "noise", "range": [1300001, 1300022], "scripts": ["169", "176"],
     "why": "176 noise keys equal 169 sealed noise keys for k = 0-2, on different scenarios; "
            "176 calls its seeds 'unused before'. Recorded, no result changed."},
]


def _overlap(a, b):
    lo, hi = max(a[0], b[0]), min(a[1], b[1])
    return (lo, hi) if lo <= hi else None


def validate_entries(entries=ENTRIES):
    """Return a list of problems with the table itself (empty list = fine)."""
    problems = []
    for e in entries:
        for field in ("scenario_seeds", "noise_keys"):
            for lo, hi in e[field]:
                if not (isinstance(lo, int) and isinstance(hi, int) and lo <= hi):
                    problems.append(f"{e['script']} {e['label']}: bad range {lo, hi} in {field}")
        if e["status"] not in ("spent", "dev_used", "reserved"):
            problems.append(f"{e['script']} {e['label']}: bad status {e['status']}")
        if (e["status"] == "reserved") != (e["purpose"] is not None):
            problems.append(f"{e['script']} {e['label']}: purpose is set only on reserved entries")
    return problems


def find_overlaps(entries=ENTRIES):
    """Overlaps between used entries (reserved entries are checked by check_free).
    Returns a list of {"kind", "range", "scripts"} sorted for stable comparison."""
    used = [e for e in entries if e["status"] != "reserved"]
    cells = {}
    for kind, field in (("scenario", "scenario_seeds"), ("noise", "noise_keys")):
        for i, a in enumerate(used):
            for b in used[i + 1:]:
                for ra in a[field]:
                    for rb in b[field]:
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


def check_free(scenario=None, noise=None, purpose=None, entries=ENTRIES):
    """Return a list of problems with using these ranges. Empty list = free to use.

    scenario / noise: inclusive (low, high) or None.
    purpose: 'stage2_dev' or 'stage3_sealed' to spend a reserved block. The request must then lie
    inside that reserved block. Without a purpose, any overlap with any entry is a problem.
    """
    problems = []
    requests = [("scenario", "scenario_seeds", scenario), ("noise", "noise_keys", noise)]
    for kind, field, req in requests:
        if req is None:
            continue
        req = (int(req[0]), int(req[1]))
        if req[0] > req[1]:
            problems.append(f"{kind}: low {req[0]} is above high {req[1]}")
            continue
        inside_block = False
        for e in entries:
            for rng in e[field]:
                ov = _overlap(req, rng)
                if not ov:
                    continue
                if e["status"] == "reserved" and e["purpose"] == purpose:
                    inside_block = inside_block or (rng[0] <= req[0] and req[1] <= rng[1])
                else:
                    problems.append(f"{kind} {ov[0]}-{ov[1]} overlaps {e['script']} {e['label']} "
                                    f"({e['status']}, {rng[0]}-{rng[1]})")
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
        print(f"{e['status']:9s} {e['script']:8s} {e['label']:48s} scenario {lo}")
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
