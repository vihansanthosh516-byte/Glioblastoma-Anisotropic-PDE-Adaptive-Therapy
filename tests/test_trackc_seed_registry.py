"""Track C seed registry (src/trackc_seed_registry.py): table is well formed, overlaps are only the
known ones, reserved blocks do not touch used seeds, and check_free refuses a spent seed.
Bookkeeping only; no simulation runs."""
from __future__ import annotations

import copy
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_seed_registry as reg  # noqa: E402


def test_table_is_well_formed():
    assert reg.validate_entries() == []


def test_only_known_overlaps_exist():
    assert reg.unexpected_overlaps() == []
    kinds = {(o["kind"], tuple(o["range"])) for o in reg.find_overlaps()}
    assert ("noise", (1300001, 1300002)) in kinds       # 169 vs 176, recorded as a flaw
    assert ("scenario", (20261031, 20261031)) in kinds   # shared 150 training draw, by design


def test_new_accidental_overlap_is_caught():
    extra = reg._entry("999", "accidental", "dev_used", [(20261300, 20261305)], [(1300000, 1300003)])
    bad = reg.unexpected_overlaps(reg.ENTRIES + [extra])
    assert bad and any("999" in o["scripts"] for o in bad)


def test_reserved_blocks_touch_no_used_seed():
    reserved = [e for e in reg.ENTRIES if e["status"] == "reserved"]
    used = [e for e in reg.ENTRIES if e["status"] != "reserved"]
    assert {e["purpose"] for e in reserved} == {"stage2_dev", "stage3_sealed"}
    for r in reserved:
        for u in used:
            for field in ("scenario_seeds", "noise_keys"):
                for a in r[field]:
                    for b in u[field]:
                        assert reg._overlap(a, b) is None, (r["label"], u["label"], a, b)


def test_check_free_refuses_spent_and_dev_seeds():
    assert reg.check_free(scenario=(20261122, 20261122))          # 167 sealed draw
    assert reg.check_free(scenario=(20261121, 20261121))          # shared tuning draw
    assert reg.check_free(scenario=(20261305, 20261305))          # inside the 169 block
    assert reg.check_free(noise=(1122001, 1122004))
    assert reg.check_free(scenario=(20262100, 20262100))          # 176 dev
    assert reg.check_free(scenario=(20262200, 20262200))          # 177 dev
    # Stage 2 dev seeds already spent (190, 192, 193, 194, 198): refused with or without the stage2_dev purpose
    for scen, noise in (((20261800, 20261800), (1800001, 1800016)), ((20261801, 20261802), (1810000, 1810000)),
                        ((20261803, 20261803), (1820000, 1820000)), ((20261804, 20261805), (1830000, 1830000)),
                        ((20261806, 20261806), (1840000, 1840000)), ((20261810, 20261811), (1850000, 1850000)),
                        ((20261820, 20261821), (1890000, 1890000)), ((20261812, 20261812), (1860000, 1860000)),
                        ((20261822, 20261822), (1890100, 1890100))):
        assert reg.check_free(scenario=scen)
        assert reg.check_free(noise=noise)
        assert reg.check_free(scenario=scen, noise=noise, purpose="stage2_dev")
    assert reg.check_free(scenario=(20261807, 20261807))                       # 199: open loop, no noise key
    assert reg.check_free(scenario=(20261807, 20261807), purpose="stage2_dev")


def test_check_free_accepts_reserved_block_only_for_its_purpose():
    reserved = next(e for e in reg.ENTRIES if e["status"] == "reserved" and e["purpose"] == "stage2_dev")
    scen = tuple(reserved["scenario_seeds"][0])
    noise = tuple(reserved["noise_keys"][0])
    assert reg.check_free(scenario=scen, noise=noise, purpose="stage2_dev") == []
    probe = (scen[0], min(scen[0] + 1, scen[1]))
    assert reg.check_free(scenario=probe, noise=(1895000, 1895000),
                          purpose="stage2_dev") == []                  # a free dev block no configuration uses
    # the dev block cannot be spent as a sealed block, and vice versa
    assert reg.check_free(scenario=scen, purpose="stage3_sealed")
    assert reg.check_free(scenario=(20262300, 20262399), purpose="stage2_dev")
    assert reg.check_free(scenario=(20262300, 20262399), noise=(2300000, 2300999),
                          purpose="stage3_sealed") == []
    # without a purpose a reserved block is not free
    assert reg.check_free(scenario=(20261905, 20262004))
    # a request that spills out of the block is refused
    assert reg.check_free(scenario=(20262050, 20262150), purpose="stage2_dev")
    assert reg.check_free(scenario=(20261900, 20261800))          # low above high


def test_remainder_helper_and_reserved_entry_follow_the_used_entries():
    assert reg._remainder((1, 10), [(3, 4), (8, 20), (-5, 0)]) == [(1, 2), (5, 7)]
    assert reg._remainder((1, 10), []) == [(1, 10)] and reg._remainder((1, 10), [(1, 10)]) == []
    reserved = [e for e in reg.ENTRIES if e["status"] == "reserved" and e["purpose"] == "stage2_dev"][0]
    used = [e for e in reg.ENTRIES if e["status"] != "reserved"]
    again = reg._reserved_stage2_dev(used)
    assert again["scenario_seeds"] == reserved["scenario_seeds"] and again["noise_keys"] == reserved["noise_keys"]
    lo, hi = reg.STAGE2_DEV_SCENARIO_BLOCK
    assert reserved["scenario_seeds"][0][0] > 20261806 and reserved["scenario_seeds"][-1][1] == hi
    assert reserved["noise_keys"][0] == [1800000, 1800000]            # the one free key before the 190 block


def test_unused_gap_seeds_are_free_without_purpose():
    assert reg.check_free(scenario=(20261600, 20261699)) == []


def test_evidence_text_is_found_in_saved_outputs():
    sample = REPO / "output" / "trackc_sealed_containment" / "frozen_protocol.json"
    if not sample.exists():
        pytest.skip("output/ not present (scratch copy)")
    failures, skipped = reg.verify_evidence()
    assert failures == [], failures
    assert len(skipped) < len([1 for e in reg.ENTRIES for _ in e["evidence"]])


def test_evidence_check_detects_a_wrong_citation(tmp_path):
    good = tmp_path / "out.json"
    good.write_text('{"seed": 20261121}', encoding="utf-8")
    e = reg._entry("x", "x", "dev_used", [(1, 2)], [], [("out.json", "20261121")])
    assert reg.verify_evidence(tmp_path, [e]) == ([], [])
    bad = copy.deepcopy(e)
    bad["evidence"] = [["out.json", "99999999"]]
    failures, skipped = reg.verify_evidence(tmp_path, [bad])
    assert failures and not skipped
    gone = copy.deepcopy(e)
    gone["evidence"] = [["missing.json", "1"]]
    assert reg.verify_evidence(tmp_path, [gone]) == ([], ["missing.json"])


def test_registry_document_is_stable_and_hashed():
    a, b = reg.registry_document(), reg.registry_document()
    assert a == b and len(a["sha256"]) == 64
    changed = copy.deepcopy(reg.ENTRIES)
    changed[0]["scenario_seeds"][0][0] += 1
    assert reg.registry_document(changed)["sha256"] != a["sha256"]
