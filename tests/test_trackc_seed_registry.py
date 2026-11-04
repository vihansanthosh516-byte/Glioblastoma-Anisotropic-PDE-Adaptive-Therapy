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


def test_check_free_accepts_reserved_block_only_for_its_purpose():
    assert reg.check_free(scenario=(20261800, 20261899), noise=(1800000, 1800999),
                          purpose="stage2_dev") == []
    # the dev block cannot be spent as a sealed block, and vice versa
    assert reg.check_free(scenario=(20261800, 20261899), purpose="stage3_sealed")
    assert reg.check_free(scenario=(20262300, 20262399), purpose="stage2_dev")
    assert reg.check_free(scenario=(20262300, 20262399), noise=(2300000, 2300999),
                          purpose="stage3_sealed") == []
    # without a purpose a reserved block is not free
    assert reg.check_free(scenario=(20261800, 20261899))
    # a request that spills out of the block is refused
    assert reg.check_free(scenario=(20262050, 20262150), purpose="stage2_dev")
    assert reg.check_free(scenario=(20261900, 20261800))          # low above high


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
