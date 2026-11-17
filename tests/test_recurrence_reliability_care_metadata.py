"""Synthetic checks for the CARE sample-level pair rule (no real data)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from recurrence_reliability.care_metadata import pair_rule  # noqa: E402


def spec(sid, rt="No", alk="No", idh="wildtype"):
    return {"Sample ID": sid, "Radiation before the surgery": rt, "Alkylate agents before the surgery": alk,
            "IDH mutation status": idh, "Diagnosis": "Glioblastoma, IDH-wildtype"}


def test_untreated_primary_and_rt_tmz_first_recurrence_is_eligible():
    status, why = pair_rule({"Primary": spec("P1T1"), "1st Recurrent": spec("P1T2", "Yes", "Yes (TMZ)")})
    assert status == "sample_level_eligible" and why == []


def test_third_timepoint_does_not_block_primary_first_recurrence_pair():
    st = {"Primary": spec("P1T1"), "1st Recurrent": spec("P1T2", "Yes", "Yes (TMZ)"), "2nd Recurrent": spec("P1T3", "Yes", "Yes (TMZ)")}
    assert pair_rule(st)[0] == "sample_level_eligible"


def test_second_recurrence_or_recurrence_only_pairs_are_excluded():
    assert pair_rule({"Primary": spec("P1T1"), "2nd Recurrent": spec("P1T2", "Yes", "Yes (TMZ)")})[0] == "excluded"
    assert pair_rule({"1st Recurrent": spec("P1T1", "Yes", "Yes (TMZ)"), "2nd Recurrent": spec("P1T2")})[0] == "excluded"


def test_non_standard_treatment_or_idh_mutant_excluded():
    for rt, alk in [("No", "No"), ("Yes", "No"), ("Yes", "Yes (TMZ and CCNU)"), ("No", "Yes (TMZ)")]:
        assert pair_rule({"Primary": spec("P1T1"), "1st Recurrent": spec("P1T2", rt, alk)})[0] == "excluded"
    assert pair_rule({"Primary": spec("P1T1", idh="mutant"), "1st Recurrent": spec("P1T2", "Yes", "Yes (TMZ)")})[0] == "excluded"
    assert pair_rule({"Primary": spec("P1T1", "Yes", "Yes (TMZ)"), "1st Recurrent": spec("P1T2", "Yes", "Yes (TMZ)")})[0] == "excluded"
