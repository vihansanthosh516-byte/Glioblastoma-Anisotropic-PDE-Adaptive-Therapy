import csv
import json

import pytest

from src.track_b_acrin6677_advanced_mri_linkage import audit_linkage, write_result


S0_HEADER = ["cn", "S0E1", "S0E2D", "S0E5", "S0E13", "S0E23"]
SC_HEADER = ["cn", "SCE1", "SCE2D", "SCE3DEID", "SCE4", "SCE20", "SCE28", "SCE36"]


def _write(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)
    return path


def _three_s0(patient, *, second_code=2):
    return [
        [patient, 1, 0, 2, 1, ""],  # no lesion; blank area still means valid visit
        [patient, second_code, 14, 2, 2, "AREA_VALUE_MUST_NOT_ESCAPE"],
        [patient, 3, 56, 2, 2, "AREA_VALUE_MUST_NOT_ESCAPE"],
    ]


def test_exact_same_day_and_visit_with_adequate_srcbv_is_counted_without_values(tmp_path):
    s0 = _write(tmp_path / "S0.csv", S0_HEADER, _three_s0("P1"))
    sc = _write(
        tmp_path / "SC.csv",
        SC_HEADER,
        [
            ["P1", 2, 14, "R1", 1, "RCBV_SECRET_1", "", ""],
            ["P1", 2, 14, "R2", 1, "", "RCBV_SECRET_2", ""],
        ],
    )

    result = audit_linkage([s0], [sc])
    encoded = json.dumps(result)

    assert result["patients_with_three_valid_s0_visits"] == 1
    assert result["patients_with_exact_sc_key_day_visit"] == 1
    assert result["patients_with_exact_adequate_sc"] == 1
    assert result["eligible_unique_patients_with_srcbv_presence"] == 1
    assert result["eligible_reader_count_distribution"] == {"2": 1}
    assert result["eligible_patient_slot_presence"] == {"sce20": 1, "sce28": 1, "sce36": 0}
    assert result["status"] == "stop_below_frozen_minimum"
    assert result["s0_area_values_parsed"] is False
    assert result["srcbv_values_parsed"] is False
    assert "P1" not in encoded
    assert "AREA_VALUE_MUST_NOT_ESCAPE" not in encoded
    assert "RCBV_SECRET" not in encoded


def test_link_requires_exact_day_and_timepoint_and_adequate_quality(tmp_path):
    s0_rows = _three_s0("P2") + _three_s0("P3") + _three_s0("P4")
    s0 = _write(tmp_path / "S0.csv", S0_HEADER, s0_rows)
    sc = _write(
        tmp_path / "SC.csv",
        SC_HEADER,
        [
            ["P2", 2, 15, "R1", 1, "1.4", "", ""],  # wrong day
            ["P3", 3, 14, "R1", 1, "1.4", "", ""],  # wrong timepoint
            ["P4", 2, 14, "R1", 2, "1.4", "", ""],  # suboptimal
        ],
    )

    result = audit_linkage([s0], [sc])

    assert result["patients_with_three_valid_s0_visits"] == 3
    assert result["patients_with_exact_sc_key_day_visit"] == 1
    assert result["patients_with_exact_adequate_sc"] == 0
    assert result["eligible_unique_patients_with_srcbv_presence"] == 0


def test_uninterpretable_s0_row_and_two_visit_patient_do_not_qualify(tmp_path):
    rows = _three_s0("P1")
    rows[1][3] = 1  # uninterpretable second visit
    rows += [["P2", 1, 0, 2, 1, ""], ["P2", 2, 14, 2, 2, "2.0"]]
    s0 = _write(tmp_path / "S0.csv", S0_HEADER, rows)
    sc = _write(tmp_path / "SC.csv", SC_HEADER, [["P1", 2, 14, "R1", 1, "1.0", "", ""]])

    result = audit_linkage([s0], [sc])

    assert result["s0_uninterpretable_rows_excluded"] == 1
    assert result["patients_with_three_valid_s0_visits"] == 0
    assert result["eligible_unique_patients_with_srcbv_presence"] == 0


def test_conflicting_visit_codes_on_second_day_prevent_link(tmp_path):
    rows = _three_s0("P1")
    rows.append(["P1", 3, 14, 2, 2, "AREA_VALUE_MUST_NOT_ESCAPE"])
    s0 = _write(tmp_path / "S0.csv", S0_HEADER, rows)
    sc = _write(tmp_path / "SC.csv", SC_HEADER, [["P1", 2, 14, "R1", 1, "1.0", "", ""]])

    result = audit_linkage([s0], [sc])

    assert result["patients_with_three_valid_s0_visits"] == 1
    assert result["patients_with_ambiguous_second_visit_code"] == 1
    assert result["patients_with_exact_sc_key_day_visit"] == 0
    assert result["eligible_unique_patients_with_srcbv_presence"] == 0


def test_non_dictionary_codes_fail_closed(tmp_path):
    s0 = _write(tmp_path / "S0.csv", S0_HEADER, [["P1", 99, 0, 2, 1, ""]])
    sc = _write(tmp_path / "SC.csv", SC_HEADER, [["P1", 2, 14, "R1", 1, "1.0", "", ""]])

    with pytest.raises(ValueError, match="unknown code in s0e1"):
        audit_linkage([s0], [sc])


def test_missing_schema_column_fails_closed(tmp_path):
    s0 = _write(tmp_path / "S0.csv", ["cn", "s0e1"], [["P1", 1]])
    sc = _write(tmp_path / "SC.csv", SC_HEADER, [])

    with pytest.raises(ValueError, match="required columns missing"):
        audit_linkage([s0], [sc])


def test_result_writer_refuses_overwrite(tmp_path):
    out = tmp_path / "audit.json"
    write_result({"count": 1}, out)
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        write_result({"count": 2}, out)
    assert json.loads(out.read_text(encoding="utf-8")) == {"count": 1}
