"""Presence-and-timing audit for ACRIN 6677 advanced MRI linkage.

This module never parses or returns area or rCBV magnitudes. It checks only
field presence, patient keys, visit codes, scan days, reader keys, and QC codes.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence


MINIMUM_PATIENTS = 25
S0_COLUMNS = ("cn", "s0e1", "s0e2d", "s0e5", "s0e13", "s0e23")
SC_COLUMNS = ("cn", "sce1", "sce2d", "sce3deid", "sce4", "sce20", "sce28", "sce36")
SRCBV_COLUMNS = ("sce20", "sce28", "sce36")
MISSING_TOKENS = frozenset({"", "na", "n/a", "nan", "null", "none"})
TIMEPOINT_CODES = frozenset(range(1, 17))


def _text(value: object) -> str:
    return "" if value is None else str(value).strip()


def _missing(value: object) -> bool:
    return _text(value).casefold() in MISSING_TOKENS


def _integer(value: object, *, field: str, allowed: frozenset[int] | None = None) -> int | None:
    if _missing(value):
        return None
    try:
        number = Decimal(_text(value))
    except InvalidOperation as exc:
        raise ValueError(f"invalid integer in {field}") from exc
    if not number.is_finite() or number != number.to_integral_value():
        raise ValueError(f"invalid integer in {field}")
    result = int(number)
    if allowed is not None and result not in allowed:
        raise ValueError(f"unknown code in {field}")
    return result


def _selected_rows(path: Path, required: Sequence[str]) -> Iterator[dict[str, str]]:
    """Read required columns only; never retain other fields from a row."""
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"empty CSV: {path.name}") from exc
        names = [name.strip().casefold() for name in header]
        if len(names) != len(set(names)):
            raise ValueError(f"duplicate normalized header: {path.name}")
        indexes = {name: index for index, name in enumerate(names)}
        missing = [name for name in required if name.casefold() not in indexes]
        if missing:
            raise ValueError(f"required columns missing in {path.name}: {', '.join(missing)}")
        selected_indexes = {name: indexes[name.casefold()] for name in required}
        for row in reader:
            if not row or not any(_text(value) for value in row):
                continue
            if len(row) != len(header):
                raise ValueError(f"malformed CSV row in {path.name}")
            yield {name: row[index] for name, index in selected_indexes.items()}


def audit_linkage(
    s0_paths: Sequence[str | Path],
    sc_paths: Sequence[str | Path],
    *,
    minimum_patients: int = MINIMUM_PATIENTS,
) -> dict[str, object]:
    """Return aggregate linkage counts without exposing row values or IDs."""
    if minimum_patients != MINIMUM_PATIENTS:
        raise ValueError("the frozen minimum is 25 patients")
    if not s0_paths or not sc_paths:
        raise ValueError("at least one S0 and one SC file are required")

    # One entry per patient and distinct day. Only validity/presence flags and
    # timepoint codes are retained; area text is never parsed or returned.
    visits: dict[str, dict[int, dict[str, object]]] = defaultdict(dict)
    s0_rows = 0
    s0_area_presence_rows = 0
    s0_no_lesion_rows = 0
    s0_uninterpretable_rows = 0
    s0_missing_qc_rows = 0
    s0_missing_key_time_rows = 0

    for raw_path in s0_paths:
        path = Path(raw_path)
        for row in _selected_rows(path, S0_COLUMNS):
            s0_rows += 1
            patient = _text(row["cn"])
            day = _integer(row["s0e2d"], field="s0e2d")
            visit = _integer(row["s0e1"], field="s0e1", allowed=TIMEPOINT_CODES)
            if not patient or day is None or visit is None:
                s0_missing_key_time_rows += 1
                continue

            quality = _integer(row["s0e5"], field="s0e5", allowed=frozenset({1, 2}))
            if quality is None:
                s0_missing_qc_rows += 1
                continue
            if quality == 1:
                s0_uninterpretable_rows += 1
                continue

            lesion_code = _integer(row["s0e13"], field="s0e13", allowed=frozenset({1, 2}))
            area_present = not _missing(row["s0e23"])
            no_lesion = lesion_code == 1
            if area_present:
                s0_area_presence_rows += 1
            if no_lesion:
                s0_no_lesion_rows += 1
            if not area_present and not no_lesion:
                continue

            visit_record = visits[patient].setdefault(
                day, {"codes": set(), "reader_rows": 0, "area_present": False, "no_lesion": False}
            )
            codes = visit_record["codes"]
            assert isinstance(codes, set)
            codes.add(visit)
            visit_record["reader_rows"] = int(visit_record["reader_rows"]) + 1
            visit_record["area_present"] = bool(visit_record["area_present"]) or area_present
            visit_record["no_lesion"] = bool(visit_record["no_lesion"]) or no_lesion

    second_visits: dict[str, tuple[int, int | None]] = {}
    for patient, patient_visits in visits.items():
        ordered_days = sorted(patient_visits)
        if len(ordered_days) < 3:
            continue
        first_three = ordered_days[:3]
        second_day = first_three[1]
        labels = patient_visits[second_day]["codes"]
        assert isinstance(labels, set)
        second_visits[patient] = (second_day, next(iter(labels)) if len(labels) == 1 else None)

    # The map contains only record keys, reader IDs, and Boolean field-presence
    # flags. Numeric sRCBV text is never converted, logged, or written.
    sc_records: dict[tuple[str, int, int], dict[str, object]] = defaultdict(
        lambda: {"any_rows": 0, "adequate_rows": 0, "readers": defaultdict(set), "slots": Counter()}
    )
    sc_rows = 0
    sc_missing_key_time_rows = 0
    sc_missing_reader_rows = 0
    sc_missing_qc_rows = 0
    sc_nonadequate_rows = 0

    for raw_path in sc_paths:
        path = Path(raw_path)
        for row in _selected_rows(path, SC_COLUMNS):
            sc_rows += 1
            patient = _text(row["cn"])
            day = _integer(row["sce2d"], field="sce2d")
            visit = _integer(row["sce1"], field="sce1", allowed=TIMEPOINT_CODES)
            if not patient or day is None or visit is None:
                sc_missing_key_time_rows += 1
                continue

            key = (patient, visit, day)
            record = sc_records[key]
            record["any_rows"] = int(record["any_rows"]) + 1
            quality = _integer(row["sce4"], field="sce4", allowed=frozenset({1, 2, 3}))
            if quality is None:
                sc_missing_qc_rows += 1
                continue
            if quality != 1:
                sc_nonadequate_rows += 1
                continue
            record["adequate_rows"] = int(record["adequate_rows"]) + 1
            reader = _text(row["sce3deid"])
            if not reader:
                sc_missing_reader_rows += 1
                continue

            present_slots = {name for name in SRCBV_COLUMNS if not _missing(row[name])}
            slot_counts = record["slots"]
            assert isinstance(slot_counts, Counter)
            for name in present_slots:
                slot_counts[name] += 1
            if present_slots:
                readers_by_id = record["readers"]
                assert isinstance(readers_by_id, defaultdict)
                readers_by_id[reader].update(present_slots)

    exact_key_ids: set[str] = set()
    exact_adequate_ids: set[str] = set()
    eligible_ids: set[str] = set()
    reader_distribution: Counter[int] = Counter()
    slot_presence_patient_counts: Counter[str] = Counter()
    ambiguous_second_visit_labels = 0

    for patient, (day, visit) in second_visits.items():
        if visit is None:
            ambiguous_second_visit_labels += 1
            continue
        key = (patient, visit, day)
        record = sc_records.get(key)
        if record is None:
            continue
        exact_key_ids.add(patient)
        adequate_rows = int(record["adequate_rows"])
        if adequate_rows:
            exact_adequate_ids.add(patient)
        readers_by_id = record["readers"]
        assert isinstance(readers_by_id, defaultdict)
        eligible_readers = {reader: slots for reader, slots in readers_by_id.items() if slots}
        if eligible_readers:
            eligible_ids.add(patient)
            reader_distribution[len(eligible_readers)] += 1
            all_slots = set().union(*eligible_readers.values())
            for slot in all_slots:
                slot_presence_patient_counts[slot] += 1

    cohort_hash = hashlib.sha256("\n".join(sorted(eligible_ids)).encode("utf-8")).hexdigest()
    eligible_n = len(eligible_ids)
    return {
        "record_type": "metadata_only_advanced_mri_linkage_audit",
        "protocol": "paper/trackB_acrin6677_advanced_mri_protocol_2026-11-11.md",
        "minimum_patients": MINIMUM_PATIENTS,
        "status": "eligible_for_next_step_not_scored" if eligible_n >= MINIMUM_PATIENTS else "stop_below_frozen_minimum",
        "s0_rows": s0_rows,
        "s0_area_field_present_rows": s0_area_presence_rows,
        "s0_no_lesion_rows": s0_no_lesion_rows,
        "s0_uninterpretable_rows_excluded": s0_uninterpretable_rows,
        "s0_missing_qc_rows_excluded": s0_missing_qc_rows,
        "s0_missing_key_or_time_rows_excluded": s0_missing_key_time_rows,
        "patients_with_three_valid_s0_visits": len(second_visits),
        "patients_with_ambiguous_second_visit_code": ambiguous_second_visit_labels,
        "sc_rows": sc_rows,
        "sc_missing_key_or_time_rows_excluded": sc_missing_key_time_rows,
        "sc_missing_qc_rows_excluded": sc_missing_qc_rows,
        "sc_nonadequate_rows_excluded": sc_nonadequate_rows,
        "sc_missing_reader_rows_excluded": sc_missing_reader_rows,
        "patients_with_exact_sc_key_day_visit": len(exact_key_ids),
        "patients_with_exact_adequate_sc": len(exact_adequate_ids),
        "eligible_unique_patients_with_srcbv_presence": eligible_n,
        "eligible_reader_count_distribution": {str(key): reader_distribution[key] for key in sorted(reader_distribution)},
        "eligible_patient_slot_presence": {name: slot_presence_patient_counts[name] for name in SRCBV_COLUMNS},
        "eligible_patient_key_sha256": cohort_hash,
        "s0_area_values_parsed": False,
        "srcbv_values_parsed": False,
        "patient_level_rows_saved": False,
        "images_opened": False,
    }


def write_result(result: Mapping[str, object], output_path: str | Path) -> None:
    path = Path(output_path)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path.name}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _main() -> int:
    parser = argparse.ArgumentParser(description="Audit ACRIN 6677 SC/S0 linkage without reading measurement magnitudes.")
    parser.add_argument("--s0", action="append", required=True, help="S0.csv path; repeat for each split")
    parser.add_argument("--sc", action="append", required=True, help="SC.csv path; repeat for each split")
    parser.add_argument("--output", required=True, help="new JSON output path; existing files are refused")
    args = parser.parse_args()
    result = audit_linkage(args.s0, args.sc)
    write_result(result, args.output)
    print(json.dumps({key: value for key, value in result.items() if key not in {"eligible_patient_key_sha256"}}, sort_keys=True))
    return 0 if result["status"] == "eligible_for_next_step_not_scored" else 2


if __name__ == "__main__":
    raise SystemExit(_main())
