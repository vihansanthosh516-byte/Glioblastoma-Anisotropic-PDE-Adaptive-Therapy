from pathlib import Path

from src.mu_glioma_loader import find_timepoints, list_patient_ids, load_patient_record


def test_find_timepoints_and_patient_ids(tmp_path: Path):
    root = tmp_path / "MU-Glioma-Post"
    patient = root / "PatientID_0001"
    tp2 = patient / "Timepoint_2"
    tp1 = patient / "Timepoint_1"
    tp2.mkdir(parents=True)
    tp1.mkdir()
    (tp2 / "x_tumorMask.nii.gz").touch()
    (tp1 / "x_tumorMask.nii.gz").touch()
    (patient / "Timepoint_3").mkdir()

    assert list_patient_ids(root) == ["PatientID_0001"]
    found = find_timepoints(patient)
    assert [item[0] for item in found] == [1, 2]
    record = load_patient_record("PatientID_0001", root, tmp_path / "missing.xlsx", include_volumes=False)
    assert record.n_timepoints == 2
    assert record.has_longitudinal_pair
    assert [item.number for item in record.timepoints] == [1, 2]
