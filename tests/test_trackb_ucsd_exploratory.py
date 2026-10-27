from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "167_trackb_ucsd_exploratory.py"
SPEC = importlib.util.spec_from_file_location("trackb_ucsd_exploratory", MODULE_PATH)
EXPLORATORY = importlib.util.module_from_spec(SPEC)
assert SPEC is not None and SPEC.loader is not None
sys.modules[SPEC.name] = EXPLORATORY
SPEC.loader.exec_module(EXPLORATORY)


class ExploratoryRunnerTests(unittest.TestCase):
    def test_origin_quality_screen_uses_only_first_scan_images(self) -> None:
        good = SimpleNamespace(patient_id="private-good", first_files=object(), second_files=object())
        bad = SimpleNamespace(patient_id="private-bad", first_files=object(), second_files=object())
        loaded = []

        def load_origin(scan):
            loaded.append(scan)
            return {"fixture": scan}

        def quality(images):
            return ["t2w"] if images["fixture"] is bad.first_files else []

        with (
            patch.object(EXPLORATORY.RUNNER, "_load_ucsd_origin_images", side_effect=load_origin),
            patch.object(EXPLORATORY.RUNNER, "_origin_intensity_quality_issues", side_effect=quality),
        ):
            eligible, modality_counts, failures = EXPLORATORY._screen_origin_pairs([good, bad])

        self.assertEqual(eligible, [(1, good)])
        self.assertEqual(loaded, [good.first_files, bad.first_files])
        self.assertEqual(modality_counts, {"t2w": 1})
        self.assertEqual(failures, [{"pair_order": 2, "modalities": ["t2w"]}])
        self.assertNotIn("private-bad", repr(failures))

    def test_patient_record_reports_voxel_counts_and_auc_without_patient_id(self) -> None:
        sample = EXPLORATORY.RUNNER.PatientSample(
            patient_key="private-patient",
            x=np.zeros((3, len(EXPLORATORY.RUNNER.FEATURES.ORIGIN_FEATURE_NAMES)), dtype=np.float32),
            prior=np.asarray([0.0, 0.3, 0.9]),
            instability=np.zeros(3, dtype=np.float32),
            target=np.asarray([1, 0, 0], dtype=np.uint8),
            risk_voxels=3,
            additions_total=2,
            additions_outside_risk=1,
            additions_outside_current_support=1,
            interval_days=120.0,
        )
        scores = {
            "full_hybrid": np.asarray([1.0, 0.5, 0.0]),
            "pde_only": np.asarray([0.9, 0.3, 0.2]),
        }

        row = EXPLORATORY._patient_record(7, sample, scores, "pde_only")

        self.assertEqual(row["positive_target_voxels_in_risk_set"], 1)
        self.assertEqual(row["negative_target_voxels_in_risk_set"], 2)
        self.assertEqual(row["scan_gap_stratum"], "91_180_days")
        self.assertEqual(row["paired_auc_difference"], 0.0)
        self.assertNotIn("patient_key", row)
        self.assertNotIn("private-patient", repr(row))


if __name__ == "__main__":
    unittest.main()
