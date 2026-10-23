"""Input-safety tests for the exploratory donor-split runner."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
RUNNER_PATH = REPO / "src" / "20_patient_split_methods.py"
SPEC = importlib.util.spec_from_file_location("patient_split_runner", RUNNER_PATH)
assert SPEC is not None and SPEC.loader is not None
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)


def write_inputs(folder: Path, split_rows: list[tuple[int, str]] | None = None) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    np.save(folder / "nn_X.npy", np.arange(18, dtype=float).reshape(6, 3))
    np.save(folder / "nn_y.npy", np.array([0, 0, 1, 1, 0, 1]))
    np.save(folder / "nn_patient_ids.npy", np.array(["a", "a", "b", "b", "c", "c"]))
    if split_rows is None:
        split_rows = [(0, "train"), (1, "train"), (2, "train"), (3, "train"),
                      (4, "test"), (5, "test")]
    pd.DataFrame(split_rows, columns=["index", "split"]).to_csv(
        folder / "nn_train_test_split_patient.tsv", sep="\t", index=False
    )


def test_loading_valid_patient_disjoint_split(tmp_path: Path) -> None:
    write_inputs(tmp_path)

    X, y, patient_ids, train_rows, test_rows = runner.load_split_inputs(tmp_path)

    assert X.shape == (6, 3)
    assert len(y) == len(patient_ids) == 6
    assert set(patient_ids[train_rows]).isdisjoint(set(patient_ids[test_rows]))
    assert train_rows.tolist() == [0, 1, 2, 3]
    assert test_rows.tolist() == [4, 5]


def test_loading_rejects_patient_overlap(tmp_path: Path) -> None:
    write_inputs(tmp_path, [(0, "train"), (1, "test"), (2, "train"), (3, "train"),
                            (4, "test"), (5, "test")])

    with pytest.raises(ValueError, match="Patient overlap"):
        runner.load_split_inputs(tmp_path)


def test_loading_rejects_duplicate_or_fractional_indices(tmp_path: Path) -> None:
    write_inputs(tmp_path, [(0, "train"), (0, "test"), (2, "train"), (3, "train")])
    with pytest.raises(ValueError, match="unique"):
        runner.load_split_inputs(tmp_path)

    write_inputs(tmp_path, [(0.5, "train"), (1, "train"), (2, "train"), (4, "test")])
    with pytest.raises(ValueError, match="finite integers"):
        runner.load_split_inputs(tmp_path)


def test_run_method_writes_outputs_from_injected_arrays(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(runner, "SEEDS", [42])
    monkeypatch.setattr(runner, "N_BOOT", 5)
    result_dir = tmp_path / "results"
    result_dir.mkdir()
    Xtr = np.array([[-2.0], [-1.0], [1.0], [2.0]])
    ytr = np.array([0, 0, 1, 1])
    Xte = np.array([[-3.0], [-0.5], [0.5], [3.0]])
    yte = np.array([0, 0, 1, 1])

    def sign_classifier(X_fit, y_fit, X_eval, hp, seed):
        assert len(X_fit) == len(y_fit)
        return (X_eval[:, 0] > 0).astype(int)

    result = runner.run_method(
        "synthetic",
        [1],
        sign_classifier,
        Xtr=Xtr,
        ytr=ytr,
        Xte=Xte,
        yte=yte,
        ptr=np.array(["a", "a", "b", "b"]),
        pte=np.array(["c", "c", "d", "d"]),
        tr_i=np.array([0, 1]),
        va_i=np.array([2, 3]),
        test_rows=np.array([10, 11, 12, 13]),
        result_dir=result_dir,
    )

    assert result["macro_f1_mean"] == 1.0
    assert (result_dir / "synthetic.json").is_file()
    assert (result_dir / "pred_seed42_synthetic.npz").is_file()


def test_import_does_not_load_project_data_or_create_output(tmp_path: Path, monkeypatch) -> None:
    calls: list[str] = []

    def record_load(*args, **kwargs):
        calls.append("load")
        raise AssertionError("import must not load data")

    monkeypatch.setattr(np, "load", record_load)
    monkeypatch.setattr(pd, "read_csv", record_load)
    isolated_src = tmp_path / "src"
    isolated_src.mkdir()
    copied_runner = isolated_src / RUNNER_PATH.name
    shutil.copy2(RUNNER_PATH, copied_runner)
    spec = importlib.util.spec_from_file_location("patient_split_runner_import_check", copied_runner)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    assert calls == []
    assert not (tmp_path / "output").exists()
