# Clean-install check (2026-11-10)

Fresh virtual environment in the session scratch folder, Python 3.14.4 (`py -3.14 -m venv`), Windows 11.

## Steps and results

1. `pip install . --dry-run --no-deps` **failed**: `pyproject.toml` named a nonexistent build backend (`setuptools.backends._legacy:_Backend`). Fixed to `setuptools.build_meta`; the dry run then reports "Would install gliblastoma-pde-adaptive-therapy-0.1.0".
2. Installed the `pyproject.toml` dependencies, the `Requirements.txt` additions (SALib, dipy, h5py), `torch`, `openpyxl`, and `pytest` from PyPI with no version pins. Exit 0. Exact resolved versions: `configs/clean_install_freeze_2026-11-10.txt` (93 packages; for example numpy 2.5.4, scipy 1.18.1, scikit-learn 1.9.1, torch 2.15.0, pandas 2.3.3).
3. `openpyxl` is required by the Track A loader (`pandas.read_excel` in `src/189_trackA_yu2020_lodo.py`) but was not declared. Added to `pyproject.toml`. `pytest` is still only a test tool and is not declared.
4. Tests in a scratch copy of the tracked tree, run with the clean environment's interpreter: 15 test files covering the Track A synthetic helpers, Track B guard/estimand/coverage/LUMIERE runner and protocol, Track C prior, driver, runner, benchmark, ambiguity MPC, config parity, GE-betting, and seed registry. **199 passed, 1 skipped, 3 expected failures (C9-C11), 123 s.**

## Limits

- Not the full test suite; legacy-ML extras (torch-geometric, gymnasium, stable-baselines3) were not installed or tested.
- Versions are not pinned in `pyproject.toml`; the freeze file records what resolved on 2026-11-10 only.
- No data-dependent analysis was run in the clean environment.
