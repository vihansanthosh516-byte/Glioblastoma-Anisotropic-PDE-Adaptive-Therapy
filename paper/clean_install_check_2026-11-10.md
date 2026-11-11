# Clean-install check (2026-11-10, revised after Codex review)

Windows 11, Python 3.14.4. Fresh virtual environments in the session scratch folder. The package was always built from a scratch copy of the tracked files, so no build files were written into the repository. An `egg-info` folder that an earlier dry run (run from the repository) wrote into `src/` was removed; its dependency list lacked `openpyxl`, which dates it to that dry run.

## A. Package-install verification (installed code)

1. First real `pip install <scratch copy>` (environment `cleanenv2`): exit 0, but only 9 files from three subpackages (`hil`, `neural_pde`, `rl`) were installed, and `run-pipeline --help` failed with `ModuleNotFoundError: No module named 'src'`. Cause: `[tool.setuptools.packages.find] where = ["src"]` packaged the subfolders of `src`, not the `src` package that the entry point `src.cli:main` and the code import.
2. Fixes: `pyproject.toml` now uses `where = ["."]`, `include = ["src", "src.*"]`. `src/cli.py` now handles `--help`/`-h` (prints usage, exits 0, starts no analysis); before, `--help` was rejected as an unknown script key (exit 1, no job started). Earlier fixes stand: build backend `setuptools.build_meta`, and `openpyxl` declared.
3. Second fresh environment (`cleanenv3`), command `python -m pip install <scratch copy>`: exit 0. `importlib.metadata` lists 363 installed `.py` files. `import src, src.cli` loads from `cleanenv3\Lib\site-packages\src`. `run-pipeline --help` exits 0 and prints usage. `--help` and `--list` start no analysis. No data job was launched.
4. Installed-code test: `tests/test_track_b_coverage_scoring.py` copied alone to a folder with no `src`, run with `PYTHONPATH` empty and `--import-mode=importlib`; the run printed `src` loaded from site-packages. **10 passed.**

Limit: an installed copy resolves `PROJECT_ROOT` to site-packages, so running numbered scripts through the installed CLI needs `GBM_PROJECT_ROOT` set to the checkout. The `--help` text says this. Many modules use bare imports (for example `import trackc_sim`) that work only with `src` on `sys.path`; they are not usable as an installed library without that path.

## B. Dependency-install tests (source tree)

Same `cleanenv3` interpreter plus `torch` and `pytest`, which the core package does not declare (torch is only in the `legacy-ml` extra, but the Track B LUMIERE runner tests need it). Resolved versions: `configs/clean_install_freeze_2026-11-10.txt` (80 packages from `pip freeze`; it includes the installed project itself).

Command, run from a scratch copy of the tracked tree with `PYTHONPATH=<copy>;<copy>\src`: pytest on 15 files (Track A synthetic helpers; Track B coverage, estimand, run guard, LUMIERE runner and protocol; Track C driver, prior support, world runner, policy benchmark, ambiguity MPC, config parity, GE-betting, GE-betting shape, seed registry). The run printed `src` loaded from the scratch source copy, **not** the installed package. **203 passed, 1 skipped, 3 expected failures (C9-C11), 121 s.**

## Limits

- Not the full test suite. Legacy-ML extras (torch-geometric, gymnasium, stable-baselines3) not installed or tested.
- Versions are not pinned in `pyproject.toml`; the freeze file records what resolved on 2026-11-10.
- No data-dependent analysis ran in either environment.
