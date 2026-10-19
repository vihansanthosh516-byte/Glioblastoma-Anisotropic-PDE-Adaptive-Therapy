> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Historical packaging changes: corrected evidence record

Corrected 2026-10-19 against source and Git history. The former application date,
2026-06-15, is author-reported and unverified. The first tracked addition is
`d0d441a` (2026-09-17). No source audit or saved execution transcript was cited
for the original account, so its “comprehensive audit” description is withdrawn.

## What the repository establishes

Commit `d0d441a` adds `src/__init__.py`, `pyproject.toml`, `src/cli.py`, and
`run_all.ps1`. It changes 61 paths under `src/`: 52 modified, seven added and two
deleted. Thus 61 is not the number of modified source files. The four named
package/runner files are tracked; the former instruction to stage them is stale.

`pyproject.toml` declares distribution `gliblastoma-pde-adaptive-therapy` (spelling
retained as existing metadata), dependencies and entry point
`run-pipeline = "src.cli:main"`. It declares Python >=3.14. This is metadata,
not a verified supported-runtime range or successful installation. The unusual
build backend and package discovery still need a separate clean build/install
check; no successful pip-install claim is made here.

The CLI currently has 56 mapping keys for 55 distinct targets, including alias
`16c` and script 53. `--list` displays historical mappings; a numeric key launches
its mapped script. Existence in this list does not establish safe execution or
validation. `--run-all` is disabled following the audit. `--help` is not
implemented. The old “53/53 succeeded” and script-53 usage transcripts lack saved
logs and are withdrawn as verification evidence.

The CLI derives its default root from `__file__` and accepts `GBM_PROJECT_ROOT`.
Four current source files mention `GBM_PROJECT_ROOT` or `GBM_DATA_DIR`:
`src/cli.py`, `src/54_cgat_build_graph.py`, `src/55_cgat_train.py`, and
`src/build_full_graph.py`. This does not establish 40+ configurable scripts.
Script 01 still uses a fixed data directory and relative output; script 02 uses
relative output. Neither docstring documents the claimed environment override.
There is no current Windows-versus-WSL hardcoded default in the CLI.

`run_all.ps1` is a historical orchestrator with a Month10 branch selecting script
45. Audit corrections repair the missing `$src` sigil, malformed condition and
root resolution. Only PowerShell parsing is checked; the scientific chain is not
run. “Fully functional,” portability and easier onboarding are not established.

## Dated inventory and remaining evidence limits

As of 2026-10-19, `CHANGELOG.md`, `CONTRIBUTING.md`, `API.md`, `usage_guide.md` and
`.github/workflows/` are absent. Source already contains type annotations (for
example `src/cli.py`). `Makefile` has a test target; the contrary claim was false.
Ten root `test_*.py` files exist; they are not claimed git-ignored without an
exact per-file ignore result. These are separate from the supported `tests/`
inventory in `../CLAUDE.md`.

Seven source files set `OMP_NUM_THREADS` to `2` (04–08 and 10–11); scripts 01–02
set it to `22`. Potential underutilization is an inference, not a measured effect.
`docs/VALIDATION_SUMMARY.md` in this archive reports four historical P0 gaps
(wave speed 31.6%, necrotic fraction 6.7%, TI_max=1.0, healthy-zone C=1.0).
Those rows describe that archived report, not current calibration results.

Current goal and commands: `../START_HERE.md` and `../CLAUDE.md`. No pipeline,
installation or scientific analysis was run to produce these corrections.
