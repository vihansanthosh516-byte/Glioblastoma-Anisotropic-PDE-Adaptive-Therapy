# Track C v2 + EARF: freeze amendment before evaluation (2026-11-19)

Written after the outcome-blind development pilot and BEFORE any evaluation world is run. After this commit, no code or config changes until the evaluation is complete. Everything here is **simulation**.

## Source of the choices
- Tufts HPC pilot job 7359807 (batch partition, 32 cores, 16 h 32 min, exit 0), clean checkout of commit `636a8cf`. Files received 2026-11-19 (`trackc_pilot_results_for_vihan.zip`: `design.json`, `pilot_7359807.log`, `tests_7359550.log`); kept outside git.
- `design.json` was written by `src/trackc_dev_pilot_design.py` (script sha256 `188662ed…`) from 40 development worlds (m = 12) plus 6 timing worlds (m = 2). Only elapsed time and the outcome-blind variance terms the design rule uses were read. No policy contrast was examined.

## The four closed gates (config `closed_gates_2026_11_19`)
1. **Cohort size:** m = 2 (`choice_m`).
2. **World count:** 317 evaluation worlds per family, 2 families = 634 worlds per cell (`choice_n_per_family`). Fits the reserved `trackc_v2_evaluation` seed block (checked with `build_schedule("evaluation", 317)`).
3. **Runtime cap:** `--max-hours 22` per Slurm job; jobs are resubmitted until the driver prints complete. The driver refuses to mix code or config hashes.
4. **Clean install:** Tufts clean checkout, Python 3.13.16, NumPy 2.5.5, SciPy 1.18.1 (uv venv): 585 passed, 1 xfailed, 0 failed, 0 skipped.

## Deviations, stated
- Package versions differ from `requirements-trackc.lock` (NumPy 2.4.4, SciPy 1.18.0, Python 3.14). Reason: no ready-made SciPy wheel for Python 3.14 on the first Arm machine tried; the Tufts run used the newest wheels. All tests passed. Evaluation must use the same Tufts environment as the pilot.
- The runbook expected 2 skipped tests; Tufts had 0 skipped, because the LFS outputs the skips depend on were downloaded. Not a failure.
- The runbook path for `design.json` lacked the `m12/` folder; the file is at `output/trackc_v2_research/dev_pilot/primary/m12/design.json`.
- Editing the config changes its hash. Pilot worlds (development only) carry the old hash; evaluation worlds are written to separate folders and all carry the new hash.

## Evaluation plan (unchanged from the runbook)
`sbatch --export=ALL,M=2,N=317 hpc/tufts_eval.sh`: primary, scan28, scan56, scan112, miss10, penalty0, penalty10, EARF evidence (locked integration test), EARF primary_anchor (descriptive). Only the primary cell (Track C) and EARF evidence (Integration) can be "positive", under the frozen rules in `paper/trackc_amendment_unblock_2026-11-12.md` and `paper/integration_earf_protocol_2026-11-14.md`. Every cell is reported, whatever the result.

## Cost (estimate from the pilot time model, not measured on evaluation)
Time model in `design.json`: 3,484 + 1,132 x m seconds per world, so about 1.6 core-hours per world at m = 2 on a Tufts core. About 1,012 core-hours per cell. 7 cells about 7,100 core-hours; the 2 EARF runs add an assumed 15-25% each, about 2,300-2,500. **Total about 9,400-9,600 core-hours**, about twice the 4,400 in the runbook (that figure assumed PC timing and fewer worlds). With 9 jobs x 30 workers: about 34 hours per job, so two 22-hour submissions each, if the queue starts them at once.

## Weak points
- The cost is an estimate; scan28 and the EARF arm are slower than primary.
- m = 2 virtual patients per cohort is small; the design rule chose it because more worlds buy more precision per core-hour than more patients.
- Simulation only; EARF is partly circular (see its protocol).
