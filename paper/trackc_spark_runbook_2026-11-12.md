# Track C v2: DGX Spark runbook (2026-11-12)

What to do on the DGX Spark, in order. Everything before step 1 is done on the PC (`paper/trackc_amendment_unblock_2026-11-12.md`). The Track C code is CPU-only NumPy and SciPy, so the Spark's GPU is not used. Its 20 Arm CPU cores do the work.

Do not look at any policy result between steps unless a step says so. Do not change any setting after a run starts.

## 0. Connect and set up (about 30 minutes)

```bash
git clone <this repository> trackc && cd trackc
python3.14 -m venv venv
venv/bin/pip install -r requirements-trackc.lock
```

If Python 3.14 is not on the Spark, use the newest Python 3 available and write the version into the run log. Use `tmux` so a lost connection does not stop a run.

## 1. Clean-install test (closes C-14)

```bash
PYTHONPATH=src venv/bin/python -m pytest -q -p no:cacheprovider tests/test_trackc_*.py
```

Expected (PC result in the commit message of the pre-Spark fix commit): all pass except two skips and one known strict xfail (the `RobustBeliefMPC` default stop rule; the driver passes the stop rule explicitly). Any other failure: stop and fix before step 2. Write the count into `paper/trackc_spark_log_<date>.md`.

## 2. Timing at two cohort sizes (closes the time model for D3)

```bash
venv/bin/python -I src/trackc_research_driver.py --phase dev_pilot --cell primary \
    --patients 2 --n-worlds 3 --workers 6 --confirm-development
```

This also gives `elapsed_seconds` per world at m = 2. These worlds are development data only.

## 3. Development pilot (decision D3)

```bash
venv/bin/python -I src/trackc_research_driver.py --phase dev_pilot --cell primary \
    --patients 12 --n-worlds 20 --workers 18 --max-hours 5 --confirm-development
```

Repeat the same command until it prints `complete: 40 of 40 worlds`. Each repeat resumes and skips finished worlds. Then:

```bash
venv/bin/python -I src/trackc_dev_pilot_design.py \
    --pilot-dir output/trackc_v2_research/dev_pilot/primary/m12 \
    --timing-dir output/trackc_v2_research/dev_pilot/primary/m02
```

It writes `design.json` with the chosen cohort size `choice_m` and world count `choice_n_per_family`. This is a design choice from development data, not a result.

## 4. Freeze before evaluation (on the PC or the Spark)

Write a dated amendment with `choice_m`, `choice_n_per_family`, total hours, and a runtime cap. Close the four open config gates (cohort size, world count, runtime cap, clean install). Commit and push. The evaluation phase stays refused until the gate list is empty.

**After this commit, change nothing.** The driver stores the code and config hashes in every world file and refuses to mix files from different code.

## 5. Evaluation run (only after step 4)

```bash
venv/bin/python -I src/trackc_research_driver.py --phase evaluation --cell primary \
    --patients <choice_m> --n-worlds <choice_n_per_family> --workers 18 --max-hours 5
```

Repeat until complete. Then run each sensitivity cell with the same command and `--cell scan28`, `scan56`, `scan112`, `miss10`, `penalty0`, `penalty10`.

## 6. Bring results back

Copy `output/trackc_v2_research/` back to the PC. Run `--verify-only` with the same arguments for every cell; each must say `complete`. Then run the frozen analysis once per cell:

```bash
venv/bin/python -I src/trackc_v2_analysis.py --cell primary --patients <choice_m> --n-worlds <choice_n_per_family>
```

Only the primary cell can be "positive". Report every cell, whatever the result.

## Size check

Rough cost (estimate; replace with step 2 and 3 numbers): see `paper/trackc_amendment_unblock_2026-11-12.md`, section "Cost".
