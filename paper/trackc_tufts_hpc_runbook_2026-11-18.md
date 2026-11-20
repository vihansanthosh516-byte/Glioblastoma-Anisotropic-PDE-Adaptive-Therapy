# Track C + EARF: Tufts HPC runbook (2026-11-18)

How to run all Track C simulations (and the EARF integration test) on the Tufts HPC cluster. It replaces the DGX Spark plan (`paper/trackc_spark_runbook_2026-11-12.md`); the science steps are the same. Everything is **simulation**, not real patients.

Repo: https://github.com/vihansanthosh516-byte/Glioblastoma-Anisotropic-PDE-Adaptive-Therapy

## Before you start

- **Check the rules.** Tufts HPC accounts are for Tufts work. Confirm with Tufts Research Technology that running this project is allowed.
- The code is CPU-only (NumPy/SciPy). Do not request GPUs.
- Do not change any code or setting after a run starts. The driver stores code and config hashes in every result file and refuses to mix runs.
- Do not open or read policy results between steps. Only the steps below say what to look at.
- Partition name `batch`, 32 cores, 32 GB and 24 h in the scripts are **guesses**. Check with `sinfo` and the Tufts docs, then edit the `#SBATCH` lines in `hpc/*.sh` if needed (more cores = faster; the scripts use cores - 2 workers automatically).

## Step 0: set up (once, on the login node)

```bash
git lfs install
GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/vihansanthosh516-byte/Glioblastoma-Anisotropic-PDE-Adaptive-Therapy.git trackc
cd trackc
git lfs pull --include="configs/*,output/trackc*"
head -c 40 configs/trackABC_integrated_v2.json     # must start with "{"
mkdir -p logs
```

`GIT_LFS_SKIP_SMUDGE=1` skips the large files; only configs and Track C outputs (about 8 MB) are needed.

Python (if `uv` is not installed):

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
source $HOME/.local/bin/env
uv venv --python 3.13 venv
uv pip install --python venv/bin/python numpy scipy pytest
```

If Tufts has a Python 3.12+ module instead, `module load python` then `python -m venv venv && venv/bin/pip install numpy scipy pytest` also works. Each job log prints the Python, NumPy and SciPy versions; keep the logs.

## Step 1: pilot (one job)

```bash
sbatch hpc/tufts_pilot.sh
squeue -u $USER          # shows if it is waiting (PD) or running (R)
```

The job does, in order:
1. Runs the Track C tests. Expected: all pass, 2 skipped, 1 xfailed. Any failure stops the job; send the log (`logs/pilot_<jobid>.log`).
2. Timing run (2 patients, 3 worlds).
3. EARF timing run.
4. Development pilot (12 patients, 20 worlds per family = 40 worlds).
5. Writes `output/trackc_v2_research/dev_pilot/primary/m12/design.json`.

If the job ends before the pilot prints `complete: 40 of 40 worlds`, run `sbatch hpc/tufts_pilot.sh` again. It resumes and skips finished worlds.

## Step 2: STOP and send results

Send these to Vihan:
- `output/trackc_v2_research/dev_pilot/primary/m12/design.json`
- the log files in `logs/`

Vihan writes the freeze amendment (chosen cohort size `choice_m`, worlds `choice_n_per_family`, runtime cap) and pushes it. The evaluation phase is refused by the code until this is done. Then:

```bash
git pull
```

After this pull, change nothing.

## Step 3: evaluation (9 jobs at once)

Replace M and N with `choice_m` and `choice_n_per_family` from `design.json`:

```bash
sbatch --export=ALL,M=<choice_m>,N=<choice_n_per_family> hpc/tufts_eval.sh
```

This starts a job array (9 jobs that run side by side):

| Array index | Run | Role |
|---|---|---|
| 0 | primary | main result (only cell that can be "positive") |
| 1-6 | scan28, scan56, scan112, miss10, penalty0, penalty10 | sensitivity checks |
| 7 | EARF evidence | locked integration test |
| 8 | EARF primary_anchor | descriptive only |

Check progress: `squeue -u $USER` and the logs `logs/eval_<jobid>_<index>.log`. If any job ends before its log says complete, rerun only that index, for example index 3:

```bash
sbatch --export=ALL,M=<choice_m>,N=<choice_n_per_family> --array=3 hpc/tufts_eval.sh
```

## Step 4: send results back

Zip and send the whole results folder plus the logs:

```bash
tar czf trackc_results.tar.gz output/trackc_v2_research logs
```

Vihan runs `--verify-only` and the frozen analysis on the PC (`src/trackc_v2_analysis.py`). Every cell is reported, whatever the result.

## Time and size (estimates, not measured)

Pilot-based estimate (`paper/trackc_freeze_amendment_2026-11-19.md`): about 1.6 core-hours per world at m = 2 on a Tufts core, about 1,012 core-hours per cell, **about 9,400-9,600 core-hours in total** for the 9 jobs (the earlier 4,400 figure used PC timing and fewer worlds). With 9 jobs x 30 workers: about 34 hours per job, two 22-hour submissions each, if the queue starts them at once.

## If something goes wrong

- Test fails with `version https://git-lfs` in the error: a file was not downloaded. Run `git lfs pull --include="<that path>"`.
- `scipy` fails to build: use Python 3.13 (Step 0), which has ready-made packages.
- Job killed for memory: raise `--mem` in the script and resubmit (it resumes).
