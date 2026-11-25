# Run Track C + Integration (EARF) on Tufts HPC

One command runs both. Simulation only. The pilot is already done (job 7359807).
Full detail: `paper/trackc_tufts_hpc_runbook_2026-11-18.md`.

## 1. Log in

```powershell
ssh YOUR_TUFTS_USERNAME@login.pax.tufts.edu
```

Or use the browser: `ondemand.pax.tufts.edu` → Clusters → Shell Access.

## 2. Go to the project and update

```bash
cd trackc
git pull
```

After this, change nothing until the run is done.

## 3. Start the run (9 jobs)

```bash
sbatch --export=ALL,M=2,N=317 hpc/tufts_eval.sh
```

| Job index | What it is |
|---|---|
| 0 | Track C main result |
| 1-6 | Track C sensitivity checks |
| 7 | Integration (EARF) locked test |
| 8 | Integration (EARF) description only |

## 4. Check on it

```bash
squeue -u $USER
```

`PD` = waiting. `R` = running. Nothing listed = all ended.
Logs: `logs/eval_<jobid>_<index>.log`.

## 5. If a job stops early (24 h limit)

Rerun only that index. Example for index 3:

```bash
sbatch --export=ALL,M=2,N=317 --array=3 hpc/tufts_eval.sh
```

It resumes. Expect 1-2 reruns per job.

## 6. When all jobs are done

```bash
tar czf trackc_results.tar.gz output/trackc_v2_research logs
```

Send `trackc_results.tar.gz` to Claude. Claude checks and scores it on the PC.

## Size (estimate)

About 9,400-9,600 core-hours total, about 34 h per job if all start at once.

## Rules

- Do not look at results while it runs.
- Do not edit code or config.
- Job killed for memory: raise `--mem` in `hpc/tufts_eval.sh`, resubmit.
