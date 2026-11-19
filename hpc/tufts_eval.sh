#!/bin/bash
# Track C evaluation on Tufts HPC (ONLY after the freeze commit). See paper/trackc_tufts_hpc_runbook_2026-11-18.md.
# Submit: sbatch --export=ALL,M=<choice_m>,N=<choice_n_per_family> hpc/tufts_eval.sh
# Array index -> run: 0 primary, 1-6 sensitivity cells, 7 EARF evidence (locked), 8 EARF primary_anchor (descriptive).
# Safe to resubmit: finished worlds are skipped.
#SBATCH --job-name=trackc_eval
#SBATCH --partition=batch
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=0-8
#SBATCH --output=logs/eval_%A_%a.log
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
: "${M:?set M=choice_m}" "${N:?set N=choice_n_per_family}"
W=$(( SLURM_CPUS_PER_TASK - 2 ))
CELLS=(primary scan28 scan56 scan112 miss10 penalty0 penalty10 primary primary)
EARF=("" "" "" "" "" "" "" "--earf evidence" "--earf primary_anchor")
i=$SLURM_ARRAY_TASK_ID
venv/bin/python -I src/trackc_research_driver.py --phase evaluation --cell "${CELLS[$i]}" \
    --patients "$M" --n-worlds "$N" --workers "$W" --max-hours 22 ${EARF[$i]}
