#!/bin/bash
# Track C development pilot on Tufts HPC. See paper/trackc_tufts_hpc_runbook_2026-11-18.md.
# Safe to resubmit: finished worlds are skipped.
#SBATCH --job-name=trackc_pilot
#SBATCH --partition=batch
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/pilot_%j.log
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
W=$(( SLURM_CPUS_PER_TASK - 2 ))
PY=venv/bin/python
venv/bin/python -c "import sys,numpy,scipy;print('versions',sys.version,numpy.__version__,scipy.__version__)"
PYTHONPATH=src $PY -m pytest -q -p no:cacheprovider tests/test_trackc_*.py
$PY -I src/trackc_research_driver.py --phase dev_pilot --cell primary --patients 2 --n-worlds 3 --workers "$W" --confirm-development
$PY -I src/trackc_research_driver.py --phase dev_pilot --cell primary --patients 2 --n-worlds 3 --workers "$W" --confirm-development --earf evidence
$PY -I src/trackc_research_driver.py --phase dev_pilot --cell primary --patients 12 --n-worlds 20 --workers "$W" --max-hours 22 --confirm-development
$PY -I src/trackc_dev_pilot_design.py --pilot-dir output/trackc_v2_research/dev_pilot/primary/m12 --timing-dir output/trackc_v2_research/dev_pilot/primary/m02
