# Project State: last updated 2026-10-04 (evening)

## What this project is
Computational glioblastoma modelling with three tracks (A single cell to clinical, B PDE and adaptive therapy, C digital twin and RL).
The paper title is fixed. The paper is NOT written. Methods is the restart point. See [[PAPER_FINDINGS_LEDGER]] for every number and [[LIMITS_CURRENT]] for every limit.

## Where things stand
- Vault is audited against output files ([[AUDIT_REPORT]] is superseded by the ledger).
- 14 paper figures are in `paper/figures/` (made by `scripts/make_paper_figures.py`, `scripts/make_framework_figure.py`).
- Track A is the weak track: the random cell-level split inflates it (scripts 86, 88). Track A does not feed Track B (inflammation score is 1.0 for all 61 patients).
- Track B and C fixes done: script 48 at literature D (holds), script 89 fair dual-drug baseline (the second drug explains the rescue), script 87 checks.
- External checks done: CGGA (script 90), LUMIERE volumes (script 91).

## Running
- Script 88 (leak-free C-GAT), WSL, CPU, 5 folds x 3 variants, about 30 min per fold. Output `output/cgat_leak_free.json` (written per fold).

## Not done
- RHUH-GBM analysis; script 81 whole-tumour target (about 12 h); E_MAX sweep; power analysis for probe-then-commit; stronger forecast baselines; Dice forecast on new cohorts (needs masks).
- References in [[PAPER_REFERENCES]] are all unchecked. "Weidner 2023" has no verified title.

## Rules
- No commit of a paper section before the user reviews it. Every number cites a source file. Unknown numbers are `[MISSING]`. Abstract last.
- One commit per verified fix; `git status --short` first; never `git add -A`; push at once.

## Where things live
- `src/` scripts 01-91; `output/` JSON evidence; `data/external/` new cohorts (git-ignored); `vault/` notes; `paper/figures/`.
