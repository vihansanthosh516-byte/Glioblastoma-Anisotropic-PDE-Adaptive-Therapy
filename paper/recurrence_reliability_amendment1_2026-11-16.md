# Recurrence reliability protocol amendment 1 (2026-11-16, before any scoring)

## What happened

The first execution of the locked run (`run_20261116_locked`, commit 642f9ac) stopped while loading specimens: `GSM5319524_SF4810: unannotated nuclei pass filters; refusing to fill labels` (`output/recurrence_reliability/run_20261116_stdout.txt`). The stop occurred in the loading loop, before the reference pools, control genes, or any per-nucleus score were computed. No program score, delta, or direction was produced or seen. The empty run directory was removed.

## Cause

Two frozen rules conflict. Phase B eligibility admits a library whose author annotation covers at least 95% of its barcodes (SF4810: 575 of 577; SF10433: 1,172 of 1,173). Protocol section 2 says any filtered nucleus without a label stops the run. The conflict is a specification error, not a data problem.

## Change

Nuclei without an author label are dropped from every analysis and are never assigned a label. Their number per library is written to `result.json` (`unlabelled_nuclei_dropped`). A library in which more than 5% of filtered nuclei lack a label still stops the run (the Phase B threshold). This affects at most 3 nuclei in two libraries, none labelled tumor, so the malignant-cell primary analysis uses exactly the same nuclei as before; only the whole-tissue secondary analysis (S1) loses at most 3 nuclei.

Everything else in the protocol and config is unchanged. New run ID: `run_20261116_locked_a1`. Code change: `Specimen.__init__` in `src/recurrence_reliability/run.py`.
