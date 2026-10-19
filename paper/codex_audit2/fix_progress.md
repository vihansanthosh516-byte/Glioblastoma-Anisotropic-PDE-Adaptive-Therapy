# Audit correction progress

Audit evidence and original progress.md remain unchanged. Status tracks corrections, not scientific validation. Existing unrelated working-tree changes are excluded. No project analysis scripts may run.

| Chunk | Status | Changed files and verification |
|---|---|---|
| 001 | DONE | vault/PAPER_FINDINGS_LEDGER.md and vault/Decisions.md: removed uncited claims, fixed script-78 count, added sources/scopes, superseded stray 10/23. Static JSON checks; TI minimum -0.978598582 confirms original rounded -0.98 despite report minimum error. |
| 002 | DONE | vault/PAPER_FINDINGS_LEDGER.md entries 2g/2h/2e/2i/2j/3. Read JSON keys, computed velocity error and cohort summaries, counted src/*.py. Report 002 also misstates PatientID_0251 swanson sparing (actual 0.887222); used artifact. No output regenerated. |
| 003 | TODO | Pending report/source review and corrections. |
| 004 | TODO | Pending report/source review and corrections. |
| 004a | TODO | Pending report/source review and corrections. |
| 004b | TODO | Pending report/source review and corrections. |
| 005 | TODO | Pending report/source review and corrections. |
| 006 | TODO | Pending report/source review and corrections. |
| 007 | TODO | Pending report/source review and corrections. |
| 008 | TODO | Pending report/source review and corrections. |
| 009 | TODO | Pending report/source review and corrections. |
| 010 | TODO | Pending report/source review and corrections. |
| 011 | TODO | Pending report/source review and corrections. |
| 012 | TODO | Pending report/source review and corrections. |
| 013 | TODO | Pending report/source review and corrections. |
| 014 | TODO | Pending report/source review and corrections. |
| 015 | TODO | Pending report/source review and corrections. |
| 016 | TODO | Pending report/source review and corrections. |
| 017 | TODO | Pending report/source review and corrections. |
| 018 | TODO | Pending report/source review and corrections. |
| 019 | TODO | Pending report/source review and corrections. |
| 020 | TODO | Pending report/source review and corrections. |
| 021 | TODO | Pending report/source review and corrections. |
| 022 | TODO | Pending report/source review and corrections. |
| 023 | TODO | Pending report/source review and corrections. |
| 024 | TODO | Pending report/source review and corrections. |
| 025 | TODO | Pending report/source review and corrections. |
| 026 | TODO | Pending report/source review and corrections. |
| 027 | TODO | Pending report/source review and corrections. |
| 028 | TODO | Pending report/source review and corrections. |
| 029 | TODO | Pending report/source review and corrections. |
| 030 | TODO | Pending report/source review and corrections. |
| 031 | TODO | Pending report/source review and corrections. |
| 032 | TODO | Pending report/source review and corrections. |
| 033 | TODO | Pending report/source review and corrections. |
| 034 | TODO | Pending report/source review and corrections. |
| 035 | DOING | vault/Decisions.md findings corrected: historical banner, simulated/post hoc scope, removed biological inference, marked absent threshold/sweep provenance unverified. Other four files remain TODO. |
| 036 | TODO | Pending report/source review and corrections. |
| 037 | TODO | Pending report/source review and corrections. |
| 038 | TODO | Pending report/source review and corrections. |
| 039 | TODO | Pending report/source review and corrections. |
| 040 | TODO | Pending report/source review and corrections. |
| 041 | TODO | Pending report/source review and corrections. |
| 042 | TODO | Pending report/source review and corrections. |
| 043 | TODO | Pending report/source review and corrections. |
| 044 | TODO | Pending report/source review and corrections. |
| 045 | TODO | Pending report/source review and corrections. |
| 046 | TODO | Pending report/source review and corrections. |
| 047 | TODO | Pending report/source review and corrections. |
| 048 | TODO | Pending report/source review and corrections. |
| 049 | TODO | Pending report/source review and corrections. |
| 050 | TODO | Pending report/source review and corrections. |
| 051 | TODO | Pending report/source review and corrections. |
| 052 | TODO | Pending report/source review and corrections. |
| 053 | TODO | Pending report/source review and corrections. |
| 054 | TODO | Pending report/source review and corrections. |
| 055 | TODO | Pending report/source review and corrections. |
| 056 | TODO | Pending report/source review and corrections. |
| 057 | TODO | Pending report/source review and corrections. |
| 058 | TODO | Pending report/source review and corrections. |
| 059 | TODO | Pending report/source review and corrections. |
| 060 | TODO | Pending report/source review and corrections. |
| 061 | TODO | Pending report/source review and corrections. |
| 062 | TODO | Pending report/source review and corrections. |
| 063 | TODO | Pending report/source review and corrections. |
| 064 | TODO | Pending report/source review and corrections. |
| 065 | TODO | Pending report/source review and corrections. |
| 066 | TODO | Pending report/source review and corrections. |
| 067 | TODO | Pending report/source review and corrections. |
| 068 | TODO | Pending report/source review and corrections. |
| 069 | TODO | Pending report/source review and corrections. |
| 070 | TODO | Pending report/source review and corrections. |
| 071 | TODO | Pending report/source review and corrections. |
| 072 | TODO | Pending report/source review and corrections. |
| 073 | TODO | Pending report/source review and corrections. |
| 074 | TODO | Pending report/source review and corrections. |

## Safety and continuation

Initial state: d719f95; paper/codex_audit2/progress.md already modified; audit reports/SUMMARY and many data/cache directories already untracked. Do not stage these unrelated files. No analyses, output writes, compilation, or tests run. Batch 001/002 verified by read-only JSON arithmetic, filesystem inventory and git diff --check. Full diff reviewed; only ledger/Decisions and this tracker belong to the batch. Documentation-only batch needs no Python compilation or pytest; full scratch-suite verification remains pending for code/test changes.
