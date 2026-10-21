# Probabilistic RANO: is the +40% progression line above the measurement noise?

Glioblastoma project for the Greater Austin Regional Science and Engineering Fair (GARSEF; tentative February 2027 deadline, official event dates/forms unverified). One author. Status date: 2026-10-19. Start with `START_HERE.md`.

## Question
RANO 2.0 calls progression at a tumour volume rise of 40% or more and partial response at a fall of 65% or more (progression against the nadir, response against baseline; the 0.5 cm3 floor used here is our own operational rule, not in the RANO 2.0 text). Volumes now come from automated segmentation. How noisy are they, and what does that noise do to single patients and to trials?

## Results so far (all exploratory or simulation; each declared before it ran)
| Test | Result | Weakest point | Source |
|---|---|---|---|
| A15: do two tools give the same consecutive-scan ratio class? | Consecutive-scan ratio class differs in 54.5% of 121 eligible LUMIERE pairs (CI 45.6% to 63.8%) | Tools are not a scan-rescan pair; 71.9% of pairs fell below the floor or lacked a volume; one centre | `output/rano_noise/rano_tool_dependence.json`, ledger 2br |
| A15: minimal detectable change vs the +40% line | Both proxy MDC95 estimates exceed line in smallest band; proxy estimates straddle line in larger bands | Same | same |
| A18: sensitivity (fine noise grid, correlated noise, death, floors, nadir) | Under the declared materiality rules, bias and power loss start near change-SD 0.10 for single-scan and 0.15 for two-scan rules; correlation reduces but does not remove them; flip rate 53% to 63% across floors, 54.4% against the nadir (68 pairs, 20 patients) | Noise level in real tools still unknown; simulation assumptions; nadir subset small | `output/rano_sensitivity/rano_sensitivity_summary.json`, ledger 2bt |
| A19: simulation anchored to a published repeat-scan study (change-SD 0.14 to 0.34) | Earlier PFS by 7.6 to 20.4 weeks and power loss 0.20 to 0.82 (N=100, 8-week scans, single-scan rule); two-scan simulated bias ranges from 3.6 to 12.3 weeks earlier | One study, 20 patients, large tumours, BraTumIA not our tools; assumed RC%-to-SD conversion | `output/anchored_sim/anchored_sim_summary.json`, ledger 2bu |
| A21: audit of the A19 conversion (three readings, change-SD 0.079 to 0.342) | Earlier modeled PFS across the seven selected simulated values (17 distinct converted values exist) (single-scan, point estimates: -5.8 to -15.0 weeks at N=100, 8-week scans); single-scan PFS-bias magnitude varies about 2.5x across these point estimates; two-scan rule keeps it under a month at the lowest reading | Scale of the published RC unknown; one study; independent normal error | `output/conversion_audit/conversion_audit_summary.json`, ledger 2bx |
| A16: trial simulation | Noise makes PFS earlier and cuts power a lot (N=100, D=8: power 0.97 with no noise, 0.29 at s=0.25, 0.05 at s=1.0). Two of four declared predictions NOT supported (H-S3, H-S4) | Effect sizes depend on the assumed noise model. Simulation, assumed parameters | `output/trial_sim/trial_sim_summary.json`, ledger 2bs |

## Tool
`python src/rano_calc.py V0_mm3 V1_mm3 [s]` gives the probability that the true change is progression, response or stable, and the minimal detectable change (`src/rano_calc.py`, tests `tests/test_rano_calc.py`, ledger 2bv). Model-based; not for clinical use. Pass the nadir volume as V0 for progression.

Viewer: `output/viewer_rano/viewer_rano.html` (open in a browser; 62 LUMIERE patients; bands, nadir, +40% line, class probabilities; display only; ledger 2bw).

## Not done
Paper, poster. Novelty is unverified (limited snippet-level searches; no retained count log). Nothing here is a clinical claim.

## Files
`START_HERE.md` (state and rules), `analysis_plan_v1.md` (declarations), `vault/PAPER_FINDINGS_LEDGER.md` (results), `paper/` (plans), `src/140_*` through `src/145_*` and `src/rano_calc.py` (current analyses and tools), `archive_old_gbm/` (old GBM documents), `CLAUDE.md` (run commands).
