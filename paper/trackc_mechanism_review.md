# Track C: biology mechanism review and dev test (2026-10-30, Claude + science team)

SIMULATION ONLY. No patient data. All mechanism rates below are ASSUMED unless a source is given.

## 1. Check of scripts 171, 173 and chart 175
- 171 (`output/trackc_sealed_tail_hybrid/sealed_results.json`): H1 FAILED (state_linked_mod CI lower -7.97); H2, H3 passed.
- 173 (`output/trackc_sealed_tail_trigger/sealed_results.json`): H1 FAILED (nominal CI lower -8.04); H2 passed. The tuned arm was the 171 hybrid again, so 173 is a new-seed replicate.
- Chart 175: all 32 points match the files. Weak points: no error bars on points (the files hold only paired-difference CIs), no -7 margin line, pass/fail not marked per suite, mixture H1 is 0.0 [0, 0] because both arms have the same worst-10% value there (not evidence). A forest plot of H1/H2 differences with the -7 line would show this better.

## 2. Mechanisms reviewed (literature agent; items marked "search only" were not read in full)
| Mechanism | Best source | What it supports | What it does not |
|---|---|---|---|
| TMZ-driven MGMT / resistance selection | Lasri et al. 2020, R Soc Open Sci, doi 10.1098/rsos.191243 (read) | Stochastic MGMT expression gives selection under TMZ; in vitro rates (one patient-derived line) | Dosing benefit; clinical dynamics |
| TMZ shifts MGMT methylation | Storey et al. 2019, JCO CCI, doi 10.1200/CCI.18.00062 (read) | Model of methylation shift; GBM schedule optimisation (prior art for schedule search) | Mechanism proof; patient outcome |
| Transient TMZ-tolerant state | Rabé et al. 2020, Cell Death Dis, doi 10.1038/s41419-019-2200-2 (search only) | Qualitative transient tolerant, slow-cycling state in U251 | Reversion time scale |
| Reversible GSC persisters | Liau et al. 2017, Cell Stem Cell, doi 10.1016/j.stem.2016.11.003 (search only) | Reversible drug-tolerant state (kinase inhibitors, not TMZ) | Day-scale rates; TMZ |
| Quiescent stem-cell reservoir | Chen et al. 2012, Nature, doi 10.1038/nature11287 (search only) | Quiescent cells regrow tumour after TMZ (mouse) | Human rates |
| Hypermutation / MMR loss | Touat et al. 2020, Nature, doi 10.1038/s41586-020-2209-9 (search only) | TMZ-linked hypermutation in recurrent glioma | Dosing time scales |
| Pseudoprogression | RANO, Wen et al. 2010 (not read) | Imaging confound | Not a growth mechanism |
| Adjacent prior art | Liu DM et al. 2025, Phys Rev E 112:024404, doi 10.1103/xjh4-mfvx (abstract read) | Alternating MTD and low-dose cycles delays resistance in a 3-type model | Not GBM; not belief-based; no missed scans |

No source found gives a day-scale rate for reversion of TMZ tolerance. No GBM fitness-cost-of-resistance paper was found.

## 3. Ranking (team view)
1. Reversible drug-induced tolerance (mu on drug, nu off drug): qualitative support, needs no new data, runs on this PC, predicts that drug breaks help. Risk: designed to favour breaks (circular); rates assumed. **Chosen for a dev test.**
2. MGMT-like recoverable efficacy loss: similar support, also favours breaks by construction; the existing drift axis (script 169) already covers the non-recovering part.
3. Quiescent reservoir: favours front-loaded dosing, so it makes H1 harder; useful only as a stress test.

## 4. Dev test (script 176, development seeds 20262100+k only; `output/trackc_dev_plasticity/dev_results.json`, commit 8f98c2c)
n 150 per suite, 2 replicates. mu, nu log-uniform: mod [0.005, 0.02] /day, strong [0.02, 0.05] /day (ASSUMED).
- Forced-holiday containment (riskh) cut the worst-10% CVaR: 96.0-164.5 d vs front-loaded 164.3-267.7 d in every suite.
- Plain containment minus front-loaded worst-10%: -1.6 (no mechanism), 0.0 (mod), -12.5 (strong).
- Best mean: containment (704.4, 773.9) or riskh max_run 3 (720.3). Best simple non-filter baseline by mean: raw 0.6 (673.4, 686.0, 741.0).
- Weak point: n 150 and only 15 tail scenarios per suite; the per-replicate values were not stored, so no CI. The small nominal gap may be noise; the strong-suite gap is large.

**Verdict:** the mechanism did not give a dosing rule room to beat front-loaded dosing in the worst 10%. Under strong plasticity it made the tail gap larger. The skeptic had set NO-GO conditions for a sealed test (third look at H1, circular mechanism); with this dev result, no sealed test was opened. No fresh test seeds were used.

## 5. What this means
- The worst-10% scenarios need early, continuous dosing. Rules that wait or pause lose there. This held across scripts 171, 173 and 176.
- The honest Track C finding is a trade-off: containment gains mean days (scripts 171/173: +16 to +60 d vs the comparator, +9 to +31 d vs front-loaded) at a small, replicated worst-10% cost (about 2-4 d). H1 failed twice.
- Claims allowed: simulation under assumed parameters; sealed pre-specified tests; trade-off. Not allowed: clinical benefit, "first", "novel", "validated", tolerance biology as measured.

## 6. Open option (not done; would be a declared third attempt at H1)
State-triggered front-loading: front-load only when the first 1-2 scans suggest fast growth, else contain. Same family as 171/173; threshold rules are established prior art; must be frozen on dev seeds and labelled as a third attempt.
