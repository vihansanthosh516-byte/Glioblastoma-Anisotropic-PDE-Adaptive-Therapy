# Track C prior-art search (first pass)

Date: 2026-10-28 (Claude). Method: 7 web searches plus one arXiv abstract read. Search terms: particle filter adaptive therapy; adaptive therapy misspecification / missing measurements / monitoring frequency; chance-constrained or risk-constrained adaptive therapy; imperfect monitoring containment (Viossat, Noble); GBM adaptive temozolomide scheduling; POMDP belief state adaptive therapy.

**Limit:** this is a first pass. I read abstracts and search summaries, not full texts. Google Scholar, PubMed and bioRxiv were not searched one by one. A full search can still find closer work. No "first" or "novel" claim is allowed from this note.

## Closest work found

| Work | What it does | Overlap with Track C |
|---|---|---|
| Browning, Crossley, Murphy, Byrne & Hamis (2026), arXiv:2608.18387, "Adaptive therapy under parametric, structural, and measurement uncertainty" | Bayesian virtual cohort calibrated to prostate PSA data; noisy, infrequent measurements; misspecification gives unreliable predictions; AT50 still gives higher TTP | High: robustness of adaptive therapy to measurement and model uncertainty. The abstract does not describe a belief-state controller, missing scans, or a progression-probability trigger (full text not read) |
| Wang, Scott & Vladimirsky (2024), PLoS Comput Biol, threshold-aware adaptive therapy | Stochastic tumor model; maximizes the probability of success within a cost budget | Medium: probability-based dosing objective. Not a filter on noisy scans |
| Fischer, Vázquez-García & Mustonen (2015), PNAS, "The value of monitoring to control evolving populations" (arXiv:1406.6957) | Finite monitoring causes loss of control; preemptive strategies partly recover it | High for the missing-scan idea: monitoring limits are known to matter |
| Gallagher et al. (2025), Bull Math Biol, "Deriving optimal treatment timing for adaptive therapy" (doi:10.1007/s11538-025-01525-y) | Discrete-visit threshold protocol; trade-off between monitoring frequency and TTP | High: visit frequency vs TTP |
| Cheke et al. (2026), single-threshold adaptive therapy (Univ. of Greenwich repository) | High-frequency monitoring not strictly required | Medium: single-threshold rule, like our comparator |
| de Saporta and coworkers (2024–2026), POMDP / Bayes-adaptive POMDP for cancer follow-up (multiple myeloma; arXiv:2401.03972) | Belief-state treatment and visit-timing decisions | High: belief-state dosing exists; not GBM, not this misspecification grid |
| Carleton University student project (2026), particle-filter belief controller for adaptive therapy | Particle filter estimates resistant fraction from burden; belief-aware controller | High in idea, but not peer-reviewed; details not verified |
| Viossat & Noble (2021), containment strategies | Containment can beat aggressive treatment | Background: our containment rule follows this family |

## What this means

- **Belief-based dosing is not new.** Particle filters, POMDPs and probability-based objectives for cancer dosing all exist.
- **Robustness to noise and infrequent monitoring is an active topic** (Browning 2026; Gallagher 2025; Fischer 2015).
- **Not found in this first pass:** a sealed, preregistered test of a belief-based containment rule in a GBM two-population model across a misspecification grid, including the "told vs not told about a missing scan" comparison. This is a gap in what I found, not proof that no such work exists.
- **The safest uniqueness statement today:** "a sealed, pre-specified robustness evaluation of belief-based containment dosing in a simulated GBM model, including missed scans". Do not call the method new.

## Next steps

1. Read the full text of Browning et al. 2026 and Gallagher et al. 2025. Check whether they test missed visits.
2. Search PubMed and bioRxiv for "missed visit" / "missing measurement" adaptive therapy.
3. Add the confirmed references to the plan bibliography.

## Sources

- https://arxiv.org/abs/2608.18387
- https://doaj.org/article/f3365eb8ae534a3f93a56a9f92a2794d
- https://ar5iv.arxiv.org/html/1406.6957
- https://link.springer.com/article/10.1007/s11538-025-01525-y
- https://gala.gre.ac.uk/id/eprint/53918/
- https://www.arxiv.org/pdf/2401.03972
- https://hal.inrae.fr/MIAT/hal-05054392v4
- https://report-hub.scs.carleton.ca/projectPage/2864/
- https://ar5iv.labs.arxiv.org/html/2211.17041

## Addendum 2026-11-03 (14-expert audit; read from abstracts and search snippets, not full texts; verify in Stage 1 of `paper/trackc_restart_plan_2026-11-03.md`)

Missed by the first pass:
- Browning et al. 2026 (arXiv:2608.18387): the abstract also introduces a metastasis-risk metric, and adaptive therapy may be disadvantageous when sustained burden is counted. Relevant to the tail framing. Full text not read.
- Corigliano et al. (bioRxiv 10.1101/2024.11.29.626082): recovery periods and sub-clinical doses can outperform continuous dosing (colorectal targeted therapy; validation not stated). Closest to forced holidays.
- Segura-Collar et al. 2022, Neuro-Oncol Adv 4(1):vdac155 (PMC9616068): glioblastoma, longer temozolomide spacing, persister compartment; fixed schedules, proof of concept. Its rates sit in Supplementary Table S2 (not read).
- Sun et al. 2025 (arXiv:2502.09392): success depends on the pause/resume threshold; no optimal strategy under some constraints.
- Sargun et al. 2026 (arXiv:2606.10376): Bayes-filter belief with an entropy-triggered measurement rule (colorectal). Closest to a belief-based trigger.
- Viossat & Noble 2021, Nat Ecol Evol 5(6):826-835 (PMC8967123, read by one expert): containment can beat maximum tolerated dose without a fitness cost; intermittent containment is slightly worse than continuous.
- Gallagher et al. 2025 preprint authors: Gallagher, Strobl, Anderson, Maini (medRxiv 10.1101/2025.04.01.25325056).

Searches for a worst-decile (CVaR) objective in adaptive therapy found none. That used 6 WebSearch queries; PubMed, Google Scholar and Europe PMC were not searched, so it is a gap in the search, not proof.
Not opened: Zhang 2017 (login redirect), Hansen & Read 2017, Strobl 2022, Fischer 2015 full text, Wang 2024 full text, Kiram 2026.
The uniqueness statement above still holds. No novelty claim is allowed.
