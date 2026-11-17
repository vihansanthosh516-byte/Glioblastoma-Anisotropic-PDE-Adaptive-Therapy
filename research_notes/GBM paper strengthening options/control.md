# Track C and integration strengthening options

## Which additions provide distinct scientific value?

### Takeaway
The strongest feasible addition is a model-falsification and decision-stability study that asks whether EARF's apparent benefit survives a generating mechanism not used to construct its controller. A measured-perturbation bridge offers greater biological value but has much larger access and translation risks. These are proposals, not verified new methods or demonstrated improvements.

### Cited Findings
- The current EARF protocol draws world truth and constructs its ambiguity box from the same two reliability anchors. It explicitly labels this partly circular and labels MRI-to-MES and HGCC associations as proxies, not patient drug-resistance measurements. It inherits the unchanged Track C precision gate. [Local protocol](../../paper/integration_earf_protocol_2026-11-14.md).
- Hiremath et al. (2025), *Identifiability and model selection frameworks for mathematical modeling of high-grade glioma*, develops a family of biology-based HGG models. Search-returned abstract establishes close prior art for model-family comparison and identifiability; full text retrieval did not finish. [Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12105797/).
- Chaudhuri et al. (2023) explicitly develops predictive digital twins for risk-aware clinical decision-making and radiotherapy optimization under uncertainty. Generic risk-aware optimization therefore is already established. Search abstract only in this retrieval; project prior-art report also discusses the paper. [Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10598726/).
- Segura-Collar et al. (2022), *On optimal temozolomide scheduling for slowly growing glioblastomas*, already investigates altered TMZ scheduling and its relationship to resistance and toxicity. Search abstract only here; local prior-art review records having opened full article earlier. [Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9616068/).

### Inferences
**Rank 1: independent-mechanism decision stability.** Hypothesis: a reliability-anchored controller retains benefit over no-signal and fixed-policy comparators when the true observation or treatment mechanism differs from the controller's model. This differs from existing parameter grids only if entire mechanisms are excluded from development and then tested, such as reversible state switching, delayed drug effect, transient treatment-associated MRI change, or signal reliability changing over time. Begin with a dated inventory of mechanisms already exposed in Track C; do not call existing family/sensitivity tests a new track. Match drug budgets, observation access, patient draws and objective across policies. Primary endpoint: paired utility difference versus no-signal controller on a prespecified held-out mechanism family; secondary worst-case regret, policy disagreement and harmful decisions. Failure: no superiority under the unchanged primary uncertainty/precision rule, or benefits disappear in the held-out mechanism. Propose only after current run gates close; no bypass of stopped evaluation. Rough effort: several weeks, estimated, dependent on software scope and runtime. Novelty would be the specific benchmark/evidence, not robust control or model selection themselves.

**Rank 2: observationally equivalent worlds, different decisions.** Hypothesis: two models that fit the same observed burden trajectories can recommend opposite schedules, even when forecast errors are similar. Construct paired synthetic models with a fixed tolerance for observation agreement, then assess whether their preferred actions diverge. Baselines: point-estimate MPC, existing robust MPC, EARF, no-signal and best prespecified simple rule. Endpoints: decision disagreement under observational equivalence and regret in each world, with equal treatment use. Failure: cannot produce indistinguishable observation histories under the locked tolerance, or proposed protection offers no prespecified regret reduction. This is scientifically useful even as a negative result: determines whether good forecasts support treatment decisions. Existing identifiability diagnostics, prior recovery, and signal/action sensitivity mean it is an extension, not a wholly new track. Novelty needs a focused review and theorem/benchmark beyond identifying uncertain parameters. Rough effort: 2-4 weeks, estimate.

**Rank 3: treatment-toxicity mechanism.** Hypothesis: controller rankings change when course count is replaced by explicit toxicity/recovery dynamics. Existing course penalties and same-dose comparisons do not necessarily represent bone-marrow recovery, delayed toxicity or treatment interruptions. Add prespecified toxicity and recovery states and compare matched cumulative dose and toxicity risk. Endpoints: tumor burden plus toxicity violation frequency and tail outcomes. Failure: no robust change or improvement across plausible parameters, or no independent parameter support. Do not label course count a measured toxicity score. This is less feasible than rank 1 because this search did not establish an open, patient-linked GBM blood-count/dose/time dataset. Novelty limited: TMZ scheduling and toxicity already have prior art. Rough effort: 4-8 weeks or longer, estimate; biological calibration remains blocked.

### Gaps
- No broad novelty claim is established. Literature coverage was focused and incomplete.
- Full-text Firecrawl scrape calls stalled and the combined call was interrupted after about 7.5 minutes. Therefore the web facts above are explicitly abstract/search-level, not newly full-text verified.
- No new simulation results, policy outcomes, patient masks, expression values, or reserved cohorts were read. No protocol changed.

## Can measured perturbation data replace assumed resistance signals?

### Takeaway
A direct pharmacotyping validation track is biologically stronger than adding an unsupported MRI-to-resistance model. It requires independent cultures/organoids with measured drug response and reliable sample origins; organoid assays alone cannot identify patient treatment benefit or longitudinal kill dynamics.

### Cited Findings
- Verduin et al. (2023), *Patient-derived glioblastoma organoids reflect tumor heterogeneity and treatment sensitivity*, studies organoids as representatives of source tumors for treatment-response/resistance investigation. Search abstract confirms topic; file-level access, treatment schedule data and assay joins were not verified here. [Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10733660/).
- Park et al. (2026), *Patient-derived organoids predict personalized drug ...*, reports 18 glioblastoma organoid lines in the search-returned abstract. It also states organoid establishment and growth to sufficient size generally required four weeks. These are abstract/search-level facts pending full-text verification. [Primary paper](https://www.sciencedirect.com/science/article/pii/S2666379126002673).
- An EGA study record describes patient-derived organoids as predictive drug-response models. EGA listing does not prove unrestricted access, independent origins, usable assay matrices or sufficient sample-level joins. [Study EGAS50000001221](https://ega-archive.org/studies/EGAS50000001221).
- Jungwirth et al. (2026) describes automated personalized glioblastoma organoid drug screening. Its existence rules out novelty based merely on individualized organoid screening. Supplement data availability, usable sample count and independence were not verified. [Primary paper](https://www.nature.com/articles/s41698-026-01552-5).
- The project already tested the frozen Ntafoulis expression score against HGCC TMZ response without association; this must remain a reported failure. [Local current state](../../START_HERE.md).

### Inferences
**Rank 4: independent perturbation transfer.** Hypothesis: a small, prespecified molecular or functional feature predicts measured TMZ response across independent biological origins better than MGMT-only and simple ridge models. First do metadata-only access/provenance screening of the organoid sources above. Require a specimen-to-culture-to-assay crosswalk, exact drug exposure and response definitions, untreated controls, origin-level independence, and enough independent origins for a power plan. Freeze one target and comparator before opening response values. Use origin-held-out development and an untouched external source. Failure: no eligible independent validation cohort, unstable mapping across assays, or no improvement over frozen simple baseline. Do not refit a panel on HGCC and call the same HGCC cohort external validation. Do not turn viability AUC into a patient-specific drug kill rate without exposure/time-course and independent model validation. Rough effort: 1-3 months if accessible; unbounded with controlled-access approvals. This would bridge A to measured response, with a separate later C transfer hypothesis; it does not create a matched MRI/expression/treatment patient join.

### Gaps
- No candidate here has passed sample-level eligibility, terms, overlap or minimum-size gates.
- The 18-line source is small; number of lines need not equal independent patients.
- Single-time viability curves may be insufficient to distinguish growth inhibition, cell killing, resistant fractions and reversible state changes.

## Which tempting additions should be deferred?

### Takeaway
More RL agents, generic scan scheduling, or another learned proxy do not directly repair current inference gaps. Choose a new track only if its endpoint and independent evidence differ from existing work.

### Cited Findings
- The project's existing review already covers risk-aware twins, MPC, RL scheduling, response-guided stopping, signal ambiguity, identifiability/recovery and value-of-information prototypes. [Local prior-art review](../../paper/trackc_prior_art_full_2026-11-04.md); [local novelty-ceiling report](../../reports/GBM%20therapy%20novelty%20ceiling.md).
- Current latest handoff states Track C still needs cohort/world count, runtime cap and clean installation, while EARF is locked and not run. [Local current state](../../START_HERE.md).

### Inferences
- Treat scan scheduling as a later extension only if measuring decision value of an extra observation, including its cost, under an independent held-out mechanism. Existing perfect-information headroom is not empirical value of a new MRI scan.
- Do not add an immune compartment solely for complexity. Need a distinct hypothesis, data-supported parameters, observable consequences and a failure test.
- Most publishable near-term package: preserve negative proxy results, evaluate existing EARF fairly, then add independent-mechanism falsification and decision identifiability. Biological perturbation work becomes a conditional next track after access gates.

### Gaps
- This note proposes scientific direction, not publication acceptance, clinical validation, a frozen protocol, or authorization to run stopped comparisons.
- Research date: 2026-11-16. Time estimates are planning judgments, not sourced empirical benchmarks.
