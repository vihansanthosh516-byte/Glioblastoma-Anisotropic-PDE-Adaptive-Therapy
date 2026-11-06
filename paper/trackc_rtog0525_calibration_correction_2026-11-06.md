# Track C RTOG 0525 calibration correction

**Status:** dated scientific correction to the draft calibration design. It does not change or delete prior outputs. No new policy-outcome simulation was run for this correction.

## What the source establishes

The RTOG 0525 primary report defines progression-free survival (PFS) from random assignment to progression or death, or censoring at the last assessment while alive and progression-free. Random assignment occurs after concomitant radiotherapy and temozolomide; adjuvant temozolomide begins four weeks after radiotherapy. The trial reports a standard-arm median PFS of 5.5 months (95% CI 4.7–6.1 months). It uses Macdonald-style progression: at least a 25% increase in tumor size or a new lesion. ([Gilbert et al., 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3816958/))

## Why the existing match is not RTOG PFS calibration

The draft simulation starts its clock at the first adjuvant course, observes on an 84-day schedule, scores a RANO-like volume increase from the prior nadir, and has no death process. Its 143–186-day screening interval therefore differs from RTOG 0525 in **time origin, event definition, scan timing, and event types**. The project code itself records that deaths and drop-outs are not simulated (`src/192_trackc_calibration.py`). The previous counts of worlds passing this screen remain historical results under the stated simulation rules; they cannot be described as reproducing clinical RTOG PFS or validating treatment policies.

It would not be adequate to repair this by adding a fixed 28 days after the fact. The treatment start, scan schedule, progression criteria, trial eligibility/risk set, censoring, and progression-or-death process must be handled together and frozen before a new run.

## Corrected analysis rule

For the current research plan, Track C will report **simulated time to modeled progression** under its prespecified baseline and RANO-like endpoints. The model will not call these PFS, will not filter worlds to match the RTOG PFS interval, and will not treat RTOG or GEINO PFS hazard ratios as directly comparable calibration gates. Trial summaries may be shown as context only.

If the project later adds a defensible trial-aligned PFS process, that work will require a dated amendment that defines the randomization clock, adjuvant-treatment delay, assessment schedule, progression rule, death/censoring process, risk set, calibration targets, and new unused seeds before policy outcomes are examined. Existing runs and seed allocations will remain unchanged and labeled as historical.

## Primary-source checks

- RTOG 0525 randomized after completion of concurrent chemoradiation and initiated adjuvant treatment four weeks after radiotherapy; its PFS event is progression or death from random assignment. ([Primary trial report](https://pmc.ncbi.nlm.nih.gov/articles/PMC3816958/))
- The RTOG 0525 progression/death analysis describes failure as death or progression and plots PFS from randomization. ([Wang et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC5654351/))
