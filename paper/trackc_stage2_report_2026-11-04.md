# Track C Stage 2 report (interim, started 2026-11-04)

**Status: INTERIM.** Everything here is a computer simulation with assumed parameter ranges, run on development seeds. The numbers are exploratory. None is a test result, none is a clinical claim, none is a patient-specific claim, and none says a digital twin is validated. Seeds are listed in `output/trackc_seed_registry.json`. Each run folder has a `manifest.json` (git commit, protocol hash, seed blocks, clean-tree flag). Plan: `paper/trackc_restart_plan_2026-11-03.md`. Design: `paper/trackc_anchored_spec_2026-11-04.md`.

Each test below answers four questions: what was asked, did it work, what is the weakest point, what comes next.

## Test A (script 190): how the worst-decile score moves with the number of noise repeats

Status: running at the time of writing; no number exists yet. It will be added when `output/trackc_stage2_A_tail_estimator/results.json` is written and read back.

## Test B (script 191): exact identities and gain-versus-burden checks

- **Asked.** Which exact statements follow from the equations, and where can a feedback rule gain at all?
- **Did it work.** Yes. Three lemmas and one proposition (`paper/trackc_lemmas_2026-11-04.md`), checked by 13 passing tests in a scratch copy (`tests/test_trackc_lemmas.py`) and by an independent read-only algebra check. With no crowding, no cost of resistance and no switching, front-loaded dosing gives a time to progression at least as long as any allowed schedule. So any gain of feedback in the full model has to come from crowding (starting burden B0), the cost c of resistance, or switching.
- **Weakest point.** It is a statement about the assumed equations. It does not show that real tumours behave this way.
- **Next.** Used to choose the regime axes in tests C to E.

## Test C (script 192): which assumed worlds reproduce the published standard arm?

Files: `output/trackc_stage2_C_calibration/{manifest.json,worlds.csv,kept_worlds.csv,results.json}`. Development seeds 20261801 (worlds) and [20261802, world] (patients). 500 worlds, 400 simulated patients each.

- **Asked.** The published standard arm (RTOG 0525, Gilbert 2013) has a median PFS of 5.5 months (95% CI 4.7 to 6.1, which is 143 to 186 days) with scans about every three months. For which settings of the parameter ranges does a simulated standard arm (6 courses, open loop, scan every 84 days, RANO-like progression call) land inside that interval?
- **Did it work.** Yes. 81 of 500 worlds (16.2%) pass (`results.json`, `base.n_pass`). The earlier ranges, unchanged, give a median of 240 days (`old_model.base.median_interp`), above the interval.
- **What it shows (simulation, exploratory).** The anchor constrains one axis strongly: the kill strength relative to growth (the m scale). The median kept world has an m-scale of 0.135 times the old range, so m is about 0.8 to 2.7 instead of 6 to 20. Among the 116 worlds with an m-scale of at least 0.5, only 2 pass; their median PFS sits at 196 to 321 days (5th to 95th percentile), all too long. The other axes are close to unconstrained: pass fractions by third are flat for the f0 shift, the cost c and the switching rate mu (0.14 to 0.19), and lower only for the highest reversion rate nu = 0.04 (0.10, against 0.18 and 0.20). Faster growth r shortens the median PFS and a higher B0 lengthens it (Spearman -0.53 and +0.45 against the world's median PFS). In the kept worlds about one third of patients are called progressed at the first scan (S(84) about 0.67).
- **Weakest points.**
  1. The set of consistent worlds depends on how progression is called. With a 1.25 factor instead of 1.4, 54 worlds pass but the Jaccard overlap with the base set is 0.03; with the old 1.2 x baseline rule, 30 pass and the overlap is 0.03. Scans every 56 or 112 days give overlaps of 0.86 and 0.68; the RANO floor (0, 0.1, 0.2 x baseline) changes nothing. The progression rule is an assumed translation of the trial rule and uses the true burden without measurement error.
  2. One published median fixes one condition. Many worlds pass and the anchor does not say which one is true.
  3. The kill strength is not anchored by any volume-change-per-cycle measurement (none was found). The weak kill the anchor implies may only mean this model class needs a weak drug to produce a short PFS.
  4. The model has no patients who never respond to TMZ and no early deaths or drop-outs; the early part of the PFS curve is not checked.
  5. 400 patients per world: a world near an edge of the interval can flip. Test D re-draws the kept worlds.
- **Next.** Test D (negative controls) on the 81 kept worlds.

## Test D (script 193): do the kept worlds also pass the two negative controls?

Files: `output/trackc_stage2_D_negative_controls/{manifest.json,controls.csv,controlled_worlds.csv,results.json}`. Development seed [20261803, world]. 3,000 fresh simulated patients for each of the 81 kept worlds.

- **Asked.** Real trials found little or no gain from more drug: dose-dense against standard TMZ (RTOG 0525, PFS HR 0.87, 95% CI 0.75 to 1.00) and 12 against 6 cycles among patients still progression-free after cycle 6 (GEINO 14-01, PFS HR 0.99, 95% CI 0.71 to 1.38). Which kept worlds show a model hazard ratio inside both intervals?
- **Did it work.** Yes. Of 81 kept worlds, 78 still land in 143 to 186 days on fresh patients (3 flipped, a small winner's-curse effect). The dose-dense control is met by 65 worlds (80 at a drug-effect multiplier of 1.5, 48 at 3.0). The 12-versus-6 control is met by 41 worlds. 39 worlds pass everything and form the **controlled set** (`controlled_worlds.csv`). The kill rule (no world passes) did not fire.
- **What it shows (simulation, exploratory).** The two controls select slower growth, weaker kill and more initial resistance: in the 41 landmark-passing worlds the mean log shift of r is -0.29 against +0.37 in the others, the mean log m-scale is -2.29 against -1.48, and the mean log f0 shift is 1.16 against 0.60. They do **not** select on the four axes where feedback dosing could matter: the mean B0 centre is 0.314 against 0.315, the cost c 0.125 against 0.106, mu 0.021 against 0.013 and nu 0.012 against 0.013 (`results.json`, `axes`). The published data constrain growth, kill and resistant fraction; they leave B0, c, mu and nu open.
- **Weakest points.**
  1. The dose-dense arm is the standard arm with the drug effect multiplied by 2.0. That ratio of drug exposure is assumed (the label doses of the dose-dense arm were not verified); the 1.5 and 3.0 runs are the sensitivity.
  2. The model never reproduces the published landmark point estimate. Its hazard ratio for 12 against 6 courses has a 5th to 95th percentile of 0.26 to 0.89 across the 81 worlds (median 0.73), so it shows some benefit from extra courses in nearly every world, while the trial point estimate is 0.99. The 41 passing worlds are inside the trial's wide interval (lower bound 0.71), not at its centre. The trial was a small phase II study.
  3. Both controls compare open-loop arms on the same simulated patients, so the Cox interval is not a valid confidence interval; the flags use the point estimate.
  4. Toxicity, drop-out and death before progression are not modelled. Real trials include them.
- **Next.** Test E on the 39 controlled worlds.

## Test E (script 194): policy contrasts across the controlled worlds

Status: written and smoke-tested; waits for test D.
