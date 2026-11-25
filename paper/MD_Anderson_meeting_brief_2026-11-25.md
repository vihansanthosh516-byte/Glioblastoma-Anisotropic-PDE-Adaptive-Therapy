# Meeting brief: MD Anderson (2026-11-25)

Short version to say out loud first:

> "I'm a high school student doing computational glioblastoma research. I have three tracks: molecular prediction from tumor tissue (A), MRI tumor growth forecasting (B), and treatment-policy simulation (C), plus an integration test that links them. I pre-register every test and run it once, so most of my results are honest negatives. My main blocker is data, and that's what I'd like to ask about."

Every number below comes from a file in the repo (path given). Weak points are listed next to each result.

---

## Track A: molecular and spatial prediction from tumor tissue

**Question.** From tumor gene expression, can a new method (1) tell tumor core from the infiltrating edge, as labelled by a pathologist, and (2) predict temozolomide (TMZ) response?

**Region part (core vs edge)**
- New design for few patients: a pathologist labels one tissue section; the method maps core vs edge on the same patient's sections from other tumor sites. Data: GBM-Space Visium (De Jong et al., Nature 2026), where one pathologist annotated H&E sections with Ivy GAP regions and two more validated them.
- Status: plan locked (`paper/trackA_region_one_section_transfer_prereg_LOCKED_2026-11-24.md`); one run in progress. Result: **[fill in when the run finishes]**.
- Weak points: only 4 patients have labels at two or more sites, so this is a within-patient claim; the PASTE2 alignment comparator was not run.

**Drug part (TMZ)**
- Mouse test: a fixed MGMT + mismatch-repair expression score ranked TMZ survival benefit in 27 Mayo Clinic orthotopic PDX beyond MGMT methylation (partial rho 0.367, one-sided p = 0.040), beating all comparators (`paper/trackA_tmz_mayo_repair_axis_result_2026-11-24.md`).
- Human test of the same frozen score: **failed** on GSE7696 (68 patients, TMZ+RT 42 vs RT 26). Interaction z = +1.53, wrong direction, one-sided p = 0.94 (`paper/trackA_tmz_gse7696_result_2026-11-24.md`).
- Earlier: a published 102-gene TMZ panel did not rank in-vitro TMZ response in 83 HGCC lines (rho -0.102) (`paper/trackA_hgcc_tmz_102_result_2026-11-12.md`).
- Weak points: the mouse pass is small and barely significant; no public dataset links TMZ response to expression in enough patients, or per cell inside a tumor.

## Track B: MRI tumor growth forecasting

**Question.** Can a new method forecast post-operative tumor volume change on MRI better than simple baselines?
- Locked test on two external cohorts (CFB + WAW; 181 patients, 295 scan pairs): the VIG method's error was 1.1508 vs 1.1538 for "no change" (difference -0.003, 97.5% CI -0.051 to 0.034). **Not confirmed**; it beat 4 of 8 comparators (`paper/trackB_vig_cfb_waw_result_2026-11-19.md`).
- Development: three more method families (level reversion, noise-aware, image features) gave no clean gain. Diagnosis: the bottleneck is data (`START_HERE.md`).
- Tumor segmentation tool checks passed: on 243 UCSD scans the median volume error was 0.093 log units (limit 0.20) (`output/track_b_vig_cfb_waw_20261115/gate_g1.json`).
- Weak points: underpowered (about 0.5 power); a power estimate says about 365 patients would be needed for the observed effect size.
- Next: more longitudinal cohorts (Burdenko via dbGaP needs a faculty sponsor; BraTS post-treatment; TCIA).

## Track C: treatment-policy simulation

**Question.** In simulation, does an uncertainty-aware adaptive treatment controller beat fixed schedules?
- Everything is frozen. An outcome-blind pilot ran on the Tufts HPC (job 7359807, 40/40 worlds) and fixed the design (m = 2, 317 worlds per family) (`paper/trackc_freeze_amendment_2026-11-19.md`).
- The main evaluation is about 9,400-9,600 core-hours (estimate) and is not yet run.
- Weak points: simulation only; no clinical benefit or patient-specific effect can be claimed.

## Integration: linking MRI to biology to treatment

- Four real-data checks of an MRI-to-biology link all **failed**: TCGA (Gate 1), Hu 2023 (Gate 1b), CPTAC + Ivy GAP confirmation (rho 0.114, p 0.229), and UPENN survival (HR 0.869, one-sided p 0.996). Diagnosis: overfitting and selection, not code (`paper/integration_failure_diagnosis_2026-11-14.md`).
- New method EARF: the simulation controller's uncertainty about the biology signal is set from these measured (weak) links. Locked; it runs with the Track C evaluation (`paper/integration_earf_protocol_2026-11-14.md`).
- Weak point: simulation, and partly circular (see protocol).

---

## Questions to ask

**Data (most important)**
1. Does MD Anderson have, or know of, glioblastoma spatial transcriptomics or multi-region samples with **pathologist labels for tumor core vs infiltrating edge** that a student could access, through a collaboration or a data use agreement?
2. Is there any dataset where **temozolomide response** is measured with gene expression in patient tumors, ideally with a no-TMZ comparison or per-cell response (for example, ex vivo drug testing)?
3. For MRI: are there **longitudinal post-operative MRI** cohorts (several scans per patient, with tumor masks) that could be shared or analyzed under supervision?
4. Would a faculty member be willing to **sponsor a dbGaP access request** (Burdenko cohort, phs004225), or advise on how students usually do this?

**Science**
5. For core vs edge: is the infiltrating edge as labelled on H&E the right target, or would you trust a different reference (for example, 5-ALA fluorescence, MRI-guided biopsies)?
6. For TMZ: is MGMT plus mismatch repair a dead end for prediction beyond MGMT methylation, in your view?
7. What would make a result like this useful to clinicians or researchers at MD Anderson?

**Career / mentorship**
8. Is there a way to get feedback or mentorship on this work, or a summer program that fits it?

---

## If asked "what's your strongest result?"

Be honest: "My strongest finding is negative but solid: several published-style predictors fail on independent data when tested properly. My best positive result is a small mouse study (p = 0.040) that did not replicate in humans. My main contribution is the testing discipline, and I need better data to go further."
