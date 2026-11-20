# Track B: how to get a fresh test for the same question (2026-11-19)

Question, target and score are unchanged: from a post-treatment scan plus scan dates, forecast y = log(core volume + 1) at the next scan; core = labels 1 + 3; score = MAE averaged per patient first.
Status: LUMIERE, CFB and WAW are spent for VIG (`paper/trackB_vig_cfb_waw_result_2026-11-19.md`). Spent cohorts may be used for development, never as a test.

## 1. What a new test needs
- Serial post-treatment MRI, T1Gd + FLAIR at least (our frozen tool `output/track_b_vig_cfb_waw_20261115/gate_g1.json` labels core from these two), with **scan dates or day counts**.
- Not in MU, UCSD, UCSF, LUMIERE, CFB, WAW, ACRIN 6677, Duke, RHUH.
- Size: the old effect (0.027) needed about 310 patients for 80% power. The measured effect on CFB+WAW was -0.003 (`result.json`), so a bigger sample alone will not rescue VIG as it is. A new test needs a changed method AND a cohort of 250+ patients (see section 4).

## 2. Candidates, in the order to try them (sources: pages read 2026-11-19)

| # | Dataset | Access | What it offers | Needs a sponsor? | Unverified |
|---|---|---|---|---|---|
| A | **BraTS 2024 / BraTS-25 adult glioma post-treatment**, Synapse `syn59059776` (also `syn64377310`) | Free Synapse account + data terms | About 2,200 labelled post-treatment cases from 7 sites. Four sites are our spent cohorts (Duke 680, UCSF 600, Missouri 400, UCSD 350). **Unspent: Heidelberg about 300, Michigan about 100, Indiana about 70.** Labels 1-4 as in MU. (BraTS 2024 wiki, "Data Contributing Sites") | No | Whether files carry site tags, patient links across time and scan dates. Earlier read of the challenge paper found no dates. Ask the data contributors for a table of days between scans. |
| B | **TCIA "On Request" glioma collections**: Brain-Tumor-Progression (20 patients), QIN GBM Treatment Response (54), REMBRANDT (130) | TCIA access request (see each collection page) | Serial MRI. We label core with the frozen tool. Dates may be in clinical tables. | Probably no | Number of post-treatment pairs with both T1Gd and FLAIR and dates. Small. |
| C | **Burdenko-GBM-Progression** (TCIA, 180 patients, 645 studies, 4,956 series, MRI + CT + RT) | NIH Controlled Data Access (TCIA page, checked 2026-11-19) | Serial MRI after radiotherapy, with GTV. We label core with the frozen tool. | **Yes.** NIH rule: the applicant is a permanent faculty or senior researcher with an eRA Commons account; an institutional Signing Official co-signs; access lasts 1 year (NIH dbGaP access steps). A student cannot apply alone. | Whether T1Gd and FLAIR exist at every study and whether dates are in the clinical table. |
| D | **Heidelberg / DKFZ multicentre cohort** (Kickingereder et al., Lancet Oncology 2019; used as the hidden test set of BraTS-Lighthouse 2025 task 11): 306 patients, about 1,024 follow-up cases, MRI T1, T1c, T2, FLAIR | Closed. Only by collaboration with the organisers (Vollmuth, Kirchhoff, Maier-Hein) | The largest closed serial cohort I found. It matches the 250+ patient need. | Yes (data-sharing agreement) | Their labels are enhancing tumor and edema, not our core; our tool would need the images. Scan dates unknown. They may refuse: it is their hidden test set. |
| E | Private cohorts in method papers (XJYY-GBM, Oslo, Heidelberg private) | On request to the corresponding author | Unknown | Usually yes | Everything |

## 3. What you can do, in order
1. **Today, no sponsor:** make a free Synapse account; accept the BraTS data terms; tell me when the files are on disk (do not download more than the post-treatment folder, 2,200 cases is large). I check site tags, patient IDs and dates. If dates are missing, send the email in section 5a to the data contributors.
2. **Today, no sponsor:** send the TCIA requests for the three "On Request" collections (B).
3. **This week:** ask one faculty member (teacher, mentor, university contact) to be the sponsor for the NIH request (C). Draft text in section 5b.
4. **Later:** write to DKFZ (D), section 5c. Expect a low chance; ask anyway.

## 4. Method work that does not need new data (I can do this now)
- Use all spent cohorts (MU, UCSD, UCSF, LUMIERE, CFB, WAW) as development data. Leave-one-cohort-out across six cohorts. Candidate list written down before any run; one frozen choice before the new test.
- Reason: on the new test the best single method was Gompertz with absolute loss (1.119 vs persistence 1.154, `result.json`); VIG with squared-loss parts was not. A hybrid fitted with absolute loss and the Gompertz form is the first thing to try. This is a development step, not a claim.
- A new test run is only valid on a cohort that is not in this development pool.

## 5. Draft texts (fill the brackets; no personal details are written here)

**5a. To BraTS data contributors (Heidelberg, Michigan, Indiana):** "I am a student studying how glioma core volume changes between post-treatment MRI scans. I use the BraTS 2024 post-treatment training data from your site. Could you share a table with an anonymised patient ID and the number of days between each case's scans (no dates of birth, no names)? I would cite the dataset and you. Contact: [NAME, SCHOOL, EMAIL]."

**5b. To a faculty sponsor:** "I need an NIH data access request for the Burdenko-GBM-Progression collection (TCIA). It needs a faculty-level applicant with an eRA Commons account and an institutional Signing Official. The project is a locked, pre-registered test of a tumor-volume forecasting method on serial post-treatment MRI. Would you be the applicant? I will write the research-use statement and do all the analysis. Contact: [NAME]."

**5c. To DKFZ (Division of Medical Image Computing):** "We built a forecasting method for next-scan glioma core volume and tested it twice on locked outside data (both reports are in our repository: [LINK]). We would like to run it once on your multicentre longitudinal cohort under any agreement you require, with no right to publish without your approval. We need T1c and FLAIR images and the scan interval in days, or we can send you the frozen code to run yourselves and return only the error numbers. Contact: [NAME]."

Option in 5c that avoids sharing data: "send us the frozen code, you return the error numbers". This works for any closed cohort and is the most likely to be accepted.

## 6. Weakest points
- Sources: pages read on 2026-11-19 only; access rules and counts can change. Each "unverified" cell above is a check still to do.
- The BraTS option A may have no usable dates; then it gives nothing.
- A fresh test of the same VIG would be underpowered if the true effect is near zero. Only a changed method with a larger effect on development data justifies a new test.
- "Do not tune on the next test" stays in force: the method is frozen before the cohort is opened.
