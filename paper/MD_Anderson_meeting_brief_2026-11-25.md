# My project: computational glioblastoma research

I'm a high school student building and testing computational methods for glioblastoma. The project has three tracks and an integration step. For every test, I write the plan down and lock it before I look at the data, then run it once.

**Track A: tumor tissue.** I use gene expression from tumor tissue to predict two things: whether a piece of tissue is tumor core or infiltrating edge (as labelled by a pathologist), and how a tumor responds to temozolomide.

**Track B: MRI.** I forecast how a tumor's volume changes after surgery from a patient's MRI scans over time.

**Track C: treatment simulation.** I simulate tumor growth and test whether an adaptive, uncertainty-aware treatment schedule does better than fixed schedules. This is simulation only.

**Integration.** I test whether MRI features carry real biological information, and use those measured links to set how uncertain the treatment simulation should be.

---

## Questions

1. Does MD Anderson have, or know of, glioblastoma spatial transcriptomics or multi-region samples with pathologist labels for tumor core vs infiltrating edge that a student could access, through a collaboration or a data use agreement?
2. Is there any dataset where temozolomide response is measured with gene expression in patient tumors, ideally with a no-TMZ comparison or per-cell response (for example, ex vivo drug testing)?
3. Are there longitudinal post-operative MRI cohorts (several scans per patient, with tumor masks) that could be shared or analyzed under supervision?
4. Would a faculty member be willing to sponsor a dbGaP access request, or advise on how students usually do this?
5. For core vs edge, is the infiltrating edge as labelled on H&E the right target, or would you trust a different reference (for example, 5-ALA fluorescence or MRI-guided biopsies)?
6. What would make work like this useful to clinicians or researchers at MD Anderson?
7. Is there a way to get feedback or mentorship on this work, or a summer program that fits it?
