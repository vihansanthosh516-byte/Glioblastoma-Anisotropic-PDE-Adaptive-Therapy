# Track A few-patient designs: prior-art and data check (2026-11-24)

User direction (2026-11-24): build methods that work with few patients, for both the region and the drug part.

## Region: "label one section, map the rest of the same tumor"

Design: a pathologist-labelled Visium section trains a per-patient model; the model labels the patient's other sections (other sampling sites); evaluation against held-out labels on those sections, repeated per patient. Claim is within-patient only.

Prior-art check (web search by a helper agent, closest paper read directly in the browser):
- Sonpatki et al., Nat Commun 2026 (s41467-026-69716-2): anatomical-feature labels outside manual annotations are assigned "based on the shortest Euclidean distance to an AF-annotated grid within the same trajectory" (read in the article Methods). That is within-section distance propagation, not cross-section transfer evaluated on held-out labels. The "27 ROIs from 12 WSIs" are nuclei-detection training data, not region labels.
- Ravikumar et al. 2025 preprint: Ivy GAP niches by marker scoring (UCell, argmax), no transfer.
- Greenwald et al., Cell 2024: Ivy GAP programs on Visium; notes Ivy GAP misses a hypoxia-adjacent layer.
- Serial-section alignment tools (PASTE2, STalign, STAligner, moscot): generic, mostly assume adjacent or overlapping sections. PASTE2 handles partial overlap.
- No paper found doing one-section-to-other-sites pathology-label transfer within a GBM patient with held-out evaluation and per-spot confidence. Not verified: full texts of Ravi 2022 and other 2025-2026 papers; SLAT, CAST, SPACEL, Tangram, TIDEST not reached. This is a search, not proof of novelty.

Comparators to run on the same sections: alignment-based transfer (PASTE2), Ivy GAP marker scoring, GBM-Space niche abundances, a plain classifier trained on the one labelled section, and the GBM-Space cell-state abundances.

## Drug: within-tumor, cell-level TMZ response

- GSE214966 is baseline scRNA-seq for 4 patients (Lee et al., Nat Med 2024). Per-cell TMZ response exists only as pharmacoscopy scores in the paper's Supplementary Table 2 (27 patients, image-based, marker cell classes GBM / immune / other), not deposited in a repository. The per-cell drug readout is not linked to per-cell expression.
- No public dataset found with within-tumor TMZ response linked to expression per cell or per region (helper search; Eyler 2020 used dasatinib; slice-culture scRNA-seq used etoposide and panobinostat; Leipzig slice study has bulk RNA only).
- Assessment: a true few-patient, cell-level TMZ test is blocked by data. The closest usable input is the patient-level pharmacoscopy TMZ score for GBM cells (27 patients) with scRNA-seq for a subset, which is a patient-level design, not within-tumor.

Not verified: format of Supplementary Table 2; GSE214966 GEO record read directly (NCBI blocked the helper).
