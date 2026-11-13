# Integration: literature on radiogenomic methods (2026-11-12)

Source: one literature subagent, 16 searches. Numbers from abstracts unless noted. UNVERIFIED items are marked.

## What the papers say

- No paper found that predicts Neftel cell states from MRI. Our target is new.
- Best-supported target: mesenchymal (MES) vs rest. Rathore 2019 (n=112, CV): MES sens 83.9 / spec 72.8. Grossmann 2016 (91 TCGA with RNA, likely our cohort): edema volume predicts subtype, AUC 0.61 only.
- Realistic CV performance at n ~100: AUC 0.6-0.75. No independent external MES AUC > 0.8 found.
- Location (SVZ contact) predicts survival (Mistry 2020, n=502) but not molecular subtype (Berendsen 2019; Mistry 2019, n=254).
- Normalization is the largest quantified pitfall: external AUC 0.19-0.45 naive vs 0.85-0.87 after N4 + z-score/WhiteStripe (Foltyn-Dumitru 2024). Use a normal-brain reference and a fixed 32-bin discretization (Carre 2020).
- Feature selection outside CV inflates results (Nuechterlein 2021).
- Piccolo & Vindigni 2026 (IvyGAP, n=28): 21 of 24 pathways R2cv <= 0. Underpowered.
- Open 3D encoders exist (MedicalNet, Models Genesis), but benchmarks show they do not beat simpler features (Falconnier 2026; Hollet 2026). Licenses UNVERIFIED.

## What this changes in our method search

1. Normalize to normal-appearing brain (median/IQR), not a whole-image z-score. Done in `src/integration_methoddev_features.py`.
2. Add texture (3D GLCM, 32 levels), location, habitats, peritumoral ring.
3. Keep all selection inside nested CV.
4. Expect at most rho ~0.2-0.35. The confirm pass bar (0.25) is near the top of what the field reports.

## Weakest point

Abstract-level reading. Cohort sizes conflict across sources for some papers.
