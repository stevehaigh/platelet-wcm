# Mazet et al. (2020) — Platelet PI-cycle ODE model; caution against mosaic data

- **Full citation:** Mazet F., Tindall M.J., Gibbins J.M., Fry M.J. "A model of the PI cycle reveals the regulating roles of lipid-binding proteins and pitfalls of using mosaic biological data." *Scientific Reports* 10:13244, 2020. DOI 10.1038/s41598-020-70215-7.
- **Type:** primary computational model (ODE, COPASI)
- **Local PDF:** source-info/calcium-papers/Mazet et al. - 2020 - A model of the PI cycle reveals the regulating roles of lipid-binding proteins and pitfalls of using.pdf
- **Used in the model for:** canonical reference for full PI-cycle (receptor → PLCβ → IP3) kinetics, intended to replace the WCM's forced-IP3 Dolan Fig. S2 approximation; cited as the warning against mosaic rate constants.

## What the paper does
Builds an ODE model of the complete phosphatidylinositol (PI) cycle in human platelets, parameterised *exclusively* from platelet-specific quantitative proteomics and time-resolved phospholipid/inositol-phosphate data taken under the same conditions (a single cell type). It is built iteratively via reduced models in COPASI (35 parameters). It then tests portability to other cell types (mouse platelet, HeLa) and shows that combining data from unrelated cell types ("mosaic" data) gives erroneous predictions.

## Key findings / values the model uses
- "Core model" focuses on GPCR → Gαq → PLCβ → IP3 (PI3K/PIP3 arm excluded as ~2 orders of magnitude lower).
- Lipid-binding proteins stabilise phospholipid levels: ~1.3×10⁶ PI45P2-binding proteins/cell predicted, matching ~1.12×10⁶ in the proteome (Supplementary Table S6). Removing them gives biologically unrealistic dynamics.
- Homeostatic levels: PI45P2 and PI4P each ~1.2–1.8×10⁶ molecules/cell; plasma-membrane PI ~6×10⁶/cell.
- GPCR number governs IP3 output (Fig. 4): RGq = 5000 (full complement, Thrombin) → peak IP3 ~8×10⁴/platelet; RGq = 1650 (TxA2+ADP) → intermediate; RGq = 150 (ADP alone) → low. Supports "strength of GPCR signalling ∝ receptor abundance."
- Saturating-ligand readouts (Fig. 4c): Ca²⁺ ADP 300 nM, TxA2 600 nM, Thrombin 1500 nM; IP3 ADP 3.2 µM, Thrombin 17.5 µM.
- Receptor copy numbers used (Fig. 2a): P2Y12 400, P2Y1 150, TPα 1500, PAR1 1276, PAR4 1100, Gαq 14800, Gαi 33700, PLCβ 5200.
- Predicts PI4K, PIP5K, OCRL1 are co-regulated (up/down) while SAC1 is regulated differently (likely spatially segregated).
- Mosaic-data demonstration (Figs. 5–6): swapping in HeLa/U2OS protein numbers produces inconsistent or non-signalling outputs.

## Relevance to the platelet WCM
The intended canonical source for an endogenous PI-cycle / IP3-production module (the WCM currently uses forced IP3). Its central methodological lesson — do not mix rate constants/concentrations from unrelated cell types — directly informs the WCM's provenance discipline.

## Caveats / notes
Excludes the PI3K/PIP3 arm and models IP3 production, not downstream Ca²⁺ release. Much detail (the 35 parameters, k1–k21) lives in Supplementary Tables not read here. Activation modelled with secondary GPCR signalling assumed immediate.
