# Sveshnikova et al. (2015) — Compartmentalised Ca²⁺ model with mitochondria; subpopulation formation

- **Full citation:** Sveshnikova A.N., Ataullakhanov F.I., Panteleev M.A. "Compartmentalized calcium signaling triggers subpopulation formation upon platelet activation through PAR1." *Molecular BioSystems* 11:1052–1060, 2015. DOI 10.1039/c4mb00667d.
- **Type:** primary computational model (COPASI; deterministic + stochastic)
- **Local PDF:** source-info/calcium-papers/Sveshnikova et al. - 2015 - Compartmentalized calcium signaling triggers subpopulation formation upon platelet activation throug.pdf
- **Used in the model for:** comparator for compartmentalised Ca²⁺ + mitochondrial uptake; source of the underlying model later reviewed in Sveshnikova & Panteleev 2025; oscillation-frequency and store-vs-extracellular-Ca²⁺ validation targets.

## What the paper does
Builds a mechanism-based ODE model of PAR1 (SFLLRN)-stimulated platelet calcium signalling across three intracellular compartments (cytosol, DTS, mitochondria) plus extracellular space, solved deterministically and stochastically in COPASI. It shows how a gradual increase in activation degree is converted into a stepped response hierarchy, ultimately splitting an initially homogeneous population into two subpopulations (procoagulant vs not) through mPTP opening.

## Key findings / values the model uses
- 30 species described as ODEs (Experimental section, p.1059).
- Three deterministic calcium patterns (Fig. 2 / S1): stationary low Ca²⁺ (<50 nM); oscillations at 0.3–0.4 Hz with amplitude 200–600 nM; sustained high Ca²⁺ prolonged ~20 s.
- Two decision-making mechanisms: (1) IP3R–DTS system sets the activation/spiking threshold and amplifies receptor-number variation; (2) mitochondrial Ca²⁺ accumulation integrates the cytosolic signal and triggers mPTP collapse in a fraction of cells.
- Sensitivity analysis (Table 1): most sensitive parameters are the SERCA2b Hill coefficient (4.9), IP3R activation constants l6/lm6 (3.8), IP3R equilibrium L5 (3.8); mitochondrial and SOCE parameters appear at the bottom.
- ~1000 PAR1 receptors/platelet; activated PAR1 numbers much smaller → stochastic decision-making.
- Extracellular-Ca²⁺ removal: cytosolic Ca²⁺ response decreased only ~2-fold, but procoagulant (mPTP) formation dramatically impaired (Fig. 5).
- Modules adapt platelet type-2 IP3R from Sneyd et al.; SERCA isotypes 2b and 3a via Hill equations; mPTP with non-linear matrix-Ca²⁺ dependence.

## Relevance to the platelet WCM
The canonical comparator for the WCM's mitochondrial-Ca²⁺ work (MCU, #76) and for compartmentalised store dynamics; oscillation frequency (0.3–0.4 Hz) and the store-stability-in-EDTA behaviour are listed validation targets in the provenance doc. Demonstrates why stochastic simulation matters at platelet copy numbers.

## Caveats / notes
Deliberately excludes PKC, PI3K, cAMP inhibition, and ADP/TxA2 positive-feedback loops (stated in Discussion). Uses external Ca²⁺ 2 mM and a large DTS volume relative to Purvis/Dolan. Stochastic results rely on COPASI's hybrid solver.
