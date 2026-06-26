# Kleppe et al. (2018) — NO/cGMP/cAMP inhibitory-signalling model in platelets

- **Full citation:** Kleppe R., Jonassen I., Døskeland S.O., Selheim F. "Mathematical Modelling of Nitric Oxide/Cyclic GMP/Cyclic AMP Signalling in Platelets." *International Journal of Molecular Sciences* 19(2):612, 2018. DOI 10.3390/ijms19020612.
- **Type:** primary computational model (steady-state ODE; one- and two-compartment)
- **Local PDF:** source-info/calcium-papers/Kleppe et al. - 2018 - Mathematical Modelling of Nitric OxideCyclic GMPCyclic AMP Signalling in Platelets.pdf
- **Used in the model for:** comparator for the v0.7 inhibitory axis (P2Y12/Gi/cAMP/PKA and the cGMP arm); source of PDE/PKG/PKA kinetics and cGMP–cAMP cross-talk mechanism.

## What the paper does
Mechanistic steady-state modelling (one- and two-compartment) of platelet cyclic-nucleotide signalling: NO → soluble guanylyl cyclase (sGC) → cGMP → PKG, and cGMP-mediated cross-talk to cAMP via phosphodiesterases PDE2 (cGMP-stimulated cAMP degradation) and PDE3 (cGMP-inhibited cAMP degradation). It quantifies how NO and PDE5 inhibition raise cGMP and indirectly cAMP/PKA, and why a two-compartment model is needed to explain NO-mediated PKA activation when bulk cAMP is unchanged.

## Key findings / values the model uses
- PKG expressed at very high levels in platelets: 7.3 µM (giving 14.6 µM cGMP binding sites), buffering free cGMP.
- Basal platelet cGMP ~0.45 µM (low micromolar); maximal NO-stimulated total cGMP ~7.6 µM; physiological NO range 0.1–5 nM.
- Dipyridamole (moderate PDE5 inhibitor, Ki = 0.7 µM, 1.0 µM dose) raised total cGMP by ~23% at maximal NO (experiment ~24%), but predicted >72% increase in PKG activation (free cGMP roughly doubled).
- PKG-mediated PDE5 phosphorylation lowers cGMP Kd at the GAF-A domain 130 → 30 nM but had little effect on the steady-state NO–cGMP curve.
- Kinetic parameters (Table 1): sGC NO binding kf = 300 µM⁻¹s⁻¹; PDE5 Vmax 39.0, Km 4.60 µM, Ki 0.70 µM (pPDE5 Vmax 117); PDE2 GAF Kd cG 2.0 µM, cA 25.0 µM; PDE3 Vmax 1.2, Km 150 nM, Ki(cGMP) 60 nM; PKG cGMP Kd 55 nM (high-affinity) / 750 nM (low-affinity); PKA R(cA)C Kd 2.90 µM / 1.00 µM.
- Global one-compartment model: total cAMP-PDE activity stays balanced across NO levels (PDE2/PDE3 offset), so it *cannot* explain NO-mediated PKA activation — only a two-compartment model with a PDE3-rich, PDE2-poor shape-change compartment (10% of platelet volume, 2× AC) reproduces it.
- Predicts weak PDE5 inhibitors (e.g. dipyridamole) selectively enhance cGMP–cAMP cross-talk — rationale for anti-platelet use.

## Relevance to the platelet WCM
The principal kinetic/quantitative comparator for the WCM's inhibitory axis (#10, v0.7 Slices): the cAMP/PKA brake and the cGMP module. PDE2/PDE3 cross-talk and the compartmentalisation requirement inform how a cAMP-raising arm (PGI2/forskolin/cilostazol) and its brakes should behave.

## Caveats / notes
Steady-state, not dynamic — focuses on equilibrium cyclic-nucleotide levels under constant NO. PKA isotype composition only partially resolved (modelled as Iα/IIα heterodimer). Some kinetic values sourced from non-platelet preparations (per Table 1 refs).
