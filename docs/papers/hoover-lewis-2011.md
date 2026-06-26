# Hoover & Lewis (2011) — STIM1:Orai1 stoichiometry for CRAC channel trapping and gating

- **Full citation:** Hoover PJ & Lewis RS. "Stoichiometric requirements for trapping and gating of Ca²⁺ release-activated Ca²⁺ (CRAC) channels by stromal interaction molecule 1 (STIM1)." *PNAS* 108(32):13299–13304, 2011.
- **Type:** experimental (single-cell patch-clamp / fluorescence) + MWC modelling
- **Local PDF:** source-info/calcium-papers/Hoover and Lewis - 2011 - Stoichiometric requirements for trapping and gating of Ca2+ release-activated Ca2+ (CRAC) channels b.pdf
- **Used in the model for:** the Monod-Wyman-Changeux (MWC) allosteric constants for STIM1-driven Orai1/CRAC channel gating in the SOCE block.

## What the paper does
The authors co-expressed mCherry-STIM1 and Orai1-GFP at varying ratios in HEK 293 cells and measured CRAC current (I_CRAC) and the STIM1:Orai1 ratio at ER–PM junctions ("puncta") after store depletion. They quantified the stoichiometry required to (a) trap the tetrameric Orai1 channel at junctions and (b) gate it open, and fit the activation data with a modified MWC model.

## Key findings / values the model uses
- Binding of **1–2 STIM1 per channel** is sufficient to trap (immobilise) a tetrameric CRAC channel at ER–PM junctions; the minimum STIM:Orai ratio plateaus near ~0.3 (Fig. 1D).
- Channel **activation is a steep, highly nonlinear (bell-shaped) function** of the STIM1:Orai1 ratio; peak I_CRAC occurs at a ratio of ~2 STIMs:Orai, suggesting maximal activity requires ~8 STIM1 bound per tetrameric channel (Fig. 2A).
- Both activation and Ca²⁺-dependent inactivation depend on STIM:Orai stoichiometry in a similarly nonlinear way (Fig. 3).
- Data are described by an MWC model: STIM1 binds Orai with **negative cooperativity** and opens channels with **positive cooperativity** by stabilising the open state. Model parameters: opening equilibrium constant L, stabilisation factor f, binding cooperativity factor a, total STIM ligand S_total, association constant Ka; four sites per channel (one per Orai subunit) (Fig. 4, SI Methods).

## Relevance to the platelet WCM
This is the upstream source (via Dolan & Diamond 2014) for the MWC channel-gating layer in the SOCE block of `reconstruction/platelet/dataclasses/process/calcium_signalling.py`. The code uses Hoover values L = 10⁻⁴, f = 14.2, a = 0.5 for the per-channel STIM-dimer occupancy → open-probability calculation; Dolan scans only n and KM, taking L, Ka, α, f from this paper.

## Caveats / notes
Measurements are in HEK 293 cells, not platelets, and use over-expressed protein; absolute Ka values are in arbitrary fluorescence units and were rescaled in the model to platelet dimer counts. The MWC fit is a phenomenological description, not a structural mechanism.
