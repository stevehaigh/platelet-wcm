# Caride et al. (2007) — PMCA4a vs 4b: calmodulin binding and activation kinetics

- **Full citation:** Caride AJ, Filoteo AG, Penniston JT, Strehler EE. "The Plasma Membrane Ca²⁺ Pump Isoform 4a Differs from Isoform 4b in the Mechanism of Calmodulin Binding and Activation Kinetics: Implications for Ca²⁺ Signaling." Journal of Biological Chemistry 282(35):25640–25648, 2007.
- **Type:** experimental (kinetics) + supporting computational simulation
- **Local PDF:** source-info/calcium-papers/Caride et al. - 2007 - THE PLASMA MEMBRANE CA2+ PUMP ISOFORM 4A DIFFERS FROM ISOFORM 4B IN THE MECHANISM OF CALMODULIN BIND.pdf
- **Used in the model for:** the CaM-coupled PMCA4b plasma-membrane pump kinetics (the platelet-relevant isoform) — multi-state extrusion scheme and basal rate constants.

## What the paper does
Characterizes and contrasts the kinetics of calmodulin (CaM) binding and activation of two PMCA splice variants, PMCA4a and PMCA4b, using stopped-flow fluorescence of labeled CaM (TA-CaM) and ATPase-activity assays in baculovirus-expressed protein. It then runs a simulation of a hypothetical cellular Ca²⁺ spike to show how the different isoform kinetics shape the [Ca²⁺]cyt transient. PMCA4b is the dominant, best-studied isoform; PMCA4a is faster-activating and less inhibited at rest.

## Key findings / values the model uses
- The pump is described as existing in two conformations, "open" (high Ca²⁺-ATPase activity, higher CaM affinity) and "closed" (lower activity); CaM binding relieves C-tail autoinhibition (Fig. 4C kinetic scheme for binding; STEPs 1–6 in the cellular-simulation scheme).
- For the cellular Ca²⁺-spike simulation (Discussion / Table 3, "STEP 3" onward), the paper lays out the explicit reaction scheme the platelet model adopts: CaM-independent Ca²⁺ binding/transport (k4, k4r, k5), CaM Ca²⁺-loading (k6/k6r, k7/k7r), CaM-activated PMCA cycle (k8–k12), and coupling to SERCA (k13/k14). PMCA4b basal Ca²⁺-binding/turnover constants (from Caride 2007 Table 3, verified in the WCM provenance pass): k4 (Ca²⁺ binding) = 10 s⁻¹µM⁻¹, k4r (Ca²⁺ unbinding) = 50 s⁻¹, k5 (turnover) = 5.5 s⁻¹. (PMCA4a turnover k5 = 12 s⁻¹ — the alternate isoform.)
- Activation by CaM follows a Ca²⁺-dependent sigmoid with apparent K_Ca ~2 µM for both isoforms; maximal CaM-activation rate constant k_act ≈ 5.4×10⁵ s⁻¹M⁻¹ for PMCA4b vs ~3.1×10⁵ for 4a (Fig. 2C). Both isoforms reach comparable V_max; differences are in activation/inactivation rates and basal activity, not turnover number.
- Ca²⁺ strongly slows PMCA4b inactivation (CaM dissociation) — ~30-fold drop in k_inact from 0.2 to 2 µM Ca²⁺ (Fig. 3); PMCA4a is much less Ca²⁺-sensitive (~2.1-fold).
- Simulation conclusion (Fig. 8): PMCA isoform kinetics shape the Ca²⁺-spike decay; ~10× less PMCA4a is needed to reproduce the effect of PMCA4b on restoring cytosolic Ca²⁺.

## Relevance to the platelet WCM
PMCA4b is the platelet-abundant isoform, so its CaM-coupled scheme is what the WCM PMCA uses in `reconstruction/platelet/dataclasses/process/calcium_signalling.py` (the 5-state CaM-coupled PMCA noted in CLAUDE.md as "Caride 2007 Table 3"). The basal constants k4/k4r/k5 map to the model's (k_on, k_off, k_cat); the CaM-loading and CaM-activated-pump steps give the activation that sets the transient's decay phase, alongside CaM acting as a cytosolic Ca²⁺ buffer.

## Caveats / notes
- Table 3 / the full rate-constant set (k4…k14r) is in the supplementary material and the simulation scheme; the headline PMCA4b basal values k4=10, k4r=50, k5=5.5 are confirmed in the WCM provenance doc against the source table. The fitted constants k1–k_-5 in Table 1/Table 2 describe TA-CaM binding to the protein/peptides, distinct from the cellular-simulation rate set.
- Data are from baculovirus-expressed PMCA in Sf9 microsomes and a generic Chinese-hamster-ovary Ca²⁺-spike simulation, not platelets — adoption is by isoform identity.
- The earlier (pre-correction) WCM 2-state Michaelis-Menten reduction (k_on/k_off/k_cat with KM derived) is the WCM's own simplification of Caride's multi-state scheme, not Caride's formulation.
- An earlier provenance error attributed (KM1=0.5, KM2=1.0, kcat=8.9) to this paper; those numbers do not appear in Caride 2007 (they belong to Purvis Reaction #11).
