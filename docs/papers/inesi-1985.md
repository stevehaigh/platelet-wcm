# Inesi (1985) — Mechanism of calcium transport (SERCA / SR Ca²⁺-ATPase)

- **Full citation:** Inesi G. "Mechanism of calcium transport." *Annual Review of Physiology* 47:573–601, 1985.
- **Type:** review (mechanistic / biochemical background)
- **Local PDF:** source-info/calcium-papers/Inesi - 1985 - Mechanism of calcium transport.pdf
- **Used in the model for:** mechanistic background for the E1/E2 SERCA enzymatic cycle used in the SERCA term.

## What the paper does
A review of the catalytic and transport mechanism of the sarcoplasmic reticulum (SR) Ca²⁺-ATPase, the canonical model for active cation transport and coupled enzyme catalysis. It synthesises structural, biochemical and kinetic evidence on how ATP hydrolysis is coupled to vectorial Ca²⁺ movement across the SR membrane.

## Key findings / values the model uses
- The SR Ca²⁺-ATPase is a single polypeptide of ~115,000 daltons with **one catalytic (phosphorylation) site per chain**.
- **2:1 stoichiometry** between Ca²⁺ transported and ATP hydrolysed (two Ca²⁺ per ATP).
- High-affinity Ca²⁺ sites used in transport have Kd ≈ 10⁻⁶ M; the ratio of Ca²⁺ sites to catalytic (phosphorylation) sites is **2:1**.
- Ca²⁺ binding is **cooperative** (Hill n_H > 3 at low proton concentrations), requiring interaction of at least four Ca²⁺ domains and involving a protein conformational change.
- A **sequential binding mechanism** is proposed: rapid Ca²⁺ binding to a first site triggers a conformational change that renders a second site available (E → E·Ca → E'·Ca → E'·Ca₂, Eq. 1). The enzyme cycles between E1 and E2 conformational states with a phosphoenzyme intermediate.

## Relevance to the platelet WCM
Provides the mechanistic justification for modelling SERCA as an E1/E2 enzymatic cycle (the implemented kinetics are the Dode 2002 / Purvis 2008 SERCA3b rate constants) rather than a lumped Michaelis–Menten pump, and confirms the 2 Ca²⁺ : 1 ATP transport stoichiometry used in the SERCA flux/ATP-cost accounting in `calcium_signalling.py`.

## Caveats / notes
This is skeletal-muscle SR ATPase, not the platelet SERCA2b/3 isoforms; it is background for the cycle topology and stoichiometry, not a source of platelet-specific rate constants (those come from Dode 2002 via Purvis 2008). A 1985 review predating the modern structural (E1/E2/E1P/E2P) crystallographic scheme.
