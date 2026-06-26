# Dunster et al. (2015) — GPVI → Syk early signalling regulated by phosphatases

- **Full citation:** Dunster J.L., Mazet F., Fry M.J., Gibbins J.M., Tindall M.J. "Regulation of Early Steps of GPVI Signal Transduction by Phosphatases: A Systems Biology Approach." *PLoS Computational Biology* 11(11):e1004589, 2015. DOI 10.1371/journal.pcbi.1004589.
- **Type:** primary computational model (ODE, data-fitted)
- **Local PDF:** source-info/calcium-papers/Dunster et al. - 2015 - Regulation of Early Steps of GPVI Signal Transduction by Phosphatases A Systems Biology Approach.pdf
- **Used in the model for:** reference for the GPVI → Syk → PLCγ2 collagen pathway — a pathway the WCM has *not yet* built; source of platelet copy numbers and kinetics if/when GPVI is added.

## What the paper does
Develops a series of nonlinear ODE models (A, B, C) of the earliest steps of collagen-receptor GPVI signalling, focusing on Syk recruitment and autophosphorylation on Y525 (a surrogate for Syk activity) and how protein-tyrosine phosphatases regulate it. Models are fitted to high-temporal-resolution quantitative phosphorylation data (22 time points over 250 s, CRP-stimulated, one donor) and ranked by Akaike information criterion.

## Key findings / values the model uses
- A simple constitutively-active phosphatase (Model A) reproduces the Syk-Y525 steady state but not the early transient peak (Fig. 4).
- A specific negative-feedback pathway via c-Cbl (ubiquitin ligase) and the phosphatase TULA-2, acting through a second Syk site Y323 (Models B/C), is needed to capture the early peak. TULA-2 dephosphorylates Syk-Y525, returning the receptor to an inactive state.
- Copy numbers (Table 1): GPVI 5000, cytosolic Syk 2763, TULA-2 7800, c-Cbl 2581 molecules/platelet.
- Platelet volume Vp = 7.4×10⁻¹⁸ m³ (7.4 fL); extracellular volume per cell Ve = 3.3×10⁻⁹ m³ (Table 2).
- Rate constants (Table 2): ligand binding k1 = 8 m³·mol⁻¹·s⁻¹; ligand dissociation k₋₁ = 3.02×10⁻² s⁻¹; remaining rates fitted per model.
- Human platelets possess ≥18 protein-tyrosine phosphatases, >52,000 copies total per platelet (Introduction).
- Results show a clear separation between healthy and GPVI-deficient states in Syk-activation dynamics.

## Relevance to the platelet WCM
The authoritative platelet-specific kinetic/proteomic reference for a future GPVI/ITAM arm (collagen → Syk → PLCγ2 → IP3), which the WCM currently lacks (its PLC arm is GPCR/PLCβ only). Establishes phosphatase-controlled negative feedback as the rate-limiting control on activation onset.

## Caveats / notes
Models only the *early* steps up to Syk activity, not downstream PLCγ2/Ca²⁺. Forward-reaction fitting data from a single donor (validated against a wider population); secondary autocrine signals (ADP, TxA2, integrin) were pharmacologically suppressed.
