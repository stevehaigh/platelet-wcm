# Purvis et al. (2008) — Molecular signaling model of platelet phosphoinositide and calcium regulation (P2Y1/ADP)

- **Full citation:** Purvis JE, Chatterjee MS, Brass LF, Diamond SL. "A molecular signaling model of platelet phosphoinositide and calcium regulation during homeostasis and P2Y1 activation." Blood 112(10):4069–4079, 2008.
- **Type:** primary computational model
- **Used in the model for:** the upstream GPCR → Gq → PLCβ → PI cycle → IP3 → Ca²⁺ cascade, resting concentrations/volumes, and the source (via its Table 1) of SERCA, IP3R, PKC and receptor rate constants.

## What the paper does
Builds a computational model of the human platelet integrating 77 reactions, 132 fixed kinetic rate constants and 70 species across five compartments (EX, cytosol, PM, DTS, DTS inner membrane). Four interlinked modules: (1) Ca²⁺ release/uptake (IP3R + SERCA), (2) phosphoinositide (PI) metabolism, (3) P2Y1 G-protein signaling (ADP → Gq), and (4) PKC regulation of PLCβ. Modules are tuned to resting concentrations then merged and fit to ADP-stimulated Ca²⁺ release data. A stochastic version reproduces asynchronous single-cell Ca²⁺ spiking.

## Key findings / values the model uses
- Cytosolic volume of the platelet ~6 fL (their ref 27); DTS estimated at 0.5–5% of cytosol (median 2%) by kinetic Monte Carlo, and 4.3% by direct glucose-6-phosphatase staining (DTS area vs cytosol).
- Resting [Ca²⁺]i constrained to 40–100 nM (homeostasis constraint, sampled at 100 ± 10 nM); model resting value ~75 nM with basal ADP. Measured resting [IP3] ~1200 molecules/platelet (~200 nM, their ref 31); the mild-response median estimate is ~1600 IP3R channels and median 750 IP3 molecules/cell.
- The model predicts a very low IP3R/SERCA ratio for functional platelets (e.g. ~1:5200 for high-response), reflecting that a single SERCA3b transports ~0.4 Ca²⁺/s vs a single type-2 IP3R conducting ~3000 Ca²⁺/s.
- Table 1 supplies the full rate-constant set the WCM draws on: SERCA cycle (Dode 2002), six-state IP3R (Sneyd & Dufour 2002, including L1, L3, L5 and the k/l constants), IP3R single-channel conductance 10 pS (Zschauer 1988) and leak 0.7 pS·m⁻², Nernst Po driving release (Eq. 1, Po=(0.9·IP3Ra+0.1·IP3Ro)⁴), P2Y1 ternary-complex receptor activation (Kinzer-Ursem 2007 / Waldo 2004), Gq exchange/GTPase cycle, PLCβ activation and PIP2 hydrolysis, IP3 degradation, and PKC activation/feedback.
- PKC negative feedback: activated PKC phosphorylates PLCβ, reducing its hydrolytic activity, providing the time-delayed (~10–15 s) shutoff of PIP2 hydrolysis after agonist stimulation (Fig. 2C, Signal attenuation).
- P2Y1 is low-abundance, ~150 copies/platelet (their ref 42).

## Relevance to the platelet WCM
Purvis 2008 Table 1 is the verified source for the bulk of the rate constants encoded in `reconstruction/platelet/dataclasses/process/calcium_signalling.py` (SERCA cycle, IP3R Sneyd ladder, receptor/Gq/PLCβ cascade, IP3 degradation, PKC feedback). It justifies the endogenous-IP3 upstream cascade (P2Y1 → Gq → PLCβ) that replaced Dolan's forced IP3 curve, and supplies the 6 fL cytosol and resting Ca²⁺/IP3 anchors. The PKC → PLCβ phosphorylation brake is the basis of the v0.6 PKC feedback.

## Caveats / notes
- Many Table 1 constants are themselves drawn from non-platelet primary sources (Dode 2002 SERCA, Sneyd & Dufour 2002 IP3R, Kinzer-Ursem/Waldo receptor kinetics) — Purvis is the secondary aggregator, verified against its own PDF.
- The provenance pass corrected a prior mis-attribution: the (KM1=0.5, KM2=1.0, kcat=8.9) triplet belongs to Purvis Reaction #11 (CDPDG synthesis), NOT PMCA — PMCA constants come from Caride 2007, not Purvis.
- Resting [IP3] here (~200 nM) is higher than Dolan's filtered 20–40 nM; the discrepancy is noted in the provenance doc.
- Purvis identifies the stochastic regime (low molecule counts) as biologically real; the WCM's deterministic ODE misses cell-to-cell variability.
