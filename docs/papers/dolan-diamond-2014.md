# Dolan & Diamond (2014) — Systems model of platelet Ca²⁺ homeostasis, IP3-mediated release, and SOCE

- **Full citation:** Dolan AT, Diamond SL. "Systems Modeling of Ca²⁺ Homeostasis and Mobilization in Platelets Mediated by IP3 and Store-Operated Ca²⁺ Entry." Biophysical Journal 106(9):2049–2060, 2014. (Plus Supporting Material: four figures, one table.)
- **Type:** primary computational model
- **Local PDF:** source-info/calcium-papers/Dolan and Diamond - 2014 - Systems Modeling of Ca²⁺ Homeostasis and Mobilization in Platelets Mediated by IP3 and Store-Operate.pdf ; source-info/calcium-papers/Dolan 2014 SUpporting Material.pdf
- **Used in the model for:** headline validation target (Fig. 4 Ca²⁺ transient ±extracellular Ca²⁺) and source for the SOCE (STIM1/Orai1 MWC), IP3R, SERCA and PMCA sub-modules plus representative resting initial conditions.

## What the paper does
Builds a 34-species, 35-reaction ODE model of platelet calcium across five compartments (extracellular EX, cytosol CYT, dense tubular system DTS, plasma membrane PM, DTS inner membrane IM). Four molecular modules balance Ca²⁺: IP3R-mediated release from DTS, SERCA refill, calmodulin-regulated PMCA efflux, and store-operated Ca²⁺ entry (SOCE) via STIM1/Orai1 puncta. Rather than fix one initial condition, the authors do Monte Carlo sampling of a 12-dimensional space (6 protein copy numbers + 3 nonprotein ICs + VIM + two unknown puncta parameters Km, n) and filter for configurations meeting four constraints: resting steady-state Ca²⁺, stability in EDTA, IP3-responsiveness, and functional SOCE.

## Key findings / values the model uses
- Resting balance (main text): [Ca²⁺]cyt 40–100 nM, [Ca²⁺]dts 100–400 µM (~250 µM measured by Fluo-5N, their ref 45), [Ca²⁺]ex fixed 1.2 mM in scans (Fig. 3A).
- Filtered resting [IP3] constrained narrowly to 20–40 nM (<200 copies/platelet) — lower than the measured ~130 nM (Fig. 3B, Results).
- IP3R modeled with the six-state Sneyd & Dufour 2002 kinetics (Fig. 2B, "Sneyd, 2002"); channel open probability Po based on number of tetrameric channels in active conformations. SERCA modeled with the SERCA3b cycle of Dode 2002 (Fig. 2A). PMCA uses the Caride 2007 PMCA4b CaM-coupled scheme (Fig. 2C).
- SOCE: STIM1 dimerization + STIM2-in-puncta Hill function (Eq. 2, α=0.2 max puncta fraction, +0.01 offset) feeding a Monod-Wyman-Changeux Orai gating scheme (Eq. 3) with allosteric constants L, Ka, a, f from Hoover & Lewis 2011. SOC current via Nernst-driven Eq. 4; membrane RC via Eq. 5.
- DTS volume 1–10% of cytosol; cytosol ~6 fL; only n and Km were scanned as unknowns. <0.06% of 2.6M sampled configs satisfied all four constraints.
- Responsive/SOCE configs require high inner-membrane potential VIM > −70 mV (sampled −100 to −60 mV; main conclusion).
- With Ca²⁺ex, [Ca²⁺]cyt peaks higher and ~15 s later (delayed by slow SOCE inactivation, Fig. S3); without Ca²⁺ex (EDTA) it peaks within ~5 s and decays over ~1 min with no DTS refilling (Fig. 4). 10× IP3 saturates the response; stores plateau in EDTA (<5 µM [Ca²⁺]dts decline over 333 s).
- Supporting Material Table S1: a single representative filtered IC giving sub-species counts (e.g. SERCA 11,892 total; IP3R 1,328; PMCA 769; CaM 20,481; STIM 4,265; Orai 1,447), used as the model's resting initial condition.
- Supporting Material IP3 forcing function: cubic Hermite interpolant of measured ADP-stimulated [IP3] (Daniel et al. 1998 data, only to 30 s, extrapolated beyond); parametric Eq. S1 `[IP3]=[IP3]0·(a·t^0.6·e^(−bt) + ct/(t+1) + 1)` for dose-response, exponent 0.6, rising ~5-fold within 5 s and settling above baseline by 200 s (Fig. S2).

## Relevance to the platelet WCM
This is the canonical reference for the calcium core in `reconstruction/platelet/dataclasses/process/calcium_signalling.py` and the `CalciumDynamics` process. Its Fig. 4 ±Ca²⁺ex transient is the "Dolan 5/5" regression target; the EDTA (`--ca-ex-mM 0`) condition is the no-extracellular-Ca²⁺ arm. Table S1 supplies the representative resting initial condition; the SOCE MWC machinery, IP3R/SERCA/PMCA sub-module choices, and VIM/membrane-potential treatment all trace here.

## Caveats / notes
- The model uses a forced exogenous [IP3] curve; the WCM later replaced this with an endogenous IP3 state driven by the GPCR cascade (Purvis-style), so the Dolan IP3 forcing is the earlier approximation.
- Resting [IP3] in the filtered population (20–40 nM) is below the measured ~130 nM; the authors attribute this to the sampled SERCA:IP3R ratio constraints.
- VIM is experimentally unmeasured; the model's central prediction (VIM > −70 mV) is inferred, not measured.
- Table S1 is one representative configuration from many valid ones — not a unique solution. Sub-species counts <1/cell are listed as ~0.
