# Sneyd & Dufour (2002) — Dynamic model of the type-2 IP3 receptor

- **Full citation:** Sneyd J, Dufour J-F. "A dynamic model of the type-2 inositol trisphosphate receptor." Proceedings of the National Academy of Sciences USA 99(4):2398–2403, 2002.
- **Type:** primary computational model (fit to experimental data)
- **Local PDF:** source-info/calcium-papers/Sneyd and Dufour - 2002 - A dynamic model of the type-2 inositol trisphosphate receptor.pdf
- **Used in the model for:** the IP3R type-2 gating kinetics (six-state ladder and its rate constants) used inside the calcium ODE.

## What the paper does
Constructs a kinetic model of the type-2 IP3 receptor and fits it to dynamic (time-course) and steady-state Ca²⁺-release data from type-2 IP3R. The model uses a six-state reduced scheme — receptor R, open O, activated A, shut S, and two inactivated states I1 and I2 — derived by assuming fast Ca²⁺-binding sub-states are in instantaneous equilibrium, which yields saturating (Michaelis-Menten-like) rather than mass-action Ca²⁺ kinetics. The receptor is treated as four independent identical subunits.

## Key findings / values the model uses
- Best-fit rate constants (Fig. 4 caption), which propagate (via Purvis 2008 Table 1) into the WCM IP3R ladder:
  - k1 = 0.64 s⁻¹µM⁻¹, k₋1 = 0.04 s⁻¹
  - k2 = 37.4 s⁻¹µM⁻¹, k₋2 = 1.4 s⁻¹
  - k3 = 0.11 s⁻¹ (printed as "0.11 s⁻¹µM⁻¹" in the caption — the µM⁻¹ is a unit typo; the φ3 formula is dimensionally consistent only with s⁻¹, and the body text gives φ3 ≈ 0.1 s⁻¹). k₋3 = 29.8 s⁻¹
  - k4 = 4 s⁻¹µM⁻¹, k₋4 = 0.54 s⁻¹
  - L1 = 0.12 µM, L3 = 0.025 µM, L5 = 54.7 µM
  - l₋2 = 1.7 s⁻¹, l4 = 1.7 s⁻¹µM⁻¹, l6 = 4707 s⁻¹, l₋2 = 0.8 s⁻¹, l₋4 = 2.5 µM⁻¹s⁻¹, l₋6 = 11.4 s⁻¹
- Channel open probability is (0.1·O + 0.9·A)⁴ (the 0.1/0.9 weights are not crucial; similar fits with nearby values).
- Mechanistic conclusions: Ca²⁺ binds with saturating not mass-action kinetics; Ca²⁺ decreases the rate of IP3 binding while increasing steady-state IP3 sensitivity; Ca²⁺-induced activation (O→A) is fast (~an order of magnitude faster than other transitions, φ4 dominant) while inactivation is slower; IP3 and Ca²⁺ binding are sequential (Ca²⁺ binds the activating site only after IP3 binds). Detailed balance enforced (K1=L1L2, K2=L3L4, K4=L5L6).
- Ca²⁺-independent inactivation by IP3 (O→S) has time constant ~10 s (φ3 ≈ 0.1 s⁻¹).

## Relevance to the platelet WCM
The IP3R in `reconstruction/platelet/dataclasses/process/calcium_signalling.py` uses this six-state type-2 scheme. Type-2 IP3R is the platelet-relevant isoform per Burkhart proteomics and Purvis/Dolan modeling choices. These rate constants are the most numerically load-bearing IP3R parameters in the model.

## Caveats / notes
- The paper was fit to hepatic-microsome / type-2 IP3R data, not platelet data; adoption into the platelet model is by isoform analogy (via Purvis/Dolan).
- The k3 unit in the Fig. 4 caption is a typo (s⁻¹µM⁻¹ → s⁻¹). The WCM provenance pass also flagged that Purvis 2008 Table 1 transcribed k3 as "11 s⁻¹" (a 100× error); the correct value is 0.11 s⁻¹, which collapses the resting "shut" state population to ≈0, matching Dolan Table S1.
- Open-probability exponent/weights are a modeling convenience, not a measured quantity.
- The authors stress steady-state data alone cannot characterize the receptor — dynamic fitting was essential.
