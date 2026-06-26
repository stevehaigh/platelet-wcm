# Sveshnikova & Panteleev (2025) — Systems-biology review of the platelet activation cascade

- **Full citation:** Sveshnikova A.N., Panteleev M.A. "Signal Transduction and Transformation by the Platelet Activation Cascade: Systems Biology Insights." *Hämostaseologie* 45(1):49–62, 2025. DOI 10.1055/a-2486-6758.
- **Type:** review (re-analysis of an existing computational model)
- **Local PDF:** source-info/calcium-papers/Sveshnikova and Panteleev - 2025 - Signal transduction and transformation by the platelet activation cascade Systems Biology Insights.pdf
- **Used in the model for:** source of compartment volumes, resting initial conditions, and the four-compartment framing; comparator for stochastic-bottleneck and signal-transformation arguments.

## What the paper does
A state-of-the-art mini-review that walks through one mechanism-based ODE model of SFLLRN (PAR1-agonist)-induced procoagulant platelet formation, tracking the signal from PAR1 binding down to mitochondrial permeability transition pore (mPTP) opening. It analyses how the cascade filters, amplifies, and transforms the signal shape, how calcium oscillations encode and are decoded, and where stochasticity dominates. It also demonstrates model-reduction techniques.

## Key findings / values the model uses
- Four compartments: cytosol (3 fL), plasmatic membrane equivalent (0.6 pL), DTS (1.5 fL), mitochondria (0.3 fL); total platelet volume ~5 fL (Table 1, p.52).
- 29 species, 27 differential equations, 2 fixed (external thrombin, external calcium) (Table 1, p.52).
- Resting initial conditions (Table 1): cytosolic Ca²⁺ 0.013 µM; IP3 0.05 µM; [Ca²⁺]dts 1000 µM; mitochondrial Ca²⁺ 0.1 µM; mitochondrial membrane potential 138.7 mV; PIP2 200 µM; resting PAR 0.006 µM; external Ca²⁺ 2000 µM.
- Stochastic bottleneck: Gq complexes (especially Gq-bound PLC) never exceed ~1 molecule/cell at any time (Fig. 3, p.55), so all downstream stochasticity originates there. PAR1 >1000 copies but active PAR1 never exceeds ~10/cell.
- Signal-shape transformation: step (thrombin) → peak (PAR1-Gq-GDP and intermediates) → oscillation (cytosolic Ca²⁺) → peak (mitochondrial Ca²⁺) → threshold (mPTP) (Fig. 2). IP3 oscillatory threshold ≈ 0.16 µM (p.57).
- Mitochondria act as an integrator: max mitochondrial Ca²⁺ rises with the time-integral of cytosolic Ca²⁺.
- Full reaction set with rate constants given in Table 2 (pp.53–54).

## Relevance to the platelet WCM
Direct provenance for compartment volumes and several resting initial conditions documented in `reports/data/calcium-data-provenance.md`; the DTS volume (1.5 fL ≈ 50% of cytosol) and [Ca²⁺]dts (1000 µM) are flagged outliers versus Purvis/Dolan. The signal-transformation and stochastic-bottleneck framing motivates the WCM's integer-count regime caveats.

## Caveats / notes
A review, not a new model — the analysed model is Sveshnikova 2015 (ref 37). Several IC/volume values diverge sharply from the Purvis/Dolan consensus the WCM otherwise uses. Table 2 rate constants were read from the PDF but not transcribed here.
