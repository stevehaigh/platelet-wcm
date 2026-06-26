# Purvis & Lahav (2013) — Conceptual review: encoding/decoding information in signalling dynamics

- **Full citation:** Purvis J.E., Lahav G. "Encoding and Decoding Cellular Information through Signaling Dynamics." *Cell* 152(5):945–956, 2013. DOI 10.1016/j.cell.2013.02.005.
- **Type:** review (conceptual)
- **Local PDF:** source-info/calcium-papers/Purvis and Lahav - 2013 - Encoding and Decoding Cellular Information through Signaling Dynamics.pdf
- **Used in the model for:** conceptual framing — why the *temporal shape* of a signal (here cytosolic Ca²⁺) carries information, and why single-cell dynamics matter rather than population averages.

## What the paper does
Reviews the emerging principle that cells transmit information through the *dynamics* (amplitude, frequency, duration, delay, cumulative level) of signalling molecules, not just their static abundance. It surveys ERK, NF-κB, p53, Msn2, calcium, and others, covering how stimulus identity/strength is encoded in temporal patterns, how network structure shapes those patterns, and how downstream effectors decode them into distinct fates.

## Key findings / values the model uses
- Distinct stimuli produce distinct dynamics of the *same* molecule: EGF → transient ERK (proliferation) vs NGF → sustained ERK (differentiation); TNFα → oscillatory NF-κB vs LPS → sustained NF-κB; γ-radiation → pulsatile p53 (arrest) vs UV → sustained p53 (apoptosis) (Figs. 2–3).
- Population averages distort single-cell dynamics (e.g. p53 pulses appear as damped oscillations in bulk) — single-cell measurement is essential (Fig. 1D).
- Platelets cited directly: thrombin produces a characteristic temporal Ca²⁺ pattern; pretreatment with ADP attenuates the thrombin-induced pattern → cellular "memory" / pathway crosstalk (Chatterjee et al. 2010, p.948).
- Calcium frequency encodes specificity: different Ca²⁺ spike frequencies selectively activate NF-κB, NFAT, Oct/OAP; low/sustained Ca²⁺ favours high-affinity NFAT, transient bursts favour low-affinity JNK/NF-κB (Dolmetsch et al. 1997/1998).
- Decoding mechanisms: effector affinity/kinetics, feedforward loops (e.g. ERK → c-Fos persistence detector), and cumulative-signal (area-under-curve) integration.

## Relevance to the platelet WCM
Provides the conceptual justification for treating the WCM's cytosolic Ca²⁺ *trace* (amplitude, oscillation, duration) — not just peak value — as the meaningful readout, and for analysing single-cell dynamics. The explicit platelet thrombin/ADP example connects directly to the WCM's agonist-crosstalk and second-wave behaviour.

## Caveats / notes
A conceptual review with no platelet-specific parameters; quantitative content is qualitative/schematic (the figures are illustrative shape diagrams, not data). Use for framing, not for numerical provenance.
