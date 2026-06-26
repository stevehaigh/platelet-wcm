# Ghatge et al. (2026) — MCU regulates Ca²⁺ dynamics, platelet function, bioenergetics and thrombosis

- **Full citation:** Ghatge M, Flora GD, Patel RB, Nayak MK, Kumskova M, Nguyen T, Usachev YM & Chauhan AK. "The mitochondrial calcium uniporter regulates calcium dynamics to drive platelet function, bioenergetics, and thrombosis." *Journal of Thrombosis and Haemostasis* 24:716–731, 2026.
- **Type:** experimental (MCU⁻/⁻ mouse + human platelet pharmacology)
- **Local PDF:** source-info/calcium-papers/Ghatge et al. - 2026 - The mitochondrial calcium uniporter regulates calcium dynamics to drive platelet function, bioenerge.pdf
- **Used in the model for:** validation reference for the MCU term / `mcu_vmax_scale` knockout knob (and the documented model-vs-experiment divergence).

## What the paper does
Compares wild-type and global MCU⁻/⁻ mice for arterial thrombosis (FeCl₃ carotid injury), measures mitochondrial and cytosolic Ca²⁺ in Rhod-2– and Fura-2–loaded platelets by fluorometry, and profiles bioenergetics with a Seahorse flux analyzer.

## Key findings / values the model uses
- Genetic ablation of MCU **inhibited agonist-induced platelet functions**: aggregation, fibrinogen binding to integrin αIIbβ3, granule secretion, and spreading on fibrinogen.
- MCU⁻/⁻ mice were **less susceptible to in vivo arterial thrombosis** with **unaltered tail bleeding time** (normal hemostasis).
- Mechanistically, MCU loss disrupted Ca²⁺ homeostasis via **reduced mitochondrial Ca²⁺ uptake, altered release of Ca²⁺ from the DTS, diminished cytosolic Ca²⁺, and impaired store-operated Ca²⁺ entry** in agonist-stimulated platelets.
- Ca²⁺-dependent GPVI signalling (PLCγ2 and PKC substrate phosphorylation) was significantly reduced in collagen-stimulated MCU⁻/⁻ platelets.
- Disrupted mitochondrial Ca²⁺ uptake **impaired mitochondrial respiration (OCR) and ATP production** in agonist-stimulated platelets.

## Relevance to the platelet WCM
Direct experimental counterpart to the model's MCU term and the `mcu_vmax_scale` = 0 knockout knob. The paper's finding that MCU loss **reduces** cytosolic Ca²⁺ and platelet function is the experimental anchor for the project's MCU work (#76) and the proposed bioenergetic MCU→ATP→SOCE coupling.

## Caveats / notes
**Model-vs-experiment divergence (known, treated as a finding not a failure):** the WCM models MCU as a buffer-only cytosolic Ca²⁺ sink, so knocking it out **raises** cytosolic Ca²⁺ — the opposite direction to Ghatge's measured **reduction**. The divergence arises because the model lacks the bioenergetic (ATP→SERCA/SOCE) coupling this paper demonstrates. Mouse genetics (global KO) and human pharmacology; effects are agonist-stimulated, not resting.
