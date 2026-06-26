# Ajanel et al. (2025) — MCU regulates ITAM-dependent platelet activation

- **Full citation:** Ajanel A, Andrianova I, Kowalczyk M, Menéndez-Pérez J, Bhatt SR, Portier I, Boone TC, Ballard-Kordeliski A, Chaudhuri D, Paul DS, Bergmeier W, Denorme F & Campbell RA. "Mitochondrial Calcium Uniporter Regulates ITAM-Dependent Platelet Activation." *Circulation Research* 137(4):474–492, 2025. DOI: 10.1161/CIRCRESAHA.125.326443.
- **Type:** experimental (platelet-specific Mcu⁻/⁻ mouse + human platelet pharmacology)
- **Local PDF:** source-info/calcium-papers/Ajanel et al. - 2025 - Mitochondrial Calcium Uniporter Regulates ITAM-Dependent Platelet Activation.pdf
- **Used in the model for:** validation reference for the MCU term / `mcu_vmax_scale` knockout knob (receptor-pathway specificity of MCU dependence).

## What the paper does
Generates platelet-specific Mcu-deficient mice (Mcu^pl⁻/⁻) versus littermate wild-type controls, and assesses mitochondrial Ca²⁺ flux and platelet activation in response to ITAM-based receptors (e.g. GPVI/CLEC-2) versus GPCRs (thrombin, ADP), plus in vivo hemostasis and thrombosis. Human platelets were also treated with MCU inhibitors.

## Key findings / values the model uses
- Mcu^pl⁻/⁻ platelets had **significantly reduced mitochondrial Ca²⁺ flux in response to ITAM-receptor activation**, whereas mitochondrial Ca²⁺ flux on **GPCR activation was unchanged**.
- Platelet **aggregation was significantly reduced by ITAM activation** but **GPCR-induced aggregation was unchanged**; MCU inhibition of human platelets gave similar (ITAM-specific) findings.
- In vivo: Mcu^pl⁻/⁻ mice had **reduced arterial thrombosis and less ischemic stroke injury**; hemostasis (tail bleeding) was only mildly altered.
- Mechanistically, mitochondrial ROS generation was **significantly reduced** after ITAM (not GPCR) activation; reduced p-Syk and p-PLCγ2. Inhibiting mitochondrial ROS reproduced reduced aggregation/ITAM signalling; inducing mitochondrial ROS (MitoParaquat) in Mcu^pl⁻/⁻ platelets **restored/increased** ITAM-dependent aggregation and signalling.
- Conclusion: mitochondrial Ca²⁺ flux regulates ITAM-dependent platelet activation **through generation of mitochondrial ROS**.

## Relevance to the platelet WCM
Experimental counterpart to the MCU term / `mcu_vmax_scale` knockout. Importantly it establishes that MCU dependence is **pathway-specific (ITAM, not GPCR)** and is mediated by **mitochondrial ROS**, neither of which the current calcium-only ODE represents.

## Caveats / notes
**Model-vs-experiment divergence (known finding):** the WCM's buffer-only MCU raises cytosolic Ca²⁺ on knockout, whereas this paper reports MCU loss **reduces** ITAM-driven platelet activation (via reduced mito-ROS). The WCM also drives primarily GPCR agonists (ADP/thrombin/ATP) — exactly the arm Ajanel finds MCU-independent — so the model has no ITAM/ROS axis to reproduce this result. Mouse platelet-specific KO + human pharmacology.
