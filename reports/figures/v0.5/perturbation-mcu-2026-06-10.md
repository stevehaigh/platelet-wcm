# Figure legend — `perturbation-mcu-2026-06-10.png`

**Mitochondrial Ca²⁺ uptake buffers the cytosol without rescuing the
store.** Single-mechanism perturbation of the MCU uniporter maximal uptake
rate (`K_MITO['V_max_MCU']`) under standard agonist stimulation with
extracellular Ca²⁺ present (1.2 mM): knockout (×0), baseline (×1), and
over-expression (×4). **Left:** cytosolic free Ca²⁺ over 400 s. MCU
progressively lowers the cytosolic peak (651 → 436 → 380 nM from knockout
to over-expression), consistent with the mitochondrion acting as a fast
cytosolic Ca²⁺ buffer during the transient; the knockout reproduces the
elevated-cytosolic-Ca²⁺ phenotype of MCU-deficient platelets
[cite: Ghatge 2026 — verify against Zotero]. **Right:** dense-tubular-system
free Ca²⁺ depletes to near zero within ~10 s in every condition, so MCU does
not refill the store — the buffering is cytosolic, not luminal.
Deterministic ODE model, one run per condition.

*Source:* `runscripts/manual/runPerturbation.py --experiment mcu --length 400`
(+Ca²⁺, `ca_ex = 1.2 mM`).
