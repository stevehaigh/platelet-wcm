# v0.7 figures — DTS realism (V_IM thermodynamic floor + γ_IP3R resting recalibration)

The v0.7 changes made the dense-tubular-system (DTS) Ca²⁺ store behave honestly:
`V_IM = 0 mV` (the passive IP₃R can no longer drain the store below cytosolic
free Ca²⁺) and `γ_IP3R 0.075 → 0.135 pS` (re-derived so an unstimulated cell
holds a stable resting fixed point). Each figure carries a full standalone
caption.

| Figure | Shows |
|--------|-------|
| `dts_thermodynamic_floor.png` | At fixed (recalibrated) γ, `V_IM = −60 mV` drains the store *below* cytosol (impossible for a passive channel); `V_IM = 0` bottoms it at the cytosolic equilibrium. Cytosolic transient unchanged. |
| `resting_homeostasis.png` | An unstimulated cell drifts under the old params (DTS → ~500 µM, cyt → ~6 nM); re-deriving γ for `V_IM = 0` restores a stable fixed point (DTS 250→255 µM, cyt 100→96 nM). |
| `dts_calibration_sweeps.png` | Why these values: `V_IM` sets the depletion *floor* only (all V_IM pass Dolan 5/5); γ trades resting stability against active residual, so a held residual is unachievable via IP₃R conductance. |
| `dose_response_depletion.png` | The autocrine ADP loop is an all-or-none commitment switch (store empties at every dose); the loops-off machinery is dose-graded (residual ~44 µM at threshold → ~0 at saturation). |
| `dolan_fig4_5of5.png` | Dolan & Diamond 2014 Fig. 4 acceptance — **5/5 criteria still pass** (peak 430 nM +Ca / 292 nM EDTA, SOCE Δ 138 nM): the DTS-realism changes don't regress the validated cytosolic response. |

**Reproduce:** `PYTHONPATH=$PWD python runscripts/manual/plotDtsRealism.py`
(figures 1–4) and `PYTHONPATH=$PWD python runscripts/manual/runPhase3.py
out/phase3_v07 --length 200` (figure 5). `plotDtsRealism.py` sets the module
constants `V_IM_V` / `GAMMA_IP3R_S` around each run (save/restore) to produce the
old-vs-new comparisons — the only place those are varied.
