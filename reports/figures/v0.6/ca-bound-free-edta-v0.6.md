# Figure legend — `ca-bound-free-edta-v0.6.png`

**Free vs. bound Ca²⁺ during an agonist-evoked transient without
extracellular Ca²⁺ (EDTA condition), v0.6 model.** Multi-panel diagnostic of
the endogenous-IP₃ model. After a 60 s resting baseline, agonist stimulation
evokes a self-limiting cytosolic Ca²⁺ transient. **Top:** free cytosolic Ca²⁺
rises from ~121 nM to a ~290 nM peak ~18 s after stimulus onset, then recovers
slowly as PMCA and NCX extrude Ca²⁺ (~193 nM by 360 s). With extracellular
Ca²⁺ removed, SOCE is inactive, so the dense tubular system empties during the
peak and does not refill. The cytosolic bound panel (CaM, gelsolin) and the
DTS multi-buffer panel (calreticulin C+P, HSP90B1 M+L, BiP, CREC) show buffer
occupancy rising and falling in step with the free pools. The bottom panel
confirms IP₃ is generated endogenously through the PI cycle (vs the Dolan 2014
reference curve), not forced.

*Source:* `runPlateletSim.py out/relock_validation --length 360
--agonist-delay 60 --ca-ex-mM 0 --seed 0`, rendered with `plotCaBoundFree.py
--ip3-stim-onset 60 --out reports/figures/v0.6/ca-bound-free-edta-v0.6.png`.
Numbers read from the simOut CalciumTrace `ca_cyt_nM` column (361 rows):
rest = c[60] = 120.7 nM, peak = max = 289.8 nM (argmax at t = 78 s),
end = c[360] = 192.6 nM.

*Provenance note:* regenerates `v0.5/ca-bound-free-edta-2026-06-10.png` under
the current v0.6 model. The EDTA transient is modestly damped vs v0.5 (peak
313 → 290 nM, −7 %; end 220 → 193 nM, −12 %; resting 106 → 121 nM, +14 %),
consistent with the biology added since v0.5 (PKC P2Y1-desensitisation +
PLCβ-phosphorylation brakes, the v0.61 downstream loops, and the v0.7
cAMP / PKA inhibitory axis braking IP3R) plus the small fixed-point shift from
the added ODE states. The qualitative result — a single store-limited
transient with no SOCE refill — is unchanged. The v0.6 EDTA peak (~290 nM,
360 s / 60 s-settle protocol) now agrees closely with the 200 s
acceptance-run EDTA peak (~296 nM, `test_acceptance::TestDolanTransient`).
