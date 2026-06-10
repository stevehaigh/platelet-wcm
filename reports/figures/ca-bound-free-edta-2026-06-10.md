# Figure legend — `ca-bound-free-edta-2026-06-10.png`

**Free vs. bound Ca²⁺ during an agonist-evoked transient without
extracellular Ca²⁺ (EDTA condition).** Multi-panel diagnostic of the current
endogenous-IP₃ model. After a 60 s resting baseline, agonist stimulation
evokes a self-limiting cytosolic Ca²⁺ transient. **Top:** free cytosolic
Ca²⁺ rises from ~106 nM to a ~313 nM peak within ~15 s, then recovers slowly
as PMCA and NCX extrude Ca²⁺ (~220 nM by 360 s). With extracellular Ca²⁺
removed, SOCE is inactive, so the dense tubular system empties during the
peak and does not refill. The cytosolic bound panel (CaM, gelsolin) and the
DTS multi-buffer panel (calreticulin C+P, HSP90B1 M+L, BiP, CREC) show buffer
occupancy rising and falling in step with the free pools — tight buffer–free
coupling. The bottom panel confirms IP₃ is generated endogenously through the
PI cycle (vs the Dolan 2014 reference curve), not forced.

*Source:* `runPlateletSim.py out/2026-06-10_cabf_edta --length 360
--agonist-delay 60 --ca-ex-mM 0 --seed 0`, rendered with
`plotCaBoundFree.py --ip3-stim-onset 60`.

*Provenance note:* replaces the earlier `ca-bound-free-v0.4.0.png`, which was
generated under the now-removed forced-IP₃ mode and showed a recover-within-
60 s + SOCE-refill profile that the endogenous-IP₃ model does not reproduce
(recovery requires EDTA → no SOCE refill; SOCE refill requires +Ca²⁺ →
sustained activation, no recovery). The old file is retained for the
historical design docs that reference it.
