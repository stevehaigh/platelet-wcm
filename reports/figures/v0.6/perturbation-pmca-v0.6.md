# Figure legend — `perturbation-pmca-v0.6.png`

**PMCA V_max rate-limits cytosolic Ca²⁺ recovery in the absence of
extracellular Ca²⁺ (v0.6 model).** Single-mechanism perturbation of the basal
PMCA turnover number (`K_PMCA['k_cat']`, baseline 5.5 s⁻¹), scanned across
×0.25–×4 of baseline in the no-extracellular-Ca²⁺ (EDTA) condition. With
extracellular Ca²⁺ removed, SOCE cannot refill the store, so the dense tubular
system releases a single bolus and the plasma-membrane pumps clear the cytosol
— a self-limiting transient whose recovery rate PMCA sets. **Left:** cytosolic
free Ca²⁺ over 300 s, colour-graded by V_max factor. All conditions peak near
295 nM (292–298 nM); the post-peak recovery tail clears faster at higher
V_max. **Right:** recovery-tail AUC — cytosolic Ca²⁺ integrated above the
resting baseline over the post-peak window — falls monotonically with PMCA
V_max (41,137 → 31,570 nM·s across the scan), identifying PMCA as the
rate-limiting clearance step in this condition. Deterministic ODE model, one
run per condition.

*Source:* `runPerturbation.py relock_pmca --experiment pmca --length 300`
(EDTA, `ca_ex = 0`); figure copied from `out/relock_pmca/pmca_traces.png`.
Numbers from `out/relock_pmca/pmca.npz`: peak per factor 291.7–297.5 nM;
recovery-tail AUC 41,137 (×0.25) → 40,221 → 38,333 → 35,812 → 31,570 (×4)
nM·s, monotonic; resting cyt 99.9 nM.

*Provenance note:* regenerates `v0.5/perturbation-pmca-2026-06-10.png` under
the current v0.6 model. The qualitative result is fully preserved (monotonic
AUC vs V_max, near-flat peak). Magnitudes shifted down vs v0.5 (peak band
313–319 → 292–298 nM; AUC endpoints 47,300 / 36,800 → 41,137 / 31,570 nM·s),
consistent with the v0.6 PKC / inhibitory additions lowering the evoked
transient.
