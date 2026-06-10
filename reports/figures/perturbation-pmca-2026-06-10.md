# Figure legend — `perturbation-pmca-2026-06-10.png`

**PMCA V_max rate-limits cytosolic Ca²⁺ recovery in the absence of
extracellular Ca²⁺.** Single-mechanism perturbation of the basal PMCA
turnover number (`K_PMCA['k_cat']`), scanned across ×0.25–×4 of baseline in
the no-extracellular-Ca²⁺ (EDTA) condition. With extracellular Ca²⁺ removed,
SOCE cannot refill the store, so the dense tubular system releases a single
bolus and the plasma-membrane pumps clear the cytosol — a self-limiting
transient whose recovery rate PMCA sets. **Left:** cytosolic free Ca²⁺ over
300 s, colour-graded by V_max factor. All conditions peak near 317 nM
(313–319 nM); the post-peak recovery tail clears faster at higher V_max.
**Right:** recovery-tail AUC — cytosolic Ca²⁺ integrated above the resting
baseline over the post-peak window — falls monotonically with PMCA V_max
(47,300 → 36,800 nM·s across the scan), identifying PMCA as the
rate-limiting clearance step in this condition. Deterministic ODE model,
one run per condition (no error bars; see Methods on statistics for
deterministic models).

*Source:* `runscripts/manual/runPerturbation.py --experiment pmca --length 300`
(EDTA, `ca_ex = 0`). Recovery tails are smooth: the period-2 commit-grid
ripple present in earlier renders was removed by carrying the continuous ODE
state across the 1 s integer commit (`calcium_signalling.py`).
