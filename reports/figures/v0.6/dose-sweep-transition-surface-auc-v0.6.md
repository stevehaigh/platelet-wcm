# Figure legend — `dose-sweep-transition-surface-auc-v0.6.png`

**The Ca²⁺ cascade is graded upstream but binary at the store-release step
(v0.6 model).** Two-panel dose-response of the platelet Ca²⁺ system across a
9×9 grid of agonist concentrations (the "transition strip": ADP 10⁻⁴–1 µM ×
thrombin 10⁻⁶–10⁻³ nM), 200 s sims, seed 0. Both panels share one
colourblind-safe colormap (viridis). **A (left):** 3-D surface of peak
cytosolic Ca²⁺. Across the supra-threshold region the peak locks at a flat
~525 nM plateau (up to a ~588 nM ridge at low thrombin / intermediate ADP),
independent of stimulus strength (supra-threshold CV ≈ 3 %); sub-threshold
cells (low ADP *and* low thrombin) fall off a cliff toward the ~177 nM resting
floor. The DTS reservoir releases as an all-or-nothing bolus through IP3R, so
peak height is set by reservoir size and buffering, not by upstream signal
magnitude. **B (right):** heatmap of Ca²⁺ AUC above resting baseline
(time-integrated cytosolic Ca²⁺ over the 200 s window). Unlike the peak, the
integrated response is graded — ~13,550 to ~68,500 nM·s, an ~5× dynamic range
— and rises monotonically with both agonists. Deterministic ODE model, one run
per grid cell.

*Source:* `runDoseSweep.py out/relock_dosesweep --grid 9 --length 200
--adp-min 0.0001 --adp-max 1 --thrombin-min 1e-06 --thrombin-max 0.001`
produced the sweep NPZ (copied to
`reports/figures/v0.6/dose-sweep-9x9-transition.npz`); rendered by
`plotDoseSweepThesisFigure.py` (its `NPZ` / `OUT` constants repointed to v0.6).
Printed: peak Ca²⁺ 177 / 588 nM, AUC 13,552 / 68,488 nM·s, ADP 1e-4–1 µM,
thrombin 1e-6–1e-3 nM. Grid axes verified identical to the v0.5 transition NPZ.

*Provenance note:* regenerates
`v0.5/dose-sweep-transition-surface-auc-2026-06-10.png` under the current v0.6
model over the same grid. The qualitative result is preserved and sharper: the
supra-threshold peak plateau is *flatter* than v0.5 (CV 3 % vs 10 %). The
plateau shifts up (~436 → ~525 nM) and the AUC band up (~8,000–62,650 →
~13,550–68,500 nM·s; dynamic range 8× → 5×), consistent with the v0.61
autocrine loops (ADP, TXA₂) amplifying the supra-threshold response.
