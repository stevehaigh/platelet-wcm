# Figure legend — `dose-sweep-transition-surface-auc-2026-06-10.png`

**The Ca²⁺ cascade is graded upstream but binary at the store-release
step.** Two-panel dose-response of the platelet Ca²⁺ system across a 9×9
grid of agonist concentrations (the "transition strip": ADP 10⁻⁴–1 µM ×
thrombin 10⁻⁶–10⁻³ nM), 200 s sims, seed 0. Both panels share one
colourblind-safe colormap (viridis). **A (left):** 3-D surface of peak
cytosolic Ca²⁺. Across the supra-threshold region the peak locks at a flat
~436 nM plateau (up to a ~461 nM ridge at intermediate thrombin),
independent of stimulus strength; sub-threshold cells (low ADP *and* low
thrombin, front-left) fall off a cliff toward the ~157 nM resting floor.
The DTS reservoir releases as an all-or-nothing bolus through IP3R, so peak
height is set by reservoir size and buffering, not by upstream signal
magnitude. **B (right):** heatmap of Ca²⁺ AUC above resting baseline
(time-integrated cytosolic Ca²⁺ over the 200 s window). Unlike the peak,
the integrated response is graded — ~8,000 to ~62,650 nM·s, an ~8× dynamic
range — and rises monotonically with both agonists, because sustained IP₃
holds the IP3R drain open longer. Thrombin and ADP are roughly
multiplicative; ADP alone drives nearly the full response above ~0.1 µM
(consistent with P2Y1's effective Kd ~0.5 µM). Deterministic ODE model, one
run per grid cell (no error bars; see Methods on statistics for
deterministic models).

*Source:* `runscripts/manual/plotDoseSweepThesisFigure.py`, which renders
both panels from the committed sweep data
`reports/figures/dose-sweep-9x9-transition.npz`. That NPZ was produced by
`runscripts/manual/runDoseSweep.py --grid 9 --length 200` over the
transition range (see `reports/lab-books/lab-book-2026-05-22-issue-44-45-dose-sweep.md`).

*Provenance note:* the committed per-observable PNGs
(`dose-sweep-9x9-transition-*.png`) carry captions baked into the image and
mix two colormaps; this thesis figure is re-rendered from the same data with
no in-image caption and a single colormap, per the figure guidance in
`reports/thesis/review-part-2.md` §8.
