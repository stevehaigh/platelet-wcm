---
title: "Lab book — 2026-05-22: ip3_forced removal + 2-D dose sweep + observable design"
---

# Lab book — 2026-05-22: issues #44 and #45 landed, GitHub Pages, observable redesign

## Context

Two-session day. Started on issue #44 (drop the leftover `ip3_forced`
flag and the Dolan Fig. S2 IP3-overlay path) which had been blocking
#45 (2-D agonist dose sweep). Both landed cleanly; the dose-sweep work
produced surprising biology that drove a small but meaningful redesign
of the harvested observables.

## v0.4.3 / #44 — drop the `ip3_forced` flag

Most of the rename to `agonist_forced` was already done in earlier
commits; the issue body was stale. The actual remaining scope after
auditing the codebase:

- Fix the import-time default-arg binding in `thrombin_nM`, `adp_uM`,
  `atp_ex_forcing_uM` (Robust Python §3 classic — `def f(x, peak=GLOBAL)`
  captures `GLOBAL` once at function-definition, so the dose-sweep
  override pattern `cs_mod.ADP_PEAK_UM = X` did nothing). Switched to
  `peak=None` sentinel pattern with live read inside the body.
- Delete dead code: `ip3_forcing_uM` and the v0.3 `gq_signal_uM` +
  `GQ_*` constants.
- After consultation, take the bigger refactor and drop the
  `agonist_forced` Boolean entirely — replace with three optional peak
  floats (`thrombin_peak_nM`, `adp_peak_uM`, `atp_ex_peak_uM`) threaded
  through `_ode_rhs` → `molecules_to_next_time_step` → `CalciumDynamics`
  → `run_platelet_sim` → CLI flags + webapp. CLI gets a `--at-rest`
  shorthand that zeros all three.
- Webapp config schema changes from `agonist_forced` to `at_rest`. The
  display-side legacy fallback in `runs.py` was dropped per discussion;
  the execution-side fallback in `jobs.py` was kept so a queued Resting
  job from before the deploy still runs at rest.

Phase 3 acceptance: 5/5 unchanged. 28/28 tests pass. mypy clean.

Committed as `8887305e v0.4.3 / #44`. Ultrareview caught three follow-ups
which I addressed in a follow-up pass (README + docs/README + CLAUDE.md
text I missed; preset-iff-rule docstring; IP3R "6-state Markov"
docstring leftover from v0.2.5).

**Memory saved**: `feedback_rename_grep_everywhere.md` — when renaming
a flag/attribute/symbol, grep the entire repo for the old name before
declaring done. I missed README and docs/README references because they
were in `git status` already and I assumed they were handled.

## v0.5.0 / #45 — 2-D dose-response sweep

Built `runscripts/manual/runDoseSweep.py` per the design in the issue
body. Service-layer architecture: a pure `run_dose_sweep()` function
that loops cells, calls `run_platelet_sim()` with per-cell peak
overrides, harvests scalar observables; `main()` is just argparse +
call + write artefacts. Per-cell sims pruned to `simOut/CalciumTrace/`
only by default to keep on-disk grid² × ~50 KB (vs ~5 MB unpruned).

Initial harvest: `peak_ca_nM`, `peak_ip3_nM`, `tpeak_ca_s`,
`tpeak_ip3_s`, `auc_ca_nMs`. Plus `--replot-only` flag for re-harvesting
+ re-plotting an existing sweep dir without rerunning sims.

### Surprise biology #1 — the cascade is binary in peak Ca²⁺

First reference sweep (5×5, issue defaults ADP 0.1–10 µM × thrombin
0.01–1 nM, 200 s): **peak Ca²⁺ essentially flat at ~436 nM across the
entire 100× × 100× grid**, total variation ~0.5 nM = integrator noise.
Tried wider ranges (1000× × 1000×, 4×4 OoM, 7-8×7 OoM); the picture
only got clearer.

The biology: once IP3 crosses the IP3R activation threshold, the DTS
reservoir (~250 µM × 0.258 fL ≈ 39 000 ions) releases as a bolus and
peak cytosolic Ca²⁺ is set by reservoir size and buffering, not by the
upstream signal magnitude. **The IP3R → DTS release is binary.** This
is the canonical store-operated cell behaviour and not a bug.

To get a graded surface I added IP3 + AUC observables, then re-harvested
the existing run via `--replot-only`. Findings:

- **Peak IP3** spans 50 → 305 nM cleanly across 5 OoM of agonist. The
  GPCR cascade itself is graded, multiplicative in (thrombin, ADP).
- **Ca²⁺ AUC** (above resting baseline, integrated over 200 s) recovers
  an 8× dynamic range, because sustained IP3 keeps the IP3R-mediated
  drain open and the elevated-Ca²⁺ plateau lasts longer.
- Story for the writeup: **the cascade is graded upstream, binary at the
  store-release step, and graded again in the integrated response**.

### Reference sweeps

Two reference sweeps committed to `reports/figures/v0.5/`:

| Name | Range | Story |
|---|---|---|
| `dose-sweep-9x9` (wide) | ADP 10⁻⁵ → 10³ µM × thrombin 10⁻⁶ → 10⁻¹ nM | Covers biological extremes; shows the saturation plateau and the sub-threshold floor |
| `dose-sweep-9x9-transition` | ADP 10⁻⁴ → 1 µM × thrombin 10⁻⁶ → 10⁻³ nM | Focused on the transition strip; resolves the cross-axis interaction (thrombin and ADP are roughly multiplicative; ADP alone can drive nearly the full response above ~0.1 µM, consistent with P2Y1's effective Kd ~ 0.5 µM) |

`v0.5.0 / #45` commit, then follow-ups for the transition sweep +
captioning + interactive Plotly preview.

### Surprise biology #2 — `argmax`-based time-to-peak is misleading

User-spotted: in the time-to-peak Ca²⁺ heatmap, high-ADP cells looked
like they had **slower** peaks than intermediate-ADP cells. Plotting
the actual traces along a constant-thrombin row revealed why: at high
ADP the Ca²⁺ rises to plateau in ~10 s (faster, as expected), then
**drifts upward by ~0.1% over the next 180 s** as SOCE/NCX equilibrate,
so `argmax(ca_cyt)` reports t ≈ 190 s. Pure metric artefact.

**Fix**: drop `tpeak_ca_s`/`tpeak_ip3_s`; replace with `t_to_thresh_ca_s`
(first time cyt Ca²⁺ ≥ 250 nM) and `t_to_thresh_ip3_s` (first time IP3
≥ 150 nM). Thresholds chosen to sit between resting baseline and
saturated plateau. Cells that never cross are right-censored at the
sim length (200 s = "did not fire"). After the swap the rise-time
heatmaps are cleanly monotonic in both axes.

**Memory saved**: `feedback_argmax_rise_time_artefact.md` — for any
future "how fast does X respond" observable in this repo, use a
threshold-crossing time, not argmax.

### 1-D dose-response line plots

Added `plot_line_along_adp()` — one curve per thrombin level, x = ADP,
y = observable, viridis colour-graded so the legend doubles as the
thrombin axis scale. Rendered for the three story-telling observables
(t_to_thresh_ca, t_to_thresh_ip3, auc_ca). The rise-time line plot
makes the multiplicative-cascade story obvious: thrombin sets the
floor (y-intercept at low ADP per curve), ADP brings the curve down to
a shared ~10 s asymptote at saturating ADP.

## Hosting / sharing — gist → GitHub Pages

Pushed the interactive Plotly HTML to a secret gist, gave the user the
htmlpreview URL, found that htmlpreview's CSP/sandbox blocks the Plotly
CDN script (no charts render). Enabled GitHub Pages on the repo (main
branch, root). Now any HTML committed to the repo is auto-served at
`https://stevehaigh.github.io/platelet-wcm/<path>` — no third-party
renderer in the chain. Gist deleted; the URL is now stable and tied to
git history.

## Observables, final list

```
peak_ca_nM         Peak cytosolic Ca²⁺ (binary — DTS-reservoir-determined)
peak_ip3_nM        Peak cytosolic IP3 (graded — cascade output)
t_to_thresh_ca_s   First time cyt Ca²⁺ ≥ 250 nM (rise time; censored at 200 s)
t_to_thresh_ip3_s  First time cyt IP3 ≥ 150 nM (rise time; censored at 200 s)
auc_ca_nMs         Ca²⁺ AUC above resting baseline (integrated response, nM·s)
```

Captions added to the peak Ca²⁺ and AUC PNGs (the two with the most-
loaded biology); the new rise-time observables get captions explaining
the censoring convention.

## Issues closed today

- #44 (Expunge `ip3_forced` + Dolan Fig. S2 IP3 overlay) — closed via
  `8887305e`. Acceptance criteria all met.
- #45 (2-D agonist dose sweep) — closed via `e12a0e23`. Reference 9×9
  figures committed; CLI driver, harvest, plot, interactive HTML, and
  smoke test all in.

## Issues filed

- #46 (Webapp Dose Response tab — Plotly 3-D surface + heatmap, live
  out/ scan). Deferred follow-up from #45; the CLI driver lays all the
  data groundwork, just needs the Dash plumbing.

## Files added / changed (high level)

- `runscripts/manual/runDoseSweep.py` (new, ~330 LoC)
- `runscripts/manual/plotDoseSweepInteractive.py` (new, ~150 LoC) —
  one-shot prototype of #46's webapp content
- `models/platelet/tests/sim/test_dose_sweep.py` (new, ~95 LoC)
- `reports/figures/v0.5/dose-sweep-9x9-*.png/.npz/.html` (15 reference
  artefacts across two sweeps)
- `reconstruction/platelet/dataclasses/process/calcium_signalling.py`
  (dead-code purge + agonist_forced → 3 peak floats; ~270 LoC net diff)
- `models/platelet/processes/calcium_dynamics.py` (class-attr swap)
- `runscripts/manual/{runPlateletSim,runPhase3,restConvergence}.py`
  (new API)
- `wholecell/webapp/{jobs.py, tabs/configure.py, tabs/runs.py}` (config
  schema + UI)
- `CLAUDE.md` (run-time conditions table; runDoseSweep listed)
- `README.md`, `docs/README.md` (post-rename cleanups caught by review)

## Memories saved

- `feedback_rename_grep_everywhere.md` — grep the whole repo on a
  rename, including non-code surfaces.
- `feedback_argmax_rise_time_artefact.md` — threshold-crossing time
  beats argmax for rise-time metrics whenever the trace is flat-topped.

## Next session — sequence

1. Webapp Dose Response tab (#46) — lift `plotDoseSweepInteractive.py`
   into a Dash callback with a sweep-dir dropdown and live refresh.
2. Data-out-of-code refactor — move rate constants, copy numbers, and
   compartment volumes out of `calcium_signalling.py` and
   `internal_state.py` into structured TSV / YAML in `raw_data/`. There
   should be an existing issue tracking this; check before scoping.
