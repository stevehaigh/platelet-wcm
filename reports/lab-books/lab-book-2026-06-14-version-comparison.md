---
title: "Lab book — 2026-06-14: v0.5 → v0.6 → v0.61 comparison doc, figure snapshots, PR #54 re-scope"
---

# Lab book — 2026-06-14: version-comparison write-up + PR/issue housekeeping

## Goal

Produce a single, figure-backed summary of how the model changed across the
three releases on the `PKC-P2Y1-desensitisation` branch (v0.5 baseline → v0.6
brakes → v0.61 amplifiers), verify the existing v0.6/v0.61 lab notes are
accurate, then push the branch and turn it into a reviewable PR with annotated
issues. No model code changed this session — this is documentation + figures +
git/issue housekeeping.

## Verification (lab notes are accurate)

Re-ran on the pinned 3.11.5 toolchain (`OPENBLAS_NUM_THREADS=1`):

- **Tests: 68 passed** — matches the 2026-06-13 lab book.
- **mypy: clean** ("no issues found in 57 source files").
- **Second-wave headline numbers reproduced exactly** via
  `runSecondWave.py --adp-uM 0.5 --length 300` (from `second_wave_summary.json`):

  | condition | peak | Ca²⁺ @300 s | IP₃ @300 s | P2Y1 desens max |
  |---|---|---|---|---|
  | v0.6 (open loop) | 313.8 nM | 117.3 nM | 51.8 nM | 0.19 |
  | aspirin (ADP loop on) | 314.1 nM | 210.1 nM | 189.6 nM | 0.66 |
  | full v0.61 | 314.1 nM | 210.1 nM | 254.9 nM | 0.67 |

  Sustained gap full − v0.6 = **+92.7 nM** (lab book said "+93 nM, ~80%"). ✓
  aspirin Ca²⁺ identical to full (210.1 nM) → the second wave is carried by the
  **autocrine ADP arm**; TXA₂ shows only in IP₃ (full > aspirin). ✓

Conclusion: the v0.6 and v0.61 lab notes are accurate; no corrections needed.

## Figures regenerated + snapshotted

All figures regenerated from the current model via the committed scripts (no
monkeypatching — they use the supported override knobs) and snapshotted for the
thesis-committed copy:

- `reports/figures/v0.6/` — `pkc_traces.png`, `plcb_traces.png`,
  `brake_effect_on_ca.png`, `why_brake_invisible.png`
  (from `runPerturbation.py` + `plotStoreLimitedFeedbacks.py`).
- `reports/figures/v0.61/` — `second_wave_traces.png`, `secretion_release.png`,
  `autocrine_adp_loop.png`, `thromboxane_loop.png`, `amplifiers_saturating.png`
  (from `runSecondWave.py` + `plotDownstreamModules.py` +
  `plotStoreLimitedFeedbacks.py`).

(Previously these scripts wrote to `out/figures/`, gitignored. v0.5 figures were
already snapshotted under `reports/figures/v0.5/`.)

Spot-checked `second_wave_traces.png` and `pkc_traces.png` render correctly —
detailed legends + takeaway captions intact, mathtext for chemical formulae.
The regenerated `pkc_traces.png` shows the P2Y1 desensitised fraction peaking
~0.8 (the v0.6 PR body's "~0.75" was approximate; ~0.8 is the current value).

## Comparison doc

`reports/design/version-comparison-v0.5-v0.6-v0.61-2026-06-14.qmd` — biology
added / model changes / results per release, embedding the snapshotted figures
with standalone legends. Added to the top of the `reports/design/README.md`
"current canonical" table (it now supersedes `model-status-2026-05-13` as the
most-current single summary; model-status remains the best pre-PKC calcium-core
reference). This doc is also the PR #54 description.

## Git / PR

The branch already had an open PR — **#54**, scoped "v0.6 only" — because the
v0.61 commits were never pushed (the 2026-06-13 lab book noted "PR #54 still
covers v0.6 only"). GitHub allows only one open PR per head→base pair, so rather
than create a duplicate I **pushed the branch (auto-updating #54) and re-scoped
#54 to cover v0.6 + v0.61**, with the comparison doc as the body.

## Issues annotated

v0.61 implements work tracked by several issues. All seven were **annotated
(commented)**; **none were closed** — the auto-mode classifier (correctly)
flagged closing as beyond the "annotate" the task asked for, so closing is left
to the owner.

- **#33** (thromboxane A₂ / TP receptor — positive-feedback Gq arm) — delivered
  by Slice B. Commented; **recommend closing** (left open).
- **#6** (GranuleRelease — Ca²⁺-triggered exocytosis) — delivered by
  `GranuleSecretion` (deterministic bulk route, not the stochastic per-granule
  state machine). Commented; **recommend closing** (left open).
- **#7** (GranuleState listener — granule inventory / cargo secretion) —
  delivered by `SecretionTrace` (released fractions/counts, not per-state
  populations). Commented; **recommend closing** (left open).
- **#8** (secretion-kinetics plot) — figure exists as
  `plotDownstreamModules.py secretion_release`; **open** (no `analysis/`-framework
  plot was added).
- **#5** (model granules as `UniqueMolecule` instances) — **deviation**: v0.61
  used bulk-count granule pools, not `UniqueMolecule`. Commented; **open** as a
  possible future refactor.
- **#10** (P2Y₁₂ ADP amplification) — the autocrine ADP amplifier landed on
  **P2Y1**, not the Gi-coupled P2Y₁₂. Commented; **open**.
- **#41** (spatial Ca²⁺ microdomains, labelled "v0.6 scope") — clarified that
  v0.6 shipped as PKC feedback; microdomains remain deferred. Commented; **open**.

## Follow-up — validation tests + Fig 8 fix (same session)

Prompted by the question "is Dolan 5/5 still valid and sufficient?". Finding:
the five implemented criteria are all **peak / SOCE-differential** over a 30 s
window — a fine regression guardrail for the calcium core, but they barely touch
the v0.6/v0.61 biology (the peak is store-limited and precedes the ~15 s PKC
delay) and they don't score the recovery phase.

Probed the model to calibrate bands (no failing/flaky assertions):
P2Y1 desens ~0.69 @120 s (ADP-only); PLCβ IP₃ base/ko end = 224/298 nM
(~25% damping); secretion serotonin 0.97 > fibrinogen 0.84 (dense leads α — but
**ADP released-frac 0.67 < fibrinogen** because ecto-NTPDase clears ADP[e], so
serotonin is the clean dense readout); +Ca plateau ~430 nM (**above** Dolan's
200–275 nM band — recovery phase not calibrated); EDTA stays ~240–280 nM (also
above Dolan's return-to-baseline).

New tests `models/platelet/tests/sim/test_validation_targets.py` (6 pass +
1 strict xfail): PKC desens ≥ 0.5; PLCβ brake holds IP₃ ≥ 15% below unbraked;
dense (serotonin) leads α (fibrinogen) + substantial + P-selectin exposed;
second-wave gap in [40,160] nM with ADP-carries-Ca / TXA₂-adds-IP₃
decomposition; +Ca plateau sustained-active (>200 nM) and SOCE-dependent
(plus ≫ EDTA). The unmet Dolan plateau band is a **strict xfail** so a future
recovery-phase fix is flagged.

Fig 8 (`thromboxane_loop.png`) had an inset axes (graph-within-a-graph) for
ΔIP₃ — promoted to a clean third panel; figure regenerated, doc caption updated.

## Next

- Recovery-phase calibration to the Dolan plateau band is now an explicit
  (xfailed) target — the clearest validation gap surfaced this session.
- Consider a v0.61 tag/checkpoint once PR #54 merges (mirror the v0.5 snapshot).
- Deferred from v0.61: a `runPerturbation.py` aspirin experiment (needs the
  runner generalised to override the scalar `COX1_FACTOR`); integrin §3
  (per-cell PAC-1 affinity state only — aggregation is inter-cellular).
