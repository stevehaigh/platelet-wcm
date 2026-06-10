---
title: "Lab book — 2026-06-10: continuous-state integration fix + golden re-baseline"
---

# Lab book — 2026-06-10: removing the period-2 commit-grid ripple

## Context

While reviewing the freshly-generated PMCA perturbation figure (issue #53,
see `lab-book-2026-06-10-perturbation-figures.md`), a fine ripple appeared in
the cytosolic recovery tails (after ~80 s in the EDTA condition), growing more
visible at higher PMCA V_max. This entry diagnoses it as a numerical artifact
and fixes it at source.

## Diagnosis

Measured from the saved traces (`out/2026-06-10_perturbation/pmca.npz`):

- **Ripple period = 2.0 s exactly** — twice the 1-s outer timestep. A
  Nyquist-locked period-2 alternation (up one step, down the next), the
  signature of a discretisation limit cycle, not biology.
- **Amplitude ~8–9 nM, roughly constant across the V_max scan** (×0.25 → ×4),
  so *not* PMCA-driven. It only *looks* PMCA-dependent because the ×4 tail
  settles lower and flatter, making a constant wobble more visually prominent.
- The "~80 s onset" is not a real onset: the ripple is present throughout at
  constant amplitude but is invisible while the cytosol is still falling
  steeply from the peak; it surfaces once the recovery tail flattens.

**Root cause.** Each 1-s outer step the calcium ODE re-seeded from
`counts.astype(float)` — the *rounded* integer bulk-molecule counts — discarding
the fractional residual, then committed a `np.round`ed delta and carried the
remainder. The residual-carry kept the per-species *mean* unbiased (its original
purpose, see `lab-book-2026-05-12-pi-cycle-design.md` §"IP3 stuck-at-205") but
the re-seed from the rounded state injected a sub-ion perturbation each second.
In the near-marginally-stable recovery tail this sustained as a period-2 ring.

## Fix

`reconstruction/platelet/dataclasses/process/calcium_signalling.py`,
`molecules_to_next_time_step`: seed the integration from the *continuous* state
instead of the rounded counts —

```python
counts_f = counts.astype(float)
y0 = counts_f + self._residual          # was: y0 = counts.astype(float)
... integrate ...
fractional_delta = y_final - counts_f   # was: y_final - y0 + self._residual
delta = np.round(fractional_delta).astype(np.int64)
self._residual = fractional_delta - delta
```

This preserves `y_final == counts + residual` exactly across the integer commit,
so the ODE continues from where it actually was while BulkMolecules stays
integer-valued. Safe because `RestingDecay` excludes the calcium ODE species
(`_CALCIUM_SET = frozenset(MOLECULE_NAMES)`), so `CalciumDynamics` is the sole
mutator of those counts between steps.

## Verification

- Tail ripple (std of detrended cyt, 80–300 s): **8–9 nM → 0.08 nM** (99% gone)
  across all PMCA factors.
- PMCA recovery-tail AUC trend intact and monotone: **47,300 → 36,800 nM·s**.
- Golden-scenario drift < **0.6 %** on every validation observable
  (cyt, IP₃, SOCE; the 8 % on DTS is relative error on a near-zero denominator
  while the store empties — biologically nil).
- **54 tests pass, mypy clean.**

## Re-baseline (intended numeric change)

The fix intentionally changes the exact output, so the two byte-identical
golden NPZs (`models/platelet/tests/sim/golden/{at_rest_30s,default_activation_30s}.npz`)
were regenerated with `REGEN_GOLDEN=1`. The "frozen v0.4.1 byte output" is
superseded by the fixed numerics; the **biology is unchanged** (no rate
constants touched).

## Figures regenerated

- `perturbation-{pmca,mcu}-2026-06-10.png` (+ `.md` legend sidecars) — smooth
  recovery tails; results unchanged (PMCA AUC monotone; MCU peaks 651/436/380 nM).
- **`ca-bound-free` reconciled** — the draft cited `ca-bound-free-v0.4.0.png`,
  generated under the now-removed **forced-IP₃** mode. That figure's caption
  (recover-within-60 s + DTS refills via SOCE) cannot be reproduced on the
  endogenous-IP₃ model: recovery needs EDTA (→ no SOCE refill), while SOCE
  refill needs +Ca²⁺ (→ sustained activation, no recovery). Re-rendered as the
  **EDTA self-limiting transient** on the current fixed model →
  `ca-bound-free-edta-2026-06-10.png` (rest 106 → peak 313 nM → slow recovery;
  buffers track the free pools; endogenous IP₃ via the PI cycle). Draft caption
  rewritten; old `v0.4.0` file retained for the historical design docs that
  cite it.

## Files touched

- `reconstruction/platelet/dataclasses/process/calcium_signalling.py` — the fix.
- `models/platelet/tests/sim/golden/*.npz` — re-baselined.
- `reports/figures/perturbation-{pmca,mcu}-2026-06-10.{png,md}`,
  `reports/figures/ca-bound-free-edta-2026-06-10.{png,md}` — figures + legends.
- `reports/thesis/draft-thesis.qmd` — validation-figure caption + reference.
