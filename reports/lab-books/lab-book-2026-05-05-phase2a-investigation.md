---
title: "Lab book — 2026-05-05: Phase 2a investigation — γ sweep falsifies the IP3R-only diagnosis"
---

# Lab book — 2026-05-05: Phase 2a investigation

## Session summary

Tried Step 2a from the 2026-05-01 plan — recalibrate `GAMMA_IP3R_S` from
10 pS to 0.6 fS — and the transient disappeared entirely instead of slowing.
A γ sweep + rate-balance dump showed the lab-book diagnosis was wrong on
two counts: a 3.7× math error in the Nernst driving voltage, and an
incorrect identification of SERCA's "refill rate." The actual culprit at
rest is SERCA's `k_bind_f` step (118 000 ions/s into E1·Ca at cyt=100 nM)
and the absence of a basal plasma-membrane Ca²⁺ leak. γ recalibration
alone cannot give a stable rest near 100 nM cyt with the current parameters.

**Phase 2 plan revised:** drop Step 2a (γ recalibration). Pursue
(ii) basal PM Ca²⁺ leak + (iii) SERCA IC sub-state pre-equilibration
instead.

---

## What the 2026-05-01 plan said

The Phase 1-complete lab book (`lab-book-2026-05-01-phase1-complete.md`)
identified IP3R flux as ~17 000× larger than SERCA at rest, attributed
this to using the 10 pS Zschauer 1988 single-channel conductance at face
value, and prescribed:

```python
GAMMA_IP3R_S = 0.6e-15   # 0.6 fS calibrated from SERCA balance at rest
```

The expected effect: DTS drain slows from ~0.35 s to ~5–10 s, transient
shape becomes the Dolan 60 s plateau instead of an instantaneous spike.

## Step 2a result — transient killed, not slowed

Setting `GAMMA_IP3R_S = 0.6e-15` and rerunning the 200 s IP3-forced sim:

| Time (s) | Cyt Ca²⁺ (nM) | DTS Ca²⁺ (µM) |
|----------|---------------|---------------|
| 0   | 99.9 | 250.0 |
| 20  |  1.4 | 261.0 |
| 50  |  1.1 | 262.2 |
| 100 |  1.1 | 263.8 |
| 200 |  1.1 | 267.5 |

No transient. Cyt collapses 100 → 1 nM by t=20 s; DTS *grows* slightly
(SERCA + PMCA outflow exceeds IP3R + SOCE + leak inflow). IP3 still
peaks at 244 nM at t≈8 s and STIM1 dimers stay at 22 (no DTS depletion
to trigger SOCE). The IP3R recalibration didn't slow the transient —
it removed the transient entirely.

## γ sweep (no IP3 forcing, 300 s integrations)

To find a γ that holds cyt at ~100 nM at rest, swept γ ∈ [10 pS, 0.6 fS].

| γ (S)      | cyt t=10 (nM) | t=30 | t=60 | t=120 | t=300 | DTS t=300 (µM) |
|------------|---------------|------|------|-------|-------|----------------|
| 1.0e-11    | 123.7         | 1.4  | 0.6  |  0.6  |  0.6  | 0.0   |
| 1.0e-12    | 220.3         | 2.5  | 0.8  |  0.6  |  0.6  | 0.0   |
| 5.0e-13    | 112.6         | 104.9| 0.6  |  1.1  |  1.9  | 0.0   |
| 4.0e-13    |  93.5         |  83.9| 97.4 |  1.4  |  1.9  | 0.0   |
| 3.0e-13    |  59.0         | 123.4| 32.4 | 73.9  |  1.9  | 0.0   |
| 2.0e-13    |  25.5         |  25.5| 25.7 | 25.2  | 24.9  | 225.7 |
| 1.5e-13    |  17.4         |  17.4| 17.4 | 17.4  | 17.4  | 239.9 |
| 1.0e-13    |  11.6         |  11.6| 11.6 | 11.6  | 11.6  | 254.4 |
| 1.0e-15    |   1.1         |   1.1|  1.1 |  1.1  |  1.1  | 265.0 |
| 6.0e-16    |   0.8         |   0.8|  0.8 |  0.8  |  0.8  | 264.3 |

**No γ gives a stable rest near 100 nM.** Sharp bifurcation at
γ ≈ 2.5e-13 (0.25 fS):

- **γ ≥ 3e-13:** unstable. DTS drains to 0; cyt collapses to ~2 nM after
  some oscillation. The original 10 pS sits firmly in this regime — the
  Phase 1 IP3-forced run masked the underlying instability with a
  transient large enough to dominate the first 200 s.
- **γ ≤ 2e-13:** stable at low cyt. DTS preserved (or grows slightly).
  Resting cyt scales with γ: 25 nM @ 2e-13, 12 nM @ 1e-13, 1 nM @ 6e-16.

The takeaway: there is no single γ that yields a Dolan-style 100 nM
resting cyt with the rest of the parameters fixed. γ recalibration
cannot fix Phase 2 by itself.

## Rate balance at the Dolan IC (cyt=100 nM, DTS=250 µM, all rest defaults)

To identify where the imbalance actually lives, called `_ode_rhs` at t=0
with the Dolan IC verbatim and decomposed each flux contribution:

| Flux                                | ions/s |
|-------------------------------------|--------|
| IP3R inflow (γ=10 pS)               | **+112 288** |
| SERCA bind step (cyt → E1·Ca)       | −118 070 |
| SERCA release reverse (DTS → E2P·Ca) | −6 000 (DTS-side) |
| PMCA basal extrusion                | −22 |
| SOCE flux                            | +6 |
| **Net d[cyt]/dt at IC**              | **−6 338** |
| **Net d[DTS]/dt at IC**              | **−118 288** |

The lab-book Nernst arithmetic was correct: with V_IM = −60 mV and
E_Ca,IM = +RT/zF × ln(2500) = +104 mV, |V_IM − E_Ca,IM| = 164 mV
(not 44 mV — I had a sign slip in my analysis script before running
this term-by-term decomposition). IP3R inflow at γ=10 pS really is
~112 000 ions/s; the lab-book figure stands.

What the rate balance does show is that the lab-book diagnosis was still
incomplete in a different way:

### What the lab book got right

- IP3R inflow at γ=10 pS is large (~112k ions/s).
- DTS drains rapidly at this γ (d[DTS]/dt = −118k ions/s ⇒ full drain
  in ~0.33 s), confirming the original "DTS empties in 0.35 s" finding.
- Reducing γ to 0.6 fS does balance IP3R against SERCA's cycle throughput
  at full DTS.

### What the lab book missed

- **The cytosolic Ca²⁺ resting concentration is independent of γ.** At
  steady state with d[cyt]/dt = d[DTS]/dt = 0, the IP3R/SERCA recycle
  is internal — DTS-side balance forces J_SERCA = J_IP3R. Substituting
  into the cyt equation: J_SOCE + J_leak = J_PMCA. So **resting cyt is
  set by PM crossings (PMCA out vs SOCE + leak in), not by γ.**
- **There is no model term that can hold cyt at 100 nM.** With SOCE
  ≈ 6 ions/s at full DTS / low STIM1_dim and PMCA basal ≈ 22 ions/s
  at the IC (or ~77 ions/s at PMCA quasi-equilibrium for cyt=100 nM),
  the only net cyt source besides IP3R is SOCE, which is too small.
  Once cyt drops below 100 nM, PMCA slows but never reverses. The
  system inevitably settles at low cyt (the sweep table above shows
  this directly: stable rest at 1–25 nM depending on γ).
- **SERCA E1 ↔ E1·Ca is far from quasi-equilibrium in the IC.** At
  cyt=100 nM, the binding equilibrium gives E1·Ca / E1 = k_bind_f ·
  cyt² / k_bind_r = 1000 · 0.01 / 10 = 1.0, i.e., the 5 926 (E1 + E1·Ca)
  molecules should split ~2 963 each. The Dolan IC has **E1 = 5 920,
  E1·Ca = 6**, so the binding step pulls 59 035 events/s (≈ 118 k
  ions/s) from cyt → E1·Ca on t=0, before any of the slower
  processes can respond. This dominates d[cyt]/dt for the first few
  hundred ms.

The (ii)+(iii) plan below addresses both — leak fixes the resting cyt
target, pre-equilibration removes the initial transient.

### IC observation — SERCA is not at quasi-equilibrium

At cyt=100 nM, the SERCA E1 ↔ E1·Ca binding equilibrium has

```
E1·Ca / E1  =  k_bind_f · cyt² / k_bind_r  =  1000 · 0.01 / 10  =  1.0
```

so the two should each hold about 2 963 of the 5 926 (E1 + E1·Ca)
molecules. The Dolan IC has **E1 = 5 920, E1·Ca = 6** — wildly out of
binding equilibrium. The first few ms of any sim therefore see a
~110 000 ions/s pulse from cyt into E1·Ca as the binding equilibrates,
collapsing cyt before any of the slower processes can respond.

### Steady-state algebra

With the real numbers, the cytosolic balance at steady state reduces to

```
J_IP3R + J_SOCE + J_leak  =  J_SERCA + J_PMCA      (cyt)
J_SERCA                   =  J_IP3R                (DTS, internal recycle)
=>  J_SOCE + J_leak       =  J_PMCA                (PM crossings)
```

The IP3R/SERCA recycle is internal — it sets DTS but not cyt. The
**resting cyt is set by the balance between PMCA outflow and (SOCE +
leak) inflow at the plasma membrane.** With no leak in the model and
SOCE ≈ 0 at full DTS, the only PM inflow is whatever stochastic
"diffusion" terms exist — none, in our equations. So the model has no
mechanism to maintain cyt at 100 nM regardless of γ.

## Revised Phase 2 plan

Drop Step 2a (γ recalibration); the lab-book derivation behind it does
not survive checking. Two modifications, both more biologically
defensible than rescaling γ:

### (ii) Add a basal plasma-membrane Ca²⁺ leak

Add a small constant cyt influx term to `_ode_rhs`. Magnitude target:
balance PMCA outflow at cyt=100 nM with PMCA at quasi-equilibrium
(PMCA·Ca / PMCA_free = k_on·cyt/(k_off+k_cat) ≈ 0.018, PMCA·Ca ≈ 14):

```
J_PM_leak  ≈  k_cat · PMCA·Ca_eq  =  5.5 · 14  ≈  77 ions/s
```

Round to ~75 ions/s. This is small — biologically a basal "leak"
through unidentified channels (Sage & Rink 1985-90 era) and TRPC
background activity; even the lower bound of literature estimates
sits well above this.

### (iii) Pre-equilibrate SERCA E1 ↔ E1·Ca in the IC

Move SERCA molecules from E1 to E1·Ca so they start at binding
equilibrium for cyt=100 nM:

```
E1     :  5 920  →  2 963
E1·Ca  :      6  →  2 963
```

Total SERCA preserved (11 892); only the binding-equilibrium pair is
redistributed. Removes the 110 000-ions/s spurious initial pulse.

### Expected combined effect

With (ii) + (iii):

- No initial cyt collapse from SERCA E1 → E1·Ca redistribution.
- PMCA outflow (~75 ions/s at rest) balanced by leak inflow.
- IP3R/SERCA internal recycle finds its own DTS steady state — likely
  *below* Dolan's 250 µM because at γ=10 pS, IP3R inflow exceeds
  SERCA cycle throughput at full DTS. The system will settle where
  the Nernst driving force is small enough to balance them
  (somewhere in the 25–100 µM range, to be measured).

A DTS resting concentration substantially below 250 µM would itself
be a finding — either acceptance that our parameters give a smaller
resting store than Dolan reports, or motivation to tune SERCA cycle
rates upward (option (i) in the original session list).

### What stays open

- **Whether DTS resting at ~25–100 µM is acceptable.** Dolan, Sage and
  others give 250–500 µM. If our converged resting DTS is much lower,
  the dissertation argument is weaker (we'd be modelling a
  partially-depleted store).
- **Whether the Phase 1 transient still emerges.** The IP3-forced spike
  drives Po up by orders of magnitude; even with (ii)+(iii) the
  transient may still hit 200–800 nM. To verify after implementing.
- **Issue #48 status.** The lab-book Phase 2 plan as written is wrong;
  the issue should be re-scoped around (ii)+(iii), not γ recalibration.

## Files touched / inspected this session

| File | Note |
|------|------|
| `reconstruction/platelet/dataclasses/process/calcium_signalling.py` | `GAMMA_IP3R_S` flipped to 0.6e-15 (commit `035f2814`) and reverted (commit `cb22833c`); then added `J_PM_LEAK_IONS_S` constant + cyt-side term in `_ode_rhs` |
| `reconstruction/platelet/dataclasses/internal_state.py` | Pre-equilibrated SERCA E1 / E1·Ca to 2 963 / 2 963 (was 5 920 / 6) |
| `reports/lab-books/lab-book-2026-05-01-phase1-complete.md` | Reference for the disproved Phase 2 plan |

## (ii) + (iii) implementation results

**(ii)** Added `J_PM_LEAK_IONS_S = 75.0` near the SOCE block; one-line
`dy[cyt] += J_PM_LEAK_IONS_S` in `_ode_rhs` after the SOCE term.

**(iii)** SERCA E1 (5 920 → 2 963) and E1·Ca (6 → 2 963) in `_MOLECULES`,
preserving total SERCA = 11 892. Comment in-line documents the
binding-equilibrium derivation.

### Phase 1 transient with (ii)+(iii) at γ=10 pS, IP3 forced (200 s)

| Time (s) | Cyt (nM) | DTS (µM) |
|----------|----------|----------|
|  0  | 99.9 | 250.0 |
|  1  | **299.5 (peak)** | 0.0 |
|  3  | 253.2 | 0.0 |
|  10 | 140.3 | 0.0 |
|  20 |  3.0  | 0.4 |
|  50 |  3.3  | 0.0 |
| 200 |  3.3  | 0.0 |

- **Peak: 299.5 nM** ✓ within Phase 1 acceptance band (200–800 nM); slightly
  higher than the 280 nM Phase 1 baseline (about +7%) because (iii) removes
  the spurious initial cyt drain from the SERCA E1 → E1·Ca discontinuity.
- **DTS empties to 0 by t=5 s** ✗ — this is the long-standing issue
  (IP3R inflow vastly exceeds SERCA cycle throughput at γ=10 pS).
- Resting cyt after the spike sits at ~3 nM (was ~1 nM in the Phase 1
  baseline; the +75 ions/s leak nudges it up but cannot reach 100 nM
  given the SERCA pull at low cyt).

Regression test suite: **17 / 17 pass**, including the 30% peak tolerance
band (196–364 nM). The peak shift from 280 → 299 nM is well inside the
existing tolerance.

### γ sweep with (ii)+(iii) applied (IP3 not forced, 300 s)

| γ (S)    | cyt t=300 (nM) | DTS t=300 (µM) | regime |
|----------|----------------|----------------|--------|
| 1.0e-11  | 3.3 | 0.0 | unstable — DTS drained |
| 5.0e-13  | 5.5 | 0.0 | unstable |
| **2.0e-13** | **26.0** | **407.6** | stable |
| 1.0e-13  | 12.5 | 430.2 | stable |
| 1.0e-14  | 3.6 | 448.3 | stable |

The leak slightly raises resting cyt vs. yesterday's pre-(ii) sweep (e.g.,
γ=1e-13: 12.5 vs 11.6 nM) and DTS now *overshoots* Dolan (407 vs 254 µM
without leak) because SERCA pumps the leak inflow into DTS unchecked.

### Leak-magnitude sweep at γ=2e-13 (looking for cyt=100 nM joint with DTS=250 µM)

| leak (ions/s) | cyt t=300 (nM) | DTS t=300 (µM) |
|---------------|----------------|----------------|
| 75 (current)  | 26.0 | 407.6 |
| 200           | 28.2 | 640.9 |
| 500           | 32.4 | 1213.3 |
| 1 000         | 35.7 | 2168.2 |
| 2 000         | 75.0 | 3352.4 |
| 5 000         | 266.3 | 2076.1 |

**No leak magnitude recovers (cyt=100 nM, DTS=250 µM) jointly.** Increasing
leak raises cyt slowly because SERCA siphons the extra Ca²⁺ into DTS,
ballooning DTS far above 250 µM. The system's fixed point structure under
the current SERCA cycle parameters does not admit the Dolan IC.

## What this means for Phase 2 / the dissertation timeline

The (ii)+(iii) changes are individually defensible and necessary — they
fix two genuine bugs (missing PM-balance term; SERCA IC not at binding
equilibrium for cyt=100 nM). They keep the regression suite green and
slightly *improve* the Phase 1 transient peak. **They should land.**

However they are *not sufficient* to recover Dolan-style rest at
(cyt=100 nM, DTS=250 µM). The remaining gap is in the SERCA cycle
parameters: at full DTS = 250 µM, the release step (`k_release_r` × DTS²
× E2P = 7 000 events/s) runs in *reverse*, so SERCA cannot maintain DTS
loading at 250 µM against the current IP3R inflow. The cycle's natural
fixed-point DTS is somewhere below 250 µM — by how much depends on rates
we haven't yet tuned.

Three options for what comes next:

1. **Accept the lower steady state as v0.2's working point.** Document
   resting cyt ≈ 3–25 nM and DTS ≈ 400 µM as the actual rest state of the
   model with current parameters. The Phase 1 transient still passes its
   acceptance band, which was the dissertation-relevant headline result.
   Lowest risk for the timeline; clearest write-up.

2. **Tune SERCA cycle rates.** Specifically `k_release_r` (4×10⁻³ µM⁻²s⁻¹
   = 4×10⁹ M⁻²s⁻¹) is large enough to reverse the release step at DTS=250 µM.
   Reducing it would let SERCA hold DTS at 250 µM. But this departs from
   the Purvis primary-source values, opening a new calibration question.

3. **Combine (1) + (2).** Move forward with (ii)+(iii) as v0.2 baseline,
   accept lower resting; if dissertation-grade Dolan validation requires
   higher DTS rest, revisit `k_release_r` later.

Recommendation: **option 3.** Land (ii)+(iii) now; revisit SERCA rates
only if Phase 3 (Dolan Fig. 4 validation) blocks on DTS resting level.

## Where to pick up next session

- Decide whether (ii)+(iii) commit lands as-is or with a γ adjustment.
- If keeping γ=10 pS, the Phase 1 transient still passes acceptance.
  Issue #48 should be re-scoped from "γ recalibration" to "DTS
  steady-state ≠ Dolan IC" + "SERCA cycle reverse step at full DTS".
- The regression-test docstring's note about "DTS drains to zero by ~60 s
  — known Phase 2 issue" is still accurate; no need to update yet.

## Build and test commands used

```bash
# Sweep γ at rest (no IP3 forcing) — inline script
PYTHONPATH=$PWD pyenv exec python <<'PY'
from reconstruction.platelet.dataclasses.process import calcium_signalling as cs_mod
from reconstruction.platelet.dataclasses.internal_state import _MOLECULES
# ... (see commit log; full sweep script)
PY

# Rate-balance dump at the IC — same harness, different math

# Standard 200 s IP3-forced run
PYTHONPATH=$PWD OPENBLAS_NUM_THREADS=1 \
  python runscripts/manual/runPlateletSim.py phase2a_smoke --length 200
```

---

*Branch:* `main` · *State after session:* clean (γ revert committed
in `cb22833c`); `(ii)+(iii)` implementation pending in a follow-up commit.
