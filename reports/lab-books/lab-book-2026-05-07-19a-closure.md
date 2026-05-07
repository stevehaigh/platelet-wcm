---
title: "Lab book — 2026-05-07: Issue 19a closure (Option C) + STIM1 dimer-count fix"
---

# Lab book — 2026-05-07: Issue 19a closure (Option C) + STIM1 dimer-count fix

## Session goal

Execute issue 19a as planned: build a `restConvergence.py` runscript,
integrate the ODE for ~600 s with no IP3 forcing, freeze the
converged 27-species state vector as the new initial conditions.
Acceptance criteria from the issue body: max |dc/dt| < 1 count/s
across all species, startup spike eliminated, Phase 3 still ≥ 4/5,
no biology change.

## Outcome: 19a is unachievable as scoped — pivot to 19b

### What the convergence run actually found

`restConvergence.py` integrates `_ode_rhs` directly with
`ip3_forced=False`, no STIM/SOCE perturbation. After 6000 s the
system reaches a true numerical fixed point (max |dy/dt| = 0.000
count/s, well below the 1 count/s acceptance threshold), but the
converged macro-concentrations are far from biological:

| Quantity | Initial (Dolan IC) | Converged | Dolan target |
|---|---|---|---|
| Ca²⁺_cyt | 100 nM | **1786 nM** | 100 nM |
| Ca²⁺_DTS | 250 µM | **0.056 µM** | 200–300 µM |
| IP3 | 50 nM | 50 nM (held; no forcing) | 50 nM |

So the model **does** have a stable fixed point — but it sits 18×
above the Dolan resting cytosolic band, and with an essentially
empty DTS. This makes the 19a acceptance criteria internally
inconsistent: "max |dc/dt| < 1 count/s" requires landing at this
fixed point, but "no biology change" requires keeping concentrations
in the Dolan band, which is not a fixed point of the model.

### What this means for the diagnosis

The Phase 0 audit (lab-book-2026-05-07-phase-0-biology-audit.md)
showed that the Dolan Table S1 IP3R sub-state populations are not
a Markov-chain equilibrium of our Sneyd-Dufour ladder. We hypothesised
that fixing this — letting the sub-states settle to their true
equilibrium — would eliminate the startup spike. It does (no
sub-state racing) but it doesn't fix the underlying problem: even
after sub-state equilibration, the model's Ca²⁺ flux balance is
broken. Specifically:

- Cytosolic Ca²⁺ wants to sit at ~1.8 µM, not 100 nM.
- DTS wants to be empty, not at 250 µM.
- The PM-leak / PMCA balance (`J_PM_LEAK_IONS_S = 75`) was
  calibrated against the Dolan IC at cyt = 100 nM — but that's not
  the fixed point of the system, so the calibration is anchored at
  the wrong operating point.

This is a structural / parameter problem, not a Markov-chain
equilibrium problem. Issue 19a was scoped on the assumption that the
sub-state IC was the only thing keeping the Dolan IC from being a
fixed point. The Phase 0 audit + this convergence run together show
that the sub-state IC explains the **startup spike** but not the
**resting-state values**.

### Decision: Option C — close 19a, route the work to 19b / #24

Three options were on the table:

- **A** — Constrain cyt = 100 nM, DTS = 250 µM, equilibrate only
  sub-states. Eliminates startup sub-state racing but the model
  still drifts toward the 1.8 µM fixed point during the run, so
  doesn't satisfy "max |dc/dt| < 1 count/s".
- **B** — Use the converged state (cyt = 1786 nM, DTS empty) as the
  new IC. Satisfies the convergence criterion mathematically but
  abandons the biological constraint that cytosolic Ca²⁺ is ~100 nM
  at rest. Phase 3 transient peak meaning is no longer comparable.
- **C** — Close 19a as "diagnosis complete; the fix needs structural
  / parameter work, not just sub-state IC re-derivation." Pivot to
  19b (#24) — likely needs MCU (#22) plus J_PM_LEAK / SOCE retuning
  to land the fixed point at cyt ≈ 100 nM, DTS ≈ 200–300 µM.

Going with **C**. The Phase 0 audit + the 6000 s convergence
together produced the diagnostic value 19a was supposed to deliver
(rule out simple Markov-chain misalignment), and the actual fix
belongs to 19b / #22.

## Side finding (and fix): STIM1 dimer count bug

The convergence script's mass-conservation check flagged STIM1
losing 18.5% of total monomers across the run (4265 → 3475). Trace:

- ODE: `dy[STIM1_dim] += v_dim` and `dy[STIM1_free] += -2·v_dim` —
  i.e. 1 dimerisation event creates 1 unit of `STIM1_dim` and
  consumes 2 free monomers. So `STIM1_dim` is in **dimer particles**.
- IC: `STIM1_dim = 22` (with the comment "monomer-equiv count;
  11 dimers"). So the IC is in **monomer-equivalents**.
- MWC and listener divide by `STIM_MONOMERS_PER_DIMER = 2` to get
  dimer count, consistent with the IC convention.
- These two conventions cancel at *t = 0* (IC stored as 22 / 2 = 11
  matches Dolan's 11 dimers), but as the cycle runs the ODE only
  adds 1 unit per reaction — so `STIM1_dim`-as-monomer-equiv is
  growing at half the rate it should, leaking 1 monomer per event.

**Fix**: switch to dimer-particle count throughout. IC: 22 → 11
(matches Dolan Table S1's "STIM1₂ (11)" exactly), recalibrate
`k_dim_f` from 1.15×10⁻⁴ to 5.73×10⁻⁵ (= 1.0 × 11 / 438²) for
detailed balance at the new IC, remove `/STIM_MONOMERS_PER_DIMER`
from the MWC and listener (they now consume dimer count directly).
`STIM_MONOMERS_PER_DIMER` is kept as a constant for total-monomer
mass-balance accounting only.

Mass conservation now passes (STIM1: 4265 ✓). All 21 platelet tests
still pass; mypy clean.

## Phase 3 impact of the dimer fix

Re-ran Phase 3 with the STIM1 dimer fix in place:

| Metric | Pre-Caride-k12 | Post-Caride-k12 | Post-dimer fix |
|---|---|---|---|
| +Ca_ex peak | 299 nM | 380 nM | **392 nM** |
| −Ca_ex peak | 298 nM | 325 nM | 325 nM |
| SOCE differential | 1 nM | 54 nM | **67 nM** |
| Criteria pass | 3/5 | 4/5 | 4/5 |

Slight further improvement of both +Ca_ex peak and SOCE
differential. Still 4/5 — the SOCE differential criterion is the
remaining gap and is the work tracked in 19b/#22. Test bands
(±30% of 380 nM baseline) still encompass 392 nM, so no
test_phase3.py update needed.

## What's committed in this session

| File | Change |
|---|---|
| `runscripts/manual/restConvergence.py` | New diagnostic runscript |
| `reports/data/rest-converged-2026-05-07.json` | Convergence output |
| `reports/lab-books/lab-book-2026-05-07-19a-closure.md` | This entry |
| `reconstruction/platelet/dataclasses/internal_state.py` | STIM1_dim IC 22 → 11; mass updated to per-dimer |
| `reconstruction/platelet/dataclasses/process/calcium_signalling.py` | k_dim_f 1.15e-4 → 5.73e-5; remove /2 in MWC; comment updates |
| `models/platelet/listeners/calcium_trace.py` | Remove /2 in stim2_p; remove unused import |
| `runscripts/manual/restConvergence.py` | Mass-balance formula uses 2·STIM1_dim |

## Next session

19b / #24 / #22 — the resting Ca²⁺ flux balance. Two strands:

1. **Add MCU** (#22) — the missing biology. Adds a slow Ca²⁺
   buffer compartment that:
   - During the spike, diverts ~10–30% of cyt Ca²⁺ inflow into the
     mitochondrial matrix (lower peak, broader transient).
   - After the spike, slowly releases Ca²⁺ over minutes via mNCX
     (provides the missing post-transient "tail" that pushes cyt
     back to ~100 nM rather than collapsing).
   - Three primary papers in `source-info/calcium-papers/`:
     Ajanel 2025, Ghatge 2026, Shehwar 2025.
2. **Resting-state flux retune** — once MCU is in, re-run
   `restConvergence.py` and see where the new fixed point sits.
   Likely still needs a small adjustment to `J_PM_LEAK_IONS_S`
   and/or basal SOCE conductance to land cyt at ~100 nM. Volume
   uncertainty (memory: `feedback_volume_uncertainty.md`) means we
   should also sensitivity-check the conversion factors.

---

*Branch:* `main` · *Status:* 19a closed; STIM1 dimer fix landed ·
*Linked issues:* closes #19; further work in #24, #22 · *Triggered
by:* user decision Option C 2026-05-07
