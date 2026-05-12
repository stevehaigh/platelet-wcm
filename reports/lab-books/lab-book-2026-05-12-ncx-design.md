---
title: "Lab book — 2026-05-12: NCX Na⁺/Ca²⁺ exchanger (v0.3.4)"
---

# Lab book — 2026-05-12: NCX Na⁺/Ca²⁺ exchanger

## Why this change

After Phase 4 + v0.3.2 + #22 (MCU), the long-recovery sim shows DTS
overshoot to ~1 060 µM at t = 3 000 s. The MCU work confirmed that
adding *intracellular* Ca²⁺ buffering doesn't extrude Ca²⁺ from the
cell — it just redistributes it (and actually makes the DTS overshoot
*worse*, because mito Ca²⁺ slowly flows back into the DTS via SERCA).

The fundamental bottleneck is **PMCA at low cyt Ca²⁺**: at cyt = 100–200
nM, PMCA Vmax is rate-limited and can't keep up with the Ca²⁺ load
imported during the transient.

NCX (Na⁺/Ca²⁺ exchanger) is the canonical second Ca²⁺ extrusion
pathway. In cardiomyocytes and neurons it does the bulk of high-Ca²⁺
extrusion. In platelets, NCX presence is documented but functional
contribution is contested:
- Burkhart 2012 platelet proteome: NCX1 (SLC8A1) and NCX3 (SLC8A3)
  present
- Sage & Rink 1985: limited NCX activity in functional assays
- Some functional studies: NCX contributes 10–30% of platelet Ca²⁺
  extrusion at peak

Dissertation framing: **NCX is biologically defensible but uncertain**;
we model it at a moderate level to test whether it can resolve the
DTS-overshoot recovery time.

## Design

### No new species

NCX is a single flux term: cyt Ca²⁺ → extracellular Ca²⁺ (forward mode
only; reverse mode would require modelling cyt Na⁺ and Vm, both out of
scope). Gated on `CA_EX_UM > 0` (same condition as SOCE / P2X1 /
PMCA): in the EDTA condition, NCX has nothing to exchange against.

### Kinetic scheme

NCX has both a *kinetic* substrate-binding step (the Hill K_m for cyt
Ca²⁺ at the transport site) and an *allosteric* regulatory Ca²⁺ site
that activates the exchanger. Real NCX is essentially silent at low
cyt Ca²⁺ (the allosteric gate keeps it off), regardless of what the
substrate kinetics would predict alone. We model this as the product
of an allosteric gate and a kinetic factor:

```
V_NCX = V_max_NCX × g_act(ca_cyt) × ca_cyt / (K_m + ca_cyt)

g_act(ca_cyt) = (ca_cyt / K_a)^h / (1 + (ca_cyt / K_a)^h)
```

This captures the qualitative biology: NCX is essentially silent at
rest (gate near zero), turns on during a transient (gate saturates),
and remains substrate-limited at the transport site.

### Parameters

| Param | Value | Source / rationale |
|---|---|---|
| `V_max_NCX` | 2 000 ions/s | Platelet-scaled estimate; ~10% of PMCA structural Vmax. Calibration anchor against DTS overshoot. |
| `K_m_NCX` | 5 µM | Mid-range for NCX cyt-Ca²⁺ affinity (literature 1–10 µM) |
| `K_a_NCX` | 0.3 µM | Allosteric activation half-point — gate ~ 0.01 at rest, ~ 0.9 at peak |
| `h_NCX` | 4 | Allosteric Hill cooperativity; switch-like activation |

Predicted V_NCX at key states:

| cyt Ca²⁺ | g_act | kinetic factor | V_NCX (ions/s) |
|---|---|---|---|
| 0.1 µM (rest) | 0.012 | 0.020 | ~ 0.5 |
| 0.2 µM | 0.16 | 0.038 | ~ 12 |
| 0.3 µM | 0.50 | 0.057 | ~ 57 |
| 0.5 µM (peak) | 0.89 | 0.091 | ~ 162 |
| 1.0 µM | 0.99 | 0.167 | ~ 331 |

So at rest, NCX is negligible (~ 0.5 ions/s — preserves resting
balance). During the transient peak, NCX adds ~ 160 ions/s extrusion.
During recovery at cyt = 150–200 nM, NCX adds ~ 10–30 ions/s.

### Expected impact

- **Resting state**: invariant. NCX flux ~ 0.5 ions/s at rest, much
  smaller than the PMCA + PM-leak balance (~ 30 + 75 = 105 ions/s).
- **Peak cyt Ca²⁺**: slight attenuation (a few nM) from extra
  extrusion during the transient. May need to bump γ_P2X1 slightly.
- **DTS overshoot at t = 3 000 s**: should improve. Total cell Ca²⁺
  loss rate during recovery rises from ~125 ions/s (PMCA – PM_leak)
  to ~135–155 ions/s (adding NCX 10–30). Modest 10–25% improvement.
- **Phase 3**: 5/5 should still hold after small γ_P2X1 retune.

### Acceptance criteria

1. Resting cyt 95–125 nM, DTS 220–290 µM, IP3 40–60 nM.
2. Phase 3 maintains 5/5 criteria.
3. **DTS at t = 3 000 s < 850 µM** (was 1 062 with MCU only) — soft
   target. Stretch: < 500 µM (we'll see).
4. **cyt at t = 3 000 s < 130 nM** (was 134) — modest improvement.
5. All 21 regression tests pass.

### Risks

- NCX in platelets is contested. We're modelling a *plausible*
  contribution. If it turns out NCX is functionally inactive in
  platelets, the dissertation framing should be "we tested adding NCX
  and found it contributes X to recovery — consistent with reports of
  variable NCX contribution."
- Calibration: same coupling pattern as before. Likely needs 1–2
  iterations.

### Files to change

| File | Change |
|---|---|
| `calcium_signalling.py` | Add `K_NCX` dict + ODE term; gate on `CA_EX_UM > 0` |
| `internal_state.py` | No change (no new species) |
| `models/platelet/tests/sim/test_regression.py` | Update peak Ca²⁺ baseline if it shifts |
| Lab book (this doc) | Results section after implementation |
| `dissertation-notes.md` | Update §7.3 NCX entry to "modelled with caveat" |

---

---

## Implementation results

### Final calibration

| Param | Design | Final | Why adjusted |
|---|---|---|---|
| `V_max` | 2 000 ions/s | **5 000 ions/s** | First iteration at 2k gave only 5 % DTS recovery improvement. Bumped to 5k for more recovery-phase contribution. |
| `K_m` | 5.0 µM | 5.0 (unchanged) | |
| `K_a` | 0.3 µM | **0.2 µM** | Lowered slightly so NCX activates earlier in the recovery phase (cyt ~ 200 nM) rather than only at the transient peak. |
| `h` | 4 | 4 (unchanged) | |

Predicted V_NCX at the final calibration:

| cyt | g_act | kinetic | V_NCX |
|---|---|---|---|
| 0.1 (rest) | 0.059 | 0.020 | ~ 6 |
| 0.13 (recovery floor) | 0.158 | 0.025 | ~ 20 |
| 0.2 | 0.500 | 0.038 | ~ 95 |
| 0.5 (peak) | 0.975 | 0.091 | ~ 444 |
| 1.0 | 0.998 | 0.167 | ~ 833 |

### Resting state — ✓

| Quantity | Target | Result |
|---|---|---|
| cyt | 95–125 nM | 107 nM ✓ |
| DTS | 220–290 µM | 228 µM ✓ (slightly lower but in band) |
| IP3 | 40–60 nM | 50 nM ✓ |

### Phase 3 — ✓ 5/5

- +Ca_ex peak: 419 nM (was 434)
- −Ca_ex peak: 314 nM (was 314)
- SOCE differential: 105 nM (was 120; criterion ≥ 100)
- All in Dolan bands

### Long recovery — modest but real improvement

| Version | cyt @ t=3000 | DTS @ t=3000 | Notes |
|---|---|---|---|
| v0.3.2 (pre-MCU) | 126 nM | 758 µM | PMCA bottlenecked |
| v0.3.3 (MCU only) | 134 nM | 1 062 µM | MCU redistribution made overshoot *worse* |
| **v0.3.4 (MCU + NCX)** | **128 nM** | **793 µM** | NCX gives ~25 % overshoot reduction vs MCU-only |

NCX's contribution at the slow-recovery phase (cyt ~ 130 nM): ~20 ions/s
additional extrusion on top of PMCA's ~150 ions/s. About 13 % extra
clearance rate. The 25 % DTS overshoot reduction reflects the bigger
contribution during the peak (NCX at 0.5 µM cyt ≈ 444 ions/s) plus the
late-recovery contribution accumulated over 2 000 s.

### Acceptance criteria

| # | Criterion | Result |
|---|---|---|
| 1 | Resting cyt / DTS / IP3 in band | ✓ |
| 2 | Phase 3 5/5 | ✓ |
| 3 | DTS at t=3000 < 850 µM | ✓ (793) |
| 4 | cyt at t=3000 < 130 nM | ◐ (128 — just at boundary) |
| 5 | 21 tests pass | ✓ |

**4.5/5** — passing the headline goals.

### What this confirms about the model

The MCU + NCX results together demonstrate a clean mechanistic story
for the dissertation:

1. **PMCA is the bottleneck** for total Ca²⁺ extrusion at low cyt.
2. **Intracellular buffers (CaM, gelsolin, CALR, mito) only redistribute** Ca²⁺ — they don't extrude it. MCU specifically *worsens* the DTS overshoot because mito-bound Ca²⁺ feeds back into the DTS via the SERCA cycle.
3. **Only a second extrusion pathway** (NCX, or higher PMCA capacity) can accelerate recovery.
4. **The remaining DTS overshoot** is real biology: strong agonist stimulation imports ~ 500 k Ca²⁺ ions via SOCE + P2X1, and total cell-Ca²⁺ extrusion is fundamentally rate-limited even with two extrusion pathways. Real platelets after thrombin show similar prolonged Ca²⁺ tails over 10–60 min.

### Open / forward-looking

- The NCX functional contribution in platelets is contested (Sage &
  Rink 1985 vs subsequent reports). Our calibration sits in the
  "moderate contribution" zone. Sensitivity analysis would test how
  robust the recovery dynamics are to NCX magnitude.
- If we wanted to push the DTS overshoot recovery faster, options
  include: higher PMCA copy number (Burkhart 2012 range 0.8–4 k vs
  our 769), bidirectional NCX (reverse mode during Vm changes), or
  reducing SOCE during stim (less Ca²⁺ imported).

---

*Branch:* `main` · *Status:* v0.3.4 complete · *Linked: NCX entry
in dissertation-notes §7.3*
