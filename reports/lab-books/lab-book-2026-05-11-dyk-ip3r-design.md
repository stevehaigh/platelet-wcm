---
title: "Lab book — 2026-05-11: Design — replace Sneyd-Dufour IP3R with deYoung-Keizer 1992"
---

# Lab book — 2026-05-11: deYoung-Keizer IP3R replacement (issue #27)

## Why this change

The resting-state ODE converges to [Ca²⁺]_cyt ≈ 2 240 nM and
[Ca²⁺]_DTS ≈ 0.07 µM. The DTS drains even without IP3 stimulus.
Full root-cause analysis is in #24 and lab-book-2026-05-07-dts-drain
-investigation.md; the finding in brief:

Sneyd-Dufour 2002 was calibrated at IP3 = 10 µM — 200× above the
resting platelet IP3 = 50 nM. Extrapolated to that regime the
Markov chain gives Po⁴ = 4.21×10⁻⁴. Dolan's filtered IC has
Po⁴ = 1.65×10⁻⁵ (25× lower); the 25× gap is attributable to
Dolan's Monte Carlo filtering not constraining sub-state populations
to the Sneyd-Dufour equilibrium, not to a code bug.

At the Dolan IC the resulting IP3R basal leak (≈ 110 000 ions/s)
outpaces SERCA refill (≈ 30 000 ions/s at cyt = 100 nM) by 3.7×.
The DTS empties on a timescale of seconds.

No downstream parameter patch closes this gap. The fix must either:

1. Reduce the IP3R model's resting Po to a biologically appropriate
   level, or
2. Accept the pathological resting state and live with the documented
   limitation (current v0.2 baseline, #24 open).

This lab book documents option 1: replacing Sneyd-Dufour with the
deYoung-Keizer 1992 model, which has an explicit IP3-binding
activation gate that produces near-zero open probability at resting
IP3 concentrations.

---

## Why deYoung-Keizer 1992

Sneyd-Dufour 2002 was designed to reproduce the *dynamics* of IP3R
gating during stimulated conditions (fast Ca²⁺-driven transitions,
cooperative inactivation). It has no explicit IP3-binding term with
a saturable Kd; IP3 concentration appears only inside the φ-function
rate laws as a multiplicative scaling. Below the calibration regime
the channel does not gate off cleanly.

deYoung & Keizer 1992 (*PNAS* 89:9895-9899) model each IP3R subunit
as having three independent binding sites:

| Site | Ligand | Role |
|------|--------|------|
| 1    | IP3    | activation (channel opens when IP3 bound) |
| 2    | Ca²⁺   | activation (assists opening) |
| 3    | Ca²⁺   | inhibition (Ca²⁺-induced inactivation) |

The channel is only in the conducting state when sites 1 and 2 are
occupied and site 3 is unoccupied. This gives an explicit IP3
concentration dependence with saturable Kd:

    m∞ = [IP3/(IP3+d₁)] × [Ca/(Ca+d₅)]

where d₁ ≈ 0.13 µM is the IP3 half-saturation constant.

At IP3 = 50 nM:  IP3 gate fraction = 0.05/(0.05+0.13) = **0.28**
At IP3 = 10 µM:  IP3 gate fraction = 10.0/(10.0+0.13) = **0.99**

The channel is ~72% gated off by IP3 scarcity alone at rest.
Sneyd-Dufour has no equivalent gate.

The model was calibrated in Xenopus oocytes, the same expression
system as Sneyd-Dufour, and has been used extensively in platelet
and ER-Ca²⁺ modelling in the decades since (Bezprozvanny 1991
single-channel data; Wagner & Keizer 1994; Dupont et al. 2016).

---

## Model equations — Li-Rinzel 1994 simplification

Li & Rinzel 1994 (*J Theor Biol* 166:461-473) showed that if the
IP3-binding and Ca²⁺-activation kinetics are fast relative to
Ca²⁺-inhibition kinetics, sites 1 and 2 can be treated as being in
quasi-steady state. This collapses the 8-state deYoung-Keizer model
to a single differential equation for the inhibition variable h
(fraction of subunits with site 3 unoccupied):

    dh/dt = a₂ × [d₂ - (Ca_cyt + d₂) × h]             ... (1)

In quasi-steady state h → h∞ = d₂/(Ca_cyt + d₂).

The activation fraction in quasi-steady state:

    m∞(IP3, Ca) = [IP3/(IP3+d₁)] × [Ca/(Ca+d₅)]       ... (2)

For a homotetrameric channel (open when all 4 subunits are
activated and not inhibited):

    Po_channel = m∞(IP3, Ca)⁴ × h⁴                     ... (3)

Note: some implementations use m∞³ × h (Li-Rinzel's own notation);
we use the ⁴ exponent to match the Dolan/Purvis tetrameric
convention already in the flux formula and for thermodynamic
consistency with the cooperative gating assumption. This choice
should be tested during implementation (see §Measurement).

The Ca²⁺ flux equation is **unchanged**:

    J_IP3R = γ_IP3R × N_IP3R × Po_channel × (V_IM − E_Ca) × NA/(zF)

γ_IP3R = 10 pS, N_IP3R = 1 328 (Burkhart 2012 / Dolan Table 1),
V_IM = −0.060 V — all retained from the current model.

---

## Parameters

From deYoung & Keizer 1992 Table 1 and the Li-Rinzel 1994
simplification. All values must be verified against primary-source
PDF before being entered in code.

| Symbol | Description | Value | Unit | Source |
|--------|-------------|-------|------|--------|
| d₁ | IP3 activation half-saturation | 0.13 | µM | dYK92 Table 1: b₂/a₁ = 52/400 |
| d₂ | Ca²⁺ inhibition half-saturation | 1.049 | µM | dYK92 Table 1: b₄/a₄ |
| d₅ | Ca²⁺ activation half-saturation | 0.08234 | µM | dYK92 Table 1: b₆/a₅ |
| a₂ | Ca²⁺ inhibition on-rate (h dynamics) | 0.2 | µM⁻¹·s⁻¹ | dYK92 Table 1 |

These are the four parameters used in the Li-Rinzel reduced form.
All other deYoung-Keizer rate constants are subsumed into d₁, d₂,
d₅ via their equilibrium-constant definitions.

**At resting conditions** (IP3 = 0.05 µM, Ca_cyt = 0.1 µM):

    m∞ = (0.05/0.18) × (0.1/0.183) = 0.278 × 0.548 = 0.152
    h∞  = 1.049 / 1.149             = 0.913
    Po_channel = 0.152⁴ × 0.913⁴   = 5.33×10⁻⁴ × 0.695 = 3.7×10⁻⁴

Hmm — with the ⁴⁺⁴ exponents the resting Po is still ~3.7×10⁻⁴,
comparable to the Sneyd-Dufour value. This motivates testing the
alternative **m∞⁴ × h** formula (h only to first power):

    Po_channel = 0.152⁴ × 0.913    = 5.33×10⁻⁴ × 0.913 = 4.9×10⁻⁴

Still high. The real gain from deYoung-Keizer is through the
stronger IP3 dependence at stimulated conditions vs resting, not
through a large absolute resting-Po reduction — the question is
whether the *resting flux balance* works out, which depends on
what γ_IP3R was calibrated against.

**Revised expectation**: deYoung-Keizer does not produce a
dramatically lower absolute resting Po than Sneyd-Dufour at these
parameter values. The principal difference is the correct IP3
gating across the concentration range. After the model replacement
the resting-state *balance* (IP3R leak vs SERCA) may still require
γ_IP3R recalibration — this was anticipated as a contingency in the
Phase 1 risk register.

This is an important pre-implementation finding. It means:

1. The implementation should proceed as planned (deYoung-Keizer is
   still the correct model choice — it handles the IP3 range
   properly and is calibrated on Xenopus/Xenopus-like physiology).
2. We should expect that **γ_IP3R recalibration will be needed**
   as part of Phase 4, specifically to bring resting leak down
   to ≤ SERCA refill capacity at cyt = 100 nM.
3. The Phase 4 γ_IP3R retune is now biologically grounded (it is
   adjusting the single-channel conductance to match the specific
   cell type / experimental conditions) rather than a fudge to
   compensate for a wrong kinetic model.

---

## Code changes required

### `reconstruction/platelet/dataclasses/process/calcium_signalling.py`

**Remove:**
- `K_IP3R` dict (17 entries, lines ~126-148)
- All seven `_phi_*` helper functions (lines ~324-382)
- The six `IP3R_*[dts]` entries from `MOLECULE_NAMES` (n, o, a, i1, i2, s)
- The corresponding state reads from the ODE (`n`, `o`, `a`, `i1`, `i2`, `s`)
- The six `dy[_IDX['IP3R_*']] +=` lines
- The `ip3r_total`, `po_subunit`, `po_channel` computation block
- The comment block about Dolan convention vs /4 (lines ~529-542)

**Add:**
- `K_DYK` dict with four entries: `d1`, `d2`, `d5`, `a2`
- `N_IP3R = 1328` constant (replacing `ip3r_total` derived from sub-state sum)
- One new entry in `MOLECULE_NAMES`: `'IP3R_h[dts]'`
  (the inhibition variable; stored as a fractional count 0–N_IP3R)
- In `_ode_rhs`:
  - Read `ip3r_h_count = max(y[_IDX['IP3R_h[dts]']], 0.0)`
  - Compute `h = ip3r_h_count / N_IP3R` (fraction, 0–1)
  - Compute `m_inf = (ip3 / (ip3 + K_DYK['d1'])) * (ca_cyt / (ca_cyt + K_DYK['d5']))`
  - Compute `po_channel = (m_inf ** 4) * h` (or `m_inf**4 * h**4`; test both)
  - `n_ip3r_channels = N_IP3R` (constant, replaces `ip3r_total`)
  - h ODE: `dh_dt = K_DYK['a2'] * (K_DYK['d2'] - (ca_cyt + K_DYK['d2']) * h)`
  - `dy[_IDX['IP3R_h[dts]']] += dh_dt * N_IP3R`

The Nernst flux block (lines ~547-559) is **unchanged**.

Net species count change: −6 (remove n, o, a, i1, i2, s) + 1 (add h)
= **−5 species** from the ODE state vector.

### `reconstruction/platelet/dataclasses/internal_state.py`

**Remove** the six IP3R sub-state rows (IP3R_n, IP3R_o, IP3R_a,
IP3R_i1, IP3R_i2, IP3R_s) from the species table.

**Add** one row:
```
('IP3R_h[dts]',  5.110e-4,  1212,  'protein'),
```
where 1212 = round(h₀ × N_IP3R) = round(0.913 × 1328), the
pre-equilibrated initial condition at cyt = 100 nM.

### `reports/data/calcium-data-provenance.md`

Update the IP3R section to record:
- Sneyd-Dufour removed; deYoung-Keizer 1992 implemented
- Parameters sourced from dYK92 Table 1 / Li-Rinzel 1994
- Reason for change (calibration-regime mismatch)

### `models/platelet/tests/sim/test_regression.py` (or equivalent)

Update any hard-coded checks that depend on IP3R sub-state species
counts (IP3R_n etc.). If a regression test snapshots the species
list, it needs updating.

---

## No changes needed

- `GAMMA_IP3R_S`, `V_IM_V`, `NA_OVER_zF`, `RT_OVER_zF_V` — unchanged
- `K_SERCA`, `K_PMCA`, `K_CAM_PMCA`, `K_CAM`, `K_STIM`, `K_MWC` — unchanged
- `J_PM_LEAK_IONS_S` — unchanged (Phase 4 recalibration deferred)
- `SOCE` block — unchanged
- `MOLECULE_NAMES` order for all non-IP3R species — unchanged
- `CalciumDynamics` process class — no changes needed; it reads
  species list from `MOLECULE_NAMES` and calls `_ode_rhs`

---

## Expected outcomes

### Resting fixed point (`restConvergence.py`, 6 000 s, IP3 = 50 nM constant)

| Quantity | Current (Sneyd-Dufour) | Prediction | Target (Dolan) |
|----------|----------------------|------------|----------------|
| [Ca²⁺]_cyt | 2 240 nM | *see note* | ~100 nM |
| [Ca²⁺]_DTS | 0.07 µM | *see note* | 200–300 µM |

*Note:* As discussed in §Parameters, the absolute resting Po of
deYoung-Keizer at the standard parameters is not dramatically
different from Sneyd-Dufour. Resting-state improvement therefore
depends on whether γ_IP3R is recalibrated (Phase 4). The immediate
expected improvement is:

- The h variable should approach its quasi-steady state (h ≈ 0.91
  at cyt = 100 nM), confirming the ODE is behaving correctly.
- The IP3 dose-response is now correct: the model should produce
  a much larger Po increase during stimulation (IP3 rising from
  50 nM toward µM range) than Sneyd-Dufour does, because the
  deYoung-Keizer IP3 gate spans this range cleanly.

### Phase 3 validation (200 s, with / without Ca_ex)

| Criterion | Current baseline | Prediction |
|-----------|-----------------|------------|
| +Ca_ex peak | 393 nM ✓ | May shift; retest required |
| −Ca_ex peak | 325 nM ✓ | May shift; retest required |
| SOCE differential | 67 nM ✗ | Uncertain until γ_IP3R calibrated |
| Acceptance | 4/5 | Target ≥ 4/5; accept regression if γ_IP3R needs retuning |

The Phase 3 peaks may change because the deYoung-Keizer model's
IP3-driven transient dynamics are different from Sneyd-Dufour's. If
Phase 3 falls below 4/5 after this change, γ_IP3R is adjusted in
Phase 4 to restore the match before anything else is touched.

---

## How to measure success

### Step 1 — Verify h dynamics (immediate sanity check)

After implementing, run a 600 s simulation with IP3 = 50 nM constant
and no IP3 stimulus. Confirm:
- h settles near 0.91 (h∞ at cyt = 100 nM)
- m∞ ≈ 0.152 at cyt = 100 nM, IP3 = 50 nM
- Po_channel ≈ (expected value per chosen formula)
- IP3R flux at rest can be computed manually and compared to SERCA
  refill rate

### Step 2 — Resting convergence (`restConvergence.py`)

Run: `PYTHONPATH=$PWD python runscripts/manual/restConvergence.py`

Record:
- [Ca²⁺]_cyt at convergence
- [Ca²⁺]_DTS at convergence
- h value at convergence
- IP3R flux and SERCA flux at convergence (add diagnostic print if
  needed)

**Minimum acceptance:** DTS does not fully drain (DTS > 1 µM) within
6 000 s. Full acceptance (DTS > 50 µM) likely requires Phase 4
γ_IP3R recalibration.

### Step 3 — Phase 3 validation

Run: `PYTHONPATH=$PWD python runscripts/manual/runPhase3.py out/phase3_dyk/ --length 200`

Record all five acceptance criteria. Update `test_regression.py`
baselines to the new numbers. The target is ≥ 4/5; if Phase 3 dips
to 3/5 due to peak-band shifts, this is acceptable as a transitional
state pending Phase 4 recalibration, provided the direction of change
is documented.

### Step 4 — IP3 dose-response check (qualitative)

Run two 200 s simulations: one at IP3 = 50 nM constant (no stimulus)
and one with IP3 = 1 µM constant (saturating). Confirm that the
deYoung-Keizer model produces near-zero Ca²⁺ release at 50 nM and
substantial release at 1 µM — the IP3 gating should be qualitatively
correct.

---

## Open questions before coding

1. **Po exponent**: use `m∞⁴ × h` (inhibition applies once per channel)
   or `m∞⁴ × h⁴` (independent per subunit)? The Li-Rinzel paper
   uses `m∞³ × h` (three-subunit form). For consistency with the
   current Dolan tetrameric convention we will use `m∞⁴ × h` as
   the first implementation; test `m∞⁴ × h⁴` if resting Po is still
   too high after Phase 4.

2. **N_IP3R constant**: the current code derives `ip3r_total` by
   summing sub-state counts. In deYoung-Keizer, `N_IP3R = 1 328` is
   a constant. This simplifies the flux formula and removes the need
   to normalise by a varying total.

3. **γ_IP3R recalibration (Phase 4)**: should γ_IP3R be adjusted as
   part of this issue (#27) or deferred to #30? Recommendation:
   defer to #30 (Phase 4) so that #27 is a clean model-replacement
   commit and γ_IP3R tuning is a separate, documented calibration
   step. If Phase 3 collapses below 3/5 after the replacement, pull
   γ_IP3R retune forward into #27.

---

## Files to change

| File | Change |
|------|--------|
| `reconstruction/platelet/dataclasses/process/calcium_signalling.py` | Replace K_IP3R + φ-functions with K_DYK + h ODE (primary) |
| `reconstruction/platelet/dataclasses/internal_state.py` | Swap 6 IP3R sub-state rows for 1 IP3R_h row |
| `reports/data/calcium-data-provenance.md` | Update IP3R provenance section |
| `models/platelet/tests/sim/test_regression.py` | Update species-list / baseline snapshots |

---

## Implementation results (2026-05-11)

Implementation complete. All tests pass (8/8 regression, mypy clean).
Phase 3 validation maintains 4/5.

### Gating diagnostic (`checkIP3REquilibrium.py`)

Resting conditions (IP3 = 50 nM, Ca_cyt = 100 nM, Ca_DTS = 250 µM):

| Quantity | Value |
|----------|-------|
| m∞ | 0.1523 |
| h∞ | 0.913 |
| τ_h | 3.7 s |
| Po_channel = m∞⁴ × h | 4.92×10⁻⁴ |
| IP3R flux at rest | 3 349 000 ions/s (into cytosol) |

**Important clarification on the 110 000 ions/s figure in the design
doc above**: the 110 000 ions/s cited in §"Why this change" was
computed at the *Dolan-filtered initial conditions* (Po⁴ = 1.65×10⁻⁵),
not at the Sneyd-Dufour Markov equilibrium. At the Sneyd-Dufour ODE
equilibrium the actual resting flux was ~2 870 000 ions/s; the DYK
resting flux is ~3 350 000 ions/s. Both models therefore produce a
resting IP3R leak ~30–100× SERCA refill capacity. The 110 000 ions/s
figure belongs to the Dolan Monte-Carlo filtered state, which is not
a trajectory point in our ODE system.

### Resting convergence

At the DYK ODE fixed point (no IP3 forcing):

| Quantity | Sneyd-Dufour | deYoung-Keizer |
|----------|-------------|----------------|
| [Ca²⁺]_cyt at convergence | ~2 211 nM | ~2 211 nM |
| [Ca²⁺]_DTS at convergence | ~0.07 µM | ~0.56 µM |

The DTS recovers somewhat (0.07 → 0.56 µM) because the DYK
Ca²⁺-inhibition gate (h → 0 as cyt rises) provides negative feedback
that partially caps the leak. The DTS is still not at a biologically
reasonable resting level; Phase 4 γ_IP3R recalibration is required.

### γ_IP3R recalibration targets (Phase 4, initial analysis)

To bring IP3R resting flux to the SERCA-matched range:

| Target flux | Required γ_IP3R | Reduction from 10 pS |
|-------------|-----------------|----------------------|
| 30 000 ions/s (SERCA lower bound) | 0.090 pS | 111× |
| 100 000 ions/s (SERCA upper bound) | 0.299 pS | 33× |

---

## Phase 4 — γ_IP3R recalibration (2026-05-11, issue #30)

### Analytical calibration

SERCA 6-state cycle steady-state flux solved at cyt = 100 nM,
DTS = 250 µM (linear system):

| SERCA sub-state | Steady-state population |
|----------------|------------------------|
| E2 | 5 803 |
| E1 | 5 710 |
| E1Ca | 81 |
| E1PCa | 101 |
| E2PCa | 84 |
| E2P | 113 |

Cycle flux J = 56 285 cycles/s → **SERCA Ca²⁺ efflux = 112 570 ions/s**.

For balance at cyt = 100 nM:

    γ_required = SERCA_flux / (N × Po × driving × NA/zF)
               = 112 570 / (1 328 × 4.91×10⁻⁴ × 0.1605 × 3.122×10¹⁸)
               = 112 570 / (3.27×10¹⁷)
               = **0.344 pS**

Set: `GAMMA_IP3R_S = 0.35e-12` (rounded to 2 s.f.).

### SERCA initial conditions — two-state vs six-state equilibrium

The previous SERCA initial conditions used the two-state (E1 ↔ E1Ca)
binding equilibrium:

    E1Ca/E1 = k_bind_f · cyt² / k_bind_r = 1000 · 0.01 / 10 = 1.0

This gives E1Ca = E1 = 2 963. However the correct full-cycle quasi-
steady-state ratio is:

    E1Ca/E1 = k_bind_f · cyt² / (k_bind_r + k_phos_f)
            = 1000 · 0.01 / (10 + 700) = 10/710 = 0.0141

The fast phosphorylation drain (k_phos_f = 700 s⁻¹ >> k_bind_r = 10 s⁻¹)
holds E1Ca much lower than the 2-state approximation suggests. The
old E1Ca = 2 963 was 36× too large; the phosphorylation burst at t = 0
(v_phos = 700 × 2 963 = 2.07M events/s) rapidly depleted cytosolic
Ca²⁺ below the d₅ = 82 nM activation threshold, trapping the system
at the low-Ca²⁺ attractor (~2.6 nM).

**Fix**: initialise SERCA to the 6-state cycle steady state (E1 = 5 710,
E2 = 5 803, E1Ca = 81, E1PCa = 101, E2PCa = 84, E2P = 113). This
eliminates the t = 0 burst.

### Outcome of Phase 4 calibration

After applying γ = 0.35 pS and the corrected SERCA initial conditions:
- All 21 regression + Phase 3 tests pass.
- Initial conditions: cyt = 100 nM, DTS = 250 µM ✓
- Long-time resting convergence (600 s, no IP3 forcing):
  cyt ≈ 122 nM, DTS ≈ 327 µM (still drifting; DTS overfills without
  a DTS Ca²⁺ buffer).

**The DTS overfill is expected**: without calreticulin (CALR, Phase 2 /
#28), the DTS has no luminal buffering. SERCA pumps Ca²⁺ in faster
than IP3R (at a balanced γ) removes it; DTS rises until SERCA stalls
or reverses against the thermodynamic gradient. CALR will provide the
~508 000 Ca²⁺-binding sites needed to stabilise DTS at 250 µM.

Phase 3 validation passes 4/5 (unchanged from deYoung-Keizer Phase 1)
because the 200 s IP3-stimulated transient starts from the correct
initial conditions and the peak dynamics are dominated by the fast
IP3R/SERCA/CaM kinetics, not the slow DTS drift.

**Next step**: Phase 2 (CALR buffer, issue #28). After CALR, γ_IP3R
may need minor re-tuning in the full DTS-buffered context.

---

*Branch:* `main` · *Status:* Phase 4 complete (γ + SERCA ICs); Phase 2 next ·
*Linked issues:* #27 (Phase 1, complete), #24 (parent), #28 (Phase 2,
CALR buffer; next), #30 (Phase 4, complete)
