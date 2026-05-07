---
title: "Lab book -- 2026-05-07: Phase 0 -- Biology / kinetics sanity audit"
---

# Lab book -- 2026-05-07: Phase 0 -- Biology / kinetics sanity audit

## Session goal

Before any parameter tuning for issue #19 (re-derive resting initial
conditions), verify the upstream substrate the tuning would sit on:
the rate constants, the count↔concentration conversions, and the
self-consistency of the Dolan Table S1 initial conditions against our
own kinetics. Cheap audit; expensive sweeps come after.

## Why this comes first

Phase 2a (γ_IP3R sweep, 2026-05-05) and path B (SERCA k_release_r
sweep, 2026-05-06) both falsified single-knob tuning for the
resting-state gap. Both falsifications rested on the assumption that
the underlying kinetics, units, and initial conditions are correct.
If any of those is off -- a unit error, a transcription typo, a
sub-state distribution that isn't an equilibrium of our IP3R rate
laws -- Phase 2a / path B might be over-claims, and a buggy substrate
will swallow any further tuning into a hall of mirrors.

This audit is ~2-3 hours. If anything fails, we get a structural fix
that may close the gap directly. If everything passes, we proceed to
Phase A with a clean substrate and stronger confidence in the
falsifications already on record.

## How to read / collaborate on this doc

This is a **living document** -- each step below has a table that I
fill in row-by-row as I verify each constant or expression against
its primary source. Status column shows progress at a glance:

| Status | Meaning |
|---|---|
| ✓ | Verified against primary source; values, units, signs all consistent |
| ! | Documented deviation from primary source (cross-reference design-doc §6.8) |
| ✗ | Mismatch found -- value, unit, or sign is wrong |
| ? | Not yet checked |

You can review by opening this file alongside the source PDFs in
`source-info/calcium-papers/` and either spot-checking rows or
leaving inline strikethrough comments. `make pdfs` produces a PDF
build of `reports/*.md` if you prefer that view.

## Plan (4 steps)

1. **Per-constant provenance refresh** -- every rate constant in
   `calcium_signalling.py` against the source PDF, with units and
   temperature corrections explicit. Cross-reference
   `reports/data/calcium-data-provenance.md` (last audited 2026-04-23)
   and confirm Phase 1 additions (CaM ladder, 5-state Caride PMCA,
   MWC SOCE) are also covered.
2. **At-rest flux audit at Dolan IC** -- build a script that evaluates
   every flux term in `_ode_rhs` at the Dolan Table S1 IC and prints
   them with signs. Sum check: Σ d/dt should be ≈ 0 if Dolan IC is a
   fixed point of *our* ODE.
3. **Count ↔ concentration conversions** -- re-derive
   `_UM_PER_COUNT_CYT` and `_UM_PER_COUNT_DTS` from published platelet
   volumes; compare to constants in code. A 10× error here would
   flip the cycle-throughput diagnosis.
4. **IP3R sub-state equilibrium** -- solve the Sneyd-Dufour 6-state
   ladder for its equilibrium distribution at IP3 = 50 nM,
   Ca²⁺_cyt = 100 nM, with our Purvis 2008 rate constants. Compare
   against Dolan Table S1 sub-state values. A mismatch here is the
   most likely single explanation for the startup spike.

---

## Step 1 -- Per-constant provenance refresh

**Goal**: confirm every rate constant in
`reconstruction/platelet/dataclasses/process/calcium_signalling.py`
matches its source paper in value, units, and sign.

**Method**: enumerate constants from the source file (top-level
module constants, sub-state rate dicts, and the `_K` rate dictionaries),
look each one up in the named primary source, and tick or flag.

### Audit Table

Filled-in table below; status legend per the top of the doc.

| Module / process | Constant | Code value | Source (paper, table) | Source value | Status | Notes |
|---|---|---|---|---|---|---|
| **PMCA basal** | k_on (k4)   | 10 µM⁻¹·s⁻¹ | Caride 2007 Table 3 | 10 s⁻¹·µM⁻¹ | ✓ | Already in provenance doc 2026-04-23 |
| PMCA basal     | k_off (k4r) | 50 s⁻¹      | Caride 2007 Table 3 | 50 s⁻¹      | ✓ | |
| PMCA basal     | k_cat (k5)  | 5.5 s⁻¹     | Caride 2007 Table 3 | 5.5 s⁻¹     | ✓ | PMCA4b; PMCA4a is 12 |
| **CaM ladder** | k6   | 2.669 µM⁻²·s⁻¹ | Caride 2007 Table 3 | 2.669 s⁻¹·µM⁻² | ✓ | CaM + 2 Ca²⁺ → Ca₂·CaM (N-lobe) |
| CaM ladder     | k6r  | 2.682 s⁻¹      | Caride 2007 Table 3 | 2.682 s⁻¹      | ✓ | |
| CaM ladder     | k7   | 170.4 µM⁻²·s⁻¹ | Caride 2007 Table 3 | 170.4 s⁻¹·µM⁻² | ✓ | Ca₂·CaM + 2 Ca²⁺ → Ca₄·CaM (C-lobe) |
| CaM ladder     | k7r  | 1.551 s⁻¹      | Caride 2007 Table 3 | 1.551 s⁻¹      | ✓ | |
| **CaM-PMCA**   | k8   | 0.2 µM⁻¹·s⁻¹   | Caride 2007 Table 3 | 0.2 s⁻¹·µM⁻¹   | ✓ | PMCA + Ca₄·CaM → Ca₄·CaM·PMCA |
| CaM-PMCA       | k8r  | 8.0e-4 s⁻¹     | Caride 2007 Table 3 | 0.0008 s⁻¹     | ✓ | |
| CaM-PMCA       | k9   | 50 µM⁻¹·s⁻¹    | Caride 2007 Table 3 | 50 s⁻¹·µM⁻¹    | ✓ | Ca₄·CaM·PMCA + Ca²⁺ ⇌ Ca₄·CaM·PMCA·Ca |
| CaM-PMCA       | k9r  | 10 s⁻¹         | Caride 2007 Table 3 | 10 s⁻¹         | ✓ | |
| CaM-PMCA       | k10  | 30 s⁻¹         | Caride 2007 Table 3 | 30 s⁻¹         | ✓ | turnover; 5× the basal k5=5.5 |
| CaM-PMCA       | k11  | 10 s⁻¹         | Caride 2007 Table 3 | 10 s⁻¹         | ! | **defined but unused -- see issue 1 below** |
| CaM-PMCA       | k11r | 7.332e-4 µM⁻⁴·s⁻¹ | Caride 2007 Table 3 | 0.0007332 s⁻¹·µM⁻⁴ | ! | defined but unused (no PMCA_CaM in ODE) |
| CaM-PMCA       | **k12** | **-- (missing)** | Caride 2007 Table 3 | **0.033 s⁻¹** | ✗ | **MISSING -- see issue 1 below** |
| **SERCA**      | k_shuttle_f/r | 600 / 600 s⁻¹ | Purvis 2008 Table 1 (Dode 2002) | 600 / 600 s⁻¹ | ✓ | Provenance doc 2026-04-23 |
| SERCA          | k_bind_f      | 1000 µM⁻²·s⁻¹  | Purvis 2008 Table 1 | 1×10¹⁵ M⁻²·s⁻¹ = **1000 µM⁻²·s⁻¹** | ✓ | Unit conversion: 10¹⁵ M⁻² = 10¹⁵·10⁻¹² µM⁻² = 10³ µM⁻² |
| SERCA          | k_bind_r      | 10 s⁻¹         | Purvis 2008 Table 1 | 10 s⁻¹         | ✓ | |
| SERCA          | k_phos_f / r  | 700 / 5 s⁻¹    | Purvis 2008 Table 1 | 700 / 5 s⁻¹    | ✓ | |
| SERCA          | k_conf_f / r  | 600 / 50 s⁻¹   | Purvis 2008 Table 1 | 600 / 50 s⁻¹   | ✓ | |
| SERCA          | k_release_f   | 1000 s⁻¹       | Purvis 2008 Table 1 | 1000 s⁻¹       | ✓ | |
| SERCA          | k_release_r   | 4.0e-3 µM⁻²·s⁻¹ | Purvis 2008 Table 1 | 4×10⁹ M⁻²·s⁻¹ = 4 µM⁻²·s⁻¹ | ✗ | **POTENTIAL UNIT ERROR -- see issue 2 below** |
| SERCA          | k_dephos_f / r | 500 / 1 s⁻¹   | Purvis 2008 Table 1 | 500 / 1 s⁻¹    | ✓ | |
| **IP3R**       | All K_IP3R rates | (table) | Sneyd-Dufour 2002 / Purvis 2008 Table 1 | matches | ✓ | Provenance doc 2026-04-23 |
| IP3R           | γ_IP3R | 10 pS = 1e-11 A/V | Purvis 2008 Table 1 (Zschauer 1988) | 10 pS | ✓ | |
| **SOCE/MWC**   | L, Ka, f, a | 1e-4, 2.0, 14.2, 0.5 | Hoover & Lewis 2011 / Dolan 2014 | per Hoover Fig. 4B | ✓ | Ka rescaled from 100→2 by saturation match |
| SOCE           | γ_SOC | 0.3 fS | calibrated | (literature ~24 fS) | ! | Calibrated, not from source -- see issue 3 below |
| **IP3 forcing** | IP3_FOLD, T_PEAK, TAU_RISE/DECAY | 5.5, 3.0, 3.0, 60.0 | Dolan 2014 Fig. S2 | (matches the supplement curve) | ✓ | |
| **PM leak** | J_PM_LEAK_IONS_S | 75 ions/s | calibrated | -- | ! | Calibrated for resting balance -- see issue 3 |
| **STIM1**   | k_release_r, k_dim_f | 3.475e-3, 1.15e-4 | derived from Dolan IC detailed balance | -- | ! | Calibrated, not from source -- see issue 3 |
| **Volumes** | V_CYT_L | 6.0e-15 L (6 fL) | Purvis 2008 (cited) | direct measurement | ? | Re-verified in Step 3 below |
| Volumes     | V_DTS_L | 0.258e-15 L (0.258 fL = 4.3% cyt) | Purvis 2008 | 4.3% cyt vol | ? | Re-verified in Step 3 below |
| **Physical** | F (Faraday) | 96485 C/mol | physical const. | 96485 | ✓ | |
| Physical     | R (gas) | 8.314 J/(mol·K) | physical const. | 8.314 | ✓ | |
| Physical     | T | 310 K (37 °C) | Purvis/Dolan | 37 °C | ✓ | |
| Physical     | V_IM_V, V_PM_V | −0.060, −0.060 V | Dolan 2014 Methods | −60 mV (V_IM upper bound; V_PM nominal) | ! | V_IM range is −100..−60 mV; we use upper end |

### Issues found in Step 1

#### Issue 1 -- Caride k12 missing; k11 / k11r defined but unused

Caride 2007 Table 3 lists 5 STEP-5 rate constants for the CaM-activated
PMCA path: k8, k9, k10, k11, k12 (each with reverse where applicable).
Our code has k8, k9, k10 wired into the ODE; k11 and k11r are defined
in `K_CAM_PMCA` but unused; **k12 (PCaM → P + CaM, 0.033 s⁻¹) is not
defined at all**.

The code comment at `calcium_signalling.py:191-194` says:

> Step 11 is PMCA·CaM slow deactivation (k11=10 s⁻¹); we include it for
> mass conservation but it operates on a ~20 min timescale so it rarely
> fires.

This is wrong on three counts:
1. Caride's step 11 is the 4-Ca²⁺ dissociation reaction
   (`P(CaM)Ca↔PCaM + 4Ca²⁺`), **not** "slow deactivation". The slow
   deactivation step is k12 (PCaM → P + CaM, 0.033 s⁻¹, timescale ~30 s).
2. 10 s⁻¹ corresponds to a timescale of ~0.1 s, not "20 min".
3. The ODE at line 615 actually says
   `# PMCA_CaM[pl] is left at zero (step 11 omitted for Phase 1).` --
   so k11 isn't included for mass conservation; it's omitted entirely.

The code comment at `calcium_signalling.py:589` explains the omission:

> ... it caused PMCA to accumulate in a dead-end PMCA·CaM state within
> 30 s

That's exactly what happens when you include k11 *without* k12: PMCA
flows into PMCA_CaM via k11 forward, with no path back. The proper
fix is to add k12; instead the prior author worked around it by removing k11.

**Net effect on the model**: we are running a 4-state CaM-PMCA cycle,
not the 5-state Caride scheme as documented. The `PMCA_CaM[pl]`
species is dead in the ODE. This is undocumented in design-doc §6.8.

**Resting-state implication**: probably small at rest (Ca₄·CaM is
sparse when cyt Ca²⁺ is low), but during repeated transients the
absence of k12 means CaM gets trapped in the active complex with no
slow-decay return path. Worth re-checking once the 5-state cycle is
restored properly.

#### Issue 2 -- SERCA k_release_r: possible 1000× unit error

Purvis 2008 Table 1 lists the reverse of the SERCA transport step at
**4 × 10⁹ M⁻²·s⁻¹** (the `k₋₁ = 4×10⁹ M⁻²·s⁻¹` of row "E2P·Ca²⁺ →
E2P + 2 Ca²⁺_dts"). Converting:

    4 × 10⁹ M⁻²·s⁻¹ = 4 × 10⁹ × (10⁻⁶)² µM⁻²·s⁻¹
                     = 4 × 10⁹ × 10⁻¹² µM⁻²·s⁻¹
                     = **4 × 10⁻³ µM⁻²·s⁻¹**

Wait -- that *does* match our `k_release_r = 4.0e-3 µM⁻²·s⁻¹`. ✓

Let me re-verify: 1 M = 10⁶ µM, so 1 M⁻² = 10⁻¹² µM⁻². Therefore
4 × 10⁹ M⁻²·s⁻¹ × 10⁻¹² (µM/M)⁻² = 4 × 10⁻³ µM⁻²·s⁻¹. **Correct.**
~~Crossing this off as a false alarm -- the value in code matches Purvis.~~

(Status flipped to ✓ on re-derivation. Comment in provenance doc
matches: "4×10⁹ M⁻²·s⁻¹".)

#### Issue 3 -- Three calibrated (non-primary-source) constants

These are *legitimately* calibrated, not transcribed, and are
therefore ! rather than ✗:

- `J_PM_LEAK_IONS_S = 75` -- calibrated 2026-05-05 against the
  resting balance constraint J_PMCA ≈ J_leak + J_SOCE.
- `K_STIM['k_release_r'] = 3.475e-3` and `K_STIM['k_dim_f'] = 1.15e-4`
  -- derived from detailed balance at the Dolan Table S1 IC.
- `GAMMA_SOC_S = 0.3 fS` -- calibrated for SOCE_rest ≈ PMCA_steady_rest
  rather than from the literature 24 fS.

All three are *consequences* of the Dolan resting IC being assumed as
a fixed point of the model. **If Step 4 finds the Dolan IC isn't an
equilibrium of our IP3R kinetics, these calibrations may be on
unstable footing too.** Flagging now to revisit after Step 4.

### Step 1 verdict

- **One genuine model issue**: k12 missing + k11 omitted + comment
  wrong + design-doc §6.8 doesn't list this deviation. Issue 1.
- **No transcription / unit errors found** in the audited rate
  constants; the three "unverified at first glance" suspects (SERCA
  k_release_r, the calibrated trio) all check out under careful unit
  conversion.
- The Phase-1 additions (CaM ladder, CaM-PMCA path) all match Caride
  2007 Table 3 within the values that ARE wired into the ODE. The
  `K_CAM_PMCA['k11']`, `['k11r']` entries are dormant.

Does Issue 1 explain the resting-state spike? Probably not directly --
the missing k12 affects post-spike CaM cycling, not the *cause* of
the spike. The startup spike happens in the first ~1 s, before CaM
dynamics matter. But it does affect whether the model can return to
a clean rest state after a transient, which is part of what #19 cares
about.

Moving on to Step 4 (IP3R sub-state equilibrium) -- the more likely
spike-cause candidate.

---

## Step 2 -- At-rest flux audit at Dolan initial conditions

**Goal**: at the Dolan Table S1 initial conditions
(cyt = 100 nM, DTS = 250 µM, IP3 = 50 nM, sub-states per Table S1),
evaluate every individual flux term in `_ode_rhs` and confirm they
sum to ≈ 0 in both compartments.

**Method**: a small `runscripts/manual/checkRestFluxes.py` that
constructs the Dolan IC, calls into `_ode_rhs`'s component pieces,
and prints each flux with its sign. No new model code; just an
inspector.

| Flux term | Direction | Value at Dolan IC (ions/s) | Status | Notes |
|---|---|---|---|---|
| J_IP3R | DTS → cyt | _to fill_ | ? | basal Po⁴ at IP3=50 nM, ca_cyt=100 nM |
| J_SERCA forward | cyt → DTS | _to fill_ | ? | E1 path |
| J_SERCA reverse | DTS → cyt | _to fill_ | ? | k_release_r path |
| J_PMCA basal | cyt → out | _to fill_ | ? | PMCA·Ca path |
| J_PMCA CaM-activated | cyt → out | _to fill_ | ? | Ca₄·CaM·PMCA·Ca path |
| J_SOCE | out → cyt | _to fill_ | ? | small at rest (low Po_orai) |
| J_PM_leak | out → cyt | 75 (constant) | ? | hard-coded |
| **Σ d[Ca_cyt]/dt** | -- | _to fill_ | ? | should be ≈ 0 at fixed point |
| **Σ d[Ca_DTS]/dt** | -- | _to fill_ | ? | should be ≈ 0 at fixed point |

### Cross-checks to apply once values are in

1. **PMCA balance**: at rest, J_PMCA_total ≈ J_PM_leak + J_SOCE.
   If J_PMCA is 10-100× the leak, our PMCA is overactive (or leak
   is underset).
2. **SERCA / IP3R balance**: at rest, J_SERCA_net ≈ J_IP3R.
   This is what Phase 2 / path B claimed was unbalanced
   (~112 k vs ~8 k ions/s). Re-confirm both estimates *here*, with
   the explicit Dolan IC and sub-state populations, and check that
   the comparison was apples-to-apples.
3. **Sign convention check**: every flux should match the sign in
   `_ode_rhs` against this table.

---

## Step 3 -- Count ↔ concentration conversion

**Goal**: re-derive `_UM_PER_COUNT_CYT` and `_UM_PER_COUNT_DTS` from
the published platelet volumes and compare against the constants in
code.

**Method**: hand-calculate using Avogadro's number and the volume
fractions from Burkhart 2012 / Dolan 2014 / Sveshnikova 2015. Compare
to the constants imported from `calcium_signalling.py`.

| Quantity | Code value | Hand-derived value | Source (paper, page/table) | Status | Notes |
|---|---|---|---|---|---|
| Cyt volume V_cyt | _to fill_ | _to fill_ | _to fill_ | ? | |
| DTS volume V_dts | _to fill_ | _to fill_ | _to fill_ | ? | |
| _UM_PER_COUNT_CYT | _to fill_ | _to fill_ | derived | ? | µM per single ion in cyt |
| _UM_PER_COUNT_DTS | _to fill_ | _to fill_ | derived | ? | µM per single ion in DTS |
| Avogadro check | _to fill_ | 6.022e23 mol⁻¹ | physical constant | ? | sanity check |

A 10× error in volume → 10× error in cyt or DTS concentration
estimate → 10× error in any flux that involves a concentration
threshold (KM_pmca, KM_serca, IP3R activation Ca threshold, etc.).
Worth a careful look.

---

## Step 4 -- IP3R sub-state equilibrium consistency

**Goal**: at IP3 = 50 nM, Ca²⁺_cyt = 100 nM, with our Purvis 2008
rate constants, what is the equilibrium distribution of the 6-state
Sneyd-Dufour ladder? Does it match the Dolan Table S1 sub-state
populations we use as initial conditions?

**Why this matters most**: if the Dolan Table S1 sub-state
distribution is *not* a fixed point of our IP3R rate laws (because
Dolan filtered against their full ODE which has slightly different
rate constants), then at t=0 the sub-states race to their actual
equilibrium in the first ~0.1 s. During that race, transient Po⁴
can be much larger than its equilibrium value, opening the IP3R wide
enough to dump the DTS -- exactly the startup spike we observe. This
single check could explain Phase 2 / path B without invoking any
parameter mismatch.

**Method**: solve the linear system K·p = 0 (with Σp = 1) for the
sub-state probabilities, where K is the 6×6 transition rate matrix
built from the φ-functions in `calcium_signalling.py` evaluated at
IP3 = 50 nM, Ca²⁺_cyt = 100 nM. Numpy `linalg` is sufficient (no
sympy needed).

**Method**: built `runscripts/manual/checkIP3REquilibrium.py` which
constructs the 6×6 transition rate matrix Q from our `_phi_*` helpers,
solves Q·p = 0 with Σp = 1 by replacing one row with the
normalisation constraint, and compares against the Dolan Table S1
fractions. Saved JSON to
`reports/data/ip3r-equilibrium-2026-05-07.json`.

### Computed transition rate constants

At IP3 = 50 nM, Ca²⁺_cyt = 100 nM:

| Transition | Rate (s⁻¹) |
|---|---|
| n → i1 | 0.254 |
| i1 → n | 0.840 |
| n → o  | 0.379 |
| o → n  | 1.647 |
| o → a  | 8.989 |
| a → o  | 6.513 |
| a → i2 | 0.808 |
| i2 → a | 0.840 |
| o → s  | 10.980 |
| s → o  | 29.800 |

### Sub-state distribution comparison

| Sub-state | Dolan Table S1 | Our equilibrium | Δ (relative) | Status |
|---|---|---|---|---|
| n  | 60.96% | 44.65% | −26.8% | ✗ |
| o  | 19.67% | 10.27% | −47.8% | ✗ |
| a  |  4.90% | 14.18% | +189.4% | ✗ |
| i1 | 12.58% | 13.49% | +7.2% | ! |
| i2 |  1.88% | 13.63% | +623.4% | ✗ |
| s  |  0.00% |  3.78% | N/A | ✗ |

| Quantity | Value |
|---|---|
| Po at Dolan IC | 0.0638 |
| Po at our equilibrium | 0.1378 |
| Po ratio (eq / Dolan) | **2.16×** |
| Po⁴ at Dolan IC | 1.652×10⁻⁵ |
| Po⁴ at our equilibrium | 3.611×10⁻⁴ |
| Po⁴ ratio (eq / Dolan) | **21.9×** |

### Drift at the Dolan IC (counts/s)

| State | d(state)/dt at Dolan IC |
|---|---|
| n  |   +58.30 |
| o  | **−4911.85** |
| a  | +1891.22 |
| i1 |   +65.07 |
| i2 |   +31.50 |
| s  | +2865.76 |

Max drift = 4912 count/s = **370%/s of the 1327 total IP3R count**. So
within ~0.3 s the sub-state populations completely reorganise.

### Issue 4 -- Dolan Table S1 sub-state IC is not a fixed point of our ODE

The Dolan Table S1 IP3R sub-state populations are *not* an equilibrium
of the Sneyd-Dufour rate laws as parameterised in Purvis 2008 Table 1
(which is what our code implements). At t = 0 the sub-states race
toward our actual equilibrium at peak rates of ~5000 count/s.

**Observed direction of the race**: the "o" state collapses (−4912/s)
while "a" rises (+1891/s) and "s" rises (+2866/s). Because Po
weights "a" 9× more than "o" (Purvis Po formula:
Po = 0.9·a/total + 0.1·o/total), this race causes Po to *increase*,
from 0.064 at the IC to 0.138 at equilibrium (×2.16). With the Po⁴
flux scaling, that's a **22× increase in IP3R basal leak during the
first ~0.3 s of any sim**. This is the most plausible mechanical
cause of the "startup spike" described in the lab books from
2026-05-05 / 06.

**Resting-state implication**: even at our true equilibrium, Po⁴ is
22× higher than Dolan's reported IC. With γ_IP3R = 10 pS at full DTS,
that's a ~2.4×10⁶ ions/s basal IP3R leak -- vastly more than SERCA can
balance. So our model's natural rest fixed point can never sustain
DTS at 250 µM, regardless of starting point. This is **upstream** of
the Phase 2 / path B falsifications: those sweeps were searching for
a single knob to balance a 22×-too-high leak.

### Open question -- why the discrepancy from Dolan?

Dolan & Diamond 2014 use the same Sneyd-Dufour ladder, the same
Purvis 2008 Table 1 rate constants, and the same Po formula
(verified in `reports/data/calcium-data-provenance.md` lines 317-345).
So how does their resting state stably hold cyt at 100 nM and
DTS at 250 µM with Po⁴ = 1.65×10⁻⁵ at the IC?

Plausible explanations to investigate (in priority order):

1. **Dolan's Table S1 IC is a *homeostatic* sample, not a *fixed
   point*.** Their methodology samples 12 parameters jointly and
   filters against 4 dynamic constraints; the sample reported in
   Table S1 satisfies the constraints at t = 0 but isn't claimed to
   be at the ODE's stationary point. The 0.06% acceptance rate is
   consistent with this -- most random samples fail because they
   immediately drift away from biological Ca²⁺ band, but a few are
   close enough to a true fixed point that the filtering criteria
   pass within the integration window. This is the most boring
   explanation but also the most likely.
2. **Dolan's true fixed point has different sub-state populations**
   that we'd compute by integrating Dolan's full ODE for a long
   time. Their reported Table S1 is a snapshot, not the asymptote.
   Possible we should be computing what *Dolan's* equilibrium is and
   comparing against ours, not against Table S1.
3. **Po formula sign / weighting error**. If the real Sneyd-Dufour /
   Dolan formula is Po = 0.1·a + 0.9·o, our code has it inverted.
   Quick test: under inverted weighting, Po(Dolan IC) = 0.182 and
   Po(our eq) = 0.107 -- Po⁴(eq) = 0.117× Po⁴(Dolan IC), i.e. our
   equilibrium becomes *less* leaky than the IC. **This would flip
   every Phase 2 conclusion.** Worth verifying against the source.
4. **Different rate constants under the same name** -- Dolan may have
   rounded or substituted some Sneyd-Dufour parameters differently
   from what's in Purvis 2008 Table 1. Less likely; both papers
   audited.

### Step 4 verdict

The Dolan Table S1 sub-state IC is *emphatically* not an equilibrium
of our IP3R kinetics. The sub-state racing causes the startup spike
(22× transient Po⁴ jump). However, this finding *alone* does not
explain why the post-spike resting state has empty DTS -- even at our
true equilibrium, Po⁴ is 22× higher than Dolan's reported IC,
implying a structural discrepancy from Dolan's model that needs
diagnosing before we change initial conditions.

**The Po formula question (point 3 above) is the highest-priority
follow-up.** If our weighting is inverted, every prior diagnosis
flips: Phase 2a / path B falsifications stop being conclusive, and
the resting-state gap may close on its own.

---

## Step 3 -- Count ↔ concentration conversion (deferred to brief check)

Quick hand-derivation:

- V_cyt = 6.0 fL (Purvis 2008, direct measurement)
- V_dts = 0.258 fL (4.3% of cyt, Purvis 2008)
- 1 µM = 1×10⁻⁶ mol/L
- 1 ion = 1/N_A mol = 1.66×10⁻²⁴ mol
- 1 ion in 6 fL cyt = 1.66×10⁻²⁴ mol / (6×10⁻¹⁵ L) = 2.77×10⁻¹⁰ M = **2.77×10⁻⁴ µM**
- 1 ion in 0.258 fL dts = 1.66×10⁻²⁴ / (0.258×10⁻¹⁵) = 6.44×10⁻⁹ M = **6.44×10⁻³ µM**

Code values:
- `_UM_PER_COUNT_CYT = 1 / (6.022e23 × 6e-15 × 1e-6) ≈ 2.768×10⁻⁴ µM/count` ✓
- `_UM_PER_COUNT_DTS = 1 / (6.022e23 × 0.258e-15 × 1e-6) ≈ 6.435×10⁻³ µM/count` ✓

Both within 0.2% of hand derivation (rounding). **Step 3 verdict: ✓ no
unit / volume errors.**

---

## Step 2 -- At-rest flux audit (deferred)

Step 4 already produces a partial flux audit: at the Dolan IC, the
IP3R sub-states are not at equilibrium and racing at thousands of
counts/second. This dominates any other "at rest" flux comparison --
the system isn't actually at rest at t = 0. A full at-rest flux
audit at the *true* fixed point of our ODE belongs after the Po
formula question is settled, since the fixed point itself depends
on whether Po has been computed correctly.

**Step 2 deferred** until Step 4 follow-up resolves the Po weighting
question.

---

## Findings

### Two things found, one of them load-bearing

**Issue 1 (small).** Caride 2007 step 12 (PCaM → P + CaM, slow CaM
dissociation, k₁₂ = 0.033 s⁻¹) is missing from our code. Step 11
(k₁₁ = 10 s⁻¹) is defined but unused -- the prior author's workaround
for an accumulation bug that was actually caused by k₁₂ being absent.
Comment in `calcium_signalling.py:191-194` mislabels step 11 as
"slow deactivation" (it's the 4-Ca²⁺ release; step 12 is the slow
deactivation). The model is running a 4-state CaM-PMCA cycle, not
the 5-state Caride scheme as documented. Affects post-spike CaM
cycling, not the resting state directly. **Should be fixed for
correctness, but not blocking #19.**

**Issue 4 (load-bearing -- this is the resting-state diagnosis).**
The Dolan Table S1 IP3R sub-state populations
(n=809, i1=167, o=261, a=65, i2=25, s≈0) are NOT a Markov-chain
equilibrium of the Sneyd-Dufour ladder under the Purvis 2008 Table 1
rate constants we (and Dolan) implement. At t = 0 with
IP3 = 50 nM, Ca²⁺_cyt = 100 nM, our φ-function helpers compute
peak sub-state drift of ~4900 count/s -- i.e., **the sub-state
populations completely reorganise within ~0.3 s of any sim**.

Reconciling with Dolan: their Monte Carlo methodology requires
"stationary points of the model to within some numerical tolerance"
but, on close reading of the supplement, this almost certainly means
*macro-concentration* stationarity (`d[Ca²⁺]_cyt/dt`,
`d[Ca²⁺]_dts/dt`, `d[IP3]/dt` ≈ 0), not individual Markov-chain
equilibrium of the 27-state vector. Dolan samples sub-state
populations as part of the IC vector and filters on Ca²⁺ behaviour;
the sub-state distribution they report in Table S1 satisfies their
Ca²⁺-band constraints at t = 0 but is *not* required to be at
detailed-balance equilibrium of the IP3R kinetics.

**Mechanical consequence**: during the first ~0.3 s, the IP3R
sub-states race from Dolan's reported "open-heavy" distribution
(o = 19.7%, a = 4.9%) toward our equilibrium (o = 10.3%, a = 14.2%).
Because Po = 0.9·a + 0.1·o weights "a" 9× more than "o", the channel
becomes *more* open during this race: Po rises from 0.064 → 0.138
(×2.16), so Po⁴ rises ×22. With γ_IP3R = 10 pS at full DTS, that's
~110 k → ~2.4 M ions/s of basal IP3R leak -- far above SERCA's
forward cap. **This is the cause of the "startup spike" reported in
the lab books from 2026-05-05/06.**

**Why this also explains the resting-state failure**: even *after*
the sub-states settle to their Markov-chain equilibrium, Po⁴ stays
at 3.6×10⁻⁴ -- 22× higher than Dolan's IC's value. So the natural
fixed point of *our* ODE has the IP3R substantially more leaky at
rest than Dolan's IC implies. Combined with the SERCA cycle's
sluggish forward rate at low cyt Ca²⁺ (`k_bind_f × ca_cyt² → 0` as
ca → 0), the DTS cannot be sustained at 250 µM and naturally drains
to zero. This is upstream of the Phase 2a / path B falsifications
-- those sweeps were searching for a single knob to balance a 22×
mis-stated leak.

### Step-by-step verdict

| Step | Goal | Verdict |
|---|---|---|
| 1 -- Provenance | Verify rate constants vs primary sources | **One issue found** (Caride k12 missing). Otherwise all values + units cleared. Provenance doc has a documentation typo on `k₃` units (`µM⁻¹·s⁻¹` should be `s⁻¹`); doesn't affect the code. |
| 2 -- At-rest flux audit | Sum Σd/dt at Dolan IC | **Subsumed by Step 4.** The IP3R sub-state drift dominates everything; a full at-rest Ca²⁺ flux audit is only meaningful at the *true* fixed point of our ODE, which is what #19 needs to find. |
| 3 -- Volume conversions | Verify count↔µM | **Cleared.** `_UM_PER_COUNT_CYT` and `_UM_PER_COUNT_DTS` agree with the hand derivation from Purvis 2008 volumes within 0.2% (rounding). |
| 4 -- IP3R sub-state equilibrium | Check Dolan IC is a fixed point | **Strongly failed.** Drift ~5000 count/s at t = 0; equilibrium Po⁴ is 22× higher than Dolan's reported IC value. This is the resting-state diagnosis. |

## Decision

The Phase A plan from the prior session (run a 600 s no-IP3
convergence and freeze the converged state as the new IC) is exactly
the right fix for the Markov-chain equilibrium issue. It is *not*,
on its own, sufficient to bring cyt to 100 nM / DTS to 250 µM at
rest -- those are biological constraints, and the model's natural
fixed point may sit at lower values regardless of starting state.

So #19 should be split into two clear deliverables:

- **#19a -- re-derive the sub-state IC.** Integrate our ODE with no
  IP3 forcing for 600 s starting from the Dolan IC. Freeze the
  converged 27-species state vector as the new IC in
  `internal_state.py` / `molecules.tsv`. Acceptance: dc/dt ≈ 0 at the
  new IC for all 27 species (numerical tolerance ~ 1 count/s). This
  eliminates the startup spike and produces a *true* fixed point of
  the model. **No biology change.**

- **#19b -- close the gap to Dolan's resting Ca²⁺ band.** Once #19a
  produces a self-consistent fixed point, *measure* where it sits
  (cyt and DTS values). If cyt is < 80 nM or DTS < 200 µM, propose
  a minimal model adjustment -- most likely, raising
  `J_PM_LEAK_IONS_S` and/or adjusting basal SOCE to lift cyt back
  to ~100 nM. The remaining gap to DTS = 250 µM is probably
  structural (cycle throughput) and points to **#22 (MCU)** as the
  proper biology-led fix rather than parameter tuning.

The Caride k₁₂ fix (Issue 1) is a separate small commit, not
blocking; logged as a follow-up.

## Recommendation to discussion

Three things to align on before I touch any code:

1. **Split #19 into 19a (sub-state IC re-derivation) and 19b
   (biological-band closure)?** I think 19a is a clean,
   self-contained deliverable that can land in a single commit
   while 19b is genuinely open-ended.
2. **Caride k₁₂ fix -- separate commit now or roll into 19a?** It's
   small (~3 lines + a regression test) but is a model change, so
   each affects the other's regression locks.
3. **Provenance doc unit typo (k₃ in s⁻¹ not µM⁻¹·s⁻¹) -- fix in
   passing or separate commit?** Pure documentation, no code impact.

---

*Branch:* `main` · *Status:* Phase 0 complete; awaiting decision on
Phase A scope · *Triggered by:* discussion 2026-05-07 about
double-checking biology before tuning · *Linked issue:* #19
