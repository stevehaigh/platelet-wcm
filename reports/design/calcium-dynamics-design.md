---
title: "v0.2 Calcium Dynamics — Design Document"
---

# v0.2 Calcium Dynamics — Design Document

**Issues:** #24 (data + dataclass), #25 (process, including 6-state IP3R — absorbs #43), #26 (listener), #27 (analysis)
**Branch:** `platelet`
**References:** `reports/calcium-data-provenance.md`, `reports/calcium-signalling-pathway-design.md`

---

## 1. Scope

v0.2 adds the first real biochemistry to the platelet model: a **Ca²⁺ transient**
driven by IP3. A Ca²⁺ transient is the characteristic pulse of elevated cytosolic
calcium that occurs when a platelet is activated — cytosolic Ca²⁺ rises sharply
from ~100 nM at rest to a peak of 300–500 nM, then decays over several minutes.
The goal is to reproduce the Dolan & Diamond (2014) Figure 4 response curve — a
sharp peak followed by a sustained plateau driven by store-operated calcium entry
(SOCE) — using the ODE parameters published in that paper.

![Schematic Ca²⁺ transient — expected output shape](/Users/steve/github/platelet-wcm/reports/figures/ca2-transient-reference.png)
*Schematic showing the expected Ca²⁺ transient shape. The IP3 forcing function
(lower panel) drives a rapid peak in cytosolic Ca²⁺, which then decays to a
SOCE-sustained plateau before returning to baseline. This is the shape we aim
to reproduce from Dolan & Diamond 2014 Fig. 4.*

This is a **Dolan-core-first** strategy — meaning we implement the validated
Ca²⁺ core from Dolan & Diamond (2014) first, using their published parameters,
before adding the upstream receptor cascade that would generate IP3 in a real
cell:

```
v0.2:  IP3 forcing → Ca^2^+ core    (IP3R + SERCA + PMCA + SOCE)
v0.3:  Upstream receptor cascade  (P2Y1 → Gq → PLCbeta → IP3 production)
v0.4:  P2Y12 modulation           (Gi → AC → cAMP → PKA → IP3R inhibition)
```

Each milestone is independently testable.

### Why Dolan 2014 and not Purvis 2008?

Purvis (2008) is the foundational modelling paper for platelet Ca²⁺
signalling and provides the kinetic parameters for the IP3R 6-state model and
many of the rate constants we use. However, Dolan & Diamond (2014) published an
updated model specifically calibrated to match experimental Ca²⁺ transients in
intact platelets, including the SOCE component (STIM1/Orai1) which was not
characterised when Purvis wrote. Dolan Table S1 provides a complete, consistent
set of initial conditions for all 22 state variables, validated against their
Figure 4 experimental curves.

In practice, both papers are used:

- **Purvis 2008** supplies IP3R kinetics (Sneyd & Dufour 2002 rate constants),
  SERCA cycle parameters, and compartment volumes (6 fL cytosol, 4.3% DTS).
- **Dolan 2014** supplies the complete ODE system, SOCE parameters, STIM1/Orai1
  copy numbers, and the primary validation target (Fig. 4).
  
Dolan is the primary reference because it is more recent, includes SOCE, and
provides the initial conditions as a self-consistent set.

### 1.1 Framework usage philosophy

The Karr/wcEcoli framework is used here as a **software architecture scaffold**,
not as a source of biological processes. We carry over:

- Process / State / Listener base classes and the time-stepping loop
- `BulkMolecules` / `UniqueMolecules` state containers and partitioning
- `TableWriter` / `TableReader` columnar I/O
- The webapp, runscripts, and analysis infrastructure

We do **not** carry over any E. coli biological processes (transcription,
translation, metabolic network, DNA replication). Platelets are anucleate —
those processes don't exist. The platelet-specific processes are entirely new:
`CalciumDynamics` (v0.2), granule secretion, integrin signalling (future).
The E. coli model is a useful architectural template, not a biological one.

---

## 2. Signal pathway

See `reports/calcium-signalling-pathway-design.md` for the full biology.

### 2.1 Trigger and IP3 generation

Platelet Ca²⁺ signalling is initiated when an agonist — typically ADP released
from a damaged vessel wall, or thromboxane A₂ produced by the platelet itself —
binds a surface GPCR. The primary receptor for ADP-driven Ca²⁺ signalling is
**P2Y1**, a Gq-coupled receptor. Ligand binding activates the Gq α-subunit,
which stimulates **phospholipase Cβ (PLCβ)**. PLCβ cleaves the membrane
phospholipid PIP₂ (phosphatidylinositol 4,5-bisphosphate) into two second
messengers: **IP3** (inositol 1,4,5-trisphosphate), which enters the cytosol,
and DAG (diacylglycerol), which remains membrane-bound and activates PKC. IP3
is the key trigger: it binds the IP3 receptor (IP3R) on the DTS membrane and
gates Ca²⁺ release.

**In v0.2, IP3 is not produced by this pathway.** It is instead driven by a
pre-programmed time curve matching the Dolan 2014 Fig. S2 shape (§3.2). The
full upstream cascade (P2Y1 → Gq → PLCβ → IP3) is scheduled for v0.3.

### 2.2 Ca²⁺ signal flow

At rest, cytosolic Ca²⁺ is ~100 nM; DTS Ca²⁺ is ~250 µM (a ~2,500-fold
gradient maintained by SERCA). On stimulation:

```
Agonist → P2Y1 → Gq → PLCbeta → IP3   [v0.3; forced time curve in v0.2]
                                  │
                                  v
IP3R [DTS membrane; 6-state Markov model]
  Ca^2^+ floods cytosol (peak ~300-500 nM)
    │
    |- SERCA [DTS membrane; E1/E2 cycle]     2 Ca^2^+/ATP; refills DTS store
    │
    |- PMCA  [plasma membrane; 5-state]      1 Ca/ATP; ejects Ca from cell
    │    `- CaM [cytosolic Ca^2^+ buffer]      Ca_4.CaM activates PMCA ~5x
    │
    `- DTS depletion
         `- STIM1 EF-hand releases DTS Ca^2^+
              `- STIM1 monomers dimerise → active sensor
                   `- STIM1 dimers translocate to ER-PM puncta
                        `- gates Orai1 [plasma membrane; MWC allosteric]
                             SOCE: extracellular Ca^2^+ enters cytosol

Basal / resting:
  Constant PM Ca^2^+ leak (~75 ions/s; TRPC / NCX-reverse / constitutive)
    ←→  PMCA basal extrusion + minimal SOCE (full DTS; few STIM1 dimers)
```

### 2.3 Implementation status

| Component | Biology | v0.2 status |
|-----------|---------|-------------|
| IP3 production | P2Y1 → Gq → PLCβ → IP3 | **Forced time curve** (§3.2); real upstream cascade in v0.3 |
| IP3R | 6-state Sneyd & Dufour Markov; biphasic Ca²⁺ activation + inhibition | **Implemented** |
| SERCA | E1/E2 enzymatic cycle; 2 Ca²⁺ / ATP | **Implemented** |
| PMCA | 5-state CaM-coupled scheme (basal steps 4–5; CaM-activated steps 8–10) | **Implemented** (Phase 1) |
| Calmodulin | Ca₂·CaM → Ca₄·CaM ladder; cytosolic Ca²⁺ buffer + PMCA activator | **Implemented** (Phase 1) |
| SOCE | STIM1 sensor cycle + Hoover/Dolan MWC + Orai1 conductance | **Implemented** (Phase 1) |
| Basal PM Ca²⁺ leak | Background permeability (TRPC / NCX-reverse / constitutive) | **Implemented** (Phase 2a) |
| P2Y12 modulation | Gi → ↓cAMP → ↓PKA → IP3R sensitisation | **v0.4** (not started) |

### 2.4 Key compartments

| Compartment | Volume | Ca²⁺ at rest | Role |
|-------------|--------|--------------|------|
| Cytosol | 6.0 fL (Purvis 2008 direct measurement) | ~100 nM (361 molecules) | Signal integration; all soluble processes |
| DTS (dense tubular system) | 0.258 fL = 4.3% of cytosol (Purvis 2008 glucose-6-phosphatase staining) | ~250 µM (38,842 molecules) | Intracellular Ca²⁺ store; ER equivalent |
| Extracellular / OCS | infinite reservoir | 1.2 mM (fixed; Dolan 2014) | SOCE Ca²⁺ source; PMCA sink |

---

## 3. ODE system

### 3.1 State variables and the integer-count problem

**Background — why integer counts matter here.** The wcEcoli framework tracks
all molecules as discrete integer counts (e.g. 361 Ca²⁺ ions) rather than as
continuous concentrations (e.g. 100 nM). This is the *integer-count problem*:
when the number of molecules of a species is very small, treating them as a
continuous number is inaccurate, and rounding to the nearest integer at each
timestep introduces noise. The point at which this noise becomes negligible and
continuous mathematics is a safe approximation is called the *continuum limit* —
typically around 1,000+ molecules of a given species.

**Why this matters for Ca²⁺.** Cytosolic Ca²⁺ at rest is ~100 nM in a ~6 fL
cell, which works out to only ~361 molecules. That is borderline. During a
transient, the count rises to ~3,600 at peak (~1 µM), which is comfortable, but
then falls back through the borderline range as the transient decays. The DTS
store (~38,000 Ca²⁺ ions) is well above the continuum limit and not a concern.

**Cytoplasmic Ca²⁺ is at the continuum limit.** At rest (~100 nM), the count
depends on which cytoplasmic volume we adopt:

| Source | V_cyt | Ca²⁺ at 100 nM | Ca²⁺ at peak (~1 µM) |
|--------|-------|----------------|----------------------|
| Purvis 2008 (direct measurement) | 6.0 fL | 361 | ~3,600 |
| Sveshnikova 2025 | 3.0 fL | 181 | ~1,800 |
| Minimum plausible | 2.0 fL | 120 | ~1,200 |

All three are borderline for deterministic treatment. For context, the rule of
thumb for Gillespie vs ODE is usually ~1,000 molecules — we are near or below
that at rest.

> **ODE** (Ordinary Differential Equation): a mathematical description of how
> concentrations change continuously over time. Standard in biochemical modelling;
> fast to solve numerically. Assumes species counts are large enough to treat as
> continuous numbers.
>
> **Gillespie algorithm**: an exact stochastic simulation method that fires
> individual reaction events one at a time, drawn from probability distributions.
> Correct at any copy number, but computationally expensive — the simulation cost
> scales with the total reaction rate, which can be enormous for fast enzymes.
>
> **Tau-leaping**: a Gillespie approximation that fires multiple reaction events
> in one leap, trading some accuracy for much better speed. Suitable for systems
> with a mix of fast and slow reactions.

**Decision: deterministic ODE sub-stepper for v0.2.**

Rationale:

- Gillespie for Ca²⁺ dynamics is computationally impractical: the SERCA cycle
  has transitions with rates up to 1,000 s⁻¹ × 11,892 enzymes — the algorithm
  would fire millions of reaction events per simulated second.
- Tau-leaping (the standard Gillespie approximation) is a viable middle ground
  but adds significant implementation complexity for a first version.
- The ODE gives the correct *mean* behaviour. The quantisation noise when
  rounding back to integers is ~0.3–0.8% per timestep at resting concentrations
  — this is below biological measurement uncertainty.
- Both Purvis 2008 and Sveshnikova 2025 use stochastic simulation, but they
  focus on cell-to-cell variability in population studies. For single-cell
  mean dynamics (our validation target), deterministic ODEs are standard.

**Flag for v0.3:** once the upstream receptor module is added, the PLC-Gq
complex (~1 molecule at any time, Sveshnikova 2025) is a genuine stochastic
bottleneck that will require tau-leaping or a hybrid approach.

**Volume decision:** Use V_cyt = 6.0 fL (Purvis 2008, direct measurement).
This gives the most conservative count (~361 at rest) and is the basis for
the Dolan 2014 parameters we are adopting. Document the Sveshnikova discrepancy
in the analysis.

The implementation currently carries **27 ODE state variables** — the 22 in
the original v0.2 design plus 6 added in Phase 1 (CaM Ca²⁺-binding ladder
and PMCA·CaM complex sub-states) — and one fixed constant (extracellular
Ca²⁺). Two initial counts differ from Dolan 2014 Table S1; both deviations are noted
in the table and explained in §6.8.

> **Code:** species ordering for the ODE state vector is
> `MOLECULE_NAMES` at
> `reconstruction/platelet/dataclasses/process/calcium_signalling.py:55`.
> Initial counts and per-molecule masses live in `_MOLECULES` at
> `reconstruction/platelet/dataclasses/internal_state.py:38`.

| Variable | Description | Compartment | Resting count | Source / note |
|:-------------------|:------------------------|:--------------|:-----------------|:-----------------------------------|
| `CA2_CYT` | Cytosolic free Ca²⁺ | `[c]` | **361** (100 nM × 6 fL) | Purvis 2008 |
| `CA2_DTS` | DTS stored Ca²⁺ | `[dts]` | 38 842 (250 µM × 0.258 fL) | Dolan 2014 (Fluo-5N) |
| `CA2_EX` | Extracellular Ca²⁺ | `[e]` | fixed (1.2 mM reservoir) | Dolan 2014 |
| `IP3` | Inositol trisphosphate | `[c]` | 181 (50 nM × 6 fL) | Sveshnikova 2025 / Dolan 2014 |
| `IP3R_n` | IP3R neutral | `[dts]` | 809 | Dolan 2014 Table S1 |
| `IP3R_o` | IP3R open | `[dts]` | 261 | Dolan 2014 Table S1 |
| `IP3R_a` | IP3R active (Ca²⁺-bound, conducting) | `[dts]` | 65 | Dolan 2014 Table S1 |
| `IP3R_i1` | IP3R inhibited-1 (Ca²⁺ at inhibitory site) | `[dts]` | 167 | Dolan 2014 Table S1 |
| `IP3R_i2` | IP3R inhibited-2 | `[dts]` | 25 | Dolan 2014 Table S1 |
| `IP3R_s` | IP3R shut | `[dts]` | 1 | Dolan 2014 Table S1 |
| `SERCA_E1` | SERCA empty, E1 (cytosol-facing) | `[dts]` | **2 963** | Dolan 2014 Table S1 was 5 920; pre-equilibrated for binding step (§6.8 deviation D5) |
| `SERCA_E2` | SERCA empty, E2 (DTS-facing) | `[dts]` | 5 927 | Dolan 2014 Table S1 |
| `SERCA_E1Ca` | SERCA·Ca²⁺ in E1 | `[dts]` | **2 963** | Dolan 2014 Table S1 was 6; pre-equilibrated (§6.8 D5) |
| `SERCA_E1PCa` | SERCA E1 phosphorylated·Ca²⁺ | `[dts]` | 7 | Dolan 2014 Table S1 |
| `SERCA_E2PCa` | SERCA E2P·Ca²⁺ | `[dts]` | 4 | Dolan 2014 Table S1 |
| `SERCA_E2P` | SERCA E2P, empty | `[dts]` | 28 | Dolan 2014 Table S1 |
| `PMCA` | PMCA free (basal path) | `[pl]` | 765 | Dolan 2014 Table S1 |
| `PMCA_Ca` | PMCA·Ca²⁺ (basal) | `[pl]` | 4 | Dolan 2014 Table S1 |
| `Ca4_CaM_PMCA` | Ca₄·CaM·PMCA (CaM-activated, empty) | `[pl]` | 0 | **Phase 1 add** — Caride 2007 Table 3 step 8 product |
| `Ca4_CaM_PMCA_Ca` | Ca₄·CaM·PMCA·Ca²⁺ | `[pl]` | 0 | **Phase 1 add** — Caride 2007 step 9 |
| `PMCA_CaM` | PMCA·CaM (deactivating; bookkeeping only — step 11 not integrated, see §3.5) | `[pl]` | 0 | **Phase 1 add** — Caride 2007 step 11 product |
| `CaM_free` | Free calmodulin (no Ca²⁺) | `[c]` | 20 062 | **Phase 1 add** — equilibrated at 100 nM cyt; Dolan total 20 481 (§6.8 D2) |
| `Ca2_CaM` | Ca₂·CaM (N-lobe loaded) | `[c]` | 200 | **Phase 1 add** — Caride 2007 step 6 product (§6.8 D2) |
| `Ca4_CaM` | Ca₄·CaM (fully loaded; activates PMCA) | `[c]` | 219 | **Phase 1 add** — Caride 2007 step 7 product (§6.8 D2) |
| `STIM1_free` | STIM1 free monomer | `[dts]` | 438 | Dolan 2014 Table S1 |
| `STIM1_Ca` | STIM1 DTS-bound (inactive) | `[dts]` | 3 805 | Dolan 2014 Table S1 |
| `STIM1_dim` | STIM1 dimer (active sensor) | `[dts]` | 22 | Dolan 2014 Table S1 |
| `ORAI1` | Orai1 monomer (4/channel; opening via MWC, §3.6) | `[pl]` | 1 447 | Dolan 2014 Table S1 |

Conversion: V_cyt = 6 fL, V_DTS = 0.258 fL, N_A = 6.022 × 10²³ mol⁻¹.

### 3.2 IP3 sourcing

In a whole-cell model, IP3 should be a state variable produced by the GPCR
cascade (P2Y1 → Gq → PLCβ → IP3) and consumed by IP3 phosphatase. In v0.2,
it is instead driven by a **pre-programmed time curve** — a mathematical formula
that specifies exactly how IP3 concentration rises and falls after stimulation.
This is a known simplification borrowed directly from the Dolan 2014 model
(their Fig. S2). It lets us validate the Ca²⁺ core without first building the
full receptor cascade.

This means **you cannot model agonist dose-response until v0.3**. The
time curve shape is taken directly from Dolan 2014 Fig S2:

```
IP3(t) = IP3_rest × [1 + (fold-1) × (1 - exp(-t/tau_rise))
                          × exp(-max(0, t-t_peak)/tau_decay)]

fold     = 5.5   (peak amplitude relative to rest)
tau_rise   = 3.0 s
t_peak   = 3.0 s
tau_decay  = 60.0 s
```

At each timestep the IP3 count is *set* from the curve, not integrated, and
made available to the ODE solver as a boundary condition on IP3 concentration.
IP3 is still declared as a `BulkMolecule` so that downstream consumers of the
state see a valid count, but in v0.2 mass is not conserved on IP3 — the curve
creates and destroys it implicitly. This is an intentional v0.2 simplification
and closes in v0.3 when the upstream P2Y1 process produces IP3 and IP3
phosphatase consumes it as normal bulk reactions.

**v0.3 upgrade path:** Replace the time curve with a proper upstream
process (`P2Y1Signalling`) that produces IP3 as a BulkMolecule. The Ca²⁺
process then reads IP3 as a normal state variable. No change to the Ca²⁺
ODE system is needed — IP3 concentration simply becomes time-varying input.

> **Code:** the time curve is `ip3_forcing_uM(t)` at
> `reconstruction/platelet/dataclasses/process/calcium_signalling.py:291`;
> the constants `IP3_FOLD`, `IP3_T_PEAK`, `IP3_TAU_RISE`, `IP3_TAU_DECAY`,
> `IP3_REST_UM` are defined at lines 93–101 in the same file. The forcing
> mode is gated by `CalciumDynamics._ip3_forced`
> (`models/platelet/processes/calcium_dynamics.py:39`); inside the ODE
> the IP3 count is driven toward the curve at
> `calcium_signalling.py:665`.

### 3.3 IP3R: 6-state Markov model (Sneyd & Dufour 2002)

The IP3 receptor transitions between six states. Ca²⁺-dependent activation
and inhibition produce the biphasic open probability required for oscillations.
A Hill function cannot reproduce this behaviour.

![IP3R 6-state Markov model](/Users/steve/github/platelet-wcm/reports/figures/ip3r-state-machine.png)
*IP3R 6-state Markov model, adapted from Sneyd & Dufour (2002) type-2 kinetics
as parameterised in Purvis & Bhatt (2008) Table 1. Green states (o, a) are
Ca²⁺-conducting. Blue (n) is the resting neutral state. Purple states (i1, i2)
are Ca²⁺-inhibited. Red (s) is the shut state. Transitions are Ca²⁺- and IP3-
dependent.*

Rate constants are the Sneyd & Dufour 2002 type-2 fit as parameterised in
Purvis 2008 Table 1 — every value was re-verified against the Purvis PDF on
2026-04-23 (provenance pass) and is recorded in
`calcium-data-provenance.md` § "IP3R dynamics". Two corrections logged
during the verification pass: the missing `L₃ = 0.025 µM` equilibrium
constant on the n→i1 φ-function, and a sign convention check on the i1↔s
edge.

| Transition | Forward | Reverse | Source |
|------------|---------|---------|--------|
| n ↔ o (IP3 binding) | k₂ = 37.4 µM⁻¹s⁻¹, l₄ = 1.7 µM⁻¹s⁻¹ | k₋₂ = 1.4 s⁻¹, l₋₄ = 2.5 µM⁻¹s⁻¹ | Purvis 2008 T1 |
| n ↔ i1 (Ca²⁺ inhibition) | k₁ = 0.64 µM⁻¹s⁻¹, L₁ = 0.12 µM, L₃ = 0.025 µM | k₋₁ = 0.04 s⁻¹ | Purvis 2008 T1 |
| o ↔ a (Ca²⁺ activation) | k₄ = 4 µM⁻¹s⁻¹ | k₋₄ = 0.54 µM⁻¹s⁻¹ | Purvis 2008 T1 |
| o ↔ s (shutting) | k₃ = 11 µM⁻¹s⁻¹, L₅ = 54.7 µM | k₋₃ = 29.8 s⁻¹ | Purvis 2008 T1 |
| i1, i2 ↔ s (φ-function) | l₂ = 1.7 s⁻¹, l₆ = 4 707 s⁻¹ | l₋₂ = 0.8 s⁻¹, l₋₆ = 11.4 s⁻¹ | Purvis 2008 T1 |

> **Code:** rate constants in `K_IP3R`
> (`reconstruction/platelet/dataclasses/process/calcium_signalling.py:126`).
> Sneyd & Dufour φ-function rate laws live in `_phi_n_i1_fwd` /
> `_phi_n_o_fwd` / `_phi_o_a_fwd` / `_phi_a_i2_fwd` / `_phi_o_s_fwd` and
> their reverse counterparts (`calcium_signalling.py:309–367`). The
> 6-state Markov ODE block is `calcium_signalling.py:489–501`.

Open probability — fourth-power tetramer cooperativity (all four IP3R
subunits must be in conducting conformation):

```
P_o = ((0.9 . IP3R_a + 0.1 . IP3R_o) / IP3R_total)^4
```

Ca²⁺ flux through IP3R (Purvis 2008 eq. 13 / Dolan 2014 eq. 4):

```
J_IP3R = gamma_IP3R × N_IP3R × P_o × (V_IM - E_Ca,IM) × N_A/(zF)

gamma_IP3R   = 10 pS = 10 × 10^-^1^2 A/V    Zschauer 1988 single-channel
                                     (cited via Purvis 2008 Table 1, row
                                     "Ca^2^+ release from DTS")
V_IM     = -60 mV                    DTS membrane potential (Dolan 2014
                                     cluster analysis; sec.5, sec.6.8 D6)
E_Ca,IM  = (RT/zF) × ln([Ca^2^+]_dts / [Ca^2^+]_cyt)   Nernst potential, z=2
```

The driving force `(V_IM − E_Ca,IM)`: at resting [Ca²⁺]_dts/[Ca²⁺]_cyt = 2 500,
E_Ca,IM ≈ +104 mV, so with V_IM = −60 mV the driving force is ≈ −164 mV — a
strong inward (DTS → cyt) gradient when the channel opens.

> **Code:** Po⁴ tetramer cooperativity at
> `calcium_signalling.py:507–512`; Nernst flux (with the Ca²⁺-empty
> guard) at `calcium_signalling.py:524–539`. Conductance constant
> `GAMMA_IP3R_S` and membrane potential `V_IM_V` at `calcium_signalling.py:149`
> and `:117` respectively. Pre-factor `NA_OVER_zF` and Nernst coefficient
> `RT_OVER_zF_V` at `:111` / `:110`.

Implementation note (Phase 2 finding): with γ_IP3R = 10 pS and the initial conditions giving
DTS = 250 µM, the resting IP3R inflow (~112 k ions/s) exceeds SERCA cycle
throughput (~6 k ions/s after the binding step equilibrates), so the model's
natural fixed-point DTS sits well below 250 µM. γ recalibration alone does
not recover the Dolan initial conditions — see §6.8 D7 and lab-book 2026-05-05 for the full
diagnosis.

### 3.4 SERCA: E1–E2 cycle

Six-state enzymatic cycle for the SERCA3b isoform expressed in platelets
(Purvis 2008 Table 1; ref. Dode 2002 for isoform-specific kinetics):

```
E2 <→ E1 <→ E1.Ca^2^+ → E1P.Ca^2^+ <→ E2P.Ca^2^+ → E2P → E2
          ^ (cytosol)                      v (DTS)
```

Each cycle transports 2 Ca²⁺ ions from cytosol to DTS at the cost of 1 ATP.

| Step | Transition | Forward | Reverse | Source |
|------|------------|---------|---------|--------|
| 1 | E2 → E1 (shuttle) | k_shuttle_f = 600 s⁻¹ | k_shuttle_r = 600 s⁻¹ | Purvis 2008 T1 (Dode 2002) |
| 2 | E1 + 2 Ca²⁺_cyt → E1·Ca | k_bind_f = 1 000 µM⁻²s⁻¹ | k_bind_r = 10 s⁻¹ | Purvis 2008 T1 |
| 3 | E1·Ca → E1P·Ca | k_phos_f = 700 s⁻¹ | k_phos_r = 5 s⁻¹ | Purvis 2008 T1 |
| 4 | E1P·Ca ⇌ E2P·Ca (conformational) | k_conf_f = 600 s⁻¹ | k_conf_r = 50 s⁻¹ | Purvis 2008 T1 |
| 5 | E2P·Ca → E2P + 2 Ca²⁺_dts | k_release_f = 1 000 s⁻¹ | k_release_r = 4 × 10⁻³ µM⁻²s⁻¹ | Purvis 2008 T1 (corrected for "+ 2 Ca²⁺_dts" stoichiometry — see provenance log) |
| 6 | E2P → E2 (dephosphorylation) | k_dephos_f = 500 s⁻¹ | k_dephos_r = 1 s⁻¹ | Purvis 2008 T1 |

> **Code:** rate constants in `K_SERCA`
> (`reconstruction/platelet/dataclasses/process/calcium_signalling.py:156`).
> The 6-step cycle ODE block (mass-action `v_shuttle`, `v_bind`, `v_phos`,
> `v_conf`, `v_release`, `v_dephos` plus the cyt/dts Ca²⁺ accounting) is
> at `calcium_signalling.py:541–560`. ATP cost integration is in
> `CalciumSignalling.molecules_to_next_time_step`
> (`calcium_signalling.py:717–730`).

Implementation notes:

- `k_bind_f = 1×10¹⁵ M⁻²s⁻¹` is the Purvis primary-source value. An earlier
  draft of the model reduced it ~470× to compensate for an unrelated IP3R
  Po error; that calibration was reverted when the Po⁴ tetramer formula and
  Nernst flux were restored. The high `k_bind_f` is fine when E1 / E1·Ca
  start at binding equilibrium (§6.8 D5).
- `k_release_r × DTS² × E2P` runs in *reverse* at full DTS = 250 µM
  (~14 k ions/s into E2P·Ca vs ~8 k ions/s forward release). This is the
  proximate cause of the open Phase 2 DTS-resting question — see §6.8 D7
  and §7.1.

### 3.5 PMCA — full 5-state CaM-coupled scheme (Phase 1)

The original v0.2 design used a 2-state Michaelis–Menten approximation
with Caride 2007 *basal* (CaM-free) constants. Phase 1 (commit `f3080c40`,
2026-04-30) replaced it with the full Caride 2007 Table 3 5-state scheme:
PMCA has both a basal extrusion path and a Ca₄·CaM-activated path that
runs ~5.5× faster. The full scheme also requires the CaM Ca²⁺-binding
ladder (§3.5.1 below). Together these closed the cytosolic Ca²⁺ runaway
diagnosed pre-Phase 1 (lab-book 2026-05-01).

#### Basal path (Caride 2007 steps 4–5)

```
PMCA + Ca^2^+_cyt <→ PMCA.Ca^2^+ → PMCA + Ca^2^+_ex
              k_on, k_off       k_cat
```

| Constant | Value | Caride 2007 Table 3 row | Source |
|----------|-------|-------------------------|--------|
| `k_on` (k₄)    | 10 µM⁻¹s⁻¹ | step 4 fwd | Caride 2007 T3 |
| `k_off` (k₄ᵣ)  | 50 s⁻¹     | step 4 rev | Caride 2007 T3 |
| `k_cat` (k₅)   | 5.5 s⁻¹    | step 5 (V_max basal) | Caride 2007 T3 |

Derived basal KM = (k_off + k_cat)/k_on = 5.55 µM.

> **Code:** basal-path constants in `K_PMCA`
> (`reconstruction/platelet/dataclasses/process/calcium_signalling.py:174`).
> Basal-path ODE block (`v_pmca_bind`, `v_pmca_cat`) at
> `calcium_signalling.py:574–581`.

#### CaM-activated path (Caride 2007 steps 8–10; step 11 omitted)

```
PMCA + Ca_4.CaM    <→ Ca_4.CaM.PMCA                   (step 8)
Ca_4.CaM.PMCA + Ca^2^+ <→ Ca_4.CaM.PMCA.Ca               (step 9)
Ca_4.CaM.PMCA.Ca → Ca_4.CaM.PMCA + Ca^2^+_ex            (step 10)
```

Step 11 (`Ca₄·CaM·PMCA → PMCA·CaM + 4 Ca²⁺`, slow CaM deactivation) is
**not integrated** in Phase 1. It operates on a τ ≈ 20 min timescale —
far longer than the 200 s transient — and including it caused PMCA to
accumulate in a dead-end `PMCA·CaM` state within ~30 s of activation
(early Phase 1 finding, lab-book 2026-05-01). The `PMCA_CaM[pl]` state
variable is retained at an initial count of 0 for mass-conservation bookkeeping and
listener output but is never written to by the ODE. The Caride `k11` /
`k11r` constants are defined in `K_CAM_PMCA` but currently unused; they
will be re-enabled if the v0.3 longer-timescale work needs them.

| Constant | Value | Step | Source |
|----------|-------|------|--------|
| `k8`  | 0.2 µM⁻¹s⁻¹  | 8 fwd  | Caride 2007 T3 |
| `k8r` | 8.0 × 10⁻⁴ s⁻¹ | 8 rev | Caride 2007 T3 |
| `k9`  | 50 µM⁻¹s⁻¹   | 9 fwd  | Caride 2007 T3 |
| `k9r` | 10 s⁻¹       | 9 rev  | Caride 2007 T3 |
| `k10` | 30 s⁻¹       | 10 (V_max CaM) | Caride 2007 T3 |
| `k11` *(unused)* | 10 s⁻¹ | 11 fwd (slow deact, τ ≈ 20 min) | Caride 2007 T3 |
| `k11r` *(unused)* | 7.332 × 10⁻⁴ µM⁻⁴s⁻¹ | 11 rev | Caride 2007 T3 |

Phase 1 finding: a step-10 product-recycling bug initially consumed PMCA
molecules each pump cycle; fixed by adding `+v_cam_pmca_cat` to
`dy[Ca4_CaM_PMCA]` (lab-book 2026-05-01). PMCA total is now conserved
across runs at 769.

> **Code:** CaM-activated rate constants in `K_CAM_PMCA`
> (`reconstruction/platelet/dataclasses/process/calcium_signalling.py:195`).
> Steps 8–10 ODE block (`v_cam_bind_pmca`, `v_cam_pmca_bind`,
> `v_cam_pmca_cat`) at `calcium_signalling.py:583–614`. Step 11 (slow
> deactivation) is presently omitted; see the comment at
> `calcium_signalling.py:583–591` for the rationale.

> **Provenance correction (2026-04-23).** The earlier draft of this section
> quoted `KM1 = 0.5 µM, KM2 = 1.0 µM, kcat = 8.9 s⁻¹` attributed to Caride
> 2007 / Purvis 2008 Table 1. Those numbers do not appear in either paper
> for PMCA — they correspond to Reaction #11 in Purvis 2008 (CDPDG synthesis,
> a phospholipid biosynthesis enzyme) and are unrelated to Ca²⁺ extrusion.
> Recorded in `calcium-data-provenance.md` § "Provenance correction".

#### 3.5.1 CaM Ca²⁺-binding ladder (Caride 2007 steps 6–7)

The CaM-activated PMCA path consumes Ca₄·CaM, which is produced by
sequential Ca²⁺ binding to free CaM. Two-lobe cooperative scheme; both
binding events transfer 2 Ca²⁺ at once (the slow N-lobe pair, then the
fast C-lobe pair, captured as a single µM⁻²·s⁻¹ rate per step).

```
CaM_free   + 2 Ca^2^+ <→ Ca_2.CaM      (step 6, slow)
Ca_2.CaM    + 2 Ca^2^+ <→ Ca_4.CaM      (step 7, fast)
```

| Constant | Value | Step | Source |
|----------|-------|------|--------|
| `k6`  | 2.669 µM⁻²s⁻¹ | 6 fwd | Caride 2007 T3 |
| `k6r` | 2.682 s⁻¹     | 6 rev | Caride 2007 T3 |
| `k7`  | 170.4 µM⁻²s⁻¹ | 7 fwd | Caride 2007 T3 |
| `k7r` | 1.551 s⁻¹     | 7 rev | Caride 2007 T3 |

Total CaM = 20 481 (Dolan 2014 Table S1). The Dolan initial conditions (CaM_free=20 465,
Ca₂·CaM=15, Ca₄·CaM=1) is *not* at equilibrium for our explicit Caride
ladder at 100 nM cyt — using it caused a ~34 k Ca²⁺ ion loading burst on
t=0. We override with the detailed-balance solution (CaM_free=20 062,
Ca₂·CaM=200, Ca₄·CaM=219); see §6.8 D2.

CaM that is bound to PMCA is held in the `Ca4_CaM_PMCA*` and `PMCA_CaM`
sub-states (§3.5 above), not double-counted in the free CaM ladder.

> **Code:** ladder rate constants in `K_CAM`
> (`reconstruction/platelet/dataclasses/process/calcium_signalling.py:183`).
> ODE block (`v_cam_bind1`, `v_cam_bind2` and the cyt-Ca buffering term)
> at `calcium_signalling.py:562–572`.

### 3.6 SOCE — STIM1 sensor cycle + Hoover/Dolan MWC + Dolan puncta entry

Phase 1 (commit `18ed7184`, 2026-04-29) replaced the early ad-hoc 3-state
mass-action SOCE model with the Hoover & Lewis 2011 MWC (Monod–Wyman–Changeux)
allosteric scheme as adopted by Dolan 2014 (issues #45/#46). Three coupled
pieces:

#### 1. STIM1 sensor cycle (mass-action, detailed-balance initial conditions)

```
STIM1.Ca^2^+_dts <→ STIM1_free + Ca^2^+_dts        Ca^2^+ release from EF-hand
2 STIM1_free   <→ STIM1_dim                    dimerisation = active sensor
```

| Constant | Value | Source |
|----------|-------|--------|
| `k_release_f` | 0.1 s⁻¹ | chosen so detailed balance at the Dolan initial conditions gives `k_release_r` |
| `k_release_r` | 3.475 × 10⁻³ µM⁻¹s⁻¹ | derived: `k_release_f × STIM1_Ca / (STIM1_free × ca_dts)` at initial conditions |
| `k_dim_f`     | 1.15 × 10⁻⁴ count⁻¹s⁻¹ | derived: `k_dim_r × STIM1_dim / STIM1_free²` at initial conditions |
| `k_dim_r`     | 1.0 s⁻¹ | choice; pairs with derived k_dim_f to land Dolan initial conditions at detailed balance |

The Dolan 2014 Table S1 initial conditions (`STIM1_Ca`=3 805, `STIM1_free`=438, `STIM1_dim`=22)
is held at exact detailed balance by this choice — STIM1 sub-states do not
drift at rest.

> **Code:** STIM1 cycle constants in `K_STIM`
> (`reconstruction/platelet/dataclasses/process/calcium_signalling.py:210`).
> Cycle ODE block (`v_stim1_release`, `v_dim`) at
> `calcium_signalling.py:616–627`.

#### 2. Hoover & Lewis 2011 MWC channel-opening model

The Orai channel is treated as an MWC tetramer that closes by default and
opens cooperatively as STIM2 (which we identify with `STIM1_dim` here)
binds. Each bound STIM2 stabilises the open state by a factor `f`, giving
the standard MWC form:

```
P_open = (L . (1 + a.Sf.Ka)^4) /
         (L . (1 + a.Sf.Ka)^4 + (1 + Sf.Ka)^4)    [closed-favoured]
       — equivalently rearranged with the f cooperativity factor --
       (Hoover 2011 Fig. 4 best-fit parameters)
```

| Constant | Value | Source / note |
|----------|-------|---------------|
| `L`  | 1.0 × 10⁻⁴ | Hoover 2011 Fig. 4 best-fit (intrinsic open/closed equilibrium without STIM) |
| `Ka` | 2.0        | **rescaled** from Hoover 100 a.u. — see §6.8 D3 |
| `f`  | 14.2       | Hoover 2011 Fig. 4 best-fit (per-STIM2 opening cooperativity) |
| `a`  | 0.5        | Hoover 2011 Fig. 4 best-fit (binding cooperativity, < 1 = negative) |

> **Code:** MWC parameters in `K_MWC`
> (`reconstruction/platelet/dataclasses/process/calcium_signalling.py:239`).
> The MWC equilibrium solver `_mwc_open_fraction(stim2_p, n_orai)` —
> Newton-style iteration over the Sf mass balance using cumulative
> cooperativity factors `a^(i(i-1)/2)` — is at `calcium_signalling.py:376–437`.

#### 3. Dolan 2014 puncta entry (eq. 2)

Not all STIM2 dimers are competent to engage Orai — only the fraction that
has translocated into ER–PM puncta near Orai clusters. Dolan eq. 2 makes
this a Hill function of cytosolic Ca²⁺:

```
qp = α . [Ca^2^+]_cyt^n / (KM^n + [Ca^2^+]_cyt^n) + baseline
Sf = qp . STIM1_dim     # effective STIM2 ligand for the MWC
```

| Constant | Value | Source |
|----------|-------|--------|
| `α`        | 0.2   | Dolan 2014 default |
| `KM_uM`    | 0.5   | Dolan-scanned mid-range; one of two homeostatically-constrained free params |
| `n`        | 4.0   | Dolan-scanned mid-range; the other free param |
| `baseline` | 0.01  | constitutive puncta fraction at zero [Ca²⁺]_cyt |

> **Code:** puncta-entry constants in `PUNCTA`
> (`reconstruction/platelet/dataclasses/process/calcium_signalling.py:250`).
> Hill / qp / Sf evaluation inside the ODE at
> `calcium_signalling.py:629–645`.

#### 4. SOCE current (Dolan 2014 eq. 4)

```
I_SOC = gamma_SOC × N_orai_channels × P_open × (V_PM - E_Ca,PM) × N_A/(zF)

gamma_SOC = 0.3 fS                      effective conductance, calibrated
                                    (literature 24 fS, see sec.6.8 D3)
V_PM   = -60 mV                     plasma membrane potential
E_Ca,PM = (RT/zF) × ln([Ca^2^+]_ex / [Ca^2^+]_cyt)    Nernst potential, z=2
N_orai_channels = ORAI1_count / 4   (4 monomers/tetramer)
```

The `γ_SOC = 0.3 fS` is calibrated against the rest-balance condition
`J_SOCE_rest ≈ J_PMCA_steady_rest ≈ 76 ions/s` at the MWC P_open value
the rescaled `Ka` produces (~1.2 × 10⁻³ at basal STIM1_dim = 22). The
Hoover face-value 24 fS would produce spurious µM/s leaks at our integer
counts — see §6.8 D3.

> **Code:** `GAMMA_SOC_S` and the `ORAI_SUBUNITS_PER_CHANNEL` /
> `STIM_MONOMERS_PER_DIMER` stoichiometry constants at
> `reconstruction/platelet/dataclasses/process/calcium_signalling.py:266`,
> `:283`, `:288`. SOCE current ODE block (Nernst driving force +
> ions/s applied to cyt) at `calcium_signalling.py:647–657`. The
> instantaneous SOCE flux is also recomputed for the listener trace at
> `models/platelet/listeners/calcium_trace.py:120–141`.

---

## 4. Architecture

### 4.1 Signal flow

The diagram below shows how the wcEcoli process/state architecture has been
adapted for platelet Ca²⁺ signalling. The E. coli biological processes
(transcription, metabolism, etc.) are replaced by the platelet-specific
`CalciumDynamics` process, but the underlying framework — state containers,
time-stepping, ATP/ADP partitioning — is carried over unchanged.

![Ca²⁺ signalling process architecture](/Users/steve/github/platelet-wcm/reports/figures/calcium-process-architecture.png)
*ATP coupling shown in red. Green compartment = cytosol, blue = DTS, purple = plasma membrane.*

### 4.2 File contents

#### `reconstruction/platelet/dataclasses/process/calcium_signalling.py`

**Class `CalciumSignalling`** (`calcium_signalling.py:673`) — builds the ODE
system at startup (called once, not each timestep). Stores the species list
and all rate constants. In v0.2 these are hardcoded; v0.3 migration swaps in
TSV parsing without touching the process code.

Key method: **`molecules_to_next_time_step(counts, dt, t_sim, ip3_forced)`**
(`calcium_signalling.py:689`) — takes the current 27-element integer count
vector, runs the Ca²⁺ ODE sub-stepper (BDF solver, one timestep duration)
with module-level volumes and Avogadro, applies the Dolan IP3 time curve
when `ip3_forced=True`, and returns integer count deltas plus the estimated
ATP cost. The ODE right-hand side is `_ode_rhs(t, y, t_sim_start, ip3_forced)`
at `calcium_signalling.py:440`.

#### `models/platelet/processes/calcium_dynamics.py`

**Class `CalciumDynamics(Process)`** (`calcium_dynamics.py:33`) — thin
wrapper around the ODE solver. Follows the standard wcEcoli process
lifecycle:

- **`initialize(sim, sim_data)`** (`calcium_dynamics.py:44`) — creates
  `BulkMoleculesView` over all 27 ODE state species
  (`self._solver.molecule_names`) plus a separate `BulkMoleculesView` over
  `ATP[c]` for the per-cycle ATP debit.
- **`calculateRequest()`** (`calcium_dynamics.py:55`) — calls
  `molecules_to_next_time_step` to estimate changes for this timestep, then
  requests the required molecules (Ca²⁺ species needed by the ODE) and the
  estimated ATP from the framework's partitioner.
- **`evolveState()`** (`calcium_dynamics.py:68`) — applies the Ca²⁺ count
  changes and performs ATP/ADP/Pi accounting using the ATP actually
  allocated (not the requested amount, which may have been reduced by the
  partitioner).

The class attribute `_ip3_forced` (default `True`,
`calcium_dynamics.py:39`) controls whether the Dolan Fig. S2 IP3 time
curve is applied or whether IP3 evolves via the ODE alone. The
runscript exposes this as the `--no-ip3-forcing` CLI flag (and the
webapp Configure tab as a checkbox); both override the class attribute
before `PlateletSimulation` is constructed. v0.3 will set this to
`False` permanently once the upstream P2Y1 process produces IP3
endogenously.

Similarly, the module-level `cs_mod.CA_EX_UM` constant in
`calcium_signalling.py` (default 1.2 mM × 1000 = 1200 µM) is
overridden by `runPlateletSim.py`'s `--ca-ex-mM` flag — set to 0 for
the Dolan Fig. 4 EDTA condition. Both the SOCE current and the basal
PM Ca²⁺ leak are gated on `CA_EX_UM > 0` (`_ode_rhs` line ~652);
under EDTA both are physically zero.

#### `models/platelet/listeners/calcium_trace.py`

**Class `CalciumTrace`** (`calcium_trace.py:45`). The per-timestep
recording is in `update()` (`calcium_trace.py:96`); the table schema and
attribute hooks are at `tableCreate()` / `tableAppend()`
(`calcium_trace.py:143`, `:149`).

Records 14 columns per timestep (driven by Phase 1's expanded state):

| Column | Description |
|--------|-------------|
| `ca_cyt_nM` | `CA2_CYT[c]` converted to nM |
| `ca_dts_uM` | `CA2_DTS[dts]` converted to µM |
| `ip3_nM`    | `IP3[c]` converted to nM |
| `cam_free` / `ca2_cam` / `ca4_cam` | CaM ladder sub-state counts |
| `pmca_free` / `pmca_ca` | PMCA basal sub-state counts |
| `ca4_cam_pmca` / `ca4_cam_pmca_ca` / `pmca_cam` | PMCA·CaM complex counts |
| `stim1_dim` | STIM1 dimer count (active sensor) |
| `soce_flux_nMs` | estimated Orai1 inflow (nM/s; sign convention: into cyt) |

---

## 5. Parameter index

Every numeric parameter in the implemented model, with source and any
deviation from the cited primary value. Detailed equations and rate-law
context are in §3; this section is the single-place lookup. Full
literature provenance — including Zotero item keys and direct PDF
quotations — is in `calcium-data-provenance.md`.

### 5.1 Compartment volumes

| Compartment | Value | Source |
|-------------|-------|--------|
| Cytosol (V_cyt) | 6.0 fL | Purvis 2008 direct measurement |
| DTS (V_dts) | 0.258 fL (4.3% × cyt) | Purvis 2008 glucose-6-phosphatase stain |
| Extracellular | infinite, fixed reservoir | Dolan 2014 |

### 5.2 Resting concentrations and total copy numbers

| Species | Concentration | Count (calc) | Source |
|---------|---------------|--------------|--------|
| Ca²⁺_cyt | 100 nM | 361 | Purvis 2008 |
| Ca²⁺_dts | 250 µM | 38 842 | Dolan 2014 Fluo-5N measurement |
| Ca²⁺_ex (`CA_EX_UM`) | 1.2 mM | fixed reservoir | Dolan 2014 |
| IP3 | 50 nM | 181 | Sveshnikova 2025 / Dolan 2014 (mid-range) |
| IP3R total (sum of 6 sub-states) | — | 1 328 | Dolan 2014 Table S1 |
| SERCA total (sum of 6 sub-states) | — | 11 892 | Dolan 2014 Table S1 |
| PMCA total (sum of 5 sub-states) | — | 769 | Dolan 2014 Table S1 |
| STIM1 total (sum of 3 sub-states) | — | 4 265 | Dolan 2014 Table S1 |
| Orai1 (monomers; 4/channel) | — | 1 447 | Dolan 2014 Table S1 |
| CaM total (3 free + bound to PMCA) | — | 20 481 | Dolan 2014 Table S1 |

Per-sub-state initial counts are in §3.1.

### 5.3 Physical constants

| Constant | Value | Used in |
|----------|-------|---------|
| N_A (Avogadro) | 6.022 × 10²³ mol⁻¹ | volume ↔ count conversions |
| F (Faraday) | 96 485 C/mol | Nernst potentials |
| R (gas constant) | 8.314 J/(mol·K) | Nernst potentials |
| T | 310 K (37 °C) | Nernst potentials (Purvis / Dolan) |
| RT / (zF) | 0.01334 V (z=2) | Nernst pre-factor |
| N_A / (zF) | 3.121 × 10¹⁸ ions/(A·s) (z=2) | flux ↔ ion-rate conversions |
| V_IM | −60 mV | DTS membrane potential — IP3R driving force (§3.3); Dolan cluster-analysis upper bound, see §6.8 D6 |
| V_PM | −60 mV | plasma membrane potential — SOCE current (§3.6) |

### 5.4 IP3R kinetics (Sneyd & Dufour 2002 type-2; Purvis 2008 T1)

See §3.3 for the full transition table. Conductance:

| Constant | Value | Source / note |
|----------|-------|---------------|
| `γ_IP3R` | 10 pS = 10⁻¹¹ A/V | Zschauer 1988 (cited via Purvis 2008 T1). Phase 2 calibration question — see §6.8 D7 |
| `N_IP3R` (total channels in Po formula) | sum of 6 sub-states (1 328 at initial conditions) | Dolan 2014 Table S1 |
| Po form | `((0.9·a + 0.1·o)/total)⁴` | Purvis 2008 T1 footnote (4-subunit cooperativity) |

### 5.5 SERCA cycle (SERCA3b; Purvis 2008 T1, Dode 2002 isoform kinetics)

See §3.4 for the full step table. ATP coupling: 1 ATP / cycle / 2 Ca²⁺
transported (Dode 2002).

### 5.6 PMCA — basal and CaM-activated (Caride 2007 T3 5-state)

See §3.5 for the full step table.

### 5.7 CaM Ca²⁺ binding (Caride 2007 T3 steps 6–7)

See §3.5.1.

### 5.8 SOCE — STIM1 cycle, MWC, puncta entry, Orai conductance

See §3.6 for the full constant tables. Headline values:

| Constant | Value | Source / note |
|----------|-------|---------------|
| `γ_SOC` (Orai effective single-channel) | 0.3 fS | calibrated against rest-balance; literature face-value 24 fS (Prakriya & Lewis 2002), see §6.8 D3 |
| MWC `L` | 10⁻⁴ | Hoover 2011 Fig. 4 |
| MWC `Ka` | 2.0 (rescaled from 100 a.u.) | Hoover 2011 Fig. 4; rescaling rationale §6.8 D3 |
| MWC `f` | 14.2 | Hoover 2011 Fig. 4 |
| MWC `a` | 0.5 | Hoover 2011 Fig. 4 |
| Puncta `α` | 0.2 | Dolan 2014 default |
| Puncta `KM` | 0.5 µM | Dolan 2014 free param (mid-range scan) |
| Puncta `n` | 4 | Dolan 2014 free param (mid-range scan) |
| Puncta `baseline` | 0.01 | Dolan 2014 |

### 5.9 IP3 forcing curve (Dolan 2014 Fig. S2 fit)

```
IP3(t) = IP3_rest × (1 + (fold - 1)
         × (1 - e^(-t/tau_rise)) × e^(-max(0, t-t_peak)/tau_decay))
```

| Constant | Value | Source |
|----------|-------|--------|
| `IP3_REST_UM` | 0.05 µM (50 nM) | Sveshnikova 2025 / Dolan 2014 baseline |
| `IP3_FOLD`    | 5.5 | Dolan 2014 Fig. S2 |
| `IP3_T_PEAK`  | 3.0 s | Dolan 2014 Fig. S2 |
| `IP3_TAU_RISE` | 3.0 s | Dolan 2014 Fig. S2 |
| `IP3_TAU_DECAY` | 60.0 s | Dolan 2014 Fig. S2 |

### 5.10 Calibrated / model-specific constants (not in primary sources)

| Constant | Value | Reason | Reference |
|----------|-------|--------|-----------|
| `J_PM_LEAK_IONS_S` | 75 ions/s | Phase 2a addition — without a basal PM Ca²⁺ leak the ODE has no PM-side cyt source large enough to balance PMCA; calibrated against the rest condition `J_SOCE + J_leak = J_PMCA_steady` at cyt = 100 nM (defined at `reconstruction/platelet/dataclasses/process/calcium_signalling.py:280`; applied as a constant `dy[CA2_CYT]` term at `:661`) | §6.8 D4 / lab-book 2026-05-05 |

---

## 6. Implementation decisions

### 6.1 ODE state vs BulkMolecules

**Decision:** ODE solver works in concentration (µM) internally. At the start
of each timestep, integer counts are converted to concentration, integrated,
then rounded back to counts. This matches the pattern in the existing
`TwoComponentSystem` process (wcEcoli's reference implementation of a
two-component bacterial signalling pathway), which demonstrates the same
count ↔ concentration ↔ count conversion we replicate here.

**Implication:** At resting cytosolic Ca²⁺ (~361 molecules), rounding
introduces ~0.3% quantisation noise per step — a ±1 molecule round-off against
a count of 361. Crucially, this error does **not** accumulate: each timestep
starts from the previous integer count, runs the ODE in concentration space,
and rounds back to integers, so round-off is bounded at ±1 per species per
step rather than drifting over the full 300-step simulation. Sveshnikova 2025
notes that stochastic effects at this scale are biologically real. We accept
the round-off and flag it in the analysis.

### 6.2 One process or many?

**One process.** The Ca²⁺ subsystem is tightly coupled — splitting IP3R,
SERCA, and SOCE across processes would require artificial partitioning at
every interface. The ODE captures this coupling naturally. This mirrors the
`TwoComponentSystem` design decision (see §6.1).

### 6.3 IP3R model: Hill vs 6-state?

**6-state from the start.** The Sneyd & Dufour model is the consensus
implementation in Purvis 2008 and Dolan 2014. Ca²⁺ oscillations — a key
validation target — require the biphasic Ca²⁺-dependent gating that only the
6-state model provides. Issue #43 (previously flagged as an upgrade) is
therefore merged into #25 and implemented in the initial version.

### 6.4 Extracellular Ca²⁺

Treated as a fixed reservoir (1.2 mM; Dolan 2014). Not stored in
BulkMolecules — just a constant in the SOCE and PMCA rate equations.
This is the standard simplification in all reference models.

### 6.5 ATP coupling

**SERCA and PMCA are ATPases. ATP must be accounted for.**

The framework explicitly tracks ATP/ADP via `BulkMolecules`. The Ca²⁺ process
must request ATP in `calculateRequest()` and return ADP + Pi in `evolveState()`.
Ignoring this breaks the whole-cell energy budget.

Stoichiometry:
- **SERCA**: 1 ATP → 2 Ca²⁺ transported (cytosol → DTS)
- **PMCA**: 1 ATP → 1 Ca²⁺ extruded (cytosol → extracellular)

At each timestep, the ODE integration gives the net Ca²⁺ transport. The
`estimate_atp_cost` method in `CalciumSignalling` calculates:
`ATP = floor(Δca_serca / 2) + Δca_pmca`. `calculateRequest()` requests this
amount from the framework's ATP pool. `evolveState()` then uses the ATP
actually allocated (which may be less if ATP is scarce) and returns the
corresponding ADP and Pi.

At rest, SERCA and passive leak are in balance — net ATP consumption is small.
During a Ca²⁺ transient, SERCA works hard to refill the DTS, and the ATP drain
will be visible in the mass listener as a metabolic signature of activation.
This is a useful emergent property: it means platelet activation shows up
correctly as an energy cost without any special-casing.

In v0.2, ATP, ADP, and Pi are all in the molecule inventory and the process
requests ATP each step. If the ATP pool is insufficient (unlikely while
nothing else consumes it yet), the framework's partitioner will allocate
less than requested. The process scales pump activity proportionally: it
advances SERCA and PMCA only to the extent that allocated ATP supports,
and the unsupported pump turnover does not happen. This keeps mass and
energy balance correct even under starvation. A warning is logged on any
shortfall so the failure mode is visible; a proper metabolic process in a
future milestone will make shortfalls biologically meaningful.

### 6.6 ODE solver

`scipy.integrate.solve_ivp` with method `'BDF'` (backward differentiation,
suitable for stiff systems). Tolerances: `atol=1e-10`, `rtol=1e-8`.
Ca²⁺ transients resolve on 1–3 s timescale; 1-second timestep is adequate.

### 6.7 Future extensibility

The hardcoded parameter approach (Milestone 1) is identical to the TSV-based
approach from the Process's point of view. The only change for Milestone 2 is
swapping the two assignment lines in the dataclass `__init__`. See
`calcium-signalling-pathway-design.md` §"Milestone 1 shortcut" for the exact
migration pattern.

### 6.8 Deviations from primary sources

The implementation departs from Dolan 2014 / Purvis 2008 / Caride 2007 in
nine places. Each deviation is either a Phase 1/Phase 2a implementation
update, an integer-count realism adjustment, or a known open question.
Numbering (D1–D9) is referenced throughout the rest of the document.

**D1 — PMCA model: basal-only → full 5-state CaM-coupled.** The original v0.2
design used a 2-state Caride 2007 *basal* (CaM-free) PMCA with k_cat = 5.5 s⁻¹
and accepted a slower decay phase as a v0.3 follow-up. Phase 1 (commit
`f3080c40`, lab-book 2026-05-01) implemented the full Caride 2007 Table 3
5-state scheme (basal steps 4–5 plus CaM-activated steps 8–11) plus the CaM
ladder needed to feed it. The model now has both paths; basal still operates
at low cyt, the CaM-activated path (~5.5× faster k_cat) clamps the peak.
**Reason for upgrade:** without CaM-coupled PMCA the cytosolic Ca²⁺ ran away
during the transient because basal PMCA could not extrude fast enough to
match SOCE-driven inflow.

**D2 — CaM initial condition: Dolan Table S1 → detailed-balance equilibrium.**
Dolan Table S1 reports CaM split as (CaM_free=20 465, Ca₂·CaM=15, Ca₄·CaM=1)
with total 20 481. These ratios reflect Dolan's original model, which did
*not* track explicit CaM Ca²⁺-binding kinetics. Our Phase 1 implementation
adds the Caride 2007 step-6 / step-7 ladder, and at cyt = 100 nM the
detailed-balance solution is (20 062, 200, 219). Using Dolan's verbatim
values would cause a ~34 k Ca²⁺ ion loading burst on t=0 as the ladder
equilibrated. **We override with the equilibrium values; CaM total is
preserved at 20 481.**

**D3 — γ_SOC: Prakriya & Lewis 24 fS → calibrated 0.3 fS; MWC Ka rescaled.**
The CRAC-channel single-channel literature value is ~24 fS (Prakriya & Lewis
2002, Vig 2006), measured at saturating Po with patch-clamp in HEK cells.
Hoover & Lewis 2011 fit `Ka = 100` in a.u. where saturating STIM expression
`Stotal = 3.2` a.u. We rescale `Ka` to platelet dimer counts so that
`Ka_platelet × Sf_saturating ≈ Hoover Ka × Stotal = 320`, giving
`Ka_platelet = 2.0` (the MWC shape is insensitive to ~2× perturbations near
the saturating end). Combined with the ~10⁻³ MWC Po at basal STIM1_dim = 22,
applying 24 fS at face value would produce spurious µM/s leaks at our integer
counts. We therefore calibrate `γ_SOC = 0.3 fS` against the rest-balance
condition `J_SOCE_rest ≈ J_PMCA_steady_rest ≈ 76 ions/s`. **f, a, and L
transfer directly from Hoover 2011.**

**D4 — Basal plasma-membrane Ca²⁺ leak: not in Dolan / Purvis → 75 ions/s.**
The Dolan 2014 ODE has no PM leak term; cyt Ca²⁺ entry happens only through
SOCE (Orai1). The lab-book 2026-05-05 Phase 2a analysis showed that without a
constant PM leak, the model's resting cyt cannot sit at 100 nM regardless
of γ_IP3R: at full DTS the SOCE flux is too small (~6 ions/s) to balance
PMCA basal outflow (~77 ions/s at PMCA quasi-equilibrium). We add a
constant `J_PM_LEAK_IONS_S = 75 ions/s` term, calibrated against the rest
balance `J_SOCE + J_leak = J_PMCA`. **Biological motivation:** unidentified
TRPC / NCX-reverse / residual constitutive permeability (Sage & Rink
1985–1990; Brandman & Liou 2010 review). **Open question:** Phase 2a sweeps
showed (ii)+(iii) at γ=10 pS still doesn't recover Dolan rest because of D7;
the leak magnitude is correct for cyt-balance but cyt sits at ~3 nM after
the transient because SERCA siphons leak flux into DTS faster than the leak
can replenish.

**D5 — SERCA E1 / E1·Ca initial conditions: Dolan Table S1 (5 920, 6) → pre-equilibrated
(2 963, 2 963).** At cyt = 100 nM the SERCA binding equilibrium gives
`E1·Ca / E1 = k_bind_f · cyt² / k_bind_r = 1000 · 0.01 / 10 = 1.0`. Dolan's
verbatim Table S1 initial conditions are far from this — using them produced a ~118 k ions/s
spurious cyt → E1·Ca pulse on t=0 (lab-book 2026-05-05 Phase 2a (iii)). We
redistribute the (E1 + E1·Ca = 5 926) total ~equally; **SERCA total is
preserved at 11 892**. The other SERCA sub-states (E2, E1P·Ca, E2P·Ca, E2P)
are kept at their Dolan values pending §7.1 stability evaluation.

**D6 — DTS membrane potential V_IM: Dolan range −100..−60 mV → −60 mV.**
Dolan 2014 Methods §"Membrane potentials" reports a sampling range; their
SOCE-active configurations cluster at V_IM > −70 mV. We use −60 mV (upper
bound) as the design choice; if Phase 3 SOCE flux is mis-sized in either
direction this is the first parameter to revisit. The ~3.7× sensitivity of
IP3R driving force `(V_IM − E_Ca,IM)` to V_IM means small changes here have
substantial effects on both IP3R and SOCE inflows.

**D7 — γ_IP3R: Zschauer 1988 10 pS face-value (kept) — open Phase 2 question.**
The 10 pS value is a single-channel patch-clamp measurement (Zschauer 1988;
cited via Purvis 2008 Table 1). At γ = 10 pS with the Dolan initial conditions, the Phase 1
simulation produces a Ca²⁺ peak in the 200–800 nM acceptance band but
empties DTS in ~0.35 s. The 2026-05-01 lab book proposed `γ_IP3R = 0.6 fS`
calibrated against a SERCA cycle throughput of ~6 600 ions/s; the 2026-05-05
Phase 2a sweep showed γ alone cannot recover Dolan rest because the cycle's
natural fixed-point DTS sits below 250 µM regardless (the SERCA `k_release_r`
reverses at full DTS — see §3.4). **Current state: γ kept at 10 pS; the
Phase 1 transient peak passes; resting and post-transient DTS levels are an
open issue, tracked in issue #48.**

**D8 — IP3R copy number: Burkhart 2012 proteomics 4 850 → Dolan 2014 1 328.**
Burkhart proteomics counts every ITPR2 monomer present in the platelet
proteome. Dolan 2014 fits a smaller "functionally gated DTS-membrane
population" (1 328) consistent with their homeostatic constraints. We adopt
the Dolan value because the Sneyd & Dufour rate constants describe a
functionally-active gated channel, not all expressed protein. (See §9 Q4.)

**D9 — V_cyt: Sveshnikova 2025 3 fL alternative → Purvis 2008 6 fL.** Two
literature volumes for "cytosol" exist. Sveshnikova 2025 reports 3 fL
(probably a non-granular sub-fraction); Purvis 2008 measures 6 fL directly.
We use Purvis because (a) it is a direct geometric measurement, (b) it is
the basis for the Dolan 2014 parameter set we adopt downstream, and (c) it
gives the more conservative resting Ca²⁺ count of ~361 (§3.1). The
sensitivity to volume choice is documented in the analysis plots.

---

## 7. Validation strategy

The validation criteria below are split into "passing", "open" and "future"
to reflect actual implementation status as of 2026-05-06.

### 7.1 Resting state stability (no stimulus) — **OPEN (D7)**

Run 300 s with `ip3_forced=False` (process attribute). Original v0.2 pass
criteria:

- Ca²⁺_cyt stays within 80–120 nM
- Ca²⁺_dts stays within 200–300 µM
- All six IP3R state fractions stay within 10% of their Dolan Table S1
  resting values

**Current state (lab-book 2026-05-05 sweep):** at γ_IP3R = 10 pS the system
is in an unstable regime and DTS drains to 0 at rest. At γ_IP3R ≤ 2 × 10⁻¹³
the system is stable but settles at cyt ≈ 3–25 nM and DTS ≈ 400–450 µM —
neither matches Dolan's initial conditions. The cyt-side balance is correct (§6.8 D4); the
DTS side is bottlenecked by SERCA `k_release_r` reversing at full DTS (§3.4,
§6.8 D7). This is the live Phase 2 issue (#48).

### 7.2 Ca²⁺ transient shape (primary validation) — **PEAK PASSES**

Run 200 s with `ip3_forced=True` (Dolan 2014 Fig. S2 fold = 5.5, τ_rise = 3 s,
t_peak = 3 s, τ_decay = 60 s).

Original v0.2 pass criteria:

- Peak Ca²⁺_cyt: 200–500 nM, reached within 15–20 s
- Partial DTS depletion: Ca²⁺_dts drops to 30–70% of resting
- Sustained plateau above baseline (SOCE-dependent)
- Return towards baseline within 300 s

**Current state (Phase 2a smoke test):**

- Peak Ca²⁺_cyt: **299.5 nM at t = 1 s** ✓ (within 200–800 nM acceptance band
  used by the regression suite; lab-book 2026-05-01 quoted 280 nM pre-Phase 2a).
- Partial DTS depletion: ✗ — DTS empties to 0 by t ≈ 5 s. Tied to D7.
- Sustained plateau: ✗ — no SOCE plateau because DTS empties before STIM1
  dimerisation can build up enough to engage Orai.
- Return towards baseline: cyt collapses to ~3 nM by t ≈ 30 s and stays there
  (the model's natural rest, not Dolan's 100 nM; tied to D4 and D7).

The headline dissertation-relevant result — the peak amplitude is in the
biological range — does pass. The shape (sustained SOCE plateau, return to
100 nM) is blocked on the open Phase 2 issue.

### 7.3 SOCE dependence (Phase 3) — **PASS-WITH-DEVIATIONS**

Phase 3 validation against Dolan 2014 Fig. 4 was implemented as a
two-condition driver (`runscripts/manual/runPhase3.py`) and a
comparison plot (`models/platelet/analysis/phase3_dolan_fig4.py`) on
2026-05-06. See lab-book `2026-05-06-phase3-results.md` for the full
write-up; figure at `reports/figures/phase3-dolan-fig4-2026-05-06.png`.

Acceptance criteria (Dolan 2014 Fig. 3B + lab-book Phase 3 framing):

| Criterion | Rule | Measured | Result |
|---|---|---|---|
| Active (+Ca_ex) | peak Ca_cyt > 200 nM | 299 nM | ✓ PASS |
| Active (−Ca_ex) | peak Ca_cyt > 200 nM | 298 nM | ✓ PASS |
| SOCE differential | \|peak(+) − peak(−)\| ≥ 100 nM | 1 nM | ✗ FAIL |
| Peak in Dolan ±30% (+Ca_ex) | 315–585 nM | 299 nM | ✗ FAIL |
| Peak in Dolan ±30% (−Ca_ex) | 192–358 nM | 298 nM | ✓ PASS |

3/5 pass. Both failures (SOCE differential and +Ca_ex peak amplitude)
trace to the same upstream cause — DTS empties before SOCE can
establish a plateau (§6.8 D7) — and are tracked separately in
issues #19 and #22. The model reproduces Dolan's *active*
filtering criterion in both conditions but does not yet reproduce the
SOCE-dependent shape that distinguishes them.

### 7.4 Analysis plot (`models/platelet/analysis/single/calcium_trace.py`)

The Phase 1 implementation produces a 5-panel figure:

| Panel | Content |
|-------|---------|
| 1 | [Ca²⁺]_cyt vs time + Dolan 2014 ±30% reference band |
| 2 | CaM sub-states (free / Ca₂·CaM / Ca₄·CaM) — stacked area |
| 3 | PMCA sub-states (basal vs CaM-activated) — stacked area |
| 4 | DTS Ca²⁺ + STIM1 dimer count (twin axis) |
| 5 | IP₃ + SOCE flux (twin axis) |

The plot is invoked by
`runscripts/manual/analysisPlatelet.py [run_dir] --plot calcium_trace`
and the rendered PNG/SVG land in `[run_dir]/.../plotOut/low_res_plots/`
where the webapp's Explore Plots tab can browse them.

> **Code:** the 5-panel figure builder is `Plot.do_plot` in
> `models/platelet/analysis/single/calcium_trace.py:74` (class definition
> at `:71`). The Dolan Fig. 4 reference curve constants and analytical
> approximation (`_dolan_reference_nM`) are at
> `models/platelet/analysis/single/calcium_trace.py:42–60`.

---

## 8. Issue tracking

| Issue | Title | Deliverables |
|-------|-------|-------------|
| **#24** | Ca²⁺ data and dataclass | `internal_state.py` (27 species — 22 v0.2 + 6 Phase 1 CaM / PMCA·CaM), `calcium_signalling.py` dataclass, `dts` compartment |
| **#25** | CalciumDynamics process | `calcium_signalling.py` process, ODE solver, 6-state IP3R, SERCA, PMCA, SOCE, IP3 forcing |
| **#26** | Ca²⁺ listener | `calcium_dynamics.py` listener recording 5 columns per timestep |
| **#27** | Ca²⁺ analysis plot | `calcium_transient.py` — transient shape, IP3R states, dual-axis, Dolan reference lines |
| **#43** | ~~IP3R upgrade~~ | Merged into #25 — 6-state model implemented from the start |

---

## 9. Open questions for review

The decisions in §3.1 and §5 are provisional until this list is cleared. Items
marked *resolved* have a decision in the body of the doc and are listed here
only for visibility; reviewers should push back if they disagree.

1. **DTS volume — resolved (§5).** 4.3% of cytosol = 0.258 fL (Purvis direct
   measurement). The 2% alternative (Purvis Monte Carlo median, 0.12 fL) would
   give a smaller Ca²⁺ store and less stable oscillations; the larger value is
   consistent with Dolan 2014's parameter set.

2. **Cytoplasmic volume — resolved (§3.1).** V_cyt = 6 fL (Purvis), giving a
   resting count of 361. The Sveshnikova 2025 figure of 3 fL (count = 181)
   likely reflects a different definition of "cytoplasm" (total cell vs
   non-granular fraction). Documented in the analysis.

3. **V_IM:** Set to −60 mV (Dolan cluster analysis upper bound). If SOCE flux is
   too large or too small, this is the first parameter to adjust.

4. **IP3R copy number:** Use Dolan 2014 Table S1 value (1,328) rather than Burkhart
   proteomics sum (4,850). Rationale: not all proteomic copies are in the DTS
   membrane and functionally gated. The Dolan value is from a filtered population
   that satisfies homeostatic constraints.

5. **CaM / PMCA activation — superseded (§3.5, §6.8 D1).** The v0.2
   simplification (Caride 2007 *basal* constants only) was upgraded in
   Phase 1 (commit `f3080c40`) to the full Caride 2007 Table 3 5-state
   CaM-coupled scheme plus the CaM Ca²⁺-binding ladder. This was earlier
   than planned (originally a v0.3 item) because the basal-only PMCA could
   not clamp the cytosolic Ca²⁺ peak. The v0.2 transient-decay caveat in
   §7.2 no longer applies.

6. **ATP / ADP / Pi starting counts — resolved.** ATP[c] = 1.084×10⁷
   molecules (3 mM × 6 fL × N_A, lower end of the platelet 3–5 mM
   cytosolic range; Holmsen 1979/1981 metabolic pool). ADP[c] = 1.084×10⁶
   molecules (0.3 mM, ATP:ADP = 10:1 typical resting ratio). Pi[c] =
   3.61×10⁵ molecules (100 µM × 6 fL). All three are in `_MOLECULES`
   (`internal_state.py`) and accounted for by the SERCA / PMCA pumping
   stoichiometry. Note: the earlier draft's figure of ~10⁹ came from
   misinterpreting bulk-assay "~10⁻¹² mol ATP" as per-platelet rather than
   per-assay aliquot — single platelets contain femtomoles, not picomoles,
   of ATP.

7. **Stochastic future (flagged, not blocking):** The PLC-Gq bottleneck (~1 molecule)
   will require a hybrid deterministic/stochastic approach in v0.3. At that stage,
   the upstream module should use tau-leaping (Gillespie approximation) while the
   Ca²⁺ ODE core remains deterministic. This is a known architectural choice, not
   a v0.2 concern.

---

*Document status: living design / as-built reference. Last revised
2026-05-06 to capture Phase 1 (CaM ladder + 5-state CaM-coupled PMCA +
MWC SOCE) and Phase 2a (basal PM Ca²⁺ leak, SERCA initial-conditions pre-equilibration).
Cross-references to lab books `lab-book-2026-05-01-phase1-complete.md` and
`lab-book-2026-05-05-phase2a-investigation.md` are authoritative for the
diagnoses behind §6.8.*
