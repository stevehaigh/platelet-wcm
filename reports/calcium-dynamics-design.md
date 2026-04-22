---
pdf-engine: xelatex
mainfont: "STIX Two Text"
monofont: "Menlo"
fontsize: 11pt
geometry: margin=2.5cm
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

![Schematic Ca²⁺ transient — expected output shape](/Users/steve/github/wcEcoli/reports/figures/ca2-transient-reference.png)
*Schematic showing the expected Ca²⁺ transient shape. The IP3 forcing function
(lower panel) drives a rapid peak in cytosolic Ca²⁺, which then decays to a
SOCE-sustained plateau before returning to baseline. This is the shape we aim
to reproduce from Dolan & Diamond 2014 Fig. 4.*

This is a **Dolan-core-first** strategy — meaning we implement the validated
Ca²⁺ core from Dolan & Diamond (2014) first, using their published parameters,
before adding the upstream receptor cascade that would generate IP3 in a real
cell:

```
v0.2:  IP3 forcing → Ca²⁺ core (IP3R + SERCA + PMCA + SOCE)
v0.3:  Upstream receptor cascade  (P2Y1 → Gq → PLCβ → IP3 production)
v0.4:  P2Y12 modulation           (Gi → AC → cAMP → PKA → IP3R inhibition)
```

Each milestone is independently testable and produces publishable results.

### Why Dolan 2014 and not Purvis 2008?

Purvis & Bhatt (2008) is the foundational modelling paper for platelet Ca²⁺
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

## 2. Signal pathway (brief)

See `reports/calcium-signalling-pathway-design.md` for the full biology.

At rest, cytosolic Ca²⁺ is ~100 nM, DTS Ca²⁺ is ~250 µM. On stimulation:

```
IP3 spike
  └─ IP3R (DTS membrane): opens → Ca²⁺ floods cytosol (peak ~400 nM)
       ├─ SERCA (DTS membrane): Ca²⁺-ATPase refills the DTS store
       ├─ PMCA (plasma membrane): Ca²⁺-ATPase ejects Ca²⁺ from cell
       └─ DTS depletion → STIM1 oligomerises → gates Orai1 → SOCE
```

Key compartments: cytosol (6 fL, Purvis 2008 direct measurement), DTS (0.258 fL = 4.3% of cytosol, Purvis 2008 glucose-6-phosphatase staining), extracellular (treated as an infinite reservoir at fixed 1.2 mM, Dolan 2014).

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

| Variable | Description | Compartment | Resting count |
|----------|-------------|-------------|---------------|
| `CA2_CYT` | Cytosolic free Ca²⁺ | `[c]` | **361** (100 nM, 6 fL cyt) |
| `CA2_DTS` | DTS stored Ca²⁺ | `[dts]` | 38,842 (250 µM, 0.258 fL DTS) |
| `CA2_EX` | Extracellular Ca²⁺ | `[e]` | fixed (infinite reservoir) |
| `IP3` | Inositol trisphosphate | `[c]` | 181 |
| `IP3R_n` | IP3R neutral state | `[dts]` | 809 |
| `IP3R_o` | IP3R open state | `[dts]` | 261 |
| `IP3R_a` | IP3R active state | `[dts]` | 65 |
| `IP3R_i1` | IP3R inhibited-1 state | `[dts]` | 167 |
| `IP3R_i2` | IP3R inhibited-2 state | `[dts]` | 25 |
| `IP3R_s` | IP3R shut state | `[dts]` | 1 |
| `SERCA_E1` | SERCA empty, E1 | `[dts]` | 5,920 |
| `SERCA_E2` | SERCA empty, E2 | `[dts]` | 5,927 |
| `SERCA_E1Ca` | SERCA·Ca²⁺ in E1 | `[dts]` | 6 |
| `SERCA_E1PCa` | SERCA phosphorylated·Ca²⁺ | `[dts]` | 7 |
| `SERCA_E2PCa` | SERCA phosphorylated E2·Ca²⁺ | `[dts]` | 4 |
| `SERCA_E2P` | SERCA phosphorylated, empty | `[dts]` | 28 |
| `PMCA` | PMCA unbound | `[m]` | 765 |
| `PMCA_Ca` | PMCA·Ca²⁺ complex | `[m]` | 4 |
| `STIM1` | STIM1 free monomer | `[dts]` | 438 |
| `STIM1_Ca` | STIM1 DTS-bound (inactive) | `[dts]` | 3,805 |
| `STIM1_dim` | STIM1 dimer (active sensor) | `[dts]` | 22 |
| `ORAI1` | Orai1 (all channels start closed; opening tracked via the `STIM1·Orai1*` complex count in §3.6, not as a separate state variable) | `[m]` | 1,447 |

Initial counts from Dolan 2014 Table S1 representative configuration.
Conversion: 6 fL cytosol, 0.258 fL DTS, N_A = 6.022×10²³.

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
IP3(t) = IP3_rest × [1 + (fold−1) × (1 − e^{−t/τ_rise}) × e^{−max(0, t−t_peak)/τ_decay}]

fold     = 5.5   (peak amplitude relative to rest)
τ_rise   = 3.0 s
t_peak   = 3.0 s
τ_decay  = 60.0 s
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

### 3.3 IP3R: 6-state Markov model (Sneyd & Dufour 2002)

The IP3 receptor transitions between six states. Ca²⁺-dependent activation
and inhibition produce the biphasic open probability required for oscillations.
A Hill function cannot reproduce this behaviour.

![IP3R 6-state Markov model](/Users/steve/github/wcEcoli/reports/figures/ip3r-state-machine.png)
*IP3R 6-state Markov model, adapted from Sneyd & Dufour (2002) type-2 kinetics
as parameterised in Purvis & Bhatt (2008) Table 1. Green states (o, a) are
Ca²⁺-conducting. Blue (n) is the resting neutral state. Purple states (i1, i2)
are Ca²⁺-inhibited. Red (s) is the shut state. Transitions are Ca²⁺- and IP3-
dependent.*

Rate constants (Purvis 2008 Table 1, Sneyd & Dufour 2002 type-2 kinetics):

| Transition | Forward | Reverse |
|------------|---------|---------|
| n ↔ o (IP3 binding) | k₂ = 37.4 µM⁻¹s⁻¹, l₄ = 1.7 µM⁻¹s⁻¹ | k₋₂ = 1.4 s⁻¹, l₋₄ = 2.5 µM⁻¹s⁻¹ |
| n ↔ i1 (Ca²⁺ inhibition) | k₁ = 0.64 µM⁻¹s⁻¹, L₁ = 0.12 µM | k₋₁ = 0.04 s⁻¹ |
| o ↔ a (Ca²⁺ activation) | k₄ = 4 µM⁻¹s⁻¹ | k₋₄ = 0.54 µM⁻¹s⁻¹ |
| o ↔ s (shutting) | k₃ = 11 µM⁻¹s⁻¹, L₅ = 54.7 µM | k₋₃ = 29.8 s⁻¹ |
| i1, i2 ↔ s | l₂ = 1.7 s⁻¹ | l₋₂ = 0.8 s⁻¹ |

Open probability:
```
P_o = ((0.9 · IP3R_a + 0.1 · IP3R_o) / IP3R_total)⁴
```

Ca²⁺ flux through IP3R (Purvis 2008, eq. 13):
```
J_IP3R = γ_IP3R × N_IP3R × P_o × (V_IM − E_Ca,IM) / (z × F)

γ_IP3R   = 10 pS                                single-channel conductance
V_IM     = −60 mV                               DTS membrane potential (see §5)
E_Ca,IM  = (RT / zF) × ln([Ca²⁺]_dts / [Ca²⁺]_cyt)   Ca²⁺ equilibrium (Nernst) potential
                                                across the DTS membrane (z=2)
```

The driving force is `(V_IM − E_Ca,IM)`, i.e. membrane potential minus the
Ca²⁺ equilibrium potential. At resting [Ca²⁺]_dts/[Ca²⁺]_cyt = 2,500, E_Ca,IM
≈ +104 mV, so with V_IM = −60 mV the driving force is ≈ −164 mV — a strong
gradient moving Ca²⁺ out of the DTS into the cytosol when the channel opens.

### 3.4 SERCA: E1–E2 cycle

Six-state enzymatic cycle (Purvis 2008, Dode 2002):

```
E2 ⇌ E1 ⇌ E1·Ca²⁺ → E1P·Ca²⁺ ⇌ E2P·Ca²⁺ → E2P → E2
          ↑ (cytosol)                      ↓ (DTS)
```

Each cycle transports 2 Ca²⁺ ions from cytosol to DTS at the cost of 1 ATP.
Rate constants: see `calcium-data-provenance.md` § "SERCA cycle".

### 3.5 PMCA

Two-state simplified model (Caride 2007 parameters, Purvis 2008 Table 1):

```
PMCA + Ca²⁺_cyt ⇌ PMCA·Ca²⁺ → PMCA + Ca²⁺_ex

KM1 = 0.5 µM,  KM2 = 1.0 µM,  kcat = 8.9 s⁻¹
```

> **TODO (before implementation):** the KM values above are shown in concentration
> units (µM) — the original draft had them as mM⁻¹, which is dimensionally wrong
> for a Michaelis constant. The numerical values here are provisional; confirm
> against Purvis 2008 Table 1 / Caride 2007 and update.

CaM-mediated activation is simplified in v0.2: PMCA treated as constitutively
active at basal rate. Full CaM kinetics deferred to v0.3.

### 3.6 SOCE (STIM1 / Orai1)

Dolan 2014 MWC allosteric model:

```
STIM1·Ca²⁺_dts ⇌ STIM1_free  (Ca²⁺ release from DTS-bound STIM1)
STIM1_free      ⇌ STIM1_dim   (dimerisation — active sensor)
STIM1_dim + Orai1 → STIM1·Orai1*  (CRAC channel opening)
Ca²⁺_ex  --Orai1*--> Ca²⁺_cyt     (Ca²⁺ entry)
```

Key parameter: DTS membrane potential V_IM = −60 mV (Dolan clustering analysis
shows SOCE-active configurations cluster at V_IM > −70 mV; use upper bound).

SOCE current:
```
I_SOC = g_SOC × P_open × (V_PM − E_Ca)
```

where E_Ca is the Ca²⁺ Nernst potential and g_SOC is set by the MWC model.

---

## 4. Architecture

### 4.1 Signal flow

The diagram below shows how the wcEcoli process/state architecture has been
adapted for platelet Ca²⁺ signalling. The E. coli biological processes
(transcription, metabolism, etc.) are replaced by the platelet-specific
`CalciumDynamics` process, but the underlying framework — state containers,
time-stepping, ATP/ADP partitioning — is carried over unchanged.

![Ca²⁺ signalling process architecture](/Users/steve/github/wcEcoli/reports/figures/calcium-process-architecture.png)
*ATP coupling shown in red. Green compartment = cytosol, blue = DTS, purple = plasma membrane.*

### 4.2 File contents

#### `reconstruction/platelet/dataclasses/process/calcium_signalling.py`

**Class `CalciumSignalling`** — builds the ODE system at startup (called once,
not each timestep). Stores the species list and all rate constants. In v0.2
these are hardcoded; v0.3 migration swaps in TSV parsing without touching the
process code.

Key method: **`molecules_to_next_time_step(counts, volume_cyt, volume_dts,
nAvogadro, dt, t_sim)`** — converts integer molecule counts to concentrations
(µM), runs the Ca²⁺ ODE sub-stepper (BDF solver, one timestep duration),
injects the IP3 time curve at `t_sim`, and returns integer count deltas and
estimated ATP cost.

#### `models/platelet/processes/calcium_signalling.py`

**Class `CalciumDynamics(Process)`** — thin wrapper around the ODE solver.
Follows the standard wcEcoli process lifecycle:

- **`initialize(sim, sim_data)`** — creates `BulkMolecule` views for all 22
  Ca²⁺ species and for ATP, ADP, and Pi.
- **`calculateRequest()`** — calls `molecules_to_next_time_step` to estimate
  changes for this timestep, then requests the required molecules (Ca²⁺ species
  needed by the ODE) and the estimated ATP from the framework's partitioner.
- **`evolveState()`** — applies the Ca²⁺ count changes and performs ATP/ADP/Pi
  accounting using the ATP actually allocated (not the requested amount, which
  may have been reduced by the partitioner).

#### `models/platelet/listeners/calcium_dynamics.py`

Records each timestep:

| Column | Description |
|--------|-------------|
| `ca_cyt_nM` | `CA2_CYT[c]` count converted to nM |
| `ca_dts_uM` | `CA2_DTS[dts]` count converted to µM |
| `ip3_nM` | IP3 count converted to nM |
| `ip3r_open_fraction` | `(IP3R_o + IP3R_a) / IP3R_total` |
| `soce_flux` | estimated Orai1 current (µM/s) |

---

## 5. Parameter and initial condition table

Full provenance in `calcium-data-provenance.md`. Summary for implementation:

### Volumes

| Compartment | Value | Source |
|-------------|-------|--------|
| Cytosol | 6.0 fL | Purvis 2008 (direct measurement) |
| DTS | 0.258 fL (4.3% × cyt) | Purvis 2008 glucose-6-phosphatase stain |
| Extracellular | infinite (fixed reservoir) | — |

### Resting concentrations

| Species | Conc | Count (calc) | Source |
|---------|------|-------------|--------|
| Ca²⁺_cyt | 100 nM | 361 | Purvis 2008 |
| Ca²⁺_dts | 250 µM | 38,842 | Dolan 2014 Fluo-5N measurement |
| Ca²⁺_ex | 1.2 mM | fixed | Dolan 2014 |
| IP3 | 50 nM | 181 | Sveshnikova 2025 / Dolan 2014 middle |
| IP3R total | — | 1,328 | Dolan 2014 Table S1 |
| SERCA total | — | 11,892 | Dolan 2014 Table S1 |
| PMCA total | — | 769 | Dolan 2014 Table S1 |
| STIM1 total | — | 4,265 | Dolan 2014 Table S1 |
| Orai1 total | — | 1,447 | Dolan 2014 Table S1 |

### Key physical constants

| Constant | Value |
|----------|-------|
| N_A | 6.022×10²³ mol⁻¹ |
| V_IM (DTS membrane potential) | −60 mV (Dolan upper bound; used in IP3R driving force (V_IM − E_Ca,IM) §3.3 and SOCE current §3.6) |
| V_PM (plasma membrane potential) | −60 mV (used in SOCE current I_SOC §3.6) |
| T | 310 K (37°C) |

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

In v0.2, ATP is added to the molecule inventory and the process requests it.
If the ATP pool is insufficient (unlikely in the first version since nothing
else consumes it yet), the framework's partitioner will allocate less than
requested. The process scales pump activity proportionally: it advances SERCA
and PMCA only to the extent that allocated ATP supports, and the unsupported
pump turnover does not happen. This keeps mass and energy balance correct
even under starvation. A warning is logged on any shortfall so the failure
mode is visible; a proper metabolic process in a future milestone will make
shortfalls biologically meaningful.

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

---

## 7. Validation strategy

### 7.1 Resting state stability (first test, no stimulus)

Run 300 s with no IP3 forcing (fold=1). Pass criteria:
- Ca²⁺_cyt stays within 80–120 nM
- Ca²⁺_dts stays within 200–300 µM
- All six IP3R state fractions (n, o, a, i1, i2, s; see §3.1 state table)
  stay within 10% of their Dolan Table S1 resting values

### 7.2 Ca²⁺ transient shape (primary validation)

Run 300 s with IP3 forcing (fold=5.5, τ_rise=3 s, τ_decay=60 s).
Pass criteria (Dolan 2014 Fig 4, +extracellular Ca²⁺ condition):
- Peak Ca²⁺_cyt: 200–500 nM, reached within 15–20 s
- Partial DTS depletion: Ca²⁺_dts drops to 30–70% of resting
- Sustained plateau above baseline (SOCE-dependent)
- Return to ~resting levels within 300 s

### 7.3 SOCE dependence

Run with Ca²⁺_ex = 0 (EDTA condition).
Pass criteria:
- Transient peak similar, but plateau absent / faster decay
- Matches Dolan 2014 Fig 4C (no-extracellular-Ca²⁺ curve)

### 7.4 Analysis plot (`calcium_transient.py`)

- Time series: Ca²⁺_cyt (nM), Ca²⁺_dts (µM) on dual axes
- Inset: IP3R open fraction over time
- Annotation: forcing function equation, Dolan 2014 reference values as dashed lines
- Second panel: IP3R state fractions stacked area chart

---

## 8. Issue tracking

| Issue | Title | Deliverables |
|-------|-------|-------------|
| **#24** | Ca²⁺ data and dataclass | `internal_state.py` (22 species), `calcium_signalling.py` dataclass, `dts` compartment |
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

5. **CaM / PMCA activation:** Simplify to constitutive PMCA in v0.2. Full CaM
   kinetics in v0.3. Agree?

6. **ATP molecule inventory:** ATP, ADP, Pi must be added to the platelet molecule
   inventory to support pump accounting. Starting count: ATP ~1–2 × 10⁷ molecules
   per platelet, estimated as cytosolic [ATP] ~3–5 mM in 6 fL (standard platelet
   range; Holmsen 1981, Gerson 2008 bulk assays rescaled per cell). The earlier
   draft's figure of ~10⁹ came from interpreting Gerson's "~10⁻¹² mol ATP" as
   per-platelet rather than per-assay aliquot, which is two orders too high —
   single platelets contain femtomoles, not picomoles, of ATP. Confirm the
   precise starting count before implementation.

7. **Stochastic future (flagged, not blocking):** The PLC-Gq bottleneck (~1 molecule)
   will require a hybrid deterministic/stochastic approach in v0.3. At that stage,
   the upstream module should use tau-leaping (Gillespie approximation) while the
   Ca²⁺ ODE core remains deterministic. This is a known architectural choice, not
   a v0.2 concern.

---

*Document status: DRAFT — for review before implementation*
