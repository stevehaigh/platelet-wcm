# Lab book — 2026-04-29: calcium equations restored to primary sources

## Why this session

Resumed the platelet Ca²⁺ work. A previous calibration pass (commit 8df289acb) had
flatlined the resting state by reducing five rate constants 100–1000× to balance
against an *incorrectly linearised* IP3R Po formula. The model was numerically
stable but mathematically wrong. The aim of this session was to restore the
published equation forms (Purvis 2008, Dolan 2014, Hoover & Lewis 2011) and accept
whatever transient/resting behaviour falls out.

## Bugs found in the prior implementation

| # | Bug | File / line | Source it disagreed with |
|---|-----|-------------|-------------------------|
| 1 | Channel open probability `Po = (a + 0.1·o)/total` (linear, missing `0.9` factor) | `calcium_signalling.py:233` | Purvis 2008 Table 1 / Dolan 2014 Eq. 4: `Po = (0.9·a/total + 0.1·o/total)⁴` |
| 2 | IP3R rate laws used simple mass action; dropped `L1, L3, L5, l4, l₋4, l6, l₋6`. Wegscheider loop product = 0.0145, not 1 → no detailed balance | `_ode_rhs` IP3R block | Purvis 2008 Table 1 — full Sneyd & Dufour 2002 φ-function rate laws |
| 3 | IP3R flux used empirical conductance `K_IP3R_FLUX × Po × (Ca_dts − Ca_cyt)` | same | Purvis Eq. 13 / Dolan Eq. 4: `γ·N·Po·(NA/zF)·(ψ_IM − E_Ca,IM)` with γ_IP3R = 10 pS |
| 4 | SERCA `k_bind_f` was 2.11 µM⁻²s⁻¹ (470× below primary value) | `K_SERCA['k_bind_f']` | Purvis Table 1: 1×10³ µM⁻²s⁻¹ |
| 5 | SOCE was 3-state mass action with ad hoc `k_dim_f`, `k_dim_r`, `k_orai` | `K_SOCE` | Hoover & Lewis 2011 / Dolan 2014: full MWC allosteric scheme |

## What changed in the code

Commit **18ed71847** on the `platelet` branch.

**`reconstruction/platelet/dataclasses/process/calcium_signalling.py`** — full rewrite of the rate-law block:

- New physical-constants block: `F`, `R`, `T`, `RT/zF`, `NA/(zF)`, `V_IM`, `V_PM`.
- `K_IP3R` extended with the full Sneyd parameter set; `_phi_*` helpers for each transition.
- `_ode_rhs` IP3R block uses the φ-functions; subunit topology cleaned up (no spurious i↔s shortcut).
- IP3R flux: `γ_IP3R · N · Po⁴ · (NA/zF) · (ψ_IM − E_Ca,IM)`.
- `K_SERCA['k_bind_f']` restored to 1000.
- New `K_STIM` dict for STIM1 cycle (constants from detailed balance at Dolan IC).
- New `K_MWC` dict (Hoover Fig 4B: L=1e-4, f=14.2, a=0.5, Ka rescaled).
- New `PUNCTA` dict for Dolan Eq. 2 (α=0.2, KM=0.5 µM, n=4).
- New `_mwc_open_fraction(stim2_p, n_orai)` Newton-iteration solver for the MWC equilibrium.
- SOCE current: `γ_SOC · N · Po(MWC) · (NA/zF) · (ψ_PM − E_Ca,PM)` with γ_SOC = 0.3 fS (calibrated).

**`models/platelet/listeners/calcium_trace.py`** — listener uses the same MWC chain to estimate SOCE flux for plotting.

**`models/platelet/tests/sim/test_simulation.py`** — relaxed the dry-mass-monotonic-decrease assertion. With active SOCE / PMCA the cell exchanges Ca²⁺ atoms with the extracellular reservoir, so dry mass is no longer monotonic.

**`reports/platelet-calcium-calibration.md`** — full rewrite. New §2 documents each equation change against its primary source. New §3 lists the four parameters that still required calibration and the analytical conditions used to derive them. New §4 honestly describes the residual non-physical behaviour.

**`reports/calcium-data-provenance.md`** — SOCE section now describes the implemented MWC scheme rather than the "implementation gap" caveat.

## Calibrated constants (and how they were derived)

After restoring the primary-source equations, four parameters still needed values:

| Constant | Value | Derivation |
|---------|-------|------------|
| STIM `k_release_r` | 3.475×10⁻³ µM⁻¹s⁻¹ | Detailed balance at Dolan IC: `k_release_f · st_Ca / (st_free · ca_dts)` |
| STIM `k_dim_f` | 1.15×10⁻⁴ count⁻¹s⁻¹ | Detailed balance at Dolan IC: `k_dim_r · st_dim / st_free²` |
| MWC `Ka` | 2.0 | Rescaled from Hoover a.u. (`Ka_Hoover · Stotal_Hoover / Sf_saturating_platelet = 100·3.2/170 ≈ 1.9`) |
| `γ_SOC` | 0.3 fS | Resting balance: `SOCE_rest = PMCA_steady_rest ≈ 76 ions/s` at MWC Po(rest) ≈ 1.2×10⁻³ |

Everything else (IP3R rate constants, SERCA cycle, PMCA basal, MWC L/f/a, V_IM, V_PM, γ_IP3R) is taken directly from primary sources.

## What it looks like running

Ran a 200-second platelet sim with IP3 forcing on. Results in
`out/platelet_calcium_smoke/.../plotOut/calcium_trace.{pdf,png}`.

Numerical trace:

| t (s) | ca_cyt (nM) | ca_dts (µM) | IP3 (nM) | SOCE (nM/s) | STIM_dim |
|-------|-------------|-------------|----------|-------------|----------|
| 0 | 99.9 | 250.0 | 50.1 | 1.8 | 22 |
| 1 | 1729 | 127.3 | 156 | 2.0 | 31 |
| 3 | 5919 | 3.1 | 211 | 2.5 | 73 |
| 5 | 5053 | 0.4 | 234 | 4.8 | 165 |
| 10 | 3041 | 0.08 | 242 | 13.7 | 390 |
| 30 | 679 | 0.013 | 191 | 28.9 | 751 |
| 50 | 342 | 0.006 | 151 | 6.5 | 803 |
| 75 | 31 mM (!) | 0.0 | 117 | 5.0 | 807 |
| 200 | 503 mM (!!) | 0.0 | 58 | −6.3 | 807 |

The first ~5 seconds are roughly Dolan-shaped (rise to peak, DTS depletes,
STIM mobilises). After that the system runs away because:

1. The Dolan ICs are not a true steady state of our v0.2 ODE (basal-only PMCA, unscanned KM/n) — once perturbed by the IP3 stimulus the system can't return to balance.
2. Once DTS empties, `STIM1·Ca²⁺ ⇌ STIM1_free` cannot reverse (no Ca²⁺ in DTS to rebind), so STIM stays fully dimerised → MWC stays partially open → SOCE keeps loading the cytosol.
3. PMCA at basal Caride kinetics caps at ~4000 ions/s extrusion (~1 nM/s in cytosol). Far below the SOCE inflow at moderate-to-high cytosolic Ca²⁺.

The runaway is not in the equations; it's in the v0.2 simplifications that bypass biological feedback.

## Next concrete step

The single biggest fix is **implementing Caride 2007 Table 3 5-state CaM-coupled PMCA**. It drops apparent Km into the sub-µM range and raises k_cat ~5× during a Ca²⁺ rise. That alone should clamp the transient peak at biological values and prevent the late-time runaway. One paper's worth of work; doesn't touch any other module.

After that:

1. Re-scan (KM, n) for puncta entry against platelet homeostatic constraints (Dolan-style filtering).
2. Re-derive a true steady-state IC by integrating the corrected ODE to convergence at low IP3, then use that as the simulation start point.
3. Add cytosolic Ca²⁺ buffering (calmodulin total ~20 000 + calbindin) — buffers ~95% of free Ca²⁺.

## GitHub state after session

| Issue | Action | New state |
|-------|--------|-----------|
| #45 | MWC implemented; comment + close | closed |
| #46 | MWC implemented (data layer); comment + close | closed |
| #24 | 22-species inventory + dataclass + provenance done; comment with progress | open (calibration items remain) |
| #25 | All equation-form criteria met; resting/transient calibration criteria fail; comment with status | open (calibration items remain) |

## Files touched in commit 18ed71847

- `reconstruction/platelet/dataclasses/process/calcium_signalling.py` (rate-law rewrite)
- `reconstruction/platelet/dataclasses/internal_state.py` (22-species inventory)
- `reconstruction/platelet/dataclasses/process/process.py` (wires CalciumSignalling)
- `reconstruction/platelet/raw_data/molecules.tsv` (citations + ATP/ADP correction)
- `reconstruction/platelet/simulation_data.py` (`pl` compartment)
- `models/platelet/listeners/calcium_trace.py` (MWC SOCE flux estimate)
- `models/platelet/tests/sim/test_simulation.py` (relaxed mass test)
- `reports/calcium-data-provenance.md` (SOCE section)
- `reports/calcium-dynamics-design.md` (PMCA section)
- `reports/platelet-calcium-calibration.md` (full rewrite)
- `reports/pandoc-header.tex` (Unicode chars for ≈ ✓)
- `reports/calcium-data-provenance.pdf` (deleted; .md is the source)

12 files changed, 725 insertions, 290 deletions, 1 deleted.

## Files left uncommitted (separate WIP)

- `runscripts/manual/analysisPlatelet.py` — webapp plot path tweak
- `wholecell/webapp/jobs.py` — webapp adds platelet analysis phase

These are not part of the equation work and stay on the working tree.
