# Platelet Ca²⁺ Dynamics: Implementation and Calibration Notes

**Branch:** `platelet`
**Relevant commits:** issues #24, #25, #26, #27, #45, #46
**Status:** Equations restored to primary-source form (Sneyd & Dufour IP3R, Po⁴ tetramer cooperativity, Nernst flux, Hoover/Dolan MWC SOCE, Purvis SERCA Vmax). Integrator runs stably. Resting state does **not** hold 100 nM at the Dolan Table S1 ICs — those ICs are not a steady state of the corrected model. See §4.

---

## 1. What Was Built (current state)

Three files implement IP3-mediated Ca²⁺ signalling in the platelet WCM:

| File | Role |
|------|------|
| `reconstruction/platelet/dataclasses/process/calcium_signalling.py` | ODE right-hand side, all kinetic parameters, MWC SOCE solver |
| `models/platelet/processes/calcium_dynamics.py` | Thin process wrapper calling the ODE each timestep |
| `models/platelet/listeners/calcium_trace.py` | Records ca_cyt, ca_dts, ip3, SOCE flux each timestep |
| `models/platelet/analysis/single/calcium_trace.py` | 3-panel Ca²⁺ transient plot vs Dolan 2014 reference |

The ODE covers 21 molecular species:

- **Ca²⁺ pools**: cytosol (`CA2_CYT`), DTS (`CA2_DTS`)
- **IP3R**: 6-state Sneyd & Dufour Markov gating (n, o, a, i1, i2, s)
- **SERCA**: 6-state Purvis 2008 / Dode 2002 cycle
- **PMCA**: 2-state Caride 2007 basal kinetics (CaM-coupled scheme deferred)
- **SOCE / STIM1 / Orai1**: 3-state STIM cycle + Hoover/Dolan MWC channel gating

IP3 is forced externally via `ip3_forcing_uM(t)` — a Dolan 2014 Fig. S2 approximation. The simulation timestep is 1 s; the ODE is solved with SciPy BDF within each step.

---

## 2. Equations restored to primary-source form

A previous calibration pass (commit 8df289acb) reduced several constants to balance against an *incorrect* IP3R flux formula and thereby achieve flat resting Ca²⁺. With the corrected equations the original primary-source values are restored. Five changes:

### 2.1 IP3R open probability — Po⁴ + 0.9 weighting

Purvis 2008 Table 1 / Dolan 2014 Eq. 4:

```
Po = (0.9 · IP3R_a/total + 0.1 · IP3R_o/total)⁴
```

The previous code used `Po = (a + 0.1·o)/total` — linear, missing the `0.9` factor and (more importantly) the `^4` exponent. The exponent encodes that all four IP3R subunits in the tetrameric channel must be in conducting (active or open) conformations simultaneously.

Effect: at the Dolan IC subunit fractions (a/total = 0.049, o/total = 0.197), linear Po = 0.069; Po⁴ = 1.65×10⁻⁵. The amplitude difference is ~3000× — but, crucially, the *response* to IP3 stimulation now scales as the 4th power of the active-subunit fraction, giving the characteristic sharp transient.

### 2.2 IP3R Sneyd & Dufour rate laws (φ-functions)

Purvis 2008 Table 1 lists the full Sneyd & Dufour 2002 type-2 rate laws as rational functions of `[Ca²⁺]_cyt`:

```
n + Ca²⁺ → i1:  v = [n] · ((k₁·L₁ + l₂)·[Ca] / (L₁ + [Ca]·(1 + L₁/L₃)))
                 − [i1] · (k₋₁ + l₋₂)
n + IP3  → o:   v = [n] · [IP3] · ((k₂·L₃ + l₄·[Ca]) / (L₃ + [Ca]·(1 + L₃/L₁)))
                 − [o] · ((k₋₂ + l₋₄·[Ca]) / (1 + [Ca]/L₅))
o + Ca²⁺ → a:   v = [o] · ((k₄·L₅ + l₆)·[Ca] / (L₅ + [Ca]))
                 − [a] · (L₁·(k₋₄ + l₋₆) / (L₁ + [Ca]))
a + Ca²⁺ → i2:  v = [a] · ((k₁·L₁ + l₂)·[Ca] / (L₁ + [Ca]))
                 − [i2] · (k₋₁ + l₋₂)
o ↔ s:          v = [o] · (k₃·L₅ / (L₅ + [Ca])) − [s] · k₋₃
```

with `L₁=0.12 µM, L₃=0.025 µM, L₅=54.7 µM, l₄=1.7 µM⁻¹s⁻¹, l₋₄=2.5 µM⁻¹s⁻¹, l₆=4707 s⁻¹, l₋₆=11.4 s⁻¹` plus the existing `k₁..k₋₄, l₂, l₋₂`.

The previous code used simple `k·X·Y − k₋·Z` mass action and dropped the L1/L3/L5/l₄/l₋₄/l₆/l₋₆ constants entirely. Mass action does not satisfy detailed balance for the Sneyd network — Wegscheider's loop product around `n→o→s→i1→n` gave 0.0145 at Dolan ICs, not 1.

There is one minor topological change versus the previous code: the spurious `i1↔s` and `i2↔s` mass-action transitions (which were not in Purvis Table 1) are removed. The Sneyd parameters `l₂, l₋₂` are subsumed into the φ-function rate laws, where they belong.

### 2.3 IP3R Ca²⁺ flux — Nernst form

Purvis Table 1 row "Ca²⁺ release from DTS" / Dolan Eq. 4:

```
I_IP3R = γ_IP3R · N_channels · Po · (NA / zF) · (ψ_IM − E_Ca,IM)
       γ_IP3R = 10 pS              (Zschauer 1988)
       ψ_IM   = −60 mV             (Dolan upper bound, V_IM > −70 mV cluster)
       E_Ca,IM = (RT/zF) · ln([Ca²⁺]_dts / [Ca²⁺]_cyt)
       RT/zF (z=2, 37 °C) = 13.4 mV
```

The previous code used an empirical conductance form `K_IP3R_FLUX × Po × ([Ca]_dts − [Ca]_cyt)` and dropped the constant. The Nernst form is logarithmic in the gradient, so as DTS depletes the flux saturates rather than scaling linearly — which matters during the transient.

### 2.4 SERCA k_bind_f restored to Purvis 1×10¹⁵ M⁻²·s⁻¹

The earlier calibration reduced k_bind_f from 1×10³ µM⁻²s⁻¹ to 2.11 µM⁻²s⁻¹ (470× smaller) so that resting SERCA throughput would balance the *incorrectly small* IP3R flux produced by the linear-Po + empirical-conductance combination.

With Po⁴ + Nernst, the corrected IP3R leak is ~1.13×10⁵ ions/s into cytosol at Dolan ICs. The Purvis SERCA Vmax (k_bind_f = 1000) gives 2 × v_bind ≈ 1.18×10⁵ ions/s extrusion at 100 nM cytosolic Ca²⁺ and Dolan enzyme states — in approximate balance.

### 2.5 SOCE — full Hoover/Dolan MWC

The SOCE block now implements:

1. **STIM1 cycle** (mass-action with detailed-balance constants):
   ```
   STIM1·Ca²⁺_dts ⇌ STIM1_free + Ca²⁺_dts    (Ca²⁺ release from EF-hand)
   2 STIM1_free   ⇌ STIM1_dim                 (diffusion-limited dimerisation)
   ```
   Constants `k_release_r = 3.475×10⁻³`, `k_dim_f = 1.15×10⁻⁴` are derived from
   detailed balance at the Dolan Table S1 IC (st_Ca=3805, st_free=438, st_dim=22, ca_dts=250 µM).

2. **Puncta entry** (Dolan 2014 Eq. 2):
   ```
   (STIM2)p = qp · STIM2_dimers
   qp = α · ([Ca²⁺]_cyt^n / (KM^n + [Ca²⁺]_cyt^n)) + 0.01
   α = 0.2,  KM = 0.5 µM,  n = 4   (KM, n picked mid-range from Dolan's scan)
   ```

3. **MWC channel gating** (Hoover & Lewis 2011 / Dolan Eq. 3):
   ```
   For i = 0..4 STIM dimers bound:
     [CSi]/[C] = (4 choose i) · a^(i(i−1)/2) · (Ka·Sf)^i
     [OSi]/[CSi] = f^i · L
     [OSi]/[C]   = L · f^i · (4 choose i) · a^(i(i−1)/2) · (Ka·Sf)^i

   Po = Σ[OSi] / Σ([CSi] + [OSi])
   ```
   with Hoover Fig 4B parameters `L=10⁻⁴, f=14.2, a=0.5`. `Ka` is rescaled from
   Hoover's a.u. (100 with Stotal=3.2) to platelet dimer counts (`Ka = 2`); see code comment.

   Sf is solved each ODE step by Newton-like iteration on the mass-balance
   `Sf + bound_count(Sf) = (STIM2)p`.

4. **SOC current** (Dolan Eq. 4):
   ```
   I_SOC = γ_SOC · N_orai_channels · Po · (NA/zF) · (ψ_PM − E_Ca,PM)
   γ_SOC = 0.3 fS    (calibrated from resting balance — see §3)
   ψ_PM = −60 mV    (Dolan)
   E_Ca,PM = (RT/zF) · ln([Ca²⁺]_ex / [Ca²⁺]_cyt)
   ```

This closes issues #45 and #46 as designed.

---

## 3. Constants still requiring calibration (and how they were derived)

After restoring the primary-source equations, four parameters cannot be drawn directly from the literature for our platelet count-units:

| Constant | Value | Derivation |
|---------|-------|-------------|
| STIM `k_release_r` | 3.475×10⁻³ µM⁻¹s⁻¹ | Detailed balance at Dolan IC: `= k_release_f × st_Ca / (st_free × ca_dts) = 0.1 × 3805 / (438 × 250)` |
| STIM `k_dim_f` | 1.15×10⁻⁴ count⁻¹s⁻¹ | Detailed balance at Dolan IC: `= k_dim_r × st_dim / st_free² = 1.0 × 22 / 438²` |
| MWC `Ka` | 2.0 | Rescaled from Hoover a.u.: `Ka_platelet = Ka_Hoover × (Stotal_Hoover / Sf_saturating_platelet) = 100 × (3.2 / 170) ≈ 1.9 → 2` |
| `γ_SOC` | 0.3 fS | Resting flux balance: `SOCE_rest = PMCA_steady_rest ≈ 76 ions/s` at Po(MWC, Sf_rest) ≈ 1.2×10⁻³ |

All other rate constants are taken directly from primary sources (Purvis 2008 Table 1 for IP3R/SERCA/IP3R conductance, Caride 2007 Table 3 for PMCA basal, Hoover 2011 Fig 4B for MWC L/f/a, Dolan 2014 Eqs. 2–4 for SOCE structure).

---

## 4. Resting state with corrected equations

### 4.1 Why the Dolan ICs are no longer a steady state

The Dolan 2014 Table S1 representative ICs are filtered against four homeostatic/dynamic constraints **for the Dolan ODE system** — full Sneyd rate laws, Po⁴, Nernst flux, Caride/CaM-coupled PMCA, MWC SOCE with their fitted (n, KM). Our v0.2 differs in two ways that prevent the ICs from being a true equilibrium:

1. **PMCA is basal-only (no CaM coupling).** Caride 2007's full 5-state CaM-coupled scheme drops the apparent KM into the sub-µM range and raises k_cat during a Ca²⁺ rise. With basal-only kinetics, PMCA cannot match the SOCE+IP3R inflow at peak Ca²⁺.
2. **(KM, n) for puncta entry are not from the same scan as the rest of the ICs.** Dolan scanned KM and n jointly with the protein copy numbers and V_IM to satisfy resting balance. Picking mid-range values without rerunning the scan means SOCE at moderate `[Ca²⁺]_cyt` is unbalanced.

### 4.2 Observed behaviour (no IP3 forcing)

Starting from the Dolan IC, the corrected ODE (see `tests/sim/test_simulation.py`) produces:

- `[Ca²⁺]_cyt`: 100 nM → ~600 nM peak at t≈0.5 s (SERCA enzyme states relax) → settles toward ~50 nM by t=100 s
- `[Ca²⁺]_dts`: 250 µM → empties to ~0 µM by t≈30 s
- STIM1_dim: grows from 22 → ~800 (full dimerisation as STIM_Ca releases)
- Total Ca²⁺ in cell: drops from ~39 000 atoms to ~hundreds (PMCA extrusion exceeds SOCE replenishment at this calibration)

This is biologically wrong but **mathematically self-consistent** with the corrected primary-source equations and the unscanned (KM, n).

### 4.3 Observed behaviour with IP3 forcing

Under the Dolan 2014 Fig. S2 IP3 curve, peak `[Ca²⁺]_cyt` reaches several µM (vs. Dolan Fig 4's ~400 nM). The DTS empties completely. The transient *shape* is qualitatively right (rise on IP3 stimulus, then decay) but amplitude and recovery dynamics are off because:

- SERCA has no upper bound from CaM/PMCA pulling Ca²⁺ out of the cytosol fast enough during peak
- DTS empties → IP3R flux saturates at zero gradient
- SOCE plateau is small (γ_SOC was calibrated for resting balance, not peak)

---

## 5. What still needs doing

| Priority | Task |
|----------|------|
| High | Re-scan (KM, n) for puncta entry within the Dolan homeostatic/dynamic constraint framework so that the IC is a genuine steady state of *our* ODE — not just Dolan's |
| High | Implement Caride 2007 full 5-state CaM-coupled PMCA — with proper CaM activation, peak PMCA can clamp the Ca²⁺ transient at biological values |
| Medium | Re-derive Dolan-style "representative IC" by integrating to steady state at low IP3, then use that IC as the simulation start — guarantees zero drift |
| Medium | Add cytosolic Ca²⁺ buffering (calmodulin total ~20 000, plus calbindin) — currently all Ca²⁺ is free |
| Low | Validate the analytical detailed-balance derivations of `k_release_r`, `k_dim_f`, `Ka`, `γ_SOC` against an independent numerical fit |

---

## 6. Methodology notes

### Unit conventions

All concentrations inside the ODE run in µM. Count ↔ µM conversion uses platelet geometry (Purvis 2008 direct measurement):

```python
_UM_PER_COUNT_CYT = 2.77e-4 µM/count   # cytosol volume 6 fL
_UM_PER_COUNT_DTS = 6.44e-3 µM/count   # DTS volume 0.258 fL (4.3 % cyt)
```

The same number of Ca²⁺ atoms crossing the DTS membrane appears as different concentration changes on each side (DTS is ~23× smaller than cytosol). Fluxes in the ODE are always in counts/s — the IP3R Nernst form returns ions/s directly via `NA/(zF) × current`.

### MWC unit rescaling

Hoover & Lewis 2011 fit their MWC scheme to ICRAC vs. Orai expression in HEK cells, where STIM and Orai are reported in a.u. proportional to GFP fluorescence. The fitted `Ka = 100` works in *those* units, not in molecule counts. Following Dolan 2014's algebraic-equation approach, we transfer Hoover's structure intact (L, f, a are dimensionless and unchanged) and rescale `Ka` so that the saturating biology (`Sf ≈ 170 platelet STIM2 dimers`) maps to Hoover's saturating regime (`Ka × Stotal ≈ 320`). This gives `Ka_platelet ≈ 2`. The MWC Po(Sf) curve has the same shape as Hoover's, just with Sf measured in dimer counts rather than a.u.

### γ_SOC calibration

CRAC channels have a literature single-channel conductance of ~24 fS (Prakriya & Lewis 2002, Vig 2006). With L = 10⁻⁴ giving ~0.04 channels open at rest in a 361-channel cell, that produces ~5×10³ ions/s of resting SOCE — far above the PMCA basal extrusion (~76 ions/s) and enough to drive Ca²⁺ runaway. The literature value reflects measurements at saturating Po, in cell types where many channels are open. For the platelet's small effective open count at rest, we use a calibrated effective γ_SOC of 0.3 fS, derived from the requirement `SOCE_rest = PMCA_steady_rest`.

This is a v0.2 simplification. The full single-channel calibration would require either (a) re-scanning Dolan's (n, KM) phase space and using a literature γ_SOC, or (b) treating SOCE flux probabilistically below the integer-count threshold.

---

*Report rewritten 2026-04-29 after restoring primary-source equations. See also: `calcium-dynamics-design.md`, `calcium-signalling-pathway-design.md`, `calcium-data-provenance.md`.*
