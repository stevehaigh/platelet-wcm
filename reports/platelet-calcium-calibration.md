# Platelet Ca²⁺ Dynamics: Implementation and Calibration Notes

**Branch:** `platelet`  
**Relevant commits:** issues #24, #25, #26, #27  
**Status:** Partially calibrated — simulation runs without runaway, but resting Ca²⁺ drifts to ~235 nM instead of the target 100 nM (see §4).

---

## 1. What Was Built (v0.2)

Three files implement IP3-mediated Ca²⁺ signalling in the platelet WCM:

| File | Role |
|------|------|
| `reconstruction/platelet/dataclasses/process/calcium_signalling.py` | ODE right-hand side, all kinetic parameters |
| `models/platelet/processes/calcium_dynamics.py` | Thin process wrapper calling the ODE each timestep |
| `models/platelet/listeners/calcium_trace.py` | Records ca_cyt, ca_dts, ip3, SOCE flux each timestep |
| `models/platelet/analysis/single/calcium_trace.py` | 3-panel Ca²⁺ transient plot vs Dolan 2014 reference |

The ODE covers 21 molecular species:

- **Ca²⁺ pools**: cytosol (ca_cyt), DTS (ca_dts)
- **IP3R**: 6-state De Young–Keizer gating (n, o, a, i1, i2, s)
- **SERCA**: 6-state Purvis 2008 cycle (E1, E2, E1Ca, E1PCa, E2PCa, E2P)
- **PMCA**: 2-state Ca²⁺/CaM-activated pump
- **SOCE/STIM1/Orai1**: 3-state simplified MWC (STIM1_Ca, STIM1_free, STIM1_dim)

IP3 is forced externally via `ip3_forcing_uM(t)` — a Dolan 2014 Fig. S2 approximation (50 nM rest → 275 nM peak at t≈3–5 s → τ=60 s decay). The simulation timestep is 1 s; the ODE is solved with SciPy BDF within each step.

---

## 2. The Five Calibration Bugs

The initial parameter set caused runaway behaviour: ca_cyt reached ~60 µM within seconds (biological target: ~400 nM peak) and DTS filled rather than depleting. Numerical audit identified five miscalibrated parameters:

### 2.1 IP3R resting flux too large (`K_IP3R_FLUX`)

| | Value | Resting flux |
|---|---|---|
| Original | 0.30 | 5143 nM/s |
| Fixed | **0.004** | 69 nM/s |

At `K_IP3R_FLUX = 0.30`, with Dolan Table S1 resting open probability (~0.069) and DTS–cytosol gradient (~250 µM), the IP3R alone would drain the entire DTS in 1.6 s. The correct value is derived by requiring the resting IP3R flux to match the SERCA pumping capacity (~69 nM/s).

### 2.2 SOCE conductance too large (`k_orai`)

| | Value | Resting SOCE |
|---|---|---|
| Original | 0.001 | 26,400 nM/s |
| Fixed | **7.99×10⁻⁷** | 21 nM/s |

The fix requires SOCE = PMCA efflux at rest:

```
k_orai = PMCA_efflux_µMs / (st_dim × (Ca_ex − ca_cyt))
       = 21.09×10⁻³ / (22 × 1200) = 7.99×10⁻⁷
```

### 2.3 STIM1 Ca-binding not at equilibrium (`k_release_r`)

| | Value | v_STIM1 at Dolan IC |
|---|---|---|
| Original | 0.001 | +271 ct/s (net release) |
| Fixed | **3.475×10⁻³** | ≈ 0 ct/s |

The fix requires detailed balance at the Dolan 2014 Table S1 resting IC (st_ca=3805, st_free=438, ca_dts=250 µM):

```
k_release_r = k_release_f × st_ca / (st_free × ca_dts)
            = 0.1 × 3805 / (438 × 250) = 3.475×10⁻³ µM⁻¹s⁻¹
```

### 2.4 SERCA binding rate too large (`k_bind_f`, `k_bind_r`)

| | k_bind_f | k_bind_r | Resting throughput |
|---|---|---|---|
| Original | 1000 µM⁻²s⁻¹ | 10 s⁻¹ | 59,140 ct/s |
| Fixed | **2.1101 µM⁻²s⁻¹** | **0.021101 s⁻¹** | 124 ct/s |

The Purvis 2008 value (1×10¹⁵ M⁻²s⁻¹) converts correctly to 1×10³ µM⁻²s⁻¹, but this is the *in-vitro* maximum rate for isolated SERCA — far too fast for a system with only 5920 SERCA molecules and 100 nM cytosolic Ca²⁺. The correct `k_bind_f` is derived by requiring SERCA throughput to balance the corrected IP3R + SOCE − PMCA at rest (~124 ct/s), while maintaining the 0.1 µM Km (so `k_bind_r = k_bind_f × Km²`).

---

## 3. Current Behaviour After Fixes

With the five fixes applied, the simulation no longer diverges. Qualitative behaviour:

- Ca_cyt starts at 100 nM, drifts upward to ~235 nM over ~5 s, then stabilises
- DTS depletes slowly (250 → 242 µM over 10 s)
- STIM1 and Orai1 states remain near their resting values
- IP3 forcing does produce an elevated Ca²⁺ response (not yet fully validated against Dolan 2014 Fig. 3)

---

## 4. Known Limitation: Residual Resting Drift

**The system does not hold 100 nM at rest.** This is a calibration limitation, not a code bug.

### Root cause

The Dolan 2014 Table S1 initial conditions for the SERCA protein states:

| State | Dolan IC | Equilibrium at new k_bind_f |
|---|---|---|
| SERCA_E1Ca | 6 | ~0.18 |
| SERCA_E1PCa | 7 | ~0.26 |
| SERCA_E2PCa | 4 | ~0.64 |
| SERCA_E2P | 28 | ~2.1 |

The Dolan ICs were derived for the fast-SERCA regime (`k_bind_f = 1000`). With the re-calibrated `k_bind_f = 2.11`, the true equilibrium for these enzyme states is ~14× smaller. As the simulation runs, the SERCA protein states relax toward their true equilibrium, which shifts the Ca²⁺ balance point.

At the false equilibrium (~235 nM), the IP3R flux (1489 ct/s) exactly balances SERCA (1384 ct/s) + PMCA (172 ct/s) at the higher cytosolic Ca²⁺. This is a self-consistent but biologically wrong steady state.

### Why this wasn't trivially fixed

The IP3R flux at rest is proportional to (ca_dts − ca_cyt) × Po. As ca_cyt rises from 100 to 235 nM, Po barely changes (~0.069) but the SERCA throughput increases steeply (∝ ca_cyt²: from 124 to 692 ct/s). The system finds a new balance at 235 nM rather than at 100 nM. To force the equilibrium back to 100 nM, `K_IP3R_FLUX` would need to be further reduced — but this would also reduce the IP3-evoked peak response.

The fundamental issue is that the three timescales — IP3R gating, SERCA cycling, and STIM1 activation — are calibrated from different papers with different cell types and scales, and no consistent fitting has been done. This is tracked as issue #45 (SOCE MWC) and issue #24 (parameter curation).

---

## 5. What Still Needs Doing

| Priority | Task |
|----------|------|
| High | Re-derive SERCA ICs at new `k_bind_f` and update `internal_state.py` — would eliminate the SERCA transient |
| High | Global parameter fitting: tune `K_IP3R_FLUX` and `k_orai` jointly so both resting and peak Ca²⁺ match Dolan targets |
| Medium | Implement Dolan 2014 MWC allosteric STIM1/Orai1 model (issue #46) |
| Medium | Add cytosolic Ca²⁺ buffering (calmodulin, calbindin) — currently all Ca²⁺ is free |
| Low | Replace Purvis 6-state SERCA with a simpler Michaelis–Menten pump consistent with the Dolan ODE level of detail |

---

## 6. Methodology Notes

### Analytical derivation approach

All five parameter fixes were derived analytically from detailed-balance / flux-balance conditions at the Dolan 2014 Table S1 resting IC, then verified numerically by integrating the ODE for 5–200 s and inspecting all flux terms.

### Unit conventions

All concentrations inside the ODE run in µM. Count ↔ µM conversion uses platelet geometry:

```python
_UM_PER_COUNT_CYT = 2.77e-4 µM/count  # cytosol volume 6 fL
_UM_PER_COUNT_DTS = 6.44e-3 µM/count  # DTS volume 0.258 fL
```

The same number of Ca²⁺ atoms crossing the DTS membrane appears as different concentration changes on each side (DTS is ~23× smaller than cytosol). Fluxes in the ODE are always in counts/s; µM values are only used for reaction rate laws.

### Why the Purvis SERCA k_bind_f was wrong by 3 orders of magnitude

The Purvis 2008 paper gives `k_bind_f = 1×10¹⁵ M⁻²s⁻¹`. This converts to `1×10³ µM⁻²s⁻¹` — which looks correct. The error is that this rate was measured for isolated detergent-solubilised SERCA at saturating ATP in vitro. In vivo at 100 nM cytosolic Ca²⁺ (below Km) with ~6000 SERCA copies, the effective throughput is determined by the Km, not the Vmax. The system needs the throughput to be ~124 ct/s, which back-calculates to `k_bind_f ≈ 2.11 µM⁻²s⁻¹`.

---

*Report written 2026-04-29. See also: `calcium-dynamics-design.md`, `calcium-signalling-pathway-design.md`, `reconstruction-platelet-design.md`.*
