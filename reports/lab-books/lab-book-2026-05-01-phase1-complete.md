---
title: "Lab book — 2026-05-01: Phase 1 complete, Phase 2 roadmap"
---

# Lab book — 2026-05-01: Phase 1 complete, Phase 2 roadmap

## Session summary

Two sessions (2026-04-30 / 2026-05-01) completed the Phase 1 work: CaM Ca²⁺ buffering and the 5-state
Caride PMCA are now implemented and debugged. The primary Ca²⁺ runaway is eliminated. Three bugs were
found and fixed, and the acceptance-criteria plot is now a comprehensive 5-panel figure.

---

## Where we started (pre-Phase 1)

| Symptom | Root cause |
|---------|-----------|
| Ca²⁺ runaway to ~6 µM | PMCA was basal-only (no CaM activation), no cytosolic buffering |
| Resting drift | IP3R Nernst flux formula uses 10 pS single-channel conductance: IP3R >> SERCA at rest |
| IP3R flux wrong form | Missing Po⁴ tetramer cooperativity (prior code used linear Po) |
| SOCE wrong model | Ad hoc 3-state mass action; replaced with Hoover/Dolan MWC in commit `18ed71847` |

---

## GitHub issues — status

| Issue | Title | Status |
|-------|-------|--------|
| #17–#22 | v0.1 scaffold (sim loop, states, listeners, analysis) | ✓ Closed |
| #23 | Webapp platelet support | ✓ Closed |
| #24 | Data + dataclass (inventory, internal state) | Open — will close after Phase 2 IC is finalised |
| #25 | CalciumDynamics process — resting stability + transient | Open — peak criterion ✓; DTS depletion pending Phase 2 |
| #26 | CalciumTrace listener | ✓ Closed (listener updated in Phase 1) |
| #27 | Analysis: calcium_trace plot | ✓ Closed (5-panel figure written) |
| #45 | SOCE: replace ad hoc model with MWC | ✓ Closed (`18ed71847`) |
| #46 | γ_SOC calibration | ✓ Closed (`18ed71847`) |
| #47 | Phase 1: CaM + Caride 5-state CaM-coupled PMCA | Effectively done — peak criterion ✓; DTS pending Phase 2 |
| #48 | Phase 2: re-derive resting IC | **Open — next** |
| #49 | Phase 3: Dolan 2014 Fig 4 validation | Open |

---

## Work done in Phase 1 (commits `f3080c40a`, `ec511f7be`, `e2ffdbb70`)

### 1. CaM Ca²⁺-binding ladder (Caride 2007 Table 3 steps 6–7)

CaM is now a state-variable cytosolic Ca²⁺ buffer with three sub-states:

$$
\text{CaM}_\text{free} \xrightarrow[k_{6r}]{k_6 [\text{Ca}^{2+}]^2} \text{Ca}_2\text{·CaM} \xrightarrow[k_{7r}]{k_7 [\text{Ca}^{2+}]^2} \text{Ca}_4\text{·CaM}
$$

| Parameter | Value | Source |
|-----------|-------|--------|
| $k_6$ | 2.669 µM⁻²·s⁻¹ | Caride 2007 Table 3 |
| $k_{6r}$ | 2.682 s⁻¹ | Caride 2007 Table 3 |
| $k_7$ | 170.4 µM⁻²·s⁻¹ | Caride 2007 Table 3 |
| $k_{7r}$ | 1.551 s⁻¹ | Caride 2007 Table 3 |

Total CaM = 20 481 (Dolan 2014 Table S1). ICs set at Caride detailed-balance
equilibrium with 100 nM Ca²⁺: CaM\_free = 20 062, Ca₂·CaM = 200, Ca₄·CaM = 219.

### 2. 5-state CaM-coupled PMCA (Caride 2007 Table 3 steps 4–5, 8–10)

The PMCA now has both a basal path and a CaM-activated path that is ~5.5× faster:

```
Basal path (steps 4–5):
  PMCA + Ca²⁺ ⇌ PMCA·Ca  →  PMCA + Ca²⁺_ex      k_cat = 5.5 s⁻¹

CaM-activated path (steps 8–10):
  PMCA + Ca₄·CaM ⇌ Ca₄·CaM·PMCA  +  Ca²⁺  ⇌  Ca₄·CaM·PMCA·Ca  →  Ca₄·CaM·PMCA + Ca²⁺_ex
                    k8 = 0.2 µM⁻¹s⁻¹
                    k9 = 50 µM⁻¹s⁻¹
                    k10 = 30 s⁻¹
```

Step 11 (slow CaM dissociation, τ ≈ 20 min) is omitted for Phase 1 — it accumulates
PMCA in a dead-end `PMCA·CaM` state on a timescale well beyond the 200 s transient.

### 3. New molecule inventory (6 additional species)

Added to `internal_state.py` and `calcium_signalling.py`:
`CaM_free[c]`, `Ca2_CaM[c]`, `Ca4_CaM[c]`,
`Ca4_CaM_PMCA[pl]`, `Ca4_CaM_PMCA_Ca[pl]`, `PMCA_CaM[pl]`

### 4. Bugs found and fixed

| Bug | Symptom | Fix |
|-----|---------|-----|
| **Phantom IP3R flux** | When `CA2_DTS` count hit 0 atoms, the Nernst fallback `e_ca_im = 0` gave −60 mV driving force from an empty store, creating Ca²⁺ from nothing | Guard: `flux_ip3r_ions_s = 0.0` when `ca_dts ≤ 0` |
| **PMCA step-10 recycling missing** | Every CaM-activated pump stroke consumed a PMCA molecule (414 molecules gone by t=30 s) | Add `+v_cam_pmca_cat` to `dy[Ca4_CaM_PMCA]` so the empty complex is returned to the queue |
| **CaM ICs not at equilibrium** | First-timestep CaM loading burst grabbed ~34 k Ca²⁺ ions (Dolan Table S1 values assumed no explicit Ca²⁺ binding kinetics) | Pre-set ICs to Caride detailed-balance solution at 100 nM |

### 5. CalciumTrace listener — additional columns

The listener now records 8 extra columns: `cam_free`, `ca2_cam`, `ca4_cam`,
`ca4_cam_pmca`, `ca4_cam_pmca_ca`, `pmca_cam`, `pmca_free`, `pmca_ca`.

### 6. 5-panel analysis figure

`models/platelet/analysis/single/calcium_trace.py` rewritten from 3 panels to 5:

| Panel | Content |
|-------|---------|
| 1 | [Ca²⁺]_cyt vs time + Dolan 2014 ±30% reference band |
| 2 | CaM sub-state stacked-area (cytosolic buffering) |
| 3 | PMCA activation-state stacked-area (basal vs CaM-activated) |
| 4 | DTS Ca²⁺ + STIM1 dimer count (twin axis) |
| 5 | IP₃ + SOCE flux (twin axis) |

---

## What the model now does

### Simulation output — 2026-05-01 run (`platelet_may1`)

Plot: [`calcium_trace_phase1_2026-05-01.pdf`](calcium_trace_phase1_2026-05-01.pdf) ·
[`.png`](calcium_trace_phase1_2026-05-01.png)

**Acceptance criteria:**

| Criterion | Value | Status |
|-----------|-------|--------|
| Peak [Ca²⁺]_cyt within first 20 s | **280 nM** | ✓ PASS (200–800 nM) |
| At t = 50 s | 1 nM | ✓ PASS (<1000 nM) |
| At t = 200 s | 1 nM | ✓ PASS (<1000 nM) |
| DTS min > 0 µM | **0 µM** | ✗ FAIL — Phase 2 required |

**What each panel shows:**

- **Panel 1:** The Ca²⁺ transient is a 280 nM spike at t = 1 s (inside the Dolan ±30% band),
  followed by a rapid collapse to ~0 nM. The peak is biologically correct; the *shape* is wrong
  (should be a 3–5 s rise with a ~60 s SOCE plateau). The shape problem is the Phase 2 issue —
  see below.

- **Panel 2:** At t = 1 s, Ca₄·CaM spikes to ~8 000 molecules as the CaM ladder loads up during
  the transient, sequestering ~32 k Ca²⁺ ions. CaM then slowly unloads as Ca²⁺ drops. CaM
  conservation is maintained throughout.

- **Panel 3:** PMCA sub-state dynamics are clean and conserved (~769 PMCA total). At peak Ca²⁺,
  CaM-activated PMCA (Ca₄·CaM·PMCA·Ca) appears but is near-zero instantaneously because
  k₁₀ = 30 s⁻¹ cycles it faster than the 1 s timestep can capture.

- **Panel 4:** DTS empties in ~1 s (vertical red drop at t ≈ 0) — the Phase 2 calibration problem.
  STIM1 dimers correctly rise from 22 to ~810 as the store depletes (the store-sensing mechanism
  is working).

- **Panel 5:** IP₃ rises to ~275 nM at t ≈ 3 s, decays with τ = 60 s. SOCE flux reaches 1–3.5 nM/s
  but insufficient to sustain elevated Ca²⁺ once the initial transient collapses.

---

## The key Phase 2 problem — IP3R flux calibration

The DTS empties in ~1 s because the Nernst formula with `GAMMA_IP3R_S = 10 pS` (the
Zschauer 1988 single-channel conductance, cited in Purvis 2008 Table 1) gives:

$$
J_{\text{IP3R,rest}} = \gamma_{\text{IP3R}} \cdot N \cdot P_o \cdot \frac{N_A}{zF} \cdot (V_{IM} - E_{\text{Ca},IM})
$$
$$
\approx 10\,\text{pS} \times 1328 \times 1.65 \times 10^{-5} \times 0.164\,\text{V} \times 3.12 \times 10^{18}
= 112\,000\ \text{ions/s}
$$

SERCA (limited by E2P·Ca occupancy in the Dolan IC) refills at only ~6 600 ions/s. Net drain
rate ≈ 105 400 ions/s → 38 842-atom DTS gone in **~0.35 s**.

The 10 pS value is a biophysical single-channel electrophysiology measurement, not a
whole-cell effective conductance. In a 6 fL model, applying it at face value gives an
IP3R flux that is ~17 000× larger than SERCA can match at rest.

**Calibrated target:** for IP3R flux = SERCA at resting $P_o$:

$$
\gamma_{\text{IP3R, calibrated}} = \frac{6600}{1328 \times 1.65 \times 10^{-5} \times 0.164 \times 3.12 \times 10^{18}} \approx 0.6\ \text{fS}
$$

---

## Next steps and expected timeline

### Phase 2 — IP3R recalibration + true resting IC (issue #48)  ·  ~1 week

**Step 2a — Recalibrate `GAMMA_IP3R_S` (1 line change, 1 hour)**

In `reconstruction/platelet/dataclasses/process/calcium_signalling.py` line 149:

```python
GAMMA_IP3R_S = 0.6e-15   # calibrated from SERCA balance; 0.6 fS vs measured 10 pS
```

This should slow DTS drain from 0.35 s to ~5–10 s, producing a gradual Ca²⁺ rise rather than
an instantaneous spike.

**Step 2b — Derive true resting IC (half-day)**

The Dolan Table S1 ICs are not at mathematical steady state for our ODE (we've demonstrated
they're out of balance by a factor of ~17 on IP3R vs SERCA). After recalibrating `GAMMA_IP3R_S`:

1. Run `ip3_forced=False` with the Dolan IC for 300–500 s (baseline IP₃ = 50 nM).
2. Extract the converged counts as the new `_MOLECULES` initial counts in `internal_state.py`.
3. Verify: at the true IC, `d[Ca²⁺]_cyt/dt ≈ 0` and DTS is stable.

**Step 2c — Rerun Phase 1 validation from true IC**

Apply the IP₃ pulse (200 s, `ip3_forced=True`) from the new IC. Verify:
- Peak [Ca²⁺]_cyt: 200–800 nM
- Ca²⁺ rises over 3–5 s (not an instantaneous spike)
- DTS depletes partially (not to zero)
- SOCE plateau ~100–200 nM for 30–60 s
- Returns to ~100 nM resting by t = 200 s

*Expected result once γ_IP3R is correctly scaled:* the Dolan shape should emerge naturally
because every other component (CaM ladder, 5-state PMCA, MWC SOCE, STIM1 cycle) is already
working correctly.

### Phase 3 — Dolan 2014 Fig. 4 validation (issue #49)  ·  ~1 week after Phase 2

Dolan 2014 Fig. 4 shows Ca²⁺ transients with and without extracellular Ca²⁺:

1. **With extracellular Ca²⁺** (normal): peak ~400 nM, SOCE plateau ~200 nM, return to ~100 nM.
2. **Without extracellular Ca²⁺** (`CA_EX_UM = 0`): peak ~350 nM, no SOCE plateau, faster return.

Steps:
1. Run both conditions from the Phase 2 IC.
2. Overlay on the Dolan Fig. 4 digitised reference (or the schematic reference curve already in
   the analysis plot).
3. Pass criteria (from `calcium-next-steps-plan.md`): peak within ±30% of Dolan reference;
   plateau within ±50%; correct qualitative response to Ca²⁺ removal.
4. If SOCE is too weak/strong: adjust `KM_uM` or `n` in the `PUNCTA` dict (Dolan scans these
   parameters as free variables; our current values `KM=0.5 µM, n=4` are mid-range).

### Expected state in ~2 weeks (by ~2026-05-15)

If Phase 2 (recalibration) and Phase 3 (Dolan Fig. 4) both land cleanly:

- The simulation produces a biologically realistic Ca²⁺ transient matching the Dolan 2014
  reference within ±30% for both conditions (with/without extracellular Ca²⁺).
- Issues #24, #25, #47, #48, #49 all closed.
- The model is a complete, validated Ca²⁺ signalling module for the platelet WCM —
  the headline result for the dissertation Methods section.
- The 5-panel plot (`calcium_trace.pdf`) is dissertation-ready with minimal formatting changes.

### Phase 4 — Optional parameter scan (issue pending)

If Phase 3 validation criteria fail by more than ±30%, scan `(KM_uM, n)` for the STIM puncta
entry Hill function within Dolan's homeostatic constraints. This is estimated ~2 days; only
trigger it if Phase 3 fails.

### Phase 5 — P2Y1 upstream cascade (v0.3, stretch)

Not required for the dissertation. Can be described as future work. The architecture is
prepared (IP₃ can switch from forced to endogenous in the `ip3_forced` flag).

---

## File index

| File | Purpose |
|------|---------|
| `reconstruction/platelet/dataclasses/process/calcium_signalling.py` | ODE + all rate constants. `GAMMA_IP3R_S` is the Phase 2 change target. |
| `reconstruction/platelet/dataclasses/internal_state.py` | Initial counts. Updated with equilibrium CaM ICs. |
| `models/platelet/processes/calcium_dynamics.py` | Process wrapper (calls the ODE dataclass). |
| `models/platelet/listeners/calcium_trace.py` | Records 14 columns each timestep. |
| `models/platelet/analysis/single/calcium_trace.py` | 5-panel validation figure. |
| `reports/calcium-next-steps-plan.md` | Full 5-phase roadmap with effort estimates. |
| `reports/caride-2007-pmca-rate-constants.md` | Every Caride 2007 rate constant with implementation notes. |
| `reports/calcium_trace_phase1_2026-05-01.pdf` | Current 5-panel figure (this session's run). |

---

## Build and test commands

```bash
# Run the simulation (200 s) and generate the plot
PYTHONPATH=. python3 runscripts/manual/runPlateletSim.py --length 200 platelet_smoke
PYTHONPATH=. python3 runscripts/manual/analysisPlatelet.py \
    platelet_smoke --plot calcium_trace

# Unit tests (9 tests, all should pass)
PYTHONPATH=. python3 -m pytest models/platelet/tests/ -v

# Quick ODE sanity check (no sim runner needed)
PYTHONPATH=. python3 -c "
from reconstruction.platelet.dataclasses.process.calcium_signalling import (
    CalciumSignalling, MOLECULE_NAMES)
import numpy as np
cs = CalciumSignalling(None)
print('ODE imports OK; N_SPECIES =', len(MOLECULE_NAMES))
"
```

---

## 2026-05-05 — repo migration, CLAUDE.md refresh, webapp cleanup

Maintenance session. No biology changes; the calcium model state is unchanged
from the 2026-05-01 entry above (Phase 1 complete; Phase 2 IP3R recalibration
still next). Three threads of housekeeping:

### 1. Repo migration finished (wcEcoli → platelet-wcm)

The fork has now been fully renamed and pruned. State as of `f83d9d8b`:

| Area | Change |
|------|--------|
| Remote | `stevehaigh/platelet-wcm` (was `stevehaigh/wcEcoli`) |
| Branches | `main`, `webapp` |
| Removed | `models/ecoli/`, `reconstruction/ecoli/`, all `.pyx` Cython, `cloud/`, ParCa runscripts (`runParca.py`, `runSim.py`, `analysisSingle.py`, …) |
| Docker | `cloud/docker/` → `docker/runtime/` + `docker/webapp/` |
| Azure | ACI URL is now `platelet-wcm.uksouth.azurecontainer.io`; resource names kept as `wcecoli-*` to avoid recreating shared infra |
| README | Rewritten as a platelet-WCM project README with framework attribution to CovertLab/wcEcoli |
| CI | `.github/workflows/ci.yml` runs pytest + mypy on every push/PR to `main` |

`make compile` is no longer required — there are no remaining `.pyx` files.

### 2. `CLAUDE.md` rewritten for platelet context

The previous CLAUDE.md still described an E. coli simulation. Rewritten to cover:

- The platelet model layout (`models/platelet/`, `reconstruction/platelet/`)
- Current process / listener / analysis inventory
- Calcium signalling overview (IP3R / SERCA / PMCA / SOCE / CaM, Dolan validation target)
- Platelet compartment table (`c, dts, dg, ag, m, pl, e`) and the rationale for `pl` over `pm`
- Build & run without ParCa (`SimulationDataPlatelet` constructed directly)
- Webapp + reports/lab-books layout

Memory entries also written under `~/.claude/projects/.../memory/` so future
sessions know the dissertation deadline (~mid-July 2026) and Phase 2 is the
next live work.

### 3. Web UI cleanup

The Dash webapp at `wholecell/webapp/` was still presenting the E. coli surface
(27 variants, 10 regulation toggles, 8 nutrient-shift presets). Stripped to a
platelet-only UI.

**Removed from `tabs/configure.py`:**

- `VARIANT_NAMES` — 27 E. coli variant names
- `TOGGLES` — all 10 (ppGpp regulation, tRNA charging, D-period division,
  translation supply, superhelical density, variable elongation ×2, mechanistic
  translation/AA transport, tRNA attenuation)
- 8 of 9 `PRESETS` (wildtype, AA shift up/down, anaerobic, rich media,
  multi-seed, ppGpp sweep, acetate)
- "Variant index range / Generations / Seeds" form fields

**New `tabs/configure.py` form:**

| Field | Notes |
|-------|-------|
| Length (seconds) | Default 200 s; integer ≥ 1 |
| Random seed | Default 0; not currently consumed by any stochastic process |
| Description | Free-text, used in the run directory name |

**New presets** (3, all platelet):

| Preset | Length | Purpose |
|--------|--------|---------|
| ⚡ Smoke test | 60 s | Engine-runs-end-to-end check |
| 🩸 IP3 Ca²⁺ transient | 200 s | Reproduces Phase 1 transient (Dolan 2014 Fig. S2 IP3 curve) |
| 🛌 Resting | 300 s | Verify steady state at current resting IC |

**Removed from `jobs.py`:**

- The entire `_run_job` E. coli branch (called `runParca.py`, `runSim.py`,
  `analysisSingle.py` — all deleted from the repo)
- The `parca` phase from `PHASES` and `PHASE_DURATIONS` (job pipeline is now
  just `queued → simulating → analyzing → done|failed`)
- Job config schema simplified from
  `{variant, first_variant_index, last_variant_index, generations, init_sims, seed, toggles, description}`
  to `{length_sec, seed, description}`

**Other surface changes:**

- `tabs/runs.py` — Run Status table columns are now "Length / Seed" (was
  "Variant / Generations × Seeds")
- `tabs/inspect_data.py` — default listener prefers `CalciumTrace`, default
  column for `CalciumTrace` is `ca_cyt_nM` (was `Mass` / `dryMass`; both still
  available)
- `tabs/explore.py` — fixed stale "run analysisSingle.py" hint
- `results.py` — removed redundant `wildtype_*` directory special-case
- `app.py`, `runscripts/manual/webapp.py`, `assets/style.css` — `wcecoli_root`
  → `repo_root`; banner / docstring renames

**Old runs deleted:** `out/run-check/` (the only remaining stale run from the
E. coli era).

### Verification

| Check | Result |
|-------|--------|
| All webapp imports | OK |
| `create_app()` constructs | 14 callbacks registered, no duplicate-Output / missing-Input errors |
| 60-second platelet sim via `runPlateletSim.py` | OK; writes 7 listener directories (BulkMolecules, CalciumTrace, Environment, EvaluationTime, Main, Mass, UniqueMolecules) |
| `analysisPlatelet.py` on the run | OK; writes `calcium_trace.png` and `scaffold_summary.png` to `low_res_plots/` |
| `find_sim_dirs` / `find_variants` / `find_cells` / `find_listeners` / `find_columns` / `find_plot_images` | All discover the run correctly |
| `pytest models/platelet/tests/` | 17 passed |

### File index of this session's changes

| File | Change |
|------|--------|
| `CLAUDE.md` | Full rewrite for platelet-wcm context (committed `f83d9d8b`) |
| `wholecell/webapp/tabs/configure.py` | Strip E. coli; new platelet form + 3 presets |
| `wholecell/webapp/jobs.py` | Remove E. coli code path; simplify schema; drop `parca` phase |
| `wholecell/webapp/tabs/runs.py` | Length / Seed columns instead of Variant / Scale |
| `wholecell/webapp/tabs/inspect_data.py` | Default to CalciumTrace / ca_cyt_nM |
| `wholecell/webapp/tabs/explore.py` | Update stale runscript hint |
| `wholecell/webapp/results.py` | Drop `wildtype_*` directory special-case |
| `wholecell/webapp/app.py` | `wcecoli_root` → `repo_root`; param + docstring |
| `runscripts/manual/webapp.py` | "wcEcoli Web UI" → "Platelet WCM Web UI" |
| `wholecell/webapp/assets/style.css` | Banner comment |

### Where to pick up next session

Phase 2 IP3R recalibration is still the next live biology task — see the
"Next steps" section at the top of this lab book. The webapp is now ready for
launching short Phase 2 / Phase 3 calibration runs without command-line
gymnastics: `make run`, click a preset, watch results in Inspect Data /
Explore Plots.

---

*Branch:* `main` · *Last commit before this session:* `f83d9d8b` (CLAUDE.md
rewrite) · *Webapp cleanup uncommitted at write time*
