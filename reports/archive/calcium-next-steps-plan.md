# Calcium Model — Next Steps Plan

**Branch:** `platelet`
**Date:** 2026-04-29
**Dissertation deadline:** mid-July 2026 (~10 weeks)
**Current state:** Equations correct (Sneyd, Po⁴, Nernst, MWC), but resting state drifts and IP3-stimulated transient peaks at ~6 µM with late-time runaway. See `reports/platelet-calcium-calibration.md` §4.

---

## 1. Goal

Reproduce Dolan & Diamond 2014 Fig. 4 in the platelet WCM:

- Resting `[Ca²⁺]_cyt` = 100 ± 10 nM, stable for 300 s with no stimulus
- IP3-stimulated peak `[Ca²⁺]_cyt` = 200–500 nM within 15–20 s
- Partial DTS depletion (30–70 % of resting)
- SOCE-sustained plateau when extracellular Ca²⁺ is present
- Faster decay / no plateau when extracellular Ca²⁺ is absent (EDTA)

This is the primary calibration target named in the design doc §7.

---

## 2. Why the current model misses these targets

Three v0.2 simplifications cause the runaway and the resting drift:

1. **PMCA is basal-only (Caride 2007 Table 3, no CaM).** k_cat = 5.5 s⁻¹ caps PMCA throughput at ~4000 ions/s. Once SOCE turns on the cell can't extrude fast enough, so cytosolic Ca²⁺ accumulates without bound.
2. **No cytosolic Ca²⁺ buffering.** All free Ca²⁺ is unbuffered. In reality ~95 % is bound to calmodulin / calbindin / etc., which both slows the rise and provides a memory mechanism for inactivation.
3. **The Dolan Table S1 IC is not a steady state of our v0.2 ODE.** Dolan's IC was filtered against four homeostatic constraints with their full ODE (CaM-coupled PMCA, scanned (KM, n)). Our subset of equations has a different fixed point.

Phase 1 below addresses (1) and partially (2) by introducing CaM. Phase 2 addresses (3) by re-deriving the IC against our actual ODE.

---

## 3. Phased plan

### Phase 1 — CaM + Caride 5-state CaM-coupled PMCA

**Effort:** 1–2 working days
**Source:** Caride et al. 2007 (in `source-info/calcium-papers/`); Dolan 2014 Table S1 for CaM IC counts
**Closes:** acceptance criterion on transient peak (issue #25); a deferred item from the calibration report §5

**Scope:**
- Add CaM as a tracked species. Dolan Table S1 lists 20,481 total CaM with sub-states `CaM_free=20465`, `Ca₂·CaM=15`, `Ca₄·CaM=1`, plus CaM bound to PMCA (negligible at rest).
- Implement Caride 2007 Table 3 full 5-state PMCA: PMCA, PMCA·Ca, PMCA·CaM, PMCA·CaM·Ca, PMCA(CaM·Ca₄)·Ca. This adds ~3 sub-states beyond the current 2.
- Add Ca²⁺ ↔ CaM cooperative binding kinetics (4 Ca²⁺ per CaM, two N-lobe and two C-lobe sites with different affinities — Caride or Sveshnikova references).
- Update `internal_state.py`, `molecules.tsv`, and the SimulationDataPlatelet sub-mass routing.
- Update `_ode_rhs` in `calcium_signalling.py` with the new PMCA cycle and CaM Ca²⁺-binding reactions.

**Two physiological consequences (the reason this fixes things):**

1. **PMCA Vmax shoots up during the transient.** As Ca²⁺ rises, more CaM binds Ca²⁺, more Ca₄·CaM activates PMCA, and apparent Km drops from ~5 µM to ~0.5 µM with k_cat doubling. Peak PMCA throughput goes from ~4000 ions/s (basal) to ~80,000 ions/s (saturated CaM). That is enough to clamp the cytosolic peak at biological values.
2. **CaM is a Ca²⁺ buffer.** With ~20,000 binding sites, CaM soaks up most of the free Ca²⁺ during the rise, slowing the apparent rate of cyt accumulation.

**Validation:** With Phase 1 alone, run a 200 s sim under IP3 forcing. Pass criterion: peak `[Ca²⁺]_cyt` in 200–800 nM range, no late-time runaway, DTS partially depletes (not to zero).

**Risk:** Caride's full kinetic scheme introduces several new rate constants. If any are wrong by an order of magnitude the system could over- or under-clamp. Mitigated by deriving constants directly from Caride Table 3 and verifying the resting state numerically.

---

### Phase 2 — Re-derive a true steady-state initial condition

**Effort:** half a day
**Source:** Output of the Phase 1 model
**Closes:** acceptance criterion on resting stability (issue #25)

**Scope:**
- Start from Dolan Table S1 IC, set IP3 = 50 nM constant, integrate the ODE for 600 s.
- The system will converge to its own steady state (different from Dolan's because we use a subset of their equations and a slightly different parameterisation).
- Save the converged state as the new IC in `internal_state.py` and `molecules.tsv`.
- Document the divergence from Dolan Table S1 — it should be small (< 20 %) on most species if Phase 1 worked.

**Validation:** A 300 s rest simulation (no IP3 forcing) holds `[Ca²⁺]_cyt` at 100 ± 10 nM and `[Ca²⁺]_dts` at 200–300 µM. Matches design-doc §7.1 pass criteria.

---

### Phase 3 — Dolan Fig. 4 validation (the dissertation result)

**Effort:** 1–2 working days
**Source:** Dolan 2014 Fig. 4

**Scope:**
- Run three matched simulations:
  - **Standard:** full IP3 forcing, extracellular Ca²⁺ = 1.2 mM (Dolan Fig 4A/B baseline)
  - **EDTA:** full IP3 forcing, `Ca_ex = 0` (Dolan Fig 4C control)
  - **Resting:** no IP3 forcing, 300 s (resting stability)
- Update `models/platelet/analysis/single/calcium_trace.py` to overlay Dolan-digitised reference traces (currently uses a schematic — replace with WebPlotDigitizer-extracted curves).
- Produce a single figure suitable for dissertation inclusion.

**Validation:** Quantitative match against Dolan Fig 4:
- Peak ratio (standard / EDTA) within 30 %
- DTS depletion fraction within 20 %
- Plateau height (standard) within 50 %
- Decay τ within a factor of 2

**Outputs for dissertation:** the validation figure plus a brief table of peak/plateau/decay numbers vs Dolan.

---

### Phase 4 — Optional: re-scan (KM, n) for puncta entry

**Effort:** 1–2 days
**Trigger:** only if Phase 3 validation criteria fail by a wide margin
**Source:** Dolan 2014 Fig 3 / SI Methods (clustering analysis)

**Scope:**
- Vary KM ∈ [0.1, 2.0] µM and n ∈ [2, 8] over a coarse grid
- For each (KM, n), run the resting + standard transient and score against Dolan Fig 4
- Pick the (KM, n) on the Pareto front
- Document the chosen values and explain why they differ from Dolan's range (if they do)

This is the platelet-specific equivalent of Dolan's "12-dimensional sampling" reduced to two dimensions because we now use Dolan's other parameters as fixed.

---

### Phase 5 — Stretch: P2Y1 upstream cascade (v0.3)

**Effort:** 1–2 weeks
**Source:** Purvis 2008 Table 1 receptor module; Sveshnikova 2025 for the stochastic bottleneck
**Closes:** issues #32 (Receptor signalling), #33 (P2Y12), #34 (LocalEnvironment timeline)

**Scope:**
- Implement P2Y1 + Gq + PLCβ as a `ReceptorSignalling` process producing IP3 endogenously
- IP3 is now a real BulkMolecule produced/consumed by mass-action; the IP3 forcing curve goes away
- Add an environment-driven ADP timeline so the simulation responds to agonist rather than a hard-coded curve
- Optionally tau-leap the PLC-Gq complex (~1 molecule, the Sveshnikova 2025 stochastic bottleneck)

**Validation:** ADP dose-response curve resembling Purvis 2008 Fig 5.

This is genuinely v0.3 territory. Realistic only if Phases 1–3 land well under 6 weeks — pencilled in but not committed.

---

## 4. Timeline against the dissertation

| Week | Phase | Deliverable |
|------|-------|-------------|
| 1 | 1 (CaM-PMCA) | Code + Phase 1 validation |
| 2 | 1 → 2 | Phase 1 finished, Phase 2 IC derived |
| 3 | 3 | Dolan Fig 4 validation figure |
| 4 | 3 → 4? | Validation passes; or re-scan if not |
| 5–6 | Writing | Dissertation calcium chapter draft |
| 7–8 | 5 (stretch) | P2Y1 if time permits |
| 9–10 | Writing | Final dissertation revisions |

Phases 1–3 land the headline result. Phase 5 is a "v0.3" extension; the dissertation can describe it as future work without implementation.

---

## 5. Risks and contingencies

**Risk 1 — Caride 2007 doesn't give all the rate constants we need.**
Likelihood: low. Caride 2007 Table 3 has the full set; we already use the basal subset. Mitigation: cross-check against Dolan 2014 SI which restates the Caride kinetic scheme.

**Risk 2 — Phase 1 fixes the runaway but the peak is still wrong.**
Likelihood: medium. Could be off by 2–3× because of CaM Ca²⁺ binding rates we haven't fully verified. Mitigation: Phase 4 (re-scan) catches this; or scan k_on/k_off for CaM-Ca²⁺ binding within literature ranges.

**Risk 3 — Phase 5 (P2Y1) reveals new architectural problems.**
Likelihood: medium-high. The P2Y1 module is large (~30 species, dozens of rate constants) and introduces stochasticity. Mitigation: kept as stretch; dissertation can stand on Phases 1–3.

**Risk 4 — IDE/tooling/refactor distractions.**
Likelihood: medium. Mitigation: keep `wholecell/webapp/` work separate from calcium work in commits and on the task list.

---

## 6. Issues to create / update

| GitHub issue | Action |
|--------------|--------|
| #25 (CalciumDynamics process) | Stays open; gets a comment when each phase passes its acceptance criterion |
| #24 (data + dataclass) | Stays open; close when Phase 2 IC is in place and CaM is added to inventory |
| **New: "Model: CaM + Caride 5-state CaM-coupled PMCA"** | Phase 1 deliverable — create before starting |
| **New: "Reconstruction: re-derive resting IC from converged ODE"** | Phase 2 deliverable |
| **New: "Validation: Dolan 2014 Fig 4 reference figure"** | Phase 3 deliverable |
| #32, #33, #34 (receptor signalling) | Stays open; only relevant for Phase 5 |

Recommendation: create the three new issues before starting Phase 1, so each phase has its own acceptance-criteria checklist.

---

*This plan supersedes "What still needs doing" in `reports/platelet-calcium-calibration.md` §5.
Lab-book entry for the equation-restoration session: `reports/lab-book-2026-04-29-calcium-equations.md`.*
