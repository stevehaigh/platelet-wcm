---
title: "Important notes for the dissertation write-up"
status: living document
last-updated: 2026-05-11
---

# Important notes for the dissertation write-up

A curated, long-lived list of model assumptions, calibration choices, and
biological gaps that must be acknowledged in the dissertation. Each entry
records:

- **What** the assumption / choice is
- **Why** it matters (what would change in v0.3+ if corrected)
- **State** of the evidence
- **Where** the full diagnosis lives (lab book reference)

Add new items at the bottom of the relevant section; do not re-order
historical entries — the dissertation needs a stable list of cited points.

---

## 1. Ca²⁺ buffering — cytosol

### 1.1 Cytosolic buffering ratio is at the low end of biology

- **State (2026-05-11)**: only calmodulin (CaM) modelled as a cytosolic
  Ca²⁺ buffer. At rest: 361 free cyt ions vs 1 276 CaM-bound → **78%
  buffered, ratio bound:free = 3.5:1**.
- **Literature (non-muscle cells)**: total cyt buffering ratio is typically
  **50:1 to 100:1** (98–99% bound). Major non-CaM cytosolic Ca²⁺ binders
  in platelets that we do not model: **gelsolin** (~250 000 copies, multi-
  site EF-hand-like binding, Kd ~0.1–1 µM), **annexins**, **ATP** (3 mM ×
  Kd_Ca-ATP ~1 mM gives non-trivial Ca-ATP), various Ca²⁺-regulated kinases
  and phosphatases.
- **Effect on dissertation claims**: the Phase 3 peak Ca²⁺ values
  (393 nM with Ca²⁺_ex, 325 nM without) are calibrated against Dolan 2014,
  who used a similarly under-buffered cytosol. The model and Dolan are
  internally consistent. But **the absolute peak heights would be lower
  in a fully-buffered cytosol** — typically by 2–10× depending on the
  buffer's koff. Transient kinetics are also affected (lower koff = slower
  fall after peak).
- **Mitigation in this work (2026-05-11)**: added a coarse-grained
  "gelsolin proxy" species (`GSN_free[c]`, `GSN_Ca[c]`) with 1:1 binding
  kinetics — but only as a **scaffold at N_GSN = 5 000** (50× below
  biological). Full-biology N_GSN = 250 000 was *tested* and found to
  crash Phase 3 peak heights from ~390 nM to ~130 nM, because the
  existing IP3R / SERCA fluxes are themselves calibrated against the
  Dolan 2014 under-buffered model. Closing the gap therefore requires
  a coupled re-tune (IP3 forcing magnitude or γ_IP3R) and is left for
  v0.3+. See `lab-book-2026-05-11-dyk-ip3r-design.md §Cytosolic
  buffering pass` for the full analysis.
- **What v0.3 should do**:
  1. Scale `N_GSN` up to ~250 000 with multi-site binding (or split
     into explicit gelsolin / annexin / Ca-ATP species).
  2. Re-tune IP3 forcing amplitude and/or γ_IP3R to restore the
     Dolan Fig. 4 measured peak heights with the new (correct)
     buffer load.
  3. Benchmark resting buffer ratio against Sage & Rink 1985.

---

## 2. Ca²⁺ buffering — dense tubular system (DTS)

### 2.1 DTS luminal buffering is severely underpowered until Phase 2

- **State (2026-05-11)**: only the STIM1 EF-hand modelled (3 805 STIM1·Ca
  out of 4 265 STIM1 at rest). DTS buffering ratio: **9% bound, ratio
  bound:free ≈ 0.1**.
- **Literature**: ER / SR luminal stores are **95–99% buffered**, primarily
  by calreticulin (CALR), calsequestrin (in muscle), HSP90B1 (GRP94),
  CALU, and other CREC-family proteins. Free luminal Ca²⁺ ≈ 100–500 µM
  with total luminal Ca²⁺ ≈ 5–25 mM.
- **Consequence**: the long-time resting state in our ODE drifts (cyt rises
  to 200+ nM, DTS overfills to >1 mM in a 6 000 s run) because SERCA
  pumps Ca²⁺ into a DTS with effectively no luminal buffer to absorb it.
  The DTS only stops rising once it reaches the SERCA thermodynamic
  reversal point. **The initial conditions (t = 0) and the short Phase 3
  transients (200 s) are the only physically meaningful states in the
  current model.**
- **Dissertation framing**: must explicitly state that the v0.2.5 freeze
  is a "transient-validated" model — Phase 3 stimulus / response dynamics
  are biologically calibrated, but the model does not have a stable
  biological resting fixed point.
- **Resolution**: Phase 2 (issue #28) adds CALR (508 100 Ca²⁺-binding
  sites at Kd ≈ 1 mM) and is the dominant fix. v0.3+ should add HSP90B1,
  CALU per #25.
- **Reference**: `lab-book-2026-05-11-dyk-ip3r-design.md`, this document
  §3.1 (γ_IP3R / SERCA coupling).

---

## 3. Flux-rate calibration

### 3.1 γ_IP3R = 0.35 pS is coupled to SERCA rate constants

- **Decision (2026-05-11, Phase 4 / commit `1699ac1f`)**: γ_IP3R reduced
  from 10 pS (Zschauer 1988 bilayer) to 0.35 pS (calibrated to balance
  SERCA at the Dolan resting state).
- **Derivation**: at cyt = 100 nM, DTS = 250 µM, the 6-state SERCA cycle
  steady-state flux is 112 570 ions/s (analytical solution of the linear
  system using Purvis 2008 / Dode 2002 rate constants). γ_required =
  112 570 / (N × Po × driving × NA / zF) = 0.344 pS, rounded to 0.35.
- **Why this matters**: **γ_IP3R is not an independently measured value**
  in our model — it is the value that *balances the chosen SERCA rate
  constants* at the chosen resting state. If SERCA constants change
  (see §3.2), γ_IP3R must be re-derived.
- **Biological plausibility**: 0.35 pS sits within the cellular IP3R
  effective Ca²⁺ conductance range reported by Bezprozvanny 1991 and
  Mak & Foskett 1997 (~0.05–0.5 pS under physiological conditions). The
  10 pS bilayer value is not transferable because Zschauer used
  symmetric high Ca²⁺, where K⁺ contributes negligibly to current.
- **Dissertation framing**: cite as a *calibration anchor*, not a
  measured parameter. Disclose the SERCA coupling explicitly.

### 3.2 SERCA cycle flux is probably 2–5× too high at rest

- **Current model**: Purvis 2008 / Dode 2002 rate constants give SERCA
  cycle rate of **4.7 cycles/s per pump at cyt = 100 nM**, or **112 570
  Ca²⁺ ions/s total** for 11 892 pumps.
- **Literature SERCA3b kinetics**:
  - Vmax ≈ 30–50 cycles/s at saturating Ca²⁺ (Inesi 1985; Nishi 1992)
  - Km(Ca²⁺) ≈ 0.7–1.1 µM (Dode 2002 — SERCA3 is *less* Ca²⁺-sensitive
    than SERCA2a/b, by design)
  - At cyt = 100 nM with n = 2 Hill: v/Vmax ≈ 2% → ~1 cycle/s per pump
  - Predicted total flux: **~23 800 ions/s**
- **Inherited from Purvis 2008**: the rate constants we use are Purvis's,
  who took them from Dode's protein expression studies. But Purvis's
  k_bind_f = 1 000 µM⁻²·s⁻¹ implies a faster pump than Dode's measured
  Vmax / Km values predict.
- **Implication**: if v0.3 re-derives SERCA constants from primary sources,
  γ_IP3R will drop to ~0.07–0.10 pS (and PMCA / PM-leak balances will
  shift). The relative dynamics (Phase 3 transients) should be largely
  preserved because they are dominated by ratio of fluxes, not absolutes.
- **Dissertation framing**: cite as a known calibration question
  inherited from Purvis 2008, scoped for v0.3+ revision. Phase 3
  validation against Dolan 2014 demonstrates that the *relative* SERCA /
  IP3R balance is correct for the platelet stimulus regime.

### 3.3 PM Ca²⁺ leak (75 ions/s) is at the upper end of biological estimates

- Sage & Rink 1985 measured platelet PM Ca²⁺ entry at ~10–40 ions/s at
  rest. Our 75 ions/s is ~2× high. Minor numerically (cf. 100 k ions/s
  IP3R / SERCA), but worth noting for completeness.

---

## 4. Methodological choices

### 4.1 IP3R Po formula: m∞⁴ × h (not m∞³ × h)

- **Choice**: Po = m∞⁴ × h for the deYoung-Keizer / Li-Rinzel IP3R.
- **Alternative in literature**: Li-Rinzel 1994 original used Po = m∞³ × h
  (three-site cooperativity). Dolan 2014 used a Po⁴ tetrameric form.
- **Our rationale**: m∞⁴ × h preserves the four-fold cooperativity
  convention used by Dolan, against whose data we calibrate Phase 3.
- **Sensitivity**: at cyt = 100 nM, IP3 = 50 nM, m∞ = 0.1523:
  - m∞³ × h = 0.0032 × 0.913 = 2.94×10⁻³
  - m∞⁴ × h = 0.000493 × 0.913 = 4.92×10⁻⁴
  - Po is 6× higher under m∞³ × h. γ_IP3R calibration would scale
    inversely (~0.06 pS instead of 0.35 pS).

### 4.2 SERCA initial conditions: 6-state cycle vs 2-state binding equilibrium

- **Pre-2026-05-11**: SERCA initial conditions used the partial 2-state
  E1 ↔ E1Ca equilibrium (E1Ca/E1 = k_bind_f·cyt²/k_bind_r = 1.0). This
  ignored the fast phosphorylation drain (k_phos_f = 700 s⁻¹ >> k_bind_r
  = 10 s⁻¹).
- **Fixed 2026-05-11 (commit `1699ac1f`)**: now uses full 6-state quasi-
  steady-state populations (E1Ca/E1 = k_bind_f·cyt²/(k_bind_r + k_phos_f)
  = 0.0141, so E1Ca = 81 vs the old 2 963).
- **Effect**: eliminates a spurious 2 M event/s phosphorylation burst at
  t = 0 that previously drained cytosolic Ca²⁺ to <5 nM and trapped the
  system at the low-Ca²⁺ attractor below the d₅ activation threshold.

---

## 5. Inherited assumptions

### 5.1 Compartment volume = 6 fL (cytoplasm)

- **Assumption**: platelet cytoplasm = 6 fL (DTS = 4.3% by volume = 0.26 fL).
- **Reality**: total platelet volume is 6–10 fL; cytoplasm (excluding DTS,
  mitochondria, granules) is roughly 4–7 fL with significant inter-individual
  variation.
- **Sensitivity**: all concentration → count conversions scale with volume.
  A 6 fL → 10 fL change would reduce all concentrations by 40%, change
  the IP3R/SERCA balance, and require recalibration of γ_IP3R.
- **Dissertation framing**: cite as a fixed assumption per Burkhart 2012
  / Dolan 2014; flag sensitivity in the limitations section.

### 5.2 SERCA isoform = SERCA3b (ATP2A3)

- Burkhart 2012 reports both SERCA2b and SERCA3 in platelets. We model
  only SERCA3b. The two isoforms have different Ca²⁺ affinities (SERCA2b
  Km ~0.4 µM; SERCA3b Km ~0.7–1.1 µM). Mixed-isoform model is a v0.3+
  candidate.

### 5.3 IP3R isoform = ITPR2

- Burkhart 2012 / Dolan 2014 convention: 1 328 ITPR2 copies. We treat all
  IP3R as ITPR2. Real platelets express all three isoforms (ITPR1, 2, 3)
  with ITPR2 dominant. The three have different IP3 sensitivities; mixed-
  isoform model is v0.3+.

---

## 6. Open questions for the writeup

- How to present the v0.2.5 resting-state runaway honestly without
  undermining the Phase 3 validation result. Suggested framing:
  *"transient-validated; long-time stability awaits the DTS buffer."*
- Whether to present the SERCA flux question (§3.2) as a known
  limitation or to attempt a v0.3-style re-derivation before the
  freeze.
- How much detail on γ_IP3R derivation belongs in the main text
  versus the appendix.
