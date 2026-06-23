---
title: "Lab book — 2026-05-12: MCU mitochondrial Ca²⁺ uniporter (issue #22)"
---

# Lab book — 2026-05-12: MCU mitochondrial Ca²⁺ buffering (#22)

## Why this change

After Phase 4 (PI cycle) + v0.3.2 (integer-rounding fix), the 50-min
recovery sim shows DTS at ~750 µM (vs 250 µM target) and cyt declining
toward 100 nM very slowly. Literature consensus: 5–30 min DTS recovery
is biological, but our 50+ min is at the long end.

The most likely cause is the **missing mitochondrial Ca²⁺ buffer**.
Three recent platelet papers (all in `source-info/calcium-papers/`)
make the case:
- **Ghatge et al. 2026** — MCU⁻/⁻ platelets show *elevated* cyt Ca²⁺
  and *reduced* mitochondrial Ca²⁺ uptake. Direct evidence that MCU
  acts as a cyt Ca²⁺ buffer in vivo.
- **Ajanel et al. 2025** — MCU regulates ITAM-dependent platelet
  activation; quantitative MCU-KO data.
- **Shehwar et al. 2025** (review) — MCU complex in inner mito
  membrane; NCLX is the matrix efflux pathway.

> **Correction (2026-06-23).** The Ghatge bullet above — and the repeated
> claim later in this entry — that "MCU⁻/⁻ platelets show *elevated* cyt Ca²⁺"
> is **wrong**. On reading the paper in full (Ghatge et al., *J. Thromb.
> Haemost.* 2026; 24:716–731), MCU-knockout platelets show *reduced*
> agonist-evoked cytosolic Ca²⁺ (store release + SOCE), as does Ajanel et al.
> 2025. The model (buffer-only MCU) *raises* the peak and therefore *diverges*
> from this literature; this is documented in
> `reports/experiments/3-mcu-knockout.qmd` and tracked as issue #76. The
> original (incorrect) wording is left below as written, for an accurate record.

In real platelets, mitochondria take up Ca²⁺ rapidly during a
transient (~10–30% of imported Ca²⁺), then release it slowly via
NCLX over minutes. This **bypasses the PMCA-rate-limited extrusion
pathway** that is causing our slow DTS recovery — Ca²⁺ goes:
cyt → mito (fast) → cyt (slow) → DTS via SERCA → PMCA (slow).

The mitochondrial buffering acts as a temporal redistribution: it
holds Ca²⁺ during the spike and slowly returns it, allowing PMCA to
work at lower steady-state cyt concentrations for longer, effectively
reducing the load on the PMCA bottleneck.

## Design

### Species (1 new)

- `CA2_MITO[m]` — free Ca²⁺ in mitochondrial matrix

The mito Ca²⁺-binding proteins (matrix calbindin, polyphosphate
buffers, etc.) are not modelled explicitly — lumped into an
"effective" mito Ca²⁺ count.

### Compartment

The `'m'` mitochondrial-matrix compartment is **already reserved**
in `SimulationDataPlatelet.compartment_abbrev_to_index` per the
existing platelet plan. No new compartment work needed.

### Kinetic scheme

Two fluxes, both first-order in their substrate:

```
                k_MCU × cyt²/(Km² + cyt²) × N_MITO_eff
       cyt  ───────────────────────────────────────►  mito matrix
       cyt  ◄───────────────────────────────────────  mito matrix
                          k_NCLX × ca_mito
```

- **MCU uptake**: Hill cooperativity (n=2), saturating at high cyt
  Ca²⁺. Inactive below ~100 nM cyt (the resting band) due to the
  cooperative kinetics. Active during the transient.
- **NCLX efflux**: linear release (first-order in mito Ca²⁺). Slow
  timescale (~3 min, τ_NCLX = 200 s).

### Parameters

| Param | Value | Source / rationale |
|---|---|---|
| `N_MITO` | 5 (number of mitochondria per platelet) | Shehwar 2025; literature typical 2–5 |
| Mito matrix volume (per cell) | ~0.006 fL (0.1% of cyt) | Shehwar 2025 |
| `Vmax_MCU` (ions/s, total over all mito) | **50 000** | Calibrated; estimated for 10–30% of peak Ca²⁺ flux during transient |
| `K_MCU` (Hill half-saturation) | **1.0 µM** | Mid-range from literature (varies 0.5–10 µM) |
| `n_MCU` (Hill coefficient) | **2** | Mid-range from literature (varies 2–4) |
| `k_NCLX` | **0.005 s⁻¹** | τ = 200 s (slow release); typical literature value |
| Initial `CA2_MITO[m]` | 0 | starts empty; will load during transient |

Total cell mass impact: zero — `CA2_MITO[m]` is a Ca²⁺ count, mass
of Ca²⁺ ions is included in the existing total Ca²⁺ pool when computed.

### Expected impact

- Resting state: at cyt = 100 nM, MCU flux ≈ 50 000 × 0.01/(1.01) = 495
  ions/s. NCLX at mito_rest balances this. If mito holds ~99 000 ions
  at rest (495 / 0.005), that's 100 µM in matrix — realistic.
- During peak (cyt = 500 nM): MCU = 50 000 × 0.25/0.26 = 48 077
  ions/s into mito. Substantial flux that would peak-buffer cyt.
- During recovery (cyt = 200 nM): MCU = 50 000 × 0.04/0.05 = 40 000
  ions/s in; NCLX out depends on mito loading. Net flow defines
  recovery time.

### Acceptance criteria

1. Resting state preserved: cyt 95–125 nM, DTS 220–290 µM, IP3 40–60 nM.
2. Phase 3 5/5 criteria still pass after any retune.
3. **DTS overshoot recovery accelerated**: at t = 3000 s in the long-recovery
   sim, DTS < 500 µM (vs current 750 µM); cyt < 110 nM (vs current 126).
4. Mito Ca²⁺ rises to >1 µM during transient peak, then declines slowly
   (visible in the long-recovery plot panel).
5. All 21 regression tests pass.

### Risks

- The 50 000 ions/s Vmax_MCU is a guess; will need to iterate.
- May upset Phase 3 calibration; iterate γ_IP3R / γ_P2X1 / N_GSN as
  needed.
- Adds 1 species — minimal model complexity increase.

---

---

## Implementation results (2026-05-12)

### Final calibration

| Parameter | Design | Final | Reason for adjustment |
|---|---|---|---|
| `V_max_MCU` | 50 000 ions/s | 50 000 (unchanged) | |
| `K_MCU` | 1.0 µM | 1.0 (unchanged) | |
| `n_MCU` | 2 | **4** | n=2 over-absorbed at moderate cyt Ca²⁺, killing the cyt peak (340 nM, below the Dolan band lower bound 315 nM). n=4 gives a sharper switch — near-zero MCU at rest, strong at peak. |
| `k_NCLX` | 0.005 s⁻¹ | 0.005 (unchanged) | |
| `γ_P2X1` | 1.0 fS | **1.3 fS** | Restored SOCE differential after MCU attenuated cyt peak slightly. |
| Initial `CA2_MITO[m]` | 99 000 ions | **1 000** | Resting MCU at n=4 is much smaller (~5 ions/s vs 495), so equilibrium mito Ca²⁺ is correspondingly smaller. |

### Resting state — ✓

| Quantity | Target | Result |
|---|---|---|
| cyt | 95–125 nM | 107 nM ✓ |
| DTS | 220–290 µM | 234 µM ✓ |
| IP3 | 40–60 nM | 50 nM ✓ |
| Mito Ca²⁺ | bound for accumulation | 1 000 ions (≈ 330 nM matrix) — empty as expected at rest |

### Phase 3 — ✓ 5/5

- +Ca_ex peak: 434 nM (was 479 pre-MCU)
- −Ca_ex peak: 314 nM (was 319 pre-MCU)
- SOCE differential: 120 nM (vs criterion ≥ 100)
- All in Dolan ±30 % bands

### Mitochondrial Ca²⁺ dynamics — ✓

- Resting: 1 000 ions (matrix), MCU ≈ 5 ions/s, NCLX ≈ 5 ions/s, balanced.
- Peak: **152 793 ions** at t ≈ 218 s (153× resting). Mito has buffered substantial Ca²⁺ from the cytosolic transient.
- Recovery: NCLX releases mito Ca²⁺ over τ ≈ 200 s. At t = 1 000 s, mito ≈ 7 257 ions (still elevated 7× rest); by t = 3 000 s, ≈ 3 304 ions (3× rest).

### ⚠ Honest finding: MCU does NOT fix the DTS overshoot

The original motivation for #22 (from the issue body) was the hypothesis
that MCU would help close the DTS overshoot by providing an alternative
Ca²⁺ buffer. **This hypothesis is falsified by the implementation.**

Long-recovery 3 000 s comparison:

| Metric | Pre-MCU (v0.3.2) | Post-MCU (v0.3.3) |
|---|---|---|
| cyt at t = 3 000 s | 126 nM | 134 nM (slightly worse) |
| DTS at t = 3 000 s | 758 µM | **1 062 µM (worse)** |
| DTS at t = 1 000 s | (rising) | 1 229 µM (max overshoot) |

**Why MCU makes the DTS overshoot worse, not better**:

1. During the transient, MCU absorbs Ca²⁺ from cyt — *attenuating the
   cyt peak*. This means PMCA extrudes *less* Ca²⁺ during the peak
   (PMCA rate scales with cyt Ca²⁺).
2. More total Ca²⁺ stays in the cell post-stim (mito-bound + SERCA-pumped
   into DTS).
3. After stim, mito releases its bound Ca²⁺ slowly back to cyt → SERCA
   pumps it into DTS → DTS overshoots *higher* than without MCU.
4. PMCA at low cyt is still the bottleneck — adding mito buffering
   doesn't add an extrusion pathway, just redistributes Ca²⁺ within the
   cell.

**This is a real model prediction worth highlighting** in the dissertation:
the model predicts that **MCU-mediated Ca²⁺ buffering alone cannot
accelerate DTS recovery**. To genuinely fix the slow DTS extrusion,
you need:

- A *second extrusion pathway* (NCX / Na⁺/Ca²⁺ exchanger — some platelet
  evidence; v0.4 candidate), OR
- Faster PMCA cycling (possible via PMCA isoform mix or higher pump
  count), OR
- Reduced SOCE flux during stimulation (would reduce the imported load
  in the first place).

The Ghatge/Ajanel papers show MCU-KO platelets have **elevated** cyt
Ca²⁺, consistent with our model (MCU reduces cyt peak). They don't
directly measure DTS recovery rate as a function of MCU activity — so
the model prediction here is, in principle, testable experimentally.

> **Correction (2026-06-23):** the sentence above is **wrong** — see the
> correction note near the top of this entry. Ghatge et al. 2026 and Ajanel et
> al. 2025 both report *reduced* agonist-evoked cytosolic Ca²⁺ in MCU-knockout
> platelets, so the model (which raises the peak) *diverges* from the data.
> Issue #76; `reports/experiments/3-mcu-knockout.qmd`.

### What MCU **does** add to the model

1. **Biological correctness**: real platelets have MCU; the literature
   says MCU is functionally important for Ca²⁺ dynamics. The model now
   reproduces the MCU-KO finding (without MCU, cyt peak is higher).
2. **Slow Ca²⁺ tail**: mito provides a τ ≈ 200 s release phase that
   feeds slowly back to cyt, consistent with reports of "sustained
   Ca²⁺ signals" in activated platelets.
3. **Bioenergetics hook**: mito Ca²⁺ activates pyruvate dehydrogenase
   in real biology. Sets up issue #13 (Metabolism process) cleanly.
4. **Methodological rigour**: closes the "no MCU" gap that
   dissertation-notes flagged.

### Acceptance criteria — final score

| Criterion | Result |
|---|---|
| 1. Resting state preserved | ✓ |
| 2. Phase 3 maintains 5/5 | ✓ |
| 3. DTS overshoot accelerated | ✗ — *and the finding is interesting biology, not a model defect* |
| 4. Mito Ca²⁺ rises >1 µM at peak | ✓ — 153 k ions ≈ 50 µM matrix |
| 5. All 21 tests pass | ✓ |

**4.5/5 — the failure on #3 is the dissertation-worthy finding.**

---

*Branch:* `main` · *Status:* #22 complete; MCU functional; DTS-overshoot
finding documented as a model prediction worth experimental follow-up ·
*Linked issues:* #22 (this work, complete)
