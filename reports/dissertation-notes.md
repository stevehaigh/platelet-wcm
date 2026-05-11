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
- **Action taken (2026-05-11, Phase 2 retune)**: added the coarse-grained
  gelsolin proxy at N_GSN = 800 000 sites (~160 000 gelsolin × 5 effective
  Ca²⁺-binding sites; within Burkhart 2012 / Yin & Stossel 1979 range),
  Kd = 1 µM, k_off = 100 s⁻¹. **Coupled with CALR addition and a 2× IP3R /
  SERCA flux reduction** to keep Phase 3 peaks in the Dolan band — see
  `lab-book-2026-05-11-dyk-ip3r-design.md §Phase 2`.
- **Resulting buffering ratio is ~200:1 at rest** — *higher* than the
  Sage & Rink ~50:1 literature value. This is the calibration penalty
  for retaining the Dolan-inherited IP3R flux levels. A v0.3+ retune of
  the IP3R rate constants (separate from the SERCA-rate question in
  §3.2) would let N_GSN drop closer to ~200 000–300 000 with κ ≈ 50.
- **What v0.3 should do**:
  1. Split the coarse-grained buffer into explicit gelsolin (1–2 high-
     affinity EF-hand-like sites), annexins, and Ca-ATP equilibrium.
  2. Re-derive IP3R rate constants from primary sources (not Dolan's
     fitted values) and rebalance γ_IP3R + N_GSN to match the
     ~50:1 resting buffering ratio measured by Sage & Rink.
  3. Validate against post-peak Ca²⁺ decay kinetics (sensitive to
     buffer k_off in a way the peak isn't).

---

## 2. Ca²⁺ buffering — dense tubular system (DTS)

### 2.1 DTS luminal buffering — closed (Phase 2 / #28)

- **State (after 2026-05-11 Phase 2)**: calreticulin (CALR) added with two
  binding modes:
  - **C-domain** (low affinity, high capacity): 508 100 sites at Kd = 1 mM,
    k_off = 1 000 s⁻¹ (fast equilibrium). ~102 k Ca²⁺ bound at rest.
  - **P-domain** (high affinity, slow release): 20 324 sites at Kd = 1 µM,
    k_off = 1 s⁻¹. ~20 k Ca²⁺ bound at rest (~99 % saturated).
- **DTS buffering ratio**: from ~9 % bound → **~73 % bound** at rest.
  Long-time resting state now stable at cyt = 109 nM, DTS = 264 µM
  (vs the previous runaway to cyt = 200 nM, DTS > 1 mM).
- **Open**: the DTS still drops to ~0 µM during the IP3 transient because
  IP3R peak flux (~6.5 M ions/s after the Phase 2 retune) drains the
  buffered DTS in ~20 ms. Real biology likely retains 50–100 µM free
  DTS during stimulation. Closing this gap requires either further IP3R
  rate-constant work (§3.2) or additional DTS buffers (HSP90B1, CALU
  per #25; CaM is *not* a DTS buffer — see §2.2).
- **Reference**: `lab-book-2026-05-11-dyk-ip3r-design.md §Phase 2`.

### 2.2 Calmodulin is cytosolic — the DTS uses a distinct buffer system

A point of confusion worth pre-empting in the write-up: CaM is *only*
cytosolic / membrane-associated. It is synthesised on free ribosomes and
has no ER/SR targeting or retention signal (no signal peptide, no KDEL).
Literature: Cyert 2001, Berridge et al. 2003 *Nat Rev Mol Cell Biol*,
Chin & Means 2000 *Trends Cell Biol*.

The DTS lumen uses an entirely distinct set of Ca²⁺-binding proteins:

| Protein | Function | In our model? |
|---|---|---|
| **CALR** (calreticulin) | Dominant ER/SR luminal Ca²⁺ buffer; high-capacity C-domain + high-affinity P-domain | ✓ (Phase 2) |
| HSP90B1 (GRP94) | Chaperone with ~20 Ca²⁺ sites/molecule, mM affinity | ✗ (v0.3+, #25) |
| CALU (calumenin) | CREC-family small acidic Ca²⁺ binder | ✗ (v0.3+, #25) |
| RCN1, RCN2 (reticulocalbins) | CREC-family, multi-EF-hand | ✗ (v0.3+, #25) |
| ERp44 | Small acidic, Ca²⁺ binding, redox-regulated | ✗ |
| Calsequestrin (CASQ) | Dominant SR buffer in *muscle* — platelets don't express it | N/A |

Our model includes the dominant CALR component, leaving the smaller CREC-
family buffers for v0.3+.

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

## 6. Morphology and spatial assumptions

### 6.1 Compartments treated as well-mixed

- **Our assumption**: cyt and DTS are each a single well-mixed volume
  with uniform [Ca²⁺]. ODE-only, no diffusion.
- **Reality**: Ca²⁺ microdomains exist around IP3R clusters (puff sites)
  and at PM-DTS junctions where STIM1-Orai1 couple. Local [Ca²⁺] near
  open IP3R clusters can be 10–100× the bulk cytosolic concentration
  before equilibration.
- **Mitigation**: in a 6 fL volume with free Ca²⁺ diffusion D ≈ 200 µm²/s,
  the mixing timescale is sub-ms (L²/D ≈ 0.02 ms for a 2 µm cell). So
  the well-mixed approximation is probably fine for *bulk* dynamics.
  Microdomain effects matter for fast Ca²⁺-activated processes
  (e.g. PMCA's CaM activation rate may be under-estimated because the
  *real* local Ca²⁺ near membrane-localised CaM is higher than bulk).

### 6.2 Surface-Connected Canalicular System (SCS) is not modelled separately

- **The SCS** is a platelet-specific invaginated PM network that
  penetrates deep into the cell, continuous with the extracellular
  space. It effectively **doubles to triples the PM surface area**
  for ion entry / extrusion vs the naive sphere-surface estimate.
- **Effect on our model**: all PM-localised fluxes (PMCA, Orai1 / SOCE,
  PM_LEAK) are calibrated *as effective bulk rates* against Dolan's
  data, so the SCS is implicitly absorbed into the rate constants. But
  it makes our γ_SOC, J_PM_LEAK, k_PMCA values cell-level rates, not
  per-µm² fluxes — anyone trying to compare to PM patch-clamp data
  needs to scale accordingly.

### 6.3 IP3R clustering and "puff" dynamics

- IP3Rs cluster in real cells (~10–100 channels per cluster). Each
  cluster fires stochastically (Ca²⁺ puffs).
- Our well-mixed continuous model treats all 1 328 IP3R as independent
  with identical Po — effectively assumes population-averaged behaviour.
- Probably fine for the *macroscopic* peak heights we validate against,
  but the early transient kinetics may be off (real puff onset is
  faster locally; our population-averaged onset is smoother).

### 6.4 Volume parameters

| Parameter | Our value | Reality |
|---|---|---|
| Cytoplasm volume | 6 fL | 4–7 fL (MPV-dependent; high inter-individual variation) |
| DTS volume fraction | 4.3 % of cell | 4–13 % (Dolan low end; some EM gives higher) |
| Total platelet volume | implicit | 6–10 fL |

All concentrations scale with these volumes; sensitivity-check any
flux estimate that crosses a biological band.

---

## 7. Missing channels / pathways

### 7.1 P2X1 — the biggest gap

**P2X1 is the dominant fast Ca²⁺ entry pathway in activated platelets**
and is not in our model. It is an extracellular-ATP-gated cation channel
that opens within milliseconds of platelet activation, contributing the
*first* Ca²⁺ spike — before IP3R and well before SOCE.

| Property | Value |
|---|---|
| Activation timescale | <10 ms after ATP exposure |
| Ca²⁺ permeability | P_Ca/P_Na ~10 |
| Desensitisation | Fast (~100 ms) |
| Source in vivo | Released ATP from dense granules; autocrine activation |

**Why this matters for the SOCE differential criterion (§ Phase 3)**:
the real platelet +Ca_ex vs −Ca_ex peak difference may be driven
substantially by P2X1, not just by Orai1/SOCE. Our SOCE differential ≈ 0
result may be telling us that the missing pathway is P2X1, not slower
STIM1 dimerisation. Worth flagging as a v0.2.6 candidate.

### 7.2 Dense granule Ca²⁺ store

Platelet dense granules store Ca²⁺ at very high concentrations (total
[Ca²⁺] in the 50–100 mM range, mostly complexed with pyrophosphate and
polyphosphate). NAADP / two-pore-channel (TPC) signalling releases this
during activation. Our model has dense granules as a mass species
(`CA2_DG[dg]` would be the natural species; not currently present in
`internal_state.py`) but no flux pathway. Add as a third Ca²⁺
compartment for v0.3.

### 7.3 Other DTS / PM channels not in model

| Channel | Why not in model | Priority |
|---|---|---|
| **RyR2** | Existence in platelets is contested; Dolan ignores it; some Lopez et al. evidence | v0.3+ |
| **TRPC1/4/6** | Tethering and store-operated entry; partially covered by lumped SOCE | v0.3+ |
| **SPCA1** (ATP2C1) | Pumps Ca²⁺ into Golgi/secretory granules; ER ≠ Golgi but secretory pathway is biologically active | v0.3+ |
| **MCU** (mitochondrial uniporter) | Already issue #22 | v0.3+ |
| **NCX** (Na⁺/Ca²⁺ exchanger) | Some evidence in platelets, contested | v0.3+ |

### 7.4 Receptor-PLC-IP3 upstream of IP3R

We force IP3 directly using the Dolan Fig S2 fit. Real upstream:

GPCR (thrombin / collagen / ADP receptors)
  → Gαq
  → PLCβ
  → cleaves PIP2 → IP3 + DAG

The PIP2 pool is finite, PLC has its own kinetics, and IP3 is degraded
by IP3-3-kinase and IP3-5-phosphatase. Forcing IP3 directly:
- Misses the kinetic shape of the rise (we have τ_rise = 3 s; real may
  be slower or have a shoulder)
- Misses the finite-pool / desensitisation effects
- Decouples Ca²⁺ release from upstream receptor signalling

v0.3 receptor-signalling work (#9, #10) would close this.

### 7.5 Cytoskeletal coupling (gelsolin's dual role)

Our "GSN" species treats gelsolin as a passive Ca²⁺ buffer. Real
gelsolin is **dual-purpose**: it's a Ca²⁺-activated actin-severing
protein, and its Ca²⁺-binding state drives major cytoskeletal
rearrangement during platelet activation. The same Ca²⁺ that we count
as "GSN-bound" in the buffer accounting is mechanistically what
*activates* gelsolin to sever actin filaments — a major activation
endpoint. For a Ca²⁺-only model, the buffer aspect is correct; for any
cytoskeletal-output model, the GSN species would need a state-machine
representation.

---

## 8. Open questions for the writeup

- How to present the post-Phase-2 model honestly: it has biologically
  realistic resting state and Phase-3 peak heights, but the SOCE
  differential is missing and the DTS empties more during transient
  than real biology. Two possible framings:
  - *"v0.2.6 captures the dominant calcium pathway with biology-grade
    buffering; SOCE differential limitation traces to missing fast
    Ca²⁺ entry (P2X1, §7.1)."*
  - *"transient-validated for peak heights and resting state; full
    transient shape calibration awaits v0.3."*
- Whether to present the SERCA flux question (§3.2) as a known
  limitation or to attempt a v0.3-style re-derivation before the
  freeze.
- How much detail on γ_IP3R derivation belongs in the main text
  versus the appendix.
- Whether to add P2X1 to v0.2.7 (small commit, big biological
  improvement) or defer to v0.3 receptor-signalling work.
