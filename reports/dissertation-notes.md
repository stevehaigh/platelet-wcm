---
title: "Important notes for the dissertation write-up"
status: living document
last-updated: 2026-05-12
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

### 2.1 DTS luminal buffering — expanded (Phase 3 / #25, 2026-05-12)

- **State (after 2026-05-12 Phase 3)**: full multi-buffer DTS:
  - CALR C-domain (low aff, fast): 508 100 sites, Kd = 1 mM, ~102 k bound
  - CALR P-domain (high aff, slow): 20 324 sites, Kd = 1 µM, ~20 k bound
  - **HSP90B1 medium-aff (slow)**: 40 000 sites, Kd = 2 µM, k_off = 1 s⁻¹, ~40 k bound
  - **HSP90B1 low-aff (fast)**: 110 000 sites, Kd = 600 µM, ~32 k bound
  - **BiP / HSPA5**: 50 000 sites, Kd = 500 µM, ~17 k bound
  - **CREC pool** (CALU + RCN1 + RCN2): 60 000 sites, Kd = 1 mM, ~12 k bound
- **DTS buffering ratio**: ~73 % bound (Phase 2) → **~85 % bound** at rest.
  Total DTS Ca²⁺ at rest: ~165 k ions (Phase 2) → **~265 k ions** (Phase 3).
- **Long-time resting state stable**: cyt 110 nM, DTS 252 µM (preserved).
- **Phase 3 still 5/5**: peaks 488 +Ca_ex / 336 −Ca_ex / SOCE diff 152 nM
  (after coupled retune γ_IP3R 0.175 → 0.075 pS, SERCA k_bind_f 500 → 210,
  N_GSN 800 k → 1 400 k, γ_P2X1 0.6 → 1.0 fS).
- **Open** (partly closed): free DTS [Ca²⁺] at peak still drops to ~0.02
  µM minimum (target ≥ 1 µM) because buffer release rate cannot match
  peak IP3R drain rate. **The bound DTS pools are clearly retained** —
  HSP90B1 medium holds ~75 % of its Ca²⁺ during transient — but the free
  pool still drains. Genuine free-DTS retention > 1 µM requires IP3R
  rate refinement (§3.2) or microdomain spatial modelling (§6.3).
- **Reference**: `lab-book-2026-05-12-dts-multibuffer-design.md`.

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

### 3.1 γ_IP3R = 0.175 pS is coupled to SERCA rate constants

- **Original Phase 4 calibration (2026-05-11, commit `1699ac1f`)**: γ_IP3R
  reduced from 10 pS (Zschauer 1988 bilayer) to 0.35 pS, calibrated
  against the analytical 6-state SERCA cycle steady-state flux at the
  Dolan resting state (cyt = 100 nM, DTS = 250 µM) = 112 570 ions/s.
- **Updated Phase 2 retune (2026-05-11, commit `7f4a9ffd`)**: γ_IP3R
  halved to **0.175 pS** alongside halving SERCA `k_bind_f` (1 000 →
  500 µM⁻²·s⁻¹), preserving resting balance at the new lower flux level
  (~57 k ions/s instead of 113 k). This was the coupled flux reduction
  needed to keep Phase 3 peaks in the Dolan band after the CALR buffer
  added a Ca²⁺ reserve to the DTS.
- **Why this matters**: **γ_IP3R is not an independently measured value**
  in our model — it is the value that *balances the chosen SERCA rate
  constants* at the chosen resting state. If SERCA constants change
  (see §3.2), γ_IP3R must be re-derived. Same applies if the buffer
  load (§1.1, §2.1) is re-tuned.
- **Biological plausibility**: 0.175 pS sits within the cellular IP3R
  effective Ca²⁺ conductance range reported by Bezprozvanny 1991 and
  Mak & Foskett 1997 (~0.05–0.5 pS under physiological conditions). The
  10 pS bilayer value is not transferable because Zschauer used
  symmetric high Ca²⁺, where K⁺ contributes negligibly to current.
- **Dissertation framing**: cite as a *calibration anchor*, not a
  measured parameter. Disclose the SERCA + buffer coupling explicitly.

### 3.2 SERCA cycle flux is probably 2–3× too high at rest (post-Phase-2)

- **Current model (after Phase 2 halving of `k_bind_f`)**: ~2.4 cycles/s
  per pump at cyt = 100 nM, **~57 k Ca²⁺ ions/s total** for 11 892 pumps.
  This is closer to biology than the pre-Phase-2 value (~113 k ions/s)
  but still ~2–3× above the literature SERCA3b prediction.
- **Literature SERCA3b kinetics**:
  - Vmax ≈ 30–50 cycles/s at saturating Ca²⁺ (Inesi 1985; Nishi 1992)
  - Km(Ca²⁺) ≈ 0.7–1.1 µM (Dode 2002 — SERCA3 is *less* Ca²⁺-sensitive
    than SERCA2a/b, by design)
  - At cyt = 100 nM with n = 2 Hill: v/Vmax ≈ 2% → ~1 cycle/s per pump
  - Predicted total flux: **~23 800 ions/s**
- **Inherited from Purvis 2008**: the rate constants we use are Purvis's,
  who took them from Dode's protein expression studies. Purvis's
  original k_bind_f = 1 000 µM⁻²·s⁻¹ implies a faster pump than Dode's
  measured Vmax / Km values predict. Phase 2 halved this to 500
  µM⁻²·s⁻¹, partially closing the gap.
- **Implication**: if v0.3 re-derives SERCA constants from primary sources
  and drops them further to the ~24 k ions/s literature prediction,
  γ_IP3R would drop in tandem to ~0.07–0.10 pS, and the cytosolic buffer
  load (currently calibrated to 200:1, vs Sage & Rink's 50:1) could come
  down to biological values. The relative dynamics (Phase 3 transients)
  should be largely preserved because they are dominated by ratio of
  fluxes, not absolutes.
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

### 4.3 PMCA CaM-dissociation rate: in-vitro vs in-vivo (v0.3.1)

- **Caride 2007 in-vitro value**: k12 = 0.033 s⁻¹ (τ = 30 s) for CaM
  dissociation from purified PMCA in proteoliposomes.
- **Our value (v0.3.1)**: k12 = **1.0 s⁻¹** (30× faster, τ = 1 s).
- **Biological justification**: in vivo, PMCA's C-terminal CaM-binding
  domain is competitively occupied by PIP₂ (Penniston & Enyedi 1998
  *J. Membr. Biol.* 165:101), dramatically accelerating CaM
  dissociation. The Caride value is an in-vitro artifact of the
  purified preparation lacking PIP₂.
- **Why this matters**: with the in-vitro value, 87 % of PMCA molecules
  end up trapped in the PMCA·CaM state during a sustained Ca²⁺
  transient (one cycle, then 30 s stuck before another cycle can
  start). Effective PMCA Vmax falls to ~26 ions/s instead of the
  structural Vmax of ~23 k ions/s.
- **Consequence of fix**: DTS overshoot during recovery reduced from
  >1 mM (still rising at t = 1 200 s) to ~660 µM peak then declining.
  cyt recovery toward 100 nM resumed (was locked at 213 nM).
- **The Mazet-Tindall-Gibbins-Fry 2020 critique applies here**: this
  is exactly the kind of mosaic / in-vitro-vs-in-vivo parameter that
  needs context-aware adjustment for whole-cell modelling.
- **Reference**: `lab-book-2026-05-12-pi-cycle-design.md §v0.3.1
  follow-up`.

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

### 7.1 P2X1 — closed (Phase 2.5 / 2026-05-11)

**Note for the writeup**: P2X1 (ionotropic, ATP-gated Ca²⁺ channel) is
*different* from P2Y1 (metabotropic GPCR; ADP → Gαq → PLCβ → IP3). The
two are sometimes confused. We added **P2X1** in Phase 2.5. P2Y1 (and
the rest of the GPCR / PLC cascade) remains v0.3 receptor-signalling
scope — see §7.4.

| Property | Value |
|---|---|
| Channel class | Ionotropic (trimeric ATP-gated cation channel; P_Ca/P_Na ≈ 10) |
| Activation timescale | < 10 ms after extracellular ATP exposure |
| Desensitisation timescale | ~100 ms |
| Recovery from desensitisation | ~30 s |
| Copy number per platelet | ~600–3 000 (Mahaut-Smith 2000/2004; Vial & Evans 2002) |
| Ligand source in vivo | Released ATP from dense granules during activation |

Implementation in v0.2.7: three-state kinetic scheme (Closed → Open →
Desensitised → Closed) with rate constants from Mahaut-Smith / Vial &
Evans. γ_P2X1_Ca = 0.6 fS calibrated against the Dolan SOCE-differential
target. Ca²⁺ flux gated on `CA_EX_UM > 0` — this is precisely what
makes the channel a +Ca_ex-specific contributor and closes the SOCE
differential criterion.

**Outcome**: closed the SOCE-differential gap; Phase 3 now passes 5/5
Dolan criteria (was 4/5). See `lab-book-2026-05-11-dyk-ip3r-design.md
§Phase 2.5`. ATP forcing is a placeholder (rises with the IP3 forcing
curve, τ = 0.5 s) — v0.3 should replace with explicit dense granule
secretion + ectonucleotidase clearance.

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
| ~~NCX~~ | **Closed v0.3.4 (2026-05-12)**: modelled with Hill kinetics + allosteric Ca²⁺-activation gate; V_max = 5 000 ions/s. Reduces DTS overshoot ~25 % vs MCU-only. NCX in platelets is *contested* (Burkhart 2012 proteome detects NCX1/NCX3 but Sage & Rink 1985 reports limited functional activity); modelled here at a moderate level. See `lab-book-2026-05-12-ncx-design.md`. | Closed |

### 7.4 PI cycle — closed (Phase 4 / #31, 2026-05-12)

**Closed in v0.3.0.** Forced IP3 curve replaced by the **Mazet,
Tindall, Gibbins & Fry 2020** *Sci. Rep.* 10:13889 platelet PI cycle
framework. IP3 is now a model output, produced by PLCβ-driven
hydrolysis of PIP2:

```
[Gq activity signal — gq_signal_uM forcing function]
  → PLCβ activation (k_act × Gq)
  → PIP2 hydrolysis (k_cat × plcb_a × PIP2)
  → IP3 + DAG  →  IP3R  →  Ca²⁺ release
```

New species: PIP2, DAG, PLCb_inactive, PLCb_active. New rate dicts:
`K_PLCB` (activation/deactivation/catalysis), `K_PI_CYCLE` (PIP2
resynthesis, IP3 degradation, DAG kinase).

**Calibration**: rates calibrated against (a) resting IP3 = 50 nM
balance, and (b) Dolan Fig. S2 peak shape. Phase 3 remains 5/5 with
the new dynamics.

**What's still simplified**: the full Mazet model has 35 parameters
covering PI ↔ PI4P ↔ PIP2 phosphorylation, IP3 → IP4 / IP2 splits,
and several lipid-binding proteins. Our reduced model uses 5 effective
rates and lumps the PI/PI4P chain into a single resynthesis term.
v0.4 work (#9, receptor signalling) will replace `gq_signal_uM` with
explicit GPCR cascades (P2Y1 / PAR1/4 / GPVI) — and at that point the
PI cycle parameters may be re-derived from Mazet's supplementary
tables more directly.

The largest methodological upgrade in the project to date: IP3 has
gone from a hand-fitted input to a model output.

See `lab-book-2026-05-12-pi-cycle-design.md` for the full design,
calibration log, and results.

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

---

## 9. Numerical methods and implementation notes

This section captures *computational* methodology choices that are
not biological assumptions but still need disclosure in the methods
chapter. Useful in dialogue with mathematical-biology readers who
will (rightly) ask about the interface between continuous ODE
dynamics and the discrete-count state of a whole-cell model.

### 9.0 Model prediction: MCU buffering does *not* accelerate DTS recovery

(Note: this is a §9 entry because it's a *model-level* prediction, not
a biology assumption. It belongs alongside the other numerical /
methodological items.)

In v0.3.3 (commit pending, 2026-05-12) we added the mitochondrial
Ca²⁺ uniporter (MCU) + NCLX efflux (issue #22) on the hypothesis
that mito Ca²⁺ buffering would help close the DTS-overshoot tail
after a transient.

**The hypothesis was falsified.** With MCU active:
- Cyt peak attenuated (479 → 434 nM) ✓
- Mito Ca²⁺ rises 153× during transient and slowly releases ✓
- BUT the DTS overshoot at t = 3000 s went from 758 µM → 1062 µM
  (got *worse*)

**Why**: MCU doesn't extrude Ca²⁺ from the cell — it just
redistributes within the cell. During peak, MCU absorbs cyt Ca²⁺ →
PMCA extrudes *less* (PMCA rate scales with cyt Ca²⁺). Post-stim,
mito slowly releases its load → SERCA pumps that into DTS → DTS
overshoot grows. PMCA at low cyt is still the actual bottleneck.

**Testable prediction**: the model predicts that *MCU-targeted
interventions alone* (e.g. pharmacological MCU activation) would
not accelerate DTS recovery in platelets — and could *prolong* it.
The Ghatge / Ajanel papers measure MCU's effect on cyt Ca²⁺ but
not directly on DTS recovery rate.

**What would actually fix the slow DTS recovery**:
- A second extrusion pathway (NCX / Na⁺/Ca²⁺ exchanger), or
- Faster PMCA cycling (higher pump count / different isoform), or
- Reduced SOCE flux during stimulation.

These are v0.4+ candidates. The interesting take-home is that the
naive intuition (more buffer = faster recovery) is wrong when the
buffer doesn't have its own extrusion pathway.

Reference: `lab-book-2026-05-12-mcu-design.md`.

### 9.1 ODE → integer-state commit: fractional residual carry-over

- **The problem**: the engine stores `BulkMolecules` as 64-bit
  integer counts. The `CalciumSignalling` solver evolves with float
  precision over each 1 s timestep using SciPy's BDF integrator,
  then commits the integer delta back to the counts at the end of
  the step (`reconstruction/platelet/dataclasses/process/calcium_signalling.py`
  in `molecules_to_next_time_step`).
- **Pre-fix behaviour** (commit history before `4b34e60d`,
  2026-05-12):

      delta = np.round(y_final - y0).astype(np.int64)

  Any species whose ODE-derived rate of change satisfies
  |dy/dt| < 0.5 ions/s gets rounded to 0 each step → the count is
  *frozen* at its current integer. The exact frozen point depends
  on where in the rounding band the count happens to land.
- **Symptom we encountered**: IP3 stranded at 205 counts (57 nM)
  instead of the true ODE equilibrium of 181 counts (50 nM),
  because dIP3/dt = production − degradation = 3.62 − 4.10 = −0.48
  ions/s → rounds to 0 each step. PIP2 was similarly stranded
  (|dPIP2/dt| = 0.003 ions/s — deeply inside the rounding band).
  PLCβ active stranded at 144 (rest is 143). Together these
  produced a *false steady state* with elevated cyt (185 nM) and
  the wrong IP3 baseline.
- **The fix** (commit `4b34e60d`): keep a per-species fractional
  residual `self._residual` that persists across timesteps. Each
  step:

      fractional_delta = y_final - y0 + self._residual
      delta = np.round(fractional_delta).astype(np.int64)
      self._residual = fractional_delta - delta

  For a species drifting at −0.48 ions/s, the residual reaches
  −0.96 by step 2 → rounds to −1 → the species correctly loses 1
  ion every ~2 seconds.
- **Why this is dissertation-worthy**: whole-cell modelling sits
  at the interface between continuous biochemistry and discrete
  molecular counts. The naive "evolve as floats, commit as
  integers" recipe is the obvious approach but silently breaks
  for slowly-equilibrating species. Worth pairing in the methods
  chapter with the **Mazet, Tindall, Gibbins & Fry 2020**
  in-vitro-vs-in-vivo critique: *model parameters need context,
  and model state representation needs care.*
- **Audit**: scanned the whole repo (2026-05-12) for the same
  pattern. Only three call sites of `np.round(...).astype(int)`
  or `int(np.round(...))` exist:
  1. `calcium_signalling.py:1161` — fixed.
  2. `wholecell/states/bulk_molecules.py:312` — uses
     **stochastic rounding** via `np.random.choice` weighted by
     remainder fractions; correct by design.
  3. `wholecell/states/bulk_molecules.py:318` — int cast of
     allocations *after* the stochastic distribution; correct.
- **Future considerations**: any new ODE-based process (v0.4
  receptor signalling, granule release kinetics, etc.) will need
  the same residual-carry-over treatment. The residual array is
  *not* persisted across simulation restarts — fine for our
  short single-run dissertation simulations, but a checkpointed
  multi-day simulation would need to save and restore it.
- **Reference**: diagnosis in
  `lab-book-2026-05-12-pi-cycle-design.md §v0.3.2 follow-up`;
  audit + Mazet-Fry framing in this section.
