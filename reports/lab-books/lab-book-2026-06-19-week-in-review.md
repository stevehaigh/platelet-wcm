---
title: "Lab book — 2026-06-19: week in review — the biology added (2026-06-13 → 06-19) and its effect on the model"
---

# Lab book — 2026-06-19: week in review (biology added 06-13 → 06-19)

## One-line summary

In one week the model went from a **Ca²⁺-core + nascent Gq cascade** (validated
against Dolan & Diamond 2014) to a **signalling whole-cell model holding both
halves of the regulatory logic — activation ⊕ inhibition — with terminal
functional outputs and a working antiplatelet drug panel** (≈ v0.6 → v0.7).

## Biology added, by subsystem

### 1. PKC negative feedback — the brakes (PR #54, 06-14)
- DAG → **PKC** coincidence detector (needs DAG *and* cytosolic Ca²⁺), closing
  the previously dead-end DAG branch of the PI cycle.
- **PKC → P2Y1 desensitisation** (`P2Y1_active → P2Y1_desensitised`) and
  **PKC → PLCβ phosphorylation** (PLCβ → a state Gq can't activate).
- *Effect:* two negative-feedback loops throttle the Gq → PLCβ → IP₃ drive
  during sustained agonist. Lowers IP₃ (~−25 %) and the active-receptor pool;
  ~invisible on the cytosolic-Ca²⁺ *peak* (store-limited) — visible on the IP₃ /
  P2Y1-desensitised / PLCβ-phospho readouts.

### 2. Autocrine amplification — the amplifiers + the second wave (PR #54; `runSecondWave.py`, 06-13)
- **Autocrine ADP**: secreted dense-granule `ADP[e]` feeds back onto **P2Y1 → Gq**.
- **Autocrine TXA₂**: synthesised `TXA2[e]` → **TP → Gq**.
- *Effect:* a *weak* primary agonist now produces a sustained **second wave** of
  Ca²⁺/activation that an open-loop model misses; the loops are the lever aspirin
  and (later) P2Y12 inhibitors act on.

### 3. Terminal functional outputs — the observables (PR #54 secretion/thromboxane; PR #58 integrin)
- **Granule secretion**: PKC×Ca-gated SNARE release of dense (ADP, 5-HT) and
  α-granule (fibrinogen, surface P-selectin) cargo; ecto-NTPDase ADP→AMP
  clearance. (Framework-integrated `granule_secretion` analysis plot landed
  today, #8/PR #66.)
- **Thromboxane**: COX-1 → TXA₂ synthesis — the aspirin target.
- **Integrin αIIbβ3 inside-out activation** → the **PAC-1** readout (PR #58).
- *Effect:* the model now emits the actual *experimental* observables —
  lumi-aggregometry (ATP/5-HT release), flow cytometry (P-selectin, PAC-1) —
  rather than only an internal Ca²⁺ trace.

### 4. Ca²⁺-store realism (PR #58 NCX/recovery; PR #59 DTS depletion)
- **NCX fix**: the Na⁺/Ca²⁺ exchanger now runs under EDTA, so the cytosolic
  transient **recovers to baseline** (Dolan Fig. 4C) — the recovery phase.
- **DTS depletion**: `V_IM = 0` thermodynamic store-floor (a passive IP₃R cannot
  drain the store below cytosolic Ca²⁺) + γ_IP3R recalibration → a **stable
  resting fixed point** (fixed the slow resting drift) and honest deep depletion.
- *Effect:* more physical store/recovery dynamics; resting state no longer drifts.

### 5. The inhibitory axis — v0.7 (today, 06-19: PR #68 = #10, PR #71 = Slices 1+4)
- **P2Y12 / Gi** (Slice 2): ADP's *second* receptor — lowers cAMP (the
  clopidogrel/cangrelor target).
- **cAMP / PKA node + VASP/PRI** readout (the clinical P2Y12 assay).
- **cAMP-raising Gs arm** (Slice 1): PGI₂/iloprost (Gs→AC), forskolin, explicit
  PDE3A (cilostazol/dipyridamole).
- **PKA brakes** on IP₃R, **integrin**, and **secretion** (Slice 4) — normalised
  to 1.0 at rest (so Dolan 5/5 preserved by construction).
- *Effect:* a **bidirectional** cAMP/PKA controller. cAMP-lowering (ADP/P2Y12) →
  more activation; cAMP-raising (PGI₂/forskolin/cilostazol) → strong suppression
  of **both** PAC-1 and secretion + a rise in VASP-P; clopidogrel blocks the
  lowering. The headline antiplatelet/vasodilator pharmacology now lives in the
  model.

## Net effect on the model

- **Logic completed:** activation (Gq / PKC) ⊕ inhibition (Gi / Gs – cAMP / PKA).
  Platelet behaviour is now an emergent **threshold** (quiescent below, committed
  above) set by their competition, with autocrine amplification carrying the
  second wave.
- **From a Ca²⁺ figure to a drug panel:** the model can now represent **aspirin**
  (COX-1), **clopidogrel/cangrelor** (P2Y12), **PGI₂/iloprost** (Gs), and
  **cilostazol** (PDE3) — each on its mechanism-specific readout.
- **Validation reframed** (`reports/design/validation-map-2026-06-19.qmd`):
  Dolan "5/5" demoted to a **regression invariant on the Ca²⁺ core**; the model
  is now judged on a *portfolio* of subsystem targets (VASP/PRI, PAC-1,
  lumi-aggregometry, drug dose-response). Dolan 5/5 stayed green throughout.
- **Engineering enabler (v0.62, PR #55):** per-run `RunConfig` replaced
  module-global monkeypatching — every experiment/drug above is a clean per-run
  knob (`p2y12_block`, `pgi2_nM`, `forskolin`, `pde3_block`, `cox1_factor`,
  agonist peaks, perturbation scales), with autocrine `[e]` species fed to the
  ODE by name.

## Known limitation surfaced this week

Cytosolic Ca²⁺ is **architecturally clamped** (buffers CaM/gelsolin/CALR + the
SOCE↔clearance balance): probes showed IP₃R/PMCA/SERCA each move it <5 nM. So
feedback and inhibition are **visible on the functional outputs** (integrin,
secretion, thromboxane, IP₃, cAMP, VASP-P), **not** on free cytosolic Ca²⁺.
This is why the PKC and PKA brakes act on the outputs. Making the axis visible on
cytosolic Ca²⁺ itself would need a Ca²⁺-core buffer/SOCE recalibration — the
deeper future item on issue #70.

## Pointers

Per-session detail: `lab-book-2026-06-13-second-wave.md`,
`-2026-06-14-recovery-phase.md`, `-2026-06-16-tui-diagram-dts-vim.md`,
`-2026-06-19-p2y12-inhibitory-axis.md`, `-2026-06-19-cyclic-nucleotide-slices.md`.
Validation: `reports/design/validation-map-2026-06-19.qmd`. Figures:
`reports/figures/v0.7/`.
