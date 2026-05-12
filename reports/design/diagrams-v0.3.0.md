---
title: "Platelet WCM — Ca²⁺ pathway diagrams (v0.3.0)"
---

# Platelet WCM — Ca²⁺ pathway diagrams (v0.3.0)

Mermaid versions of the ASCII diagrams in `biology-overview-2026-05-07.md`.
GitHub renders these natively; locally try the VS Code Mermaid plugin.

Updated 2026-05-12 (post Phase 4 / #31 — PI cycle replaces forced IP3).

---

## 1. Full Ca²⁺ pathway — top-level architecture

```mermaid
flowchart TB
    %% External signal
    GqSig(["Gq activity signal<br/>(stand-in for receptor stimulation;<br/>v0.4 replaces with explicit GPCR cascade)"]):::external

    %% Extracellular
    subgraph EXT["Extracellular space"]
        CaEx["Ca²⁺_ex<br/>1.2 mM"]:::ion
        ATPex["ATP_ex<br/>(released from<br/>dense granules)"]:::ion
    end

    %% Plasma membrane
    subgraph PM["Plasma Membrane"]
        P2X1["P2X1<br/>ATP-gated<br/>(fast entry, +Ca_ex only)"]:::channel
        Orai1["Orai1<br/>CRAC (SOCE)"]:::channel
        PMCA["PMCA<br/>(basal + CaM-activated)"]:::pump
        PMleak["PM leak"]:::leak
    end

    %% Cytosol — PI cycle and buffers
    subgraph CYT["Cytosol (6 fL, ~100 nM free Ca²⁺ at rest)"]
        PIP2["PIP₂<br/>~112 k"]:::substrate
        PLCb["PLCβ<br/>(active/inactive)"]:::enzyme
        IP3["IP₃<br/>~50 nM (rest)<br/>model output"]:::messenger
        DAG["DAG"]:::messenger
        CaCyt["Ca²⁺_cyt<br/>~100 nM free<br/>~85% buffered"]:::ion
        CaM["CaM ladder<br/>(20 481 molecules)"]:::buffer
        GSN["Gelsolin proxy<br/>(1.4 M sites)"]:::buffer
    end

    %% DTS membrane
    subgraph DTSmem["DTS Membrane"]
        IP3R["IP3R<br/>deYoung-Keizer<br/>(1 328 channels)"]:::channel
        SERCA["SERCA<br/>(6-state E1/E2)"]:::pump
        STIM1["STIM1 sensor<br/>(monomer ⇌ dimer)"]:::sensor
    end

    %% DTS lumen
    subgraph DTS["DTS Lumen (0.26 fL, ~250 µM free at rest, ~85% buffered)"]
        CaDTS["Ca²⁺_DTS free<br/>~250 µM"]:::ion
        CALR["CALR<br/>C-domain + P-domain<br/>(528 k sites)"]:::buffer
        HSP90B1["HSP90B1<br/>M + L sites<br/>(150 k)"]:::buffer
        BiP["BiP / HSPA5<br/>(50 k sites)"]:::buffer
        CREC["CREC pool<br/>CALU + RCN1/2<br/>(60 k sites)"]:::buffer
    end

    %% Signal flow — PI cycle
    GqSig -.->|activates| PLCb
    PLCb -.->|cat. hydrolysis| PIP2
    PIP2 -->|"→ IP₃"| IP3
    PIP2 -->|"→ DAG"| DAG

    %% IP3 → IP3R → Ca²⁺ release
    IP3 -.->|activates| IP3R
    CaDTS -->|"flux through IP3R"| CaCyt

    %% SERCA pumps Ca²⁺ back into DTS
    CaCyt -->|"SERCA pumps"| CaDTS

    %% Buffers in cyt
    CaCyt <-.->|buffer| CaM
    CaCyt <-.->|buffer| GSN

    %% Buffers in DTS (CALR is dominant)
    CaDTS <-.->|buffer| CALR
    CaDTS <-.->|buffer| HSP90B1
    CaDTS <-.->|buffer| BiP
    CaDTS <-.->|buffer| CREC

    %% SOCE
    CaDTS -.->|senses store depletion| STIM1
    STIM1 -.->|activates| Orai1
    CaEx -->|"SOCE refill"| CaCyt

    %% P2X1 fast entry
    ATPex -.->|gates| P2X1
    CaEx -->|"P2X1 fast"| CaCyt

    %% PMCA + PM leak
    CaCyt -->|extrudes| PMCA
    PMCA --> CaEx
    CaEx -->|leak| CaCyt
    CaM -.->|"Ca₄·CaM activates"| PMCA

    classDef ion fill:#fce4a6,stroke:#c8881e,stroke-width:2px,color:#000
    classDef channel fill:#a6d8fc,stroke:#1e6ac8,stroke-width:2px,color:#000
    classDef pump fill:#d4a6fc,stroke:#7e1ec8,stroke-width:2px,color:#000
    classDef buffer fill:#a6fcc4,stroke:#1ec851,stroke-width:2px,color:#000
    classDef sensor fill:#fca6a6,stroke:#c81e1e,stroke-width:2px,color:#000
    classDef enzyme fill:#fcd4a6,stroke:#c87a1e,stroke-width:2px,color:#000
    classDef substrate fill:#e6e6e6,stroke:#666,stroke-width:1px,color:#000
    classDef messenger fill:#fff8a6,stroke:#c8b21e,stroke-width:2px,color:#000
    classDef leak fill:#cccccc,stroke:#666,stroke-width:1px,stroke-dasharray:3 3,color:#000
    classDef external fill:#ffffff,stroke:#000,stroke-width:1px,stroke-dasharray:4 2,color:#000
```

**Reading the diagram:**
- Solid arrows = Ca²⁺ flux
- Dashed double arrows = reversible binding (buffers)
- Dotted arrows = regulatory / activation (not Ca²⁺ flow)
- Eight coupled mechanisms — PI cycle generates IP₃ → IP3R releases DTS Ca²⁺ → cyt buffers absorb → SERCA and PMCA extrude → SOCE / P2X1 refill from outside.

---

## 2. PI cycle (Mazet 2020 framework) — replaces forced IP3

```mermaid
flowchart LR
    Gq(["Gq signal<br/>(GQ_REST=0.1 µM<br/>GQ_PEAK=2 µM)"]):::external

    PLCbi["PLCβ_inactive<br/>~857 at rest"]:::enzyme
    PLCba["PLCβ_active<br/>~143 at rest"]:::active
    PIP2["PIP₂<br/>~112 k"]:::substrate
    PIsource(["PI / PI4P pool<br/>(lumped, not tracked)"]):::external
    IP3["IP₃<br/>~50 nM rest<br/>~275 nM peak"]:::messenger
    DAG["DAG"]:::messenger
    IPout(["IP₂ / IP₄<br/>(out of model)"]):::external
    PAout(["PA<br/>(out of model)"]):::external

    Gq -.->|"k_act × Gq"| PLCbi
    PLCbi -->|activated| PLCba
    PLCba -->|"k_inact (τ ~3s)"| PLCbi

    PIsource -->|"k_resynth (3.62/s)"| PIP2
    PIP2 -->|"k_cat × PLCb_a × PIP2"| IP3
    PIP2 -->|"k_cat × PLCb_a × PIP2"| DAG
    PLCba -.->|catalyses| PIP2

    IP3 -->|"k_ip3_deg (τ ~50s)"| IPout
    DAG -->|"k_dag_deg (τ ~20s)"| PAout

    classDef external fill:#ffffff,stroke:#000,stroke-width:1px,stroke-dasharray:4 2,color:#000
    classDef enzyme fill:#fcd4a6,stroke:#c87a1e,stroke-width:2px,color:#000
    classDef active fill:#fc8c2e,stroke:#c84e1e,stroke-width:2px,color:#000
    classDef substrate fill:#e6e6e6,stroke:#666,stroke-width:1px,color:#000
    classDef messenger fill:#fff8a6,stroke:#c8b21e,stroke-width:2px,color:#000
```

**Key idea**: the Gq signal (currently a forced curve, replaced by
explicit receptor cascade in v0.4) activates PLCβ; active PLCβ
catalyses PIP₂ → IP₃ + DAG; IP₃ feeds the IP3R; DAG is currently
unused but kept for future PKC modelling.

---

## 3. IP3R deYoung-Keizer / Li-Rinzel model

```mermaid
flowchart LR
    Closed["Channel closed<br/>(h = 1, Po = 0 when IP₃ low)"]:::off
    Open["Channel open<br/>Po = m∞⁴ × h"]:::on
    Inhibited["Ca²⁺-inhibited<br/>(h ↓)"]:::desens

    Closed -->|"IP₃ binds (d₁=0.13 µM)<br/>+ Ca²⁺ activates (d₅=0.082 µM)"| Open
    Open -->|"Ca²⁺ at site 3 (d₂=1.05 µM)<br/>a₂=0.2 µM⁻¹·s⁻¹<br/>τ_h ≈ 3.7 s"| Inhibited
    Inhibited -->|"Ca²⁺ falls<br/>slow recovery"| Closed
    Open -->|"if Ca²⁺ rises<br/>(slow inactivation)"| Inhibited

    classDef off fill:#cccccc,stroke:#666,stroke-width:2px,color:#000
    classDef on fill:#a6fcc4,stroke:#1ec851,stroke-width:2px,color:#000
    classDef desens fill:#fca6a6,stroke:#c81e1e,stroke-width:2px,color:#000
```

`h` is the slow inactivation variable; m∞ is quasi-steady (fast).
Tetrameric cooperativity → Po = m∞⁴ × h. Flux through the open
channel is Nernst-based: I = γ × N × Po × (ψ_IM − E_Ca,IM) × NA/(zF),
with γ = 0.075 pS (calibration anchor).

---

## 4. P2X1 — 3-state ATP-gated channel

```mermaid
flowchart LR
    C["Closed<br/>(1 000 at rest)"]:::off
    O["Open<br/>(conducts Ca²⁺)"]:::on
    D["Desensitised<br/>(non-conducting)"]:::desens

    C -->|"k_act × ATP<br/>(30 µM⁻¹·s⁻¹)"| O
    O -->|"k_close<br/>(5 s⁻¹)"| C
    O -->|"k_des<br/>(10 s⁻¹, τ ~100 ms)"| D
    D -->|"k_rec<br/>(0.03 s⁻¹, τ ~30 s)"| C

    classDef off fill:#cccccc,stroke:#666,stroke-width:2px,color:#000
    classDef on fill:#a6fcc4,stroke:#1ec851,stroke-width:2px,color:#000
    classDef desens fill:#fca6a6,stroke:#c81e1e,stroke-width:2px,color:#000
```

Ca²⁺ flux through the Open state is **gated on extracellular Ca²⁺
availability** — in the EDTA (−Ca_ex) condition, the channel still
cycles through C → O → D but contributes no Ca²⁺. This is what
produces the SOCE differential in the Phase 3 ±Ca_ex comparison.

---

## 5. SERCA E1/E2 cycle (6 states, 2 Ca²⁺ per cycle)

```mermaid
flowchart LR
    E2["E2<br/>(DTS-facing, empty)"]:::state
    E1["E1<br/>(cyt-facing, empty)"]:::state
    E1Ca["E1·Ca²⁺<br/>(2 Ca bound from cyt)"]:::state
    E1PCa["E1P·Ca²⁺<br/>(phosphorylated)"]:::state
    E2PCa["E2P·Ca²⁺<br/>(DTS-facing, occluded)"]:::state
    E2P["E2P<br/>(empty, phosphorylated)"]:::state

    E2 -->|"k_shuttle (600 s⁻¹)"| E1
    E1 -->|"+ 2 Ca²⁺_cyt<br/>k_bind_f × Ca²"| E1Ca
    E1Ca -->|"k_phos (700 s⁻¹)"| E1PCa
    E1PCa -->|"k_conf (600 s⁻¹)"| E2PCa
    E2PCa -->|"k_release<br/>→ 2 Ca²⁺_DTS"| E2P
    E2P -->|"k_dephos (500 s⁻¹)"| E2

    E1 -->|"k_shuttle reverse"| E2
    E1Ca -->|"k_bind reverse"| E1
    E1PCa -.->|reverse| E1Ca
    E2PCa -.->|reverse| E1PCa
    E2P -.->|reverse| E2PCa
    E2 -.->|reverse| E2P

    classDef state fill:#d4a6fc,stroke:#7e1ec8,stroke-width:2px,color:#000
```

Solid arrows = forward cycle (cyt → DTS pumping); dashed = reverse.
At rest the cycle runs at ~14 k cycles/s × 2 Ca²⁺ = ~28 k Ca²⁺/s.

---

## 6. PMCA — Caride 2007 5-state CaM-coupled scheme

```mermaid
flowchart TB
    PMCA["PMCA (free)"]:::state
    PMCACa["PMCA·Ca²⁺<br/>(basal path)"]:::state
    CaMPMCA["Ca₄·CaM·PMCA<br/>(CaM-bound, no Ca²⁺ pumped)"]:::state
    CaMPMCACa["Ca₄·CaM·PMCA·Ca²⁺<br/>(CaM-activated, ready to pump)"]:::state
    PMCACaM["PMCA·CaM<br/>(deactivating; slow CaM release)"]:::state

    PMCA -->|"+ Ca²⁺<br/>basal"| PMCACa
    PMCACa -->|"k_cat → Ca²⁺_ex"| PMCA
    PMCA -->|"+ Ca₄·CaM<br/>(CaM activation)"| CaMPMCA
    CaMPMCA -->|"+ Ca²⁺"| CaMPMCACa
    CaMPMCACa -->|"k_cat (5× basal)<br/>→ Ca²⁺_ex"| CaMPMCA
    CaMPMCA -->|"k₁₁: → 4 Ca²⁺_cyt"| PMCACaM
    PMCACaM -->|"k₁₂ (slow τ ~30 s)<br/>+ Ca²⁺ → CaM"| PMCA

    classDef state fill:#d4a6fc,stroke:#7e1ec8,stroke-width:2px,color:#000
```

Two pumping paths: basal (slow, runs at all Ca²⁺) and CaM-activated
(5× faster, requires Ca₄·CaM binding). Step 11 releases 4 Ca²⁺ from
the CaM back to cyt during deactivation; step 12 is the slow CaM
unbinding.

---

## 7. STIM1 / Orai1 SOCE

```mermaid
flowchart LR
    STIMCa["STIM1·Ca²⁺_DTS<br/>(inactive, ~3 800 at rest)"]:::sensor
    STIMfree["STIM1_free<br/>(active monomer, ~440 at rest)"]:::sensor
    STIMdim["STIM1 dimer<br/>(~11 at rest;<br/>thousands when DTS empty)"]:::active

    Orai1["Orai1 tetramer<br/>(360 channels)"]:::channel
    OraiOpen["Orai1 open<br/>(MWC allosteric)"]:::on

    CaEx["Ca²⁺_ex 1.2 mM"]:::ion
    CaCyt["Ca²⁺_cyt"]:::ion

    STIMCa -->|"Ca²⁺ unbinds<br/>(DTS depletion)"| STIMfree
    STIMfree -->|"+ STIM1_free<br/>(dimerisation)"| STIMdim
    STIMdim -.->|"binds Orai1 tetramer<br/>(MWC cooperativity)"| Orai1
    Orai1 --> OraiOpen
    CaEx -->|"flux through Orai1"| CaCyt

    classDef sensor fill:#fca6a6,stroke:#c81e1e,stroke-width:2px,color:#000
    classDef active fill:#fc6a6a,stroke:#c81e1e,stroke-width:2px,color:#000
    classDef channel fill:#a6d8fc,stroke:#1e6ac8,stroke-width:2px,color:#000
    classDef on fill:#a6fcc4,stroke:#1ec851,stroke-width:2px,color:#000
    classDef ion fill:#fce4a6,stroke:#c8881e,stroke-width:2px,color:#000
```

Store-operated entry. As DTS [Ca²⁺] drops during the transient, STIM1
releases its Ca²⁺ → dimerises → binds Orai1 → CRAC channel opens →
Ca²⁺ flows in from outside. Slow (>10 s timescale) — too slow to
contribute to the early cyt peak, which is what P2X1 covers.

---

## BioRender assembly guide

These diagrams are good for the design/code documentation. For
**dissertation / publication figures**, BioRender is the standard.
Suggested asset/template search terms (icon name, where it would
slot into the diagrams above):

### Top-level pathway (diagram 1)
- "platelet" — template / cell outline
- "endoplasmic reticulum" / "DTS tubule" — for the DTS compartment (DTS is the platelet equivalent of the smooth ER; can use ER imagery)
- "plasma membrane" — bilayer
- "IP3 receptor" / "IP3R" — for the DTS membrane channel
- "SERCA" / "Ca²⁺ ATPase" — for the SERCA pump
- "PMCA" — plasma membrane Ca²⁺ ATPase
- "Orai1" / "CRAC channel" — for the SOCE channel
- "STIM1" / "stromal interaction molecule" — DTS Ca²⁺ sensor
- "P2X receptor" — ATP-gated cation channel
- "calmodulin" — CaM (often shown as dumbbell shape with two Ca²⁺-binding lobes)
- "gelsolin" — Ca²⁺-regulated actin-severing protein
- "calreticulin" — ER luminal buffer
- "calcium ion" / "Ca²⁺" — Ca²⁺ symbol

### PI cycle (diagram 2)
- "PLC beta" / "phospholipase C" — PLCβ
- "PIP2" / "phosphatidylinositol 4,5-bisphosphate" — substrate
- "IP3" / "inositol trisphosphate" — second messenger
- "DAG" / "diacylglycerol" — co-product
- "G alpha q" / "Gαq" — receptor-activated G-protein

### Suggested layout
1. Draw the platelet outline; place SCS (surface-connected canalicular system, distinct from DTS — visible as PM invaginations) and DTS as separate internal structures.
2. PM populated with: P2X1, Orai1, PMCA, and (implicit) the basal leak. Show the SCS as PM extensions.
3. DTS membrane: IP3R, SERCA, STIM1.
4. DTS lumen: free Ca²⁺ + buffer proteins (CALR with C/P domains, HSP90B1, BiP, CREC).
5. Cytosol: free Ca²⁺ + CaM (3 states) + gelsolin + PI cycle (PIP2 → PLCβ → IP3 + DAG).
6. Arrows in Ca²⁺-flux convention: brown/orange for Ca²⁺ ions; dashed for regulatory.

The BioRender search tools (see code-overview for accessing the
Anthropic MCP integration) can help find these assets.
