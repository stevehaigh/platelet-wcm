# Model diagram — update checklist (pre-v0.5 → v0.63)

Working tick-list for manually updating the BioRender schematic
(`haigh.bio/decks/lab-meeting-2026-05-14/img/image4.png`) to reflect the model
as of **v0.63**. Intentionally `.md` (not `.qmd`) so the checkboxes render and
tick interactively on GitHub while you work through it.

The original figure already has the full Ca²⁺ machinery (P2X1, Orai1/STIM1,
P2Y1/PAR1/PAR4 → Gq → PLC → PIP₂ → IP₃/DAG, DTS with IP₃R + SERCA, luminal
buffers, CaM, mito MCU/NCLX, PMCA, **NCX**, SCS). Four things were greyed out as
"planned": **PKC**, **granules**, **TXA₂R**, **P2Y₁₂**. Three of those are now
built; this list un-greys them, adds the brand-new nodes, and re-points the
greyed region at the *next* frontier (v0.7 inhibitory axis).

> NCX is already drawn and needs no change — the v0.63 fix only made it run in
> the EDTA / no-Ca²⁺ condition too.

---

## Proposed layout

Keep the existing canvas; add into the empty right / lower-right and tighten the
DAG branch into a PKC hub. Rough map (relative to the current figure):

```
                        ATP        ADP            Thrombin        PGI2  NO ┄►v0.7
                         │          │              │   │
══[P2X1]═══[Orai1]═══[P2Y1]══[PAR1/PAR4]══[TP]══[P2Y12•grey]═══════════════════   « MEMBRANE »
                          Gq receptors:  P2Y1 · PAR1 · PAR4 · TP                 ┌─── v0.7 · INHIBITORY AXIS (planned — greyed) ───┐
                                      │  (converge on shared Gq)                 │                                                 │
                                      ▼                                          │ ADP ──► P2Y12 ──► Gi ──┐                        │
                                     Gαq                                         │                        ▼                        │
┌─ DTS (Ca²⁺ store) ──┐               │  activates                               │ PGI2 ► IP ► Gs ► AC ──► cAMP ──► PKA            │
│ IP₃R        SERCA   │               ▼                                          │                         │                       │
│ STIM1  (Ca sensor)  │              PLCβ ─► PIP₂ ─► IP₃ ─┬─► DAG                │ NO ► sGC ► cGMP ► PKG    ▼   PDE3A ┘            │
│ buffers: CALR,      │   ◄─(opens IP₃R on DTS)───────────┘    │                 │                       VASP-P (readout)          │
│ HSP90B1, BiP, CREC  │                                        ▼                 │                                                 │
└─────────────────────┘                                        │  (DAG+Ca²⁺)     └─────────────────────────────────────────────────┘
                                                      (α/β/δ) PKC   ◄══ THE HUB
                          feedback ↑:  ⊣ P2Y1 (brake①)    ⊣ PLCβ (brake②)
                          outputs  ↓:
                                ┌──────────────┬───────────────┐
                                ▼              ▼               ▼
┌─────── mito ────────┐     SECRETION      INTEGRIN        TXA₂ SYNTHESIS
│ MCU        NCLX     │     dense: ADP,    αIIbβ3          cPLA₂►COX1►synthase
│ uptake     efflux   │       5-HT ►[e]    resting ⇌       (aspirin ⊣ COX1)
└─────────────────────┘     α: FGA ►[e],     active             │
                              SELP►surface  PAC-1 (Ab)          ▼
                                                          TXA₂ ►[e]► TP ─► Gq
                          AUTOCRINE:  ADP[e] ─► P2Y1    |    TXA₂[e] ─► TP
══[NCX]═══════════════════════════════════════════════════════════[PMCA]═══════   « Ca²⁺ extrusion »
```

*(Arrows: `─►` forward · `◄─` feedback · `⊣` inhibits/brake · `⇌` reversible.
`(α/β/δ)` = the lumped PKC isoforms (conventional α/β + novel δ). The PKC hub fans
out to 2 brakes ↑, 2 autocrine amplifiers ↺, and 3 terminal outputs ↓.)*

---

## A. PKC — promote stub → central hub

- [ ] Replace the greyed *"(PKC not yet in model)"* label on the DAG branch with
      a real **PKC** node (lumped conventional α/β + novel δ).
- [ ] Draw the activation input: **DAG + Ca²⁺ → PKC** (PKC is a coincidence
      detector; both inputs required).
- [ ] Brake arrow ①: **PKC ⊣ P2Y1** — P2Y1 desensitisation
      (`P2Y1_active → P2Y1_desensitised`; Mundell 2006 / Nicholas 2023).
- [ ] Brake arrow ②: **PKC ⊣ PLCβ** — PLCβ phosphorylation out of the
      Gq-activatable pool (`PLCβ → PLCβ-P`; Purvis 2008).
- [ ] Style PKC as the hub: 4 feedback loops (2 brakes above + 2 autocrine
      amplifiers below) **plus** 3 terminal outputs (secretion, thromboxane,
      integrin). Fan the arrows out so this reads as the new centre of gravity.

## B. Granules + secretion — un-grey

- [ ] Un-grey the *"(Granules not yet in model)"* circles; split into a
      **dense granule** (cargo: ADP, 5-HT) and an **α-granule** (cargo:
      fibrinogen/FGA, P-selectin/SELP).
- [ ] Add the **PKC × Ca²⁺ coincidence gate** driving secretion (same gate as
      thromboxane + integrin; zero at rest).
- [ ] Secretion arrows to SCS/`[e]`: dense **ADP → [e]**, **5-HT → [e]**;
      α **FGA → [e]**.
- [ ] **SELP → surface P-selectin** (`SELP_surface[pl]`) — the degranulation /
      activation marker (note it stays on the membrane, not released).

## C. Thromboxane — un-grey TXA₂R + add the synthesis arm

- [ ] Add a **TXA₂ synthesis** node: lumped **cPLA₂ → COX-1 → TXA₂-synthase**,
      Ca²⁺ × PKC-gated.
- [ ] Output **TXA₂ → [e]**; show decay **TXA₂ → TXB₂** (t½ ≈ 30 s, the stable
      ELISA metabolite).
- [ ] Un-grey **TXA₂R → TP receptor** (reversible activation by TXA₂;
      ~1000 copies).
- [ ] Arrow **TP → Gq** (active TP joins the shared `total_active_R` pool).
- [ ] (Optional) Mark the **aspirin** lever: aspirin ⊣ COX-1 (`cox1_factor = 0`
      abolishes the whole loop).

## D. Integrin αIIbβ3 — NEW terminal output

- [ ] Add **αIIbβ3** on the plasma membrane as a **resting ⇌ active**
      conformational switch (`aIIbb3_resting` 80 000 ⇌ `aIIbb3_active`).
- [ ] Input arrow: **PKC + Ca²⁺ gate** → activation (lumps CalDAG-GEFI → Rap1b →
      talin/kindlin inside-out signalling).
- [ ] Readout label: **active fraction = PAC-1** (activation-specific antibody,
      flow cytometry).
- [ ] (Optional) Mark the **Glanzmann / αIIbβ3-antagonist** knockout lever
      (`integrin_act_scale = 0`).

## E. Autocrine loops — the arrows that make PKC a hub

- [ ] **ADP loop**: secreted **ADP[e] → P2Y1** (closes PKC → secretion → ADP →
      P2Y1). Optionally show self-limit: ecto-NTPDase **ADP[e] → AMP[e]**.
- [ ] **TXA₂ loop**: secreted **TXA₂[e] → TP → Gq** (autocrine amplifier).

## F. v0.7 frontier — repaint the greyed region (planned, not yet built)

- [ ] Keep **P2Y₁₂** greyed, but reframe it as the entry to a new greyed
      **inhibitory axis** block: ADP → **P2Y₁₂ → Gi → ↓ adenylate cyclase**.
- [ ] Greyed: **PGI₂ → IP receptor → Gs → adenylate cyclase → cAMP → PKA**.
- [ ] Greyed: **NO → sGC → cGMP → PKG**.
- [ ] Greyed: **PDE3A** degrading cAMP; **VASP-P** as the activation-state
      readout.

## G. Housekeeping

- [ ] Update the figure legend / caption for v0.63 (PKC hub, 3 outputs, autocrine
      loops; NCX now active in both Ca²⁺ conditions).
- [ ] Add a small version tag (e.g. "model v0.63") so the deck slide is
      self-dating.

---

## BioRender icon reference

Grounded in a library search. ✅ = on-target & placeable; ⚠️ = stand-in or needs
a direct re-search in the BioRender app (the exact molecule has no current
placeable icon). Match the **existing receptors'** GPCR style for any new GPCR so
the figure stays visually consistent.

| Node | Suggested icon (search term) | Status |
|---|---|---|
| TP / P2Y₁₂ / IP receptors | **"GPCR with G protein (schematic)"**, **"GPCR (folded, editable, with ligand and G protein)"**, **"GPCR (editable)"** | ✅ |
| Gq / Gi / Gs | **"G protein (heterotrimeric, schematic)"**; active/inactive Gα variants | ✅ |
| PKC | re-search **"PKC"** in-app for a 2021+ icon; else stand-in **"Lck (active, schematic)"** or the oval from **"PKA bound to cAMP"** | ⚠️ no exact placeable |
| Dense / α-granule | **"Vesicle with particles"**, **"Synaptic vesicle with neurotransmitters (small)"**; or re-search **"dense granule"** | ⚠️ generic stand-in |
| Secretion event | **"Exocytosis (linear membrane, medium vesicle)"**, **"Exocytosis (x-large circle membrane, large vesicle)"** | ✅ |
| TXA₂ synthesis enzymes | **"Generic enzyme 6a (schematic)"**, **"Generic enzyme 9a (schematic)"**, or **"Phospholipase A2 (pathway)"** (no COX-1 icon exists) | ⚠️ generic stand-in |
| Aspirin ⊣ COX-1 | **"Competitive inhibition"**, **"Ibuprofen"**, or re-search **"aspirin"** | ⚠️ |
| αIIbβ3 integrin | **"VLA-4 (extended, with ligand)"** (only placeable heterodimer); or compose **"Integrin (alpha, monomer, extended 1, schematic)"** + **"Integrin (beta, monomer, bent, schematic)"**; or re-search **"integrin heterodimer 2021"** | ⚠️ αIIbβ3 schematic exists but unplaceable |
| PAC-1 antibody | **"Antibody (conjugated, editable)"** (with fluorophore — fits flow cytometry), **"Antibody"**, **"Antibody IgG (round style, editable)"** | ✅ |
| Adenylate cyclase (v0.7) | **"Adenylyl Cyclase (schematic)"**, **"Adenylyl cyclase (activated, schematic)"** | ✅ |
| PDE3A (v0.7) | **"Phosphodiesterase 3 (PDE3, pathway)"** | ✅ |
| NO (v0.7) | **"Nitric oxide (molecule)"** | ✅ |
| cAMP (v0.7) | re-search **"cAMP (simple structure)"** in-app | ⚠️ no placeable |
| sGC (v0.7) | re-search **"soluble guanylyl cyclase"** in-app | ⚠️ no placeable |
| PGI₂ (v0.7) | re-search **"PGI2"** / **"prostaglandin I2"** in-app | ⚠️ no icon found |
| VASP-P (v0.7) | **"Phosphorylated protein"** as a readout stand-in | ⚠️ no VASP icon |
| (Aggregation note) | **"Platelet (active)"**, **"Platelet (resting)"** | ✅ |

**Caveats:** PKC, COX-1, cAMP, sGC, PGI₂ and VASP have no clean placeable icon in
the API results — the table gives either a stand-in or the exact term to type into
BioRender directly (newer 2021+ assets often surface in-app but not via search).
The αIIbβ3 closed/open heterodimer schematics exist but are pre-2021 (unplaceable);
the VLA-4 heterodimer or composed α+β monomers are the placeable route.
