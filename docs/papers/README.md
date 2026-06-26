# Source papers — summaries

One Markdown summary per primary source held in `source-info/calcium-papers/`,
written for AI assistants and the author to grasp **what each paper contributes
to the model** without re-reading the PDF. Each summary states the paper's own
findings (read from the PDF) plus a "Relevance to the platelet WCM" link to the
code/parameters that use it.

> **Provenance vs summary.** For *value-by-value* provenance (which exact number
> came from which table, with verification status), see
> [`../../reports/data/calcium-data-provenance.md`](../../reports/data/calcium-data-provenance.md)
> and [`../../reports/data/platelet-proteome-data-sources.md`](../../reports/data/platelet-proteome-data-sources.md).
> These per-paper files are the narrative companion to that.

## Core model-providing papers

The papers the model's equations and parameters are built on.

| Paper | Contribution |
|-------|--------------|
| [dolan-diamond-2014](dolan-diamond-2014.md) | **The headline validation target** (Fig. 4 Ca²⁺ transient ±Ca_ex); IP3R+SERCA+PMCA+SOCE model; Table S1 representative ICs; IP3 forcing function |
| [purvis-2008](purvis-2008.md) | Upstream GPCR → Gq → PLCβ → IP₃ → Ca²⁺ cascade; resting concentrations; 6 fL cytosol; SERCA cycle + IP3R conductance |
| [sneyd-dufour-2002](sneyd-dufour-2002.md) | Type-2 IP₃R six-state gating kinetics (the Li-Rinzel reduction's parent scheme) |
| [caride-2007](caride-2007.md) | PMCA4b CaM-coupled pump kinetics (basal constants k4/k4r/k5) |
| [hoover-lewis-2011](hoover-lewis-2011.md) | STIM1–Orai1 MWC allosteric SOCE constants (L, Ka, α, f) |

## Pumps, channels, store & mitochondria

| Paper | Contribution |
|-------|--------------|
| [inesi-1985](inesi-1985.md) | SERCA / SR Ca²⁺-ATPase mechanism (2:1 Ca:ATP, cooperativity) |
| [landolfi-1998](landolfi-1998.md) | Mitochondria–ER store-refill functional interactions |
| [ghatge-2026](ghatge-2026.md) | MCU⁻/⁻: MCU loss **reduces** Ca²⁺/function/thrombosis (model diverges — buffer-only MCU) |
| [ajanel-2025](ajanel-2025.md) | MCU in ITAM-dependent activation; mito-ROS mechanism |
| [shehwar-2025](shehwar-2025.md) | Review: platelet mitochondrial Ca²⁺ machinery |

## Systems-biology comparators

Other platelet/cell ODE models — comparators and sources for pathways the WCM has
built or may build.

| Paper | Contribution |
|-------|--------------|
| [sveshnikova-panteleev-2025](sveshnikova-panteleev-2025.md) | 29-species, 4-compartment procoagulant model; volumes/ICs; ~1-molecule Gq-PLC stochastic bottleneck |
| [sveshnikova-2015](sveshnikova-2015.md) | Compartmentalised Ca²⁺ + mitochondria (COPASI); oscillations; subpopulation formation |
| [dunster-2015](dunster-2015.md) | GPVI → Syk → PLCγ2 early signalling (a **not-yet-built** collagen arm) |
| [mazet-2020](mazet-2020.md) | Full platelet PI-cycle ODE model; the "mosaic rate constants" caution |
| [kleppe-2018](kleppe-2018.md) | NO/cGMP/cAMP inhibitory-axis model (relevant to the v0.7 inhibitory axis) |
| [purvis-lahav-2013](purvis-lahav-2013.md) | Review: dynamic encoding/decoding of cell-signalling information |

## Proteome & luminal Ca²⁺ buffers

| Paper | Contribution |
|-------|--------------|
| [burkhart-2012](burkhart-2012.md) | **THE platelet proteome** — per-platelet copy numbers for the species inventory |
| [huang-2021](huang-2021.md) | Classified platelet proteome from genome-wide transcripts (cross-check) |
| [argon-simen-1999](argon-simen-1999.md) | GRP94/HSP90B1 ER chaperone — luminal Ca²⁺ buffer |
| [lievremont-1997](lievremont-1997.md) | BiP/HSPA5 — luminal Ca²⁺ buffer (~25% of ER store) |
| [honore-vorum-2000](honore-vorum-2000.md) | CREC family review (CALU/RCN1/RCN2) — low-affinity ER buffers |
| [vorum-1998](vorum-1998.md) | Calumenin cloning + Ca²⁺-binding (7 EF-hands, low affinity) |

## Integrin / clinical

| Paper | Contribution |
|-------|--------------|
| [coller-1995](coller-1995.md) | GPIIb/IIIa (αIIbβ3) as antithrombotic target; abciximab; Glanzmann thrombasthenia → the integrin module + its knockout |

---

## Key references NOT held in this local directory

These matter to the model but are not PDFs in `source-info/calcium-papers/`, so
they have no summary file here. **Verify in Zotero before citing** — only
first-author + approximate year are safe to recall; do not state a finding's
direction without reading. (See the memory note on citation verification.)

- **Karr et al. 2012** *Cell* 150:389 & **Ahn-Horst et al. 2022** *npj Syst Biol
  Appl* 8:30 — the CovertLab wcEcoli whole-cell model the engine derives from
  (cite per the README attribution).
- **De Young & Keizer 1992** / **Li & Rinzel 1994** — the IP₃R model the code's
  reduced one-ODE-for-h gating is based on (parent of the Sneyd-Dufour scheme).
- **Dode 2002** — SERCA3b cycle kinetics (used via Purvis 2008 Table 1).
- **Mundell et al. 2006** — PKC-dependent P2Y1 desensitisation (the PKC-feedback
  validation target).
- **Gachet 2012** *Purinergic Signal.* — P2Y12 / Gi biology (basis for the
  VASP/PRI inhibitory-axis readout).
- **Zou et al. 2022** / **Stolla et al. 2011** — Rap1b / integrin inside-out
  signalling (the PI3K/Akt→Rap1b arm, #73).
- **Bunne et al. 2024** *Cell* — the "AI Virtual Cell" perspective (thesis
  framing: mechanistic-vs-ML tension; anucleate-platelet counterpoint).

If any of these is later added as a PDF, give it a summary file here following the
same template and link it from the relevant section above.
