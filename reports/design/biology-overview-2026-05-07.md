# Platelet WCM — biology overview (v0.2)

A single-cell, deterministic model of intracellular Ca²⁺ dynamics in a
resting / activating human platelet. Reproduces the IP3-mediated
calcium transient that initiates platelet activation, calibrated
against Dolan & Diamond 2014 *Biophys. J.* 106:2049–2060 Fig. 4.

## The biological story we tell

A platelet at rest holds cytosolic Ca²⁺ at ~100 nM against an
extracellular concentration of ~1.2 mM and a cytoplasmic store
concentration of ~250 µM, by balancing pumps and channels across two
membranes. On stimulation (collagen, thrombin, ADP), receptors generate
inositol-1,4,5-trisphosphate (IP3); IP3 opens IP3 receptors on the
dense tubular system (DTS, the platelet's ER equivalent) and Ca²⁺
floods the cytosol. Cytosolic Ca²⁺ then drives granule secretion,
integrin activation and shape change. The transient resolves over
seconds-to-minutes via SERCA reuptake, plasma-membrane extrusion, and
store-operated entry refilling the DTS.

## The pathway, as wired in v0.2

```
                     +-- SOCE: STIM1 (DTS) --> Orai1 (PM) -- + 1.2 mM Ca²⁺_ex
                     |                                       |
[IP3 spike] --> IP3R (DTS membrane) --> Ca²⁺_cyt rises --> +-- PMCA (PM) --> Ca²⁺_ex
                                            |              |
                                            +-- CaM ladder (cytosolic buffer; activates PMCA)
                                            |
                                            +-- SERCA (DTS) --> refill DTS Ca²⁺ store
```

Five mechanisms coupled through cytosolic Ca²⁺. ATP is consumed by
both ATPase pumps (SERCA, PMCA).

## Components, with the published kinetic model adopted for each

| Mechanism | Kinetic model | Source |
|---|---|---|
| **IP3R** | 6-state Markov chain (n / o / a / i1 / i2 / s); Po = ((0.9·a + 0.1·o)/total)⁴ tetramer cooperativity; Nernst flux | Sneyd & Dufour 2002 (rates), Purvis & Bhatt 2008 (parameterisation), Zschauer 1988 (single-channel γ) |
| **SERCA** | 6-state E1/E2 enzymatic cycle (E2 ⇌ E1 ⇌ E1·Ca → E1P·Ca ⇌ E2P·Ca → E2P → E2); 1 ATP / 2 Ca²⁺ | Dode 2002 (isoform 3b kinetics), Purvis 2008 (rate constants) |
| **PMCA** | 5-state CaM-coupled scheme: basal path (Ca²⁺ binding → extrusion, V_max = 5.5 s⁻¹) **plus** Ca₄·CaM-activated path (V_max = 30 s⁻¹, ~5.5× faster); 1 ATP / Ca²⁺ | Caride et al. 2007 Table 3 |
| **Calmodulin** | Two-lobe Ca²⁺-binding ladder: CaM_free ⇌ Ca₂·CaM ⇌ Ca₄·CaM; ~20 481 molecules; acts as cytosolic Ca²⁺ buffer + PMCA activator | Caride et al. 2007 steps 6–7 |
| **STIM1 sensor cycle** | DTS-bound (Ca-loaded, inactive) ⇌ free monomer ⇌ dimer (active sensor); detailed-balance rate constants | Dolan 2014 + Hoover & Lewis 2011 |
| **Orai1 / SOCE** | Monod–Wyman–Changeux allosteric model: STIM2 dimers translocate into puncta (Hill function on cytosolic Ca²⁺), bind Orai1 tetramers cooperatively; channel opening as fraction of bound STIM2 | Hoover & Lewis 2011, Dolan 2014 puncta entry (eq. 2 + eq. 4) |
| **IP3 production** | v0.2 placeholder: pre-programmed time curve fitted to Dolan Fig. S2 (5.5× rise, τ_rise=3 s, τ_decay=60 s). v0.3 will replace with P2Y1 / Gq / PLCβ cascade. | Dolan 2014 Fig. S2 |

## Compartments and copy numbers (Dolan 2014 Table S1, except where noted)

| Compartment | Volume | Key species & counts |
|---|---|---|
| Cytosol | 6 fL (Purvis 2008 direct) | 361 Ca²⁺ (100 nM), 181 IP3 (50 nM), 20 481 CaM, 10.8 M ATP |
| DTS (cytoplasmic store) | 0.258 fL (4.3 % of cyt) | 38 842 Ca²⁺ (250 µM), 1 328 IP3R (across 6 sub-states), 11 892 SERCA (across 6 sub-states), 4 265 STIM1 (across 3 sub-states) |
| Plasma membrane | (surface) | 769 PMCA (5 sub-states), 1 447 Orai1 monomers (≈ 360 tetrameric channels) |
| Extracellular | infinite reservoir | fixed 1.2 mM Ca²⁺ |

## What we can ask the model

Three biologically-distinct run conditions, exposed in the webapp and
on the CLI:

| Condition | What it tests |
|---|---|
| 🩸 **IP3 transient (+Ca²⁺_ex)** | Canonical activation. IP3R-driven release, SERCA reuptake, SOCE refilling. Phase 1 acceptance criterion: peak in 200–800 nM band. |
| 🧪 **EDTA transient (no Ca²⁺_ex)** | Isolates IP3R contribution by removing extracellular Ca²⁺. SOCE inactive (no source); PM leak inactive. Compares against the +Ca_ex condition to test SOCE dependence (Dolan 2014 Fig. 4). |
| 🛌 **Resting (no stimulus)** | IP3 stays at 50 nM baseline; no transient driven. Inspects the model's rest behaviour. |

Plus a Phase 3 driver that runs the +Ca_ex and EDTA conditions
back-to-back and produces the Dolan Fig. 4 comparison figure.

## What v0.2 reproduces

Validated against Dolan 2014 Fig. 4 + Fig. 3B filtering criteria, run
2026-05-06 (figure: `reports/figures/phase3-dolan-fig4-2026-05-06.png`):

- ✓ Peak Ca²⁺_cyt > 200 nM with extracellular Ca²⁺ (299 nM measured)
- ✓ Peak Ca²⁺_cyt > 200 nM under EDTA (298 nM measured)
- ✓ Peak under EDTA within Dolan's ±30 % band (192–358 nM)
- ✓ SOCE current correctly zero under EDTA (no Ca²⁺ source)
- ✓ STIM1 dimers rise from 22 to ~810 on store depletion (sensing works)
- ✗ SOCE peak differential ≥ 100 nM between conditions (measured 1 nM)
- ✗ Peak with extracellular Ca²⁺ within Dolan's ±30 % band (measured 299 nM, band 315–585 nM)

3 of 5 Dolan acceptance criteria pass. The two failures share one
upstream cause — the DTS empties before SOCE can build a sustained
plateau — documented as a v0.2 known limitation (design doc §6.8 D7).

## What v0.2 does *not* yet model

| Biology | Status | Tracked as |
|---|---|---|
| Mitochondrial Ca²⁺ uniporter (MCU + mNCX) — captures Ca²⁺ during the spike, slowly releases over minutes | Not implemented; three platelet MCU papers in `source-info/calcium-papers/` (Ajanel 2025, Ghatge 2026, Shehwar 2025) | Issue **#22** |
| P2Y1 / Gq / PLCβ upstream cascade — endogenous IP3 production from ADP signalling | v0.2 uses a forced curve | Phase v0.3 |
| P2Y12 / Gi / cAMP / PKA modulation — negative regulation of activation | Not implemented | Phase v0.4 |
| GPVI signalling cascade — collagen-mediated activation pathway | Not implemented | (Dunster 2015 reference in pile) |
| Granule release (dense and α-granules) | Not implemented (cargo molecules in inventory but no exocytosis kinetics) | Phase v0.3 |
| Integrin αIIbβ3 inside-out signalling | Not implemented | Phase v0.6 (optional) |
| cAMP / cGMP / NO suppression at rest | Not implemented (Kleppe 2018 reference in pile) | Future |

## Design philosophy

Reuse a validated published model where one exists; deviate only when
the published model is genuinely incompatible with the framework or
the biology. Nine deviations from primary sources are catalogued in
design doc §6.8, each with the deviation's value, the reason, and a
pointer to the lab-book entry that diagnosed it. Examples: SERCA
E1 / E1·Ca initial counts pre-equilibrated for binding (D5);
γ_SOC calibrated against rest balance vs the Hoover face-value
single-channel conductance (D3); a basal plasma-membrane Ca²⁺ leak
(D4) added because Dolan 2014 has no PM leak term and the rest
balance requires one. The dissertation lead is *"reproduce Dolan
where possible, deviate honestly where necessary, document every
deviation."*

## Primary references

- **Dolan & Diamond 2014** *Biophys. J.* 106:2049–2060 — the validation target
- **Purvis & Bhatt 2008** *Plos Comp Biol* 4:e1000050 — kinetic parameterisation foundation
- **Caride et al. 2007** *J Biol Chem* — PMCA isoform 4b 5-state scheme
- **Sneyd & Dufour 2002** *PNAS* — IP3R 6-state Markov model
- **Hoover & Lewis 2011** *PNAS* — Orai/STIM CRAC channel MWC framework
- **Burkhart et al. 2012** *Blood* — platelet proteome reference (copy numbers)
- All in `source-info/calcium-papers/`; per-value provenance in `reports/data/calcium-data-provenance.md`.
