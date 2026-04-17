# Platelet WCM Design Notes — Summary for Lab Book

**Date:** 17 April 2026. **Branch:** `plan/calcium-signalling-pathway`.

## Overview

A set of design notes in `reports/` now defines the architecture for adapting
the wcEcoli whole-cell model framework to simulate the human platelet. The
goal is a first-generation non-spatial WCM of a resting/activated platelet,
built around signalling, calcium dynamics, and granule secretion — reusing
the wcEcoli engine without inheriting its bacterial reconstruction content.

## The design notes

### Background and orientation

- **`wcecoli-notes.md`** — reading notes on the wcEcoli codebase: three-layer
  mental model (framework / E. coli model / reconstruction), first files to
  open, simulation loop, major processes (transcription, translation, FBA
  metabolism, replication, division), output format, analysis plots. Acts as
  the orientation reference for anyone joining the work.
- **`parca-module-notes.md`** — focused explainer for the E. coli parameter
  calculator: what ParCa does, its 8-stage fitting pipeline, the
  `SimulationDataEcoli` object, and a clean split between generic (reusable)
  and E. coli-specific (not reusable) elements.

### Platelet design

- **`reconstruction-platelet-design.md`** — concrete proposal for a
  `reconstruction/platelet/` namespace. Reuses the ParCa staging pattern but
  reframes the question from *"viable growing E. coli"* to *"resting or
  activated platelet under signalling, metabolic, and compartmental
  conditions"*. Includes directory layout, flat-file inputs, dataclass
  structure, a 9-stage platelet ParCa pipeline, and a 4-phase implementation
  roadmap.
- **`calcium-signalling-pathway-design.md`** — concrete 3-layer design for
  the ADP → Ca²⁺ signalling process: raw TSV reactions, a sympy/scipy ODE
  dataclass (pattern borrowed from `TwoComponentSystem`), and a runtime
  `Process` class. Includes the full species list for the P2Y1/Gq, P2Y12/Gi,
  and SOCE pathways, across compartments `[c]`, `[m]`, `[dts]`, `[e]`.
  Validation target: reproducing the Dolan & Diamond 2014 Ca²⁺ transient
  shape.

### Literature and data

- **`platelet-proteome-data-sources.md`** — review of 14 quantitative
  proteomics datasets for initial inventory: Burkhart 2012 as the backbone,
  Zeiler 2014 for recalibration, Aslan 2021 for compartment assignments,
  granule sub-proteomes (Maynard, Hernandez-Ruiz, Melo), phosphoproteomics
  (Zahedi, Unsworth, Babur). Resolves known conflicts (αIIbβ3 copy number,
  species differences) and defines a tier structure for initial-condition
  assembly.
- **`zotero-literature-analysis.md`** — synthesis of ~90 papers across four
  Zotero collections, organised by platelet subsystem. High-priority
  additions: Kleppe 2018 (NO/cGMP/cAMP inhibitory pathway ODE), Dunster 2015
  (GPVI→Syk→PLCγ2 ODE), Mazet 2020 (PI cycle ODE), Thomas 2014
  (constraint-based metabolic reconstruction, iAT-PLT-636), Ghatge 2026
  (mitochondrial Ca²⁺), Flora 2023 (aerobic glycolysis switch).

### Feasibility study

- **`platelet-wcm-feasibility-study.pdf` / `.tex`** — formal synthesis
  document combining rationale, framework analysis, and implementation plan
  into a single citable report.

## Key technical decisions

- **Reuse wcEcoli framework, fork reconstruction and model namespaces.**
  Create `reconstruction/platelet/` and `models/platelet/` rather than
  deforming the E. coli tree.
- **Organising axis shifts from growth/regulation to activation/
  compartments.** Replaces `growth_rate_dependent_parameters` with
  `activation_parameters`; replaces TF-condition fitting with
  agonist/activation-state fitting.
- **ODE kinetics (not FBA) for signalling.** Calcium and receptor cascades
  modelled as coupled ODE systems using the same sympy/`solve_ivp` pattern as
  `TwoComponentSystem`; metabolism formalism deferred (Thomas 2014
  constraint-based reconstruction is a candidate).
- **Single ODE system for the whole Ca²⁺ signalling cascade.** Tight coupling
  (IP3 → Ca²⁺ → SERCA → STIM1) makes splitting across processes artificial.
- **Deterministic/stochastic boundary at ~1000 copies per platelet.** ODE
  treatment for species above (covers all major signalling effectors);
  stochastic treatment for rarer species below ~500 copies.
- **Literature feeds structured flat TSVs, not inline constants.** Keeps the
  reconstruction inspectable and reproducible; provenance tracked at the
  table level.

## Current status

- **Engine work**: two fixes already merged — safe time-limit termination for
  non-dividing cells (#14), ATP partition robustness under missing molecules
  (#16/17); `SimulationDataPlatelet` skeleton landed (#17).
- **Active branch**: `plan/calcium-signalling-pathway` — drafting the calcium
  module design.
- **Next step**: start implementing the minimal `reconstruction/platelet/`
  skeleton described in the design doc (Phase 1 of the roadmap).

## Phased implementation plan

1. **Engine contract** — make wcEcoli engine work for a non-dividing cell
   (largely done)
2. **Platelet stub** — minimal `models/platelet/` + `reconstruction/platelet/`
   running through the existing web UI
3. **Calcium dynamics** — IP3R / SERCA / SOCE; validate Ca²⁺ transient
   against Dolan & Diamond 2014
4. **Receptor signalling** — P2Y1 / P2Y12 / GPVI upstream of calcium
5. **Granule exocytosis** — Ca²⁺-threshold-driven, using `UniqueMolecule`
   instances
6. **Metabolism and integrin activation** — stretch goals
