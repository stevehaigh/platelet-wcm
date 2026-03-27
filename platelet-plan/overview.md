# Platelet Whole-Cell Model: Plan Overview

## 1. Big-picture feasibility

wcEcoli gives:
- A full modelling stack: reconstruction -> parameter calculation (ParCa) -> simulation -> analysis.
- A modular process/store architecture: metabolism, transcription, translation, division, etc. as separable processes.
- A tested simulation engine: time stepping, bookkeeping, output, plotting, workflows.

For platelets, the biology is very different, but the software abstractions (processes, stores, variants, analysis) are reusable. So "clone, run, prune, then rebuild" is conceptually sound -- as long as you're ruthless about what you keep.

---

## 2. Key biological mismatches to keep in mind

If you treat this as a literal "E. coli -> platelet" port, you'll fight the model. Better to treat wcEcoli as a framework template and design a new platelet model that just happens to live inside its infrastructure.

Some big differences to work on early:

- **Anucleate, no de novo transcription:**
  - No genome replication, no transcriptional regulation, no cell cycle.
  - mRNA and protein dynamics are mostly decay + translation from pre-loaded mRNA, if you include any translation at all.

- **Finite lifespan, no division:**
  - Platelets have a ~days-scale lifespan, no growth to division.
  - Time axis is more about resting -> primed -> activated -> spent than exponential growth.

- **Compartmentalisation:**
  - Plasma membrane, open canalicular system, dense granules, alpha-granules, mitochondria, cytosol, maybe microdomains.
  - You'll want explicit compartments and trafficking between them.

- **Activation and signalling:**
  - Receptors (GPVI, PARs, P2Y12, integrins, etc.), signalling cascades, Ca2+ dynamics.
  - This is where the "whole-cell" flavour really lives for platelets.

- **Granule release:**
  - Dense and alpha-granules as discrete, countable objects with state (loaded, fusing, released).
  - Exocytosis kinetics, content release into extracellular space, feedback on activation.

So: reuse the engine, not the E. coli ontology.

---

## 3. "Clone -> run -> prune" plan

### Phase 1: Treat wcEcoli as a black-box engine

1. Clone and run a minimal sim:
   - Get a single standard E. coli simulation running end-to-end with the manual runscripts.
   - Confirm you understand where:
     - Parameters live (reconstruction/ParCa).
     - Processes are defined (models/ecoli/processes).
     - Stores/state are defined (bulk vs unique molecules, etc.).

2. Map the architecture:
   - Sketch a quick diagram of:
     - Processes (what updates what).
     - Stores (what state exists where).
     - Top-level simulation object (how processes and stores are wired).

**Deliverable:** a 1-2 page "mental model" of wcEcoli as a generic whole-cell engine.

### Phase 2: Create a "blank cell" inside the same framework

Instead of hacking E. coli into a platelet:

1. Fork a new organism namespace:
   - e.g. `models/platelet/...` mirroring the E. coli layout but initially almost empty.
   - Minimal reconstruction: just enough to define:
     - Compartments.
     - A few molecule types (e.g. Ca2+, ADP, fibrinogen, a generic "granule cargo").

2. Define a minimal platelet sim object:
   - One or two toy processes, e.g.:
     - Resting turnover: slow leakage/decay of some species.
     - Simple activation switch: at time t_0, flip a "activated" flag and increase Ca2+.
   - Wire them into the existing simulation engine (time stepping, logging, output).

3. Strip dependencies:
   - Make sure this platelet sim can run without any E. coli reconstruction or ParCa.
   - You want a clean separation: E. coli is just another "organism" using the same engine.

**Deliverable:** a simple but runnable "platelet stub model" that proves the framework is reusable.

### Phase 3: Add platelet-specific biology incrementally

Now you can layer in the interesting stuff:

1. Compartment and state design:
   - Define stores for:
     - Compartments: membrane, cytosol, granules, extracellular.
     - Granules: as unique objects with attributes (type, cargo, location, state).
     - Signalling state: receptor occupancy, Ca2+, activation flags.

2. Granule release as a first serious process:
   - Start with a simple rule-based process:
     - If activated and Ca2+ above threshold -> granule fusion events at a given rate.
     - On fusion: move cargo from granule to extracellular store, change granule state.
   - This gives a clear, mechanistic, platelet-specific behaviour that's easy to visualise.

3. Activation pathways:
   - Add one or two minimal receptor -> Ca2+ pathways (even if heavily coarse-grained).
   - You can always refine later with more detailed signalling modules.

4. Optional: metabolism and energy:
   - If you care about ATP usage, glycolysis, etc., you can port or rewrite a very reduced metabolic module, but keep this late.

**Deliverable:** a first-pass platelet model where you can simulate activation and granule release over time and inspect trajectories.

---

## 4. Why start with granule release?

Given your comparative-modelling brain, granule release is a nice "anchor":
- **Discrete, countable entities:** fits well with unique-molecule stores and causal tracing.
- **Clear experimental observables:** secretion profiles, surface markers, etc.
- **Bridges scales:** receptor signalling -> Ca2+ -> exocytosis -> extracellular milieu.

You can then hang other modules off it: cytoskeleton reorganisation, integrin activation, procoagulant surface, etc.

---

## 5. Review notes and refinements

### Things to watch out for

1. **Cell division is tightly coupled to E. coli.** `wholecell/sim/divide_cell.py` hardcodes chromosome domain division, RNA-RNAP coupling, ribosome-mRNA coupling. Since platelets don't divide, you can simply not wire up a division function -- but verify the engine doesn't assume division will eventually happen.

2. **UniqueMolecules state has E. coli division modes baked in.** The division logic for unique molecules assumes domains, chromosomal segments, RNA, ribosomes. May need to subclass or simplify `UniqueMolecules` for platelet-specific unique objects (granules, receptor complexes).

3. **The reconstruction layer is the hardest part.** wcEcoli's `SimulationDataEcoli` is massive. You'll need a `SimulationDataPlatelet` equivalent, even if much smaller. Plan for this explicitly in Phase 2; define what data sources you'll use (literature values, Reactome/UniProt for platelet proteome, etc.).

4. **Time axis needs rethinking.** wcEcoli's simulation loop is built around exponential growth toward division. For platelets (resting -> primed -> activated -> spent over hours/days), modify how the time axis and "end condition" work.

### Suggested additions

- **Phase 1.5: Identify engine modifications needed.** Simulation termination condition, time-step adaptation logic, and the division hook. Minor changes but worth cataloging early.

- **Phase 2: Be explicit about the data pipeline.** Where do platelet parameters come from? Even a toy model needs Ca2+ resting concentrations, granule counts (~40-80 dense, ~40-80 alpha per platelet), receptor densities, etc. A small `reconstruction/platelet/` with curated literature values prevents magic numbers.

- **Phase 3: Consider Ca2+ dynamics as your second process** (right after or alongside granule release). Ca2+ is the central integrator of platelet activation -- it connects receptors to granule release, integrin activation, and shape change.

### v0.1 target

Define a concrete "v0.1 platelet model": simulate thrombin-induced Ca2+ rise and dense granule release over 5 minutes, plot secretion kinetics. A falsifiable deliverable before expanding scope.
