# Design for `reconstruction/platelet/`

This note proposes a concrete design for a future `reconstruction/platelet/`
namespace.

The goal is not to clone `reconstruction/ecoli/` file-for-file. The goal is to
keep the useful reconstruction architecture from wcEcoli while replacing the
parts that are tightly tied to a growing, dividing bacterium.

The design here assumes the broader platelet plan already captured elsewhere in
the repository:

- `reports/lab-book-2026-03-27.md`
- `reports/zotero-literature-analysis.md`
- `reports/parca-module-notes.md`

---

## Design goals

`reconstruction/platelet/` should do four things well:

1. turn literature-derived platelet data into structured, versioned flat inputs
2. compile those inputs into a `SimulationDataPlatelet` object
3. support a staged parameter-calculation workflow like ParCa
4. stay small and explicit enough that the first working platelet model can be
   built incrementally

The first version should optimise for **clarity and evolvability**, not for
biological completeness.

---

## The framing question

The current E. coli ParCa answers:

> "What fitted parameter set defines a viable E. coli cell under these growth
> and regulatory conditions?"

The platelet equivalent should answer:

> "What fitted parameter set defines a resting or activated platelet under
> these signalling, metabolic, and compartmental conditions?"

That single shift — from growth-and-regulation to activation-and-compartments
— drives most of the design decisions below.

---

## High-level role in the repository

The platelet reconstruction layer should occupy the same architectural role that
`reconstruction/ecoli/` does today:

- `reconstruction/platelet/` defines the platelet knowledge base and fitted
  simulation data
- `models/platelet/` defines the organism-specific runtime processes and
  listeners
- `wholecell/` remains the generic simulation engine

So the responsibility split would be:

- **reconstruction**: "what does a platelet contain and how is it parameterised?"
- **model/runtime**: "how does that platelet evolve over time in simulation?"

---

## Proposed directory layout

The first useful target is:

```text
reconstruction/platelet/
├── __init__.py
├── knowledge_base_raw.py
├── simulation_data.py
├── fit_sim_data_1.py
├── initialization.py
├── flat/
│   ├── README.md
│   ├── compartments.tsv
│   ├── proteins.tsv
│   ├── metabolites.tsv
│   ├── receptors.tsv
│   ├── receptor_binding_parameters.tsv
│   ├── calcium_kinetics.tsv
│   ├── phosphoproteins.tsv
│   ├── phosphorylation_states.tsv
│   ├── activation_states.tsv
│   ├── granules.tsv
│   ├── granule_cargo.tsv
│   ├── metabolic_reactions.tsv
│   ├── metabolite_concentrations.tsv
│   ├── metabolism_kinetics.tsv
│   ├── molecule_groups.tsv
│   ├── molecule_ids.tsv
│   ├── constants.tsv
│   ├── parameters.tsv
│   └── adjustments/
│       ├── protein_count_adjustments.tsv
│       ├── metabolite_adjustments.tsv
│       ├── receptor_adjustments.tsv
│       └── calcium_adjustments.tsv
├── dataclasses/
│   ├── __init__.py
│   ├── constants.py
│   ├── common_names.py
│   ├── getter_functions.py
│   ├── molecule_groups.py
│   ├── molecule_ids.py
│   ├── adjustments.py
│   ├── activation_parameters.py
│   ├── relation.py
│   ├── process/
│   │   ├── __init__.py
│   │   ├── process.py
│   │   ├── calcium_signaling.py
│   │   ├── receptor_signaling.py
│   │   ├── granule_exocytosis.py
│   │   ├── metabolism.py
│   │   ├── membrane_transport.py
│   │   └── integrin_activation.py
│   └── state/
│       ├── __init__.py
│       ├── internal_state.py
│       ├── external_state.py
│       ├── bulkMolecules.py
│       ├── uniqueMolecules.py
│       └── stateFunctions.py
└── scripts/
    ├── extract_literature_tables.py
    ├── build_proteome_inputs.py
    ├── build_calcium_inputs.py
    └── build_granule_inputs.py
```

This is deliberately similar to `reconstruction/ecoli/`, but smaller and more
focused.

---

## What each top-level file should do

### `knowledge_base_raw.py`

This should be the platelet equivalent of
`reconstruction/ecoli/knowledge_base_raw.py`.

Its job should be to:

- load all TSV-backed platelet inputs from `flat/`
- expose them as structured attributes on a raw knowledge-base object
- keep the loading logic simple and deterministic
- avoid performing major fitting or biological inference

This file should answer:

> "What curated platelet facts do we have on disk?"

not:

> "What fitted platelet model should we run?"

### `simulation_data.py`

This should define `SimulationDataPlatelet`.

It should play the same role as `SimulationDataEcoli`, but without bacterial
concepts like doubling time or growth-rate parameterisation as the organising
axis.

Its `initialize()` method should probably assemble:

- constants
- adjustments
- molecule IDs and groups
- getter functions
- activation parameters
- process dataclasses
- internal and external state dataclasses
- relations between those subsystems

The key difference from E. coli is that this object should be organised around
**activation state and compartmental state**, not growth state.

### `fit_sim_data_1.py`

This should be the platelet ParCa entry point.

It should keep the same staged, checkpointable style as the E. coli ParCa, but
with platelet-oriented stages. A first pass could look like:

1. `initialize()`
2. `input_adjustments()`
3. `resting_state_specs()`
4. `activation_state_specs()`
5. `fit_receptor_signaling()`
6. `fit_calcium_state()`
7. `initialize_granules()`
8. `set_conditions()`
9. `final_adjustments()`

The design intent is:

- preserve the useful workflow shape
- remove promoter-centric and growth-law-centric logic
- let more biology be added later without rearranging the whole pipeline

### `initialization.py`

This should hold reusable functions that convert `SimulationDataPlatelet` into
initial bulk and unique molecule containers.

It is a good place for:

- bulk molecule initialisation
- unique granule object initialisation
- compartment mass and volume initialisation
- consistent initial activation-state loading

This is one of the places where reuse from `reconstruction/ecoli/initialization.py`
is likely to be valuable.

---

## Proposed `flat/` inputs

The flat directory should be intentionally smaller than the E. coli one at
first. It should start from the minimum data needed for a stub platelet model
and then grow by module.

### Core identity tables

- `compartments.tsv`
- `proteins.tsv`
- `metabolites.tsv`
- `molecule_groups.tsv`
- `molecule_ids.tsv`
- `constants.tsv`
- `parameters.tsv`

These establish what exists in the platelet model at all.

### Signalling tables

- `receptors.tsv`
- `receptor_binding_parameters.tsv`
- `calcium_kinetics.tsv`
- `phosphoproteins.tsv`
- `phosphorylation_states.tsv`
- `activation_states.tsv`

These replace a large fraction of the bacterial transcription/regulation inputs.

### Secretory and compartment tables

- `granules.tsv`
- `granule_cargo.tsv`

These are especially important because platelet granules are a major biological
feature that does not exist in the E. coli model.

### Metabolism tables

- `metabolic_reactions.tsv`
- `metabolite_concentrations.tsv`
- `metabolism_kinetics.tsv`

Whether metabolism is handled by ODEs, a reduced network, or a constraint-based
formalism can be decided later. The flat-file interface should not assume the
answer too early.

### Adjustment tables

Keep an `adjustments/` subdirectory from the start.

This preserves a useful wcEcoli pattern: separate the imported literature value
from the value the model currently uses after curation and calibration.

---

## Proposed dataclass structure

The platelet dataclasses should mirror the useful parts of
`reconstruction/ecoli/dataclasses/` while dropping the growth-law stack.

### Keep or adapt directly

These concepts are still useful:

- `constants.py`
- `common_names.py`
- `getter_functions.py`
- `molecule_groups.py`
- `molecule_ids.py`
- `adjustments.py`
- `relation.py`
- `state/`

The reason is simple: every cell model still needs names, IDs, grouped molecule
sets, lookup helpers, state containers, and cross-links.

### Replace with platelet-specific logic

Instead of `growth_rate_dependent_parameters.py`, use something like:

- `activation_parameters.py`

This module should own quantities that vary with platelet state, such as:

- resting versus activated parameter sets
- receptor occupancy presets
- calcium pool initial conditions
- phosphorylation-state presets
- granule release thresholds

This change matters because it moves the organising axis from:

- **doubling time / media / TF state**

to:

- **agonist / activation / compartment state**

### Proposed `process/` dataclasses

Start with dataclasses for the processes you actually expect to simulate early:

- `calcium_signaling.py`
- `receptor_signaling.py`
- `granule_exocytosis.py`
- `metabolism.py`
- `membrane_transport.py`
- `integrin_activation.py`

Not all of these need rich logic at the beginning. Some can start as thin data
holders while the matching runtime processes are still stubs.

### Proposed `state/` dataclasses

The state dataclass layout can remain close to E. coli:

- `internal_state.py`
- `external_state.py`
- `bulkMolecules.py`
- `uniqueMolecules.py`
- `stateFunctions.py`

That is useful because the generic engine already understands these concepts.

The main platelet-specific change is the content:

- bulk molecule pools should include cytosolic signalling species and metabolites
- unique molecules should probably include granules and possibly receptor
  clusters or membrane microdomains later
- compartment definitions should include at least cytosol, membrane, dense
  granules, alpha granules, open canalicular system, and mitochondria if the
  biology requires them

---

## Recommended platelet ParCa stages

The platelet `fit_sim_data_1.py` should not try to preserve E. coli stage names
that no longer make sense.

A good first design is:

### 1. `initialize()`

- instantiate `SimulationDataPlatelet`
- load raw data
- build IDs, groups, constants, and compartments

### 2. `input_adjustments()`

- apply curated overrides
- normalise units
- resolve conflicts across literature sources

### 3. `resting_state_specs()`

- compute the baseline non-activated platelet
- set resting molecule counts and compartment occupancies

### 4. `activation_state_specs()`

- derive named activation-state parameter sets
- set phosphorylation fractions, receptor states, and calcium-pool presets

### 5. `fit_receptor_signaling()`

- compile receptor abundances and binding parameters
- derive a consistent parameterisation for the first signalling layer

### 6. `fit_calcium_state()`

- compile calcium-handling parameters
- derive initial pool sizes, exchange rates, and buffering assumptions

### 7. `initialize_granules()`

- convert granule metadata and cargo tables into bulk and unique state
  initialisation data

### 8. `set_conditions()`

- consolidate resting and activation-state data into runtime-facing structures

### 9. `final_adjustments()`

- run final consistency checks
- freeze derived lookup structures
- write any metrics useful for debugging reconstruction quality

This is close enough to wcEcoli ParCa to feel familiar, but different enough to
fit platelet biology honestly.

---

## Suggested relationship to the literature workflow

The Zotero-driven literature work should feed `reconstruction/platelet/` in two
steps:

1. **curation into flat tables**
2. **compilation by platelet ParCa**

That means the primary workflow should be:

1. identify a paper via Zotero
2. extract quantitative content into a tracked flat TSV or generated table
3. document provenance in the table or its build script
4. let `knowledge_base_raw.py` and `fit_sim_data_1.py` consume the structured
   result

This is better than hard-coding paper-derived constants inside Python files,
because it keeps the reconstruction layer inspectable and reproducible.

---

## How to stage implementation

The directory should be built in phases.

### Phase 1: minimal runnable skeleton

Create:

- `reconstruction/platelet/__init__.py`
- `reconstruction/platelet/knowledge_base_raw.py`
- `reconstruction/platelet/simulation_data.py`
- `reconstruction/platelet/fit_sim_data_1.py`
- `reconstruction/platelet/initialization.py`
- `reconstruction/platelet/flat/` with a very small set of TSVs
- `reconstruction/platelet/dataclasses/` with only the modules needed for a stub

Goal:

- produce a valid `sim_data` object for a minimal platelet namespace

### Phase 2: calcium-first biology

Add:

- `calcium_kinetics.tsv`
- `activation_states.tsv`
- `dataclasses/process/calcium_signaling.py`

Goal:

- support a single calcium-focused platelet process and its initial state

### Phase 3: secretion and activation expansion

Add:

- `granules.tsv`
- `granule_cargo.tsv`
- receptor signalling tables
- granule and receptor process dataclasses

Goal:

- move from a signalling stub to a recognisably platelet-specific model

### Phase 4: metabolism and richer condition handling

Add:

- metabolic tables
- metabolism dataclass
- more explicit agonist and inhibitor condition modelling

Goal:

- support energy constraints and richer activation scenarios

---

## What should be copied, generalized, or left behind

### Copy or adapt

- checkpointing pattern from `reconstruction/ecoli/fit_sim_data_1.py`
- raw-data loading pattern from `knowledge_base_raw.py`
- simulation-data assembly pattern from `simulation_data.py`
- bulk/unique initialisation patterns from `initialization.py`
- state dataclass structure from `dataclasses/state/`

### Generalize carefully

- `getter_functions.py`
- `molecule_ids.py`
- `molecule_groups.py`
- `relation.py`

These should be generalized only where there is clear shared logic, not because
the files happen to have the same names in E. coli.

### Leave behind

- growth-rate-dependent parameters
- TF condition fitting
- promoter occupancy fitting
- ppGpp and transcription attenuation logic
- replication and cell-cycle reconstruction logic

Those are not "temporary gaps" for a platelet model; they are the wrong
organising assumptions.

---

## Naming and modelling conventions

A few conventions would make this namespace easier to grow.

### Use state-oriented names

Prefer names like:

- `activation_states.tsv`
- `activation_parameters.py`
- `resting_state_specs()`

over bacterial legacy terms that imply gene-regulatory fitting.

### Keep provenance close to data

Every flat table or build script should clearly indicate its literature source.
That matters more for platelet reconstruction than for code elegance, because the
real difficulty will be data trustworthiness and interpretation.

### Avoid overcommitting to one mathematical formalism

The flat inputs and dataclasses should not force everything into one modelling
style too early.

Platelet subsystems may eventually mix:

- ODE-based calcium and signalling kinetics
- threshold or event-based granule release logic
- reduced metabolic dynamics

The reconstruction layer should support that pluralism.

---

## A minimal first tree

If the full proposed tree feels too large, this is the smallest version worth
creating first:

```text
reconstruction/platelet/
├── __init__.py
├── knowledge_base_raw.py
├── simulation_data.py
├── fit_sim_data_1.py
├── initialization.py
├── flat/
│   ├── compartments.tsv
│   ├── proteins.tsv
│   ├── metabolites.tsv
│   ├── activation_states.tsv
│   ├── calcium_kinetics.tsv
│   ├── granules.tsv
│   └── parameters.tsv
└── dataclasses/
    ├── __init__.py
    ├── constants.py
    ├── getter_functions.py
    ├── molecule_ids.py
    ├── molecule_groups.py
    ├── activation_parameters.py
    ├── process/
    │   ├── __init__.py
    │   ├── process.py
    │   └── calcium_signaling.py
    └── state/
        ├── __init__.py
        ├── internal_state.py
        ├── external_state.py
        ├── bulkMolecules.py
        └── uniqueMolecules.py
```

That is enough to support a first stub plus a calcium-focused early milestone.

---

## Bottom line

`reconstruction/platelet/` should be designed as a **platelet-specific
reconstruction namespace** that reuses wcEcoli's reconstruction architecture
without inheriting its bacterial assumptions.

The best design is:

- structurally similar to `reconstruction/ecoli/`
- smaller at the start
- centred on activation state, signalling, calcium, granules, and compartments
- fed by flat literature-derived inputs
- able to grow incrementally from a stub into a richer platelet model

That gives the platelet work a clean home in the repository and keeps the
boundary between framework, reconstruction, and runtime model clear.
