# ParCa Module Notes

Reading notes on the `ParCa` (parameter calculator) module in the wcEcoli
repository. These notes document what ParCa does in the current E. coli
codebase; the forward-looking platelet adaptation lives in
`reconstruction-platelet-design.md`.

---

## Short version

`ParCa` is the parameter calculator. Its job is to turn raw biological inputs
into a fitted `sim_data` object that the simulator can run. Concretely:

- load the raw E. coli knowledge base
- build condition-dependent physiological assumptions
- fit regulatory and capacity-related quantities
- assemble process/state dataclasses into a single `SimulationDataEcoli`
- serialise the result

ParCa is the bridge between curated biology on disk and a runnable whole-cell
model.

---

## Where ParCa lives

Main entry points:

- `reconstruction/ecoli/fit_sim_data_1.py`
- `reconstruction/ecoli/simulation_data.py`
- `reconstruction/ecoli/knowledge_base_raw.py`
- `wholecell/fireworks/firetasks/parca.py`

Responsibilities split cleanly:

- `knowledge_base_raw.py` loads raw structured inputs from
  `reconstruction/ecoli/flat/`
- `fit_sim_data_1.py` orchestrates fitting and transformation
- `simulation_data.py` defines the assembled `SimulationDataEcoli` object
- `parca.py` wraps the pipeline in a task that writes the fitted outputs

---

## The fitting pipeline

`fitSimData_1()` runs a fixed sequence of stages, each wrapped in `@save_state`
so intermediate results can be checkpointed and reloaded during development:

1. **`initialize()`** — instantiates `SimulationDataEcoli` and wires in
   condition, compartment, and molecular-weight metadata, helper classes
   (`Constants`, `MoleculeIds`, `MoleculeGroups`, `GetterFunctions`),
   growth-rate parameters, and the process/state/relation dataclasses.
2. **`input_adjustments()`** — curated adjustments to translation
   efficiencies, RNA expression, and RNA/protein degradation rates before
   heavier fitting. This is where the code encodes the difference between
   "raw imported data" and "values the model should actually use".
3. **`basal_specs()`** — baseline cell specification for the basal condition,
   connecting to ppGpp regulation and capacity fitting for ribosomes and RNA
   polymerase.
4. **`tf_condition_specs()`** — expands fitting from basal to multiple
   transcription-factor conditions. ParCa's core role here is constructing a
   *family of internally consistent model states* across biological conditions.
5. **`fit_condition()`** — reconciles condition-specific specs with global
   constraints, with strong coupling to metabolism and growth-related resource
   allocation.
6. **`promoter_binding()`** and **`adjust_promoters()`** — fit and refine
   transcriptional regulation using fold-change and regulatory data to derive
   promoter binding probabilities.
7. **`set_conditions()`** — consolidates fitted condition-dependent
   quantities back into the simulation-data structure.
8. **`final_adjustments()`** — final corrections and E. coli-specific
   regulatory details (attenuation, ppGpp).

---

## Inputs and outputs

### Inputs

Defined by `reconstruction/ecoli/knowledge_base_raw.py`, which loads TSVs from
`reconstruction/ecoli/flat/`:

- genes, RNAs, proteins, metabolites, compartments
- metabolic reactions and concentrations
- transcription-factor and fold-change data
- growth-rate-dependent parameters
- condition and media definitions
- degradation rates and manual adjustments

ParCa is much more than a parser: it turns those heterogeneous inputs into a
single coherent model state.

### Outputs

A serialised `sim_data` object, written through
`wholecell/fireworks/firetasks/parca.py`. The simulator consumes this compiled
result directly, not the raw flat tables.

---

## The `SimulationDataEcoli` class

The best place to see what ParCa assembles is
`reconstruction/ecoli/simulation_data.py`. `SimulationDataEcoli.initialize()`
builds an object with:

- `constants`
- `adjustments`
- `molecule_ids`
- `molecule_groups`
- `getter`
- `growth_rate_parameters`
- `mass`
- `external_state`
- `process`
- `internal_state`
- `relation`

ParCa is not just fitting scalar parameters. It prepares the entire
organism-specific configuration the simulation engine depends on: which
molecules exist, which compartments exist, which groups and aliases exist, what
each process needs to know, and the cross-links between subsystems.

The key insight for a newcomer: **ParCa is the reconstruction layer for the
organism namespace.**

---

## Generic vs E. coli-specific

The ParCa **architecture** is reusable:

1. raw-data to simulation-data pipeline
2. intermediate stage checkpoints via `@save_state`
3. structured simulation-data assembly via dataclasses
4. clean separation: loader / fitter / consumer

Strongly E. coli-specific pieces:

1. growth-rate parameterisation
2. condition logic based on media and TF state
3. transcription-centric regulation (promoter fitting, fold-change,
   attenuation, ppGpp)
4. coupling to bacterial metabolism and biosynthesis
5. DNA/RNA/protein synthesis assumptions that do not map onto an anucleate cell

The implication for the platelet work is that a new
`reconstruction/platelet/` namespace should reuse the **pattern** (staged
pipeline, checkpointing, structured `SimulationData` assembly) without
inheriting the bacterial content. The concrete design for that namespace is in
`reconstruction-platelet-design.md`.
