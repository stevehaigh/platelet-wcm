# Platelet runtime scaffold: what changed and why

## Summary

This change adds the first runnable `models/platelet/` scaffold to the
repository.

The goal is not to add real platelet biology yet. The goal is to prove that the
existing `wholecell/` simulation engine can host a platelet-specific runtime
namespace with its own:

- simulation class
- process namespace
- initial-conditions path
- tests

In practice, this is the smallest vertical slice that lets later calcium and
receptor-signalling work land in a clear place.

---

## What was added

### `models/platelet/`

New platelet runtime modules:

- `models/platelet/sim/simulation.py`
- `models/platelet/sim/initial_conditions.py`
- `models/platelet/processes/platelet_stub.py`
- empty package scaffolds for `listeners/` and `analysis/single/`

The key new object is `PlateletSimulation`, a `Simulation` subclass that uses
the generic engine states:

- `BulkMolecules`
- `UniqueMolecules`
- `LocalEnvironment`

It wires in one minimal process, `PlateletStub`, whose only job is to exercise
the normal request / partition / evolve / merge lifecycle without changing any
molecule counts.

### `reconstruction/platelet/`

New platelet reconstruction/runtime wiring:

- `reconstruction/platelet/dataclasses/process/`
- `reconstruction/platelet/initialization.py`

`SimulationDataPlatelet` now exposes a `process` namespace, so runtime code can
read platelet-specific process data the same way E. coli runtime code reads
`sim_data.process.*`.

The internal-state stub now also exposes `initial_counts` in addition to the
unit-aware `bulk_data` structure. That gives platelet initialisation code a
clear place to read starting counts from.

### Tests

New tests in `models/platelet/tests/sim/test_simulation.py` cover three things:

1. `SimulationDataPlatelet` exposes the new process namespace
2. `PlateletSimulation` initializes the platelet stub process successfully
3. a real one-step run preserves the placeholder counts through the full
   lifecycle

These tests are intentionally small and readable. They are there to prove the
runtime contract, not to validate biology.

---

## What is still placeholder code

Several parts of this scaffold are intentionally fake for now.

### Placeholder molecule inventory

`reconstruction/platelet/dataclasses/internal_state.py` still contains:

- `DUMMY_PROTEIN[c]`
- `DUMMY_LIPID[e]`

These are temporary placeholders. They exist only to prove that the platelet
simulation can initialize non-zero bulk state and step cleanly.

They do **not** represent the eventual platelet molecule inventory. That work is
deferred to the Burkhart/Purvis curation issues.

### Placeholder process

`models/platelet/processes/platelet_stub.py` is a no-op process.

It intentionally:

- builds a bulk-molecule view
- participates in `calculateRequest()`
- participates in `evolveState()`
- leaves counts unchanged

This is useful because it proves the runtime plumbing with almost no biological
complexity. It is a scaffold, not a model.

### Placeholder initial conditions

`models/platelet/sim/initial_conditions.py` and
`reconstruction/platelet/initialization.py` now load counts from
`sim_data.internal_state.bulk_molecules.initial_counts`.

That is the right interface shape, but the values are still placeholder values.
Real resting Ca2+, ATP/ADP, receptor counts, and granule cargo remain future
work.

---

## Where this replaces or avoids E. coli-specific code

The main design decision in this change is not "copy E. coli and rename it".
It is "reuse the generic engine, but stop inheriting E. coli assumptions where
they are misleading."

### Reused directly from the generic engine

The platelet scaffold currently reuses these `wholecell/` components unchanged:

- `wholecell.sim.simulation.Simulation`
- `wholecell.states.bulk_molecules.BulkMolecules`
- `wholecell.states.unique_molecules.UniqueMolecules`
- `wholecell.states.local_environment.LocalEnvironment`
- `wholecell.processes.process.Process`

This is the point of the scaffold: prove that the generic engine is already
usable for platelet runtime work.

### Replaced with platelet-specific code

The scaffold adds platelet-specific equivalents for the first runtime-facing
surfaces:

| E. coli surface | Platelet scaffold |
| --- | --- |
| `models/ecoli/sim/simulation.py` | `models/platelet/sim/simulation.py` |
| `models/ecoli/sim/initial_conditions.py` | `models/platelet/sim/initial_conditions.py` |
| `reconstruction/ecoli/dataclasses/process/process.py` | `reconstruction/platelet/dataclasses/process/process.py` |

This is the first step toward a separate platelet namespace rather than a set of
special cases inside the E. coli model.

### Deliberately not reused yet

The scaffold does **not** reuse `models/ecoli/listeners/mass.py`.

That listener assumes E. coli-specific:

- submass categories
- compartment IDs
- growth/division-oriented reporting

Using it unchanged would create confusing coupling at exactly the point where the
platelet scaffold is supposed to stay easy to reason about. For now,
`PlateletSimulation` keeps `_listenerClasses = ()` and relies only on the
default `EvaluationTime` listener that the base engine adds automatically.

This keeps the current change honest: it proves the platelet runtime hook points
without pretending that a platelet mass-reporting surface already exists.

---

## Why this structure is useful for the calcium-signalling work

The calcium-signalling design needs a place to land a real platelet runtime
process. This scaffold creates that place.

The next calcium-focused step can now add:

- a platelet process dataclass for calcium state
- a platelet runtime process that reads it
- small, targeted tests around that process

without first having to solve namespace layout, simulation construction, or
initialisation wiring again.

In short: this change does not solve platelet biology. It removes the framework
ambiguity that would otherwise make every later biology PR harder to review.
