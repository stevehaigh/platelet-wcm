# wcEcoli Architecture Map

A detailed map of the wcEcoli simulation framework, annotated with reusability
assessment for building a platelet whole-cell model.

---

## 1. High-Level Architecture

```
                     +----------------------------+
                     |   Reconstruction / ParCa   |
                     |  (reconstruction/ecoli/)    |
                     |                             |
                     |  raw_data -> SimulationData |
                     +-------------+--------------+
                                   |
                                   | sim_data
                                   v
+----------+       +--------------+---------------+       +----------+
|          |       |      Simulation Engine        |       |          |
| Loggers  |<----->|   (wholecell/sim/simulation)  |<----->|  Hooks   |
|          |       |                               |       |          |
+----------+       +--+--------+--------+--------++       +----------+
                      |        |        |        |
                      v        v        v        v
                  +------+ +------+ +------+ +------+
                  |State1| |State2| |Proc.1| |Proc.N|
                  |      | |      | |      | |      |
                  +------+ +------+ +------+ +------+
                      ^        ^        |        |
                      |        |        v        v
                      +--------+-----+----------+
                               |
                               v
                          +---------+
                          |Listeners|
                          +---------+
```

The framework has four layers:
1. **Reconstruction** -- builds `sim_data` (the knowledge base) from raw data
2. **Simulation engine** -- generic time-stepping loop
3. **Organism model** -- processes, states, listeners, initial conditions (E. coli-specific)
4. **Analysis** -- post-simulation plotting and data extraction

---

## 2. Simulation Engine Core

### File: `wholecell/sim/simulation.py` (534 lines)

The `Simulation` base class is **fully generic**. A subclass must define:

```python
class Simulation:
    _definedBySubclass = (
        "_internalStateClasses",    # tuple of InternalState subclasses
        "_externalStateClasses",    # tuple of ExternalState subclasses
        "_processClasses",          # nested tuple: ((ProcA, ProcB), (ProcC,), ...)
        "_initialConditionsFunction",  # callable(sim_data) -> sets initial state
    )

    # Optional overrides:
    _listenerClasses = ()
    _hookClasses = ()
    _divideCellFunction = None  # called on cell cycle completion
    _shellColumnHeaders = ("Time (s)",)
```

### Simulation Loop (per time step)

```
_pre_evolve_state()
    |-- _adjustTimeStep()           # binary search for valid dt
    |-- hook.preEvolveState()
    |-- state.reset_process_mass_diffs()
    |
    for each process_group in _processClasses:
        _evolveState(process_group)
            |-- state.updateQueries()       # refresh views
            |-- process.calculateRequest()  # each process declares needs
            |-- state.partition(processes)   # allocate molecules to processes
            |-- process.evolveState()       # execute biology
            |-- process.wasTimeStepShortEnough()  # validate
            |-- state.merge(processes)      # reconcile state
            |-- external_state.update()     # media timeline
    |
_post_evolve_state()
    |-- state.calculateMass()
    |-- listener.update()
    |-- hook.postEvolveState()
    |-- logger.append()
```

### Key design decisions:
- **Process groups execute sequentially**; processes within a group share the same partition step. This enforces ordering dependencies (e.g., TfUnbinding before TfBinding).
- **Time step adaptation**: each process implements `isTimeStepShortEnough()`. The engine binary-searches for the largest valid dt across all processes.
- **Cell cycle termination**: `cellCycleComplete()` sets a flag; `finalize()` calls `_divideCellFunction()`. For platelets, you can simply never set this flag and use `_lengthSec` as the end condition.

### Reusability: FULLY GENERIC
No organism-specific logic. Can be subclassed directly for a platelet model.

---

## 3. State Classes

### 3.1 Base Classes (Generic)

**`wholecell/states/internal_state.py`** -- `InternalState`
- Tracks mass arrays indexed by `sim_data.submass_name_to_index` and `sim_data.compartment_abbrev_to_index`
- Defines interface: `initialize()`, `allocate()`, `partition()`, `merge()`, `calculateMass()`
- Maintains per-process mass diffs for bookkeeping

**`wholecell/states/external_state.py`** -- `ExternalState`
- Minimal interface: `initialize(sim, sim_data, timeline)`, `update()`

### 3.2 BulkMolecules (`wholecell/states/bulk_molecules.py`)

Tracks integer copy numbers of molecules in a flat array.

**How it works:**
1. `initialize()`: reads molecule IDs and masses from `sim_data.internal_state.bulk_molecules.bulk_data`
2. Creates a `BulkObjectsContainer` (name-indexed numpy array)
3. Each process gets molecules via `BulkMoleculesView` / `BulkMoleculeView`
4. Partitioning: priority-based fractional allocation when demand exceeds supply
5. Division modes loaded from `sim_data.molecule_groups`

**sim_data interface required:**
```python
sim_data.internal_state.bulk_molecules.bulk_data['id']    # array of molecule ID strings
sim_data.internal_state.bulk_molecules.bulk_data['mass']   # mass per molecule (with units)
sim_data.constants.n_avogadro                              # Avogadro's number
sim_data.submass_name_to_index                             # e.g. {'protein': 0, 'RNA': 1, ...}
sim_data.compartment_abbrev_to_index                       # e.g. {'c': 0, 'p': 1, 'e': 2, ...}
sim_data.molecule_groups.bulk_molecules_binomial_division   # list of molecule IDs
sim_data.molecule_groups.bulk_molecules_equal_division      # list of molecule IDs
```

**Partitioning algorithm** (`calculatePartition()`):
- Processes declare requests via views
- If total request <= supply: everyone gets what they asked for
- If total request > supply: proportional allocation with stochastic rounding
- Priority levels allow some processes to be served first

**Reusability: HIGH**
The state itself is generic. You just need to provide the right `sim_data` structure. The `tableAppend` method has one E. coli-specific hardcoded reference (`ATP[c]` index) that would need to be overridden or generalized.

### 3.3 UniqueMolecules (`wholecell/states/unique_molecules.py`)

Tracks individual molecules with per-instance attributes (structured numpy arrays).

**How it works:**
1. `initialize()`: reads molecule definitions from `sim_data.internal_state.unique_molecule.unique_molecule_definitions`
   - This is a dict: `{molecule_name: {attr_name: dtype, ...}, ...}`
   - E. coli example: `{'active_RNAP': {'domain_index': np.int64, 'coordinates': np.int64, ...}}`
2. Adds `massDiff_*` attributes automatically for mass tracking
3. Creates a `UniqueObjectsContainer` with these definitions
4. Processes interact via `UniqueMoleculesView` with explicit access control (read/edit/delete)

**Division modes** (E. coli-specific, loaded in `initialize()`):
```python
self.division_mode['domain_index']          # molecules divided by chromosome domain
self.division_mode['RNA']                   # full mRNAs binomial, partial follow RNAP
self.division_mode['active_ribosome']       # follow bound mRNA
self.division_mode['chromosomal_segment']   # follow chromosome domain
```

**sim_data interface required:**
```python
sim_data.internal_state.unique_molecule.unique_molecule_definitions  # {name: {attr: dtype}}
sim_data.internal_state.unique_molecule.unique_molecule_masses       # structured array
sim_data.submass_name_to_index
sim_data.compartment_abbrev_to_index
sim_data.molecule_groups.unique_molecules_*_division   # lists for each division mode
```

**Reusability: MEDIUM**
- The container infrastructure is fully generic
- Division mode loading (lines 96-99) is hardcoded for E. coli's 4 modes
- For platelets: subclass and define your own division modes (or none, since platelets don't divide)
- `calculateMass()` hardcodes compartment 'c' for all unique molecules -- would need generalization for multi-compartment platelet model

### 3.4 LocalEnvironment (`wholecell/states/local_environment.py`)

Tracks external molecule concentrations with media timeline support.

**How it works:**
- Stores concentrations (not counts) in a `BulkObjectsContainer`
- Supports time-series media shifts (e.g., nutrient depletion experiments)
- Provides exchange data for metabolism (import constraints)

**sim_data interface required:**
```python
sim_data.external_state.make_media              # media construction helper
sim_data.external_state.current_timeline_id     # optional preset timeline
sim_data.external_state.saved_timelines         # {id: [(time, media_id), ...]}
sim_data.external_state.saved_media             # {media_id: {mol_id: concentration}}
sim_data.external_state.exchange_data_from_concentrations  # callable
sim_data.external_state.exchange_to_env_map     # internal -> env molecule mapping
sim_data.external_state.import_constraint_threshold
```

**Reusability: MEDIUM-HIGH**
- The media timeline concept maps well to platelet environments (e.g., agonist addition at time t)
- Exchange/import constraint logic is metabolism-specific; may not be needed initially

---

## 4. Process Architecture

### Base Class: `wholecell/processes/process.py` (185 lines)

**Fully generic.** Defines the interface every process must implement:

```python
class Process:
    def initialize(self, sim, sim_data):   # setup, create views
        pass
    def calculateRequest(self):            # declare molecule needs
        pass
    def evolveState(self):                 # execute biology for one dt
        pass
    def isTimeStepShortEnough(self, dt, safety):  # pre-check
        return True
    def wasTimeStepShortEnough(self):      # post-check
        return True
```

**View constructors** (called during `initialize()`):
```python
self.bulkMoleculesView(moleculeIDs)   # returns BulkMoleculesView
self.bulkMoleculeView(moleculeID)     # returns BulkMoleculeView (single)
self.uniqueMoleculesView(molName)     # returns UniqueMoleculesView
self.environmentView(moleculeIDs)     # returns EnvironmentView
```

**Listener communication:**
```python
self.writeToListener(listenerName, attrName, value)
self.readFromListener(listenerName, attrName)
```

### E. coli Process Inventory (15 processes, 8 groups)

```
Group 1: TfUnbinding
Group 2: Equilibrium, TwoComponentSystem, RnaMaturation
Group 3: TfBinding
Group 4: TranscriptInitiation, PolypeptideInitiation, ChromosomeReplication,
         ProteinDegradation, RnaDegradation, Complexation
Group 5: TranscriptElongation, PolypeptideElongation
Group 6: ChromosomeStructure
Group 7: Metabolism
Group 8: CellDivision
```

**None of these are reusable for platelets.** They all assume E. coli-specific molecules, reactions, and biology. They serve only as examples of how to write a process.

---

## 5. View System

### File: `wholecell/views/view.py`

Views are the interface between processes and states. A View:
1. Holds a reference to a state and a process index
2. Queries the state for current totals (`_updateQuery()`)
3. Accepts requests from the process (`requestIs()`, `requestAll()`)
4. Provides allocated counts after partitioning (`counts()`)
5. Allows the process to modify its allocation (`countsInc()`, `countsDec()`)

**Flow per time step:**
```
state.updateQueries()  -->  view._updateQuery()  -->  view._totalIs(current_count)
process.calculateRequest()  -->  view.requestIs(n)
state.partition()  -->  allocates from _countsRequested to _countsAllocatedFinal
process.evolveState()  -->  view.counts() / view.countsInc() / view.countsDec()
state.merge()  -->  reconciles _countsAllocatedFinal back into state containers
```

**Reusability: FULLY GENERIC**

---

## 6. Container Infrastructure

### `wholecell/containers/bulk_objects_container.py`
- Name-indexed numpy array for molecule counts
- Supports: `counts()`, `countsIs()`, `countInc()`, `countDec()`, `emptyLike()`
- Serialization: `tableCreate()`, `tableAppend()`
- Can use int64 (default for counts) or float64 (for concentrations)

### `wholecell/containers/unique_objects_container.py`
- Manages collections of unique objects as structured numpy arrays
- Each object type is a "collection" with defined attributes
- Supports: `objectsNew()`, `objectsInCollection()`, `attr()`, `attrIs()`, `delByIndexes()`
- Request-based mutation system with deferred merging (thread-safe-ish)
- Access control: `Access.EDIT`, `Access.DELETE`

**Reusability: FULLY GENERIC**

---

## 7. Listener System

### Base Class: `wholecell/listeners/listener.py`

```python
class Listener:
    def initialize(self, sim, sim_data): pass
    def allocate(self): pass
    def initialUpdate(self): pass
    def update(self): pass              # called every time step
    def tableCreate(self, tableWriter): pass
    def tableAppend(self, tableWriter): pass
```

Listeners are read-only observers. They:
- Compute derived quantities (mass, growth rate, etc.)
- Log to shell (via `registerLoggedQuantity()`)
- Write to disk (via TableWriter)
- Can be read by processes (via `readFromListener()`)

**E. coli listeners** (16 total, all organism-specific):
Mass, ReplicationData, RibosomeData, UniqueMoleculeCounts, FBAResults,
RnaDegradationListener, TranscriptElongationListener, RnapData, EnzymeKinetics,
GrowthLimits, RnaSynthProb, MonomerCounts, RNACounts, ComplexationListener,
EquilibriumListener, DnaSupercoiling, RnaMaturationListener

**Default listener** (always included): `EvaluationTime` -- tracks process timing

**Reusability: Framework GENERIC, individual listeners NOT reusable**

---

## 8. Cell Division

### File: `wholecell/sim/divide_cell.py` (564 lines)

**Entirely E. coli-specific.** The logic:
1. Check for full chromosomes (dead cell if none at end of sim)
2. `chromosomeDivision()` -- split chromosome domains into daughter trees
3. `divideBulkMolecules()` -- binomial split of bulk counts
4. `divideUniqueMolecules()` -- domain_index, RNA, active_ribosome modes
5. Save daughter state via serialization

**For platelets: NOT NEEDED.** Platelets don't divide. The engine supports this:
- Don't set `_divideCellFunction`
- The sim runs for `_lengthSec` and finalizes without division
- `finalize()` checks `self._cellCycleComplete` before calling division

---

## 9. Reconstruction / sim_data Interface

### File: `reconstruction/ecoli/simulation_data.py`

`SimulationDataEcoli` is the knowledge base. It assembles:

```
SimulationDataEcoli
    .constants          # physical constants (Avogadro, cell density)
    .molecule_ids       # canonical molecule IDs
    .molecule_groups    # groups for division, regulation, etc.
    .getter             # helper functions for lookups
    .mass               # mass parameters
    .growth_rate_parameters
    .external_state     # media definitions, exchange constraints
    .internal_state     # bulk molecule data, unique molecule definitions
    .process            # per-process parameter bundles
    .relation           # cross-process relationships
    .submass_name_to_index        # e.g. {'protein': 0, 'rRNA': 1, ...}
    .compartment_abbrev_to_index  # e.g. {'c': 0, 'p': 1, 'e': 2}
```

### Minimum sim_data interface for a new organism

Based on what `BulkMolecules`, `UniqueMolecules`, `InternalState`, and the `Simulation`
engine actually read, a minimal `SimulationDataPlatelet` must provide:

```python
class SimulationDataPlatelet:
    # Required by InternalState base:
    submass_name_to_index = {}          # e.g. {'protein': 0, 'lipid': 1, 'small_molecule': 2}
    compartment_abbrev_to_index = {}    # e.g. {'c': 0, 'm': 1, 'dg': 2, 'ag': 3, 'e': 4}

    # Required by BulkMolecules:
    constants.n_avogadro                # Avogadro's number
    internal_state.bulk_molecules.bulk_data['id']    # molecule ID array
    internal_state.bulk_molecules.bulk_data['mass']  # mass array with units
    molecule_groups.bulk_molecules_binomial_division  # list of mol IDs
    molecule_groups.bulk_molecules_equal_division     # list of mol IDs

    # Required by UniqueMolecules:
    internal_state.unique_molecule.unique_molecule_definitions  # {name: {attr: dtype}}
    internal_state.unique_molecule.unique_molecule_masses       # structured array
    molecule_groups.unique_molecules_*_division                 # division mode lists

    # Required by LocalEnvironment:
    external_state.make_media
    external_state.current_timeline_id
    external_state.saved_timelines
    external_state.saved_media
    external_state.exchange_data_from_concentrations
    external_state.exchange_to_env_map
    external_state.import_constraint_threshold
```

---

## 10. E. coli Simulation Wiring (Example)

### File: `models/ecoli/sim/simulation.py`

This is the template for what `models/platelet/sim/simulation.py` would look like:

```python
class EcoliSimulation(Simulation):
    _internalStateClasses = (BulkMolecules, UniqueMolecules)
    _externalStateClasses = (LocalEnvironment,)

    _processClasses = (
        (TfUnbinding,),
        (Equilibrium, TwoComponentSystem, RnaMaturation),
        (TfBinding,),
        (TranscriptInitiation, PolypeptideInitiation, ChromosomeReplication,
         ProteinDegradation, RnaDegradation, Complexation),
        (TranscriptElongation, PolypeptideElongation),
        (ChromosomeStructure,),
        (Metabolism,),
        (CellDivision,),
    )

    _listenerClasses = (Mass, ReplicationData, ...)
    _initialConditionsFunction = calcInitialConditions
    _divideCellFunction = divide_cell
```

---

## 11. Reusability Summary

| Component | File(s) | Generic? | Platelet action |
|-----------|---------|----------|-----------------|
| Simulation engine | `wholecell/sim/simulation.py` | YES | Use as-is, subclass |
| Process base class | `wholecell/processes/process.py` | YES | Use as-is |
| View system | `wholecell/views/view.py` | YES | Use as-is |
| InternalState base | `wholecell/states/internal_state.py` | YES | Use as-is |
| ExternalState base | `wholecell/states/external_state.py` | YES | Use as-is |
| BulkMolecules | `wholecell/states/bulk_molecules.py` | MOSTLY | Override `tableAppend` (ATP[c] hardcode) |
| UniqueMolecules | `wholecell/states/unique_molecules.py` | MOSTLY | Subclass: custom division modes, multi-compartment mass |
| LocalEnvironment | `wholecell/states/local_environment.py` | MOSTLY | Simplify exchange logic for initial work |
| Containers | `wholecell/containers/*.py` | YES | Use as-is |
| Listener base | `wholecell/listeners/listener.py` | YES | Use as-is |
| EvaluationTime | `wholecell/listeners/evaluation_time.py` | YES | Use as-is (default listener) |
| Loggers (Shell, Disk) | `wholecell/loggers/*.py` | YES | Use as-is |
| Cell division | `wholecell/sim/divide_cell.py` | NO | Not needed (platelets don't divide) |
| E. coli processes | `models/ecoli/processes/*.py` | NO | Write new platelet processes |
| E. coli listeners | `models/ecoli/listeners/*.py` | NO | Write new platelet listeners |
| E. coli reconstruction | `reconstruction/ecoli/` | NO | Write new `reconstruction/platelet/` |

---

## 12. Proposed Platelet File Structure

```
models/platelet/
    sim/
        simulation.py           # PlateletSimulation(Simulation)
        initial_conditions.py   # set initial molecule counts
    processes/
        __init__.py
        calcium_dynamics.py     # Ca2+ store release, influx, buffering
        granule_release.py      # dense/alpha granule exocytosis
        receptor_signalling.py  # GPVI, PAR, P2Y12 -> Ca2+
        protein_decay.py        # slow protein/mRNA turnover
        # later: integrin_activation.py, shape_change.py, metabolism.py
    listeners/
        __init__.py
        mass.py                 # total mass tracking
        activation_state.py     # resting/primed/activated/spent
        granule_counts.py       # granule inventory over time
        calcium_trace.py        # Ca2+ concentration trace

reconstruction/platelet/
    simulation_data.py          # SimulationDataPlatelet
    knowledge_base_raw.py       # curated literature values
    dataclasses/
        constants.py
        molecule_groups.py
        state/
            internal_state.py   # bulk + unique molecule definitions
            external_state.py   # agonist environment
```

---

## 13. Key Engine Behaviors to Understand for Platelet Work

### Time step adaptation
- Default max time step: 1.0 s (`MAX_TIME_STEP`)
- Default sim length: 3 hours (`lengthSec = 3*60*60`)
- Binary search for valid dt: each process can veto a too-large dt
- For platelet activation (seconds to minutes), the defaults are reasonable
- For platelet lifespan (days), you'd increase `lengthSec` significantly

### Partitioning (resource allocation)
- Only applies to `BulkMolecules` (not `UniqueMolecules`)
- When multiple processes compete for the same molecule, the engine allocates proportionally
- Priority levels allow critical processes to be served first
- This is useful for platelet ATP/Ca2+ competition between processes

### Simulation termination
- `run()` runs for `_lengthSec` seconds of simulated time
- If `_cellCycleComplete` is set (by a process like CellDivision), sim ends early
- If `_raise_on_time_limit` is True, an exception is raised if time runs out without division
- For platelets: set `_raise_on_time_limit = False`, let sim run to completion
- Or: write a `PlateletDeath` process that sets `sim._isDead = True` when the platelet is spent

### DEFAULT_SIMULATION_KWARGS
Many kwargs are E. coli-specific (`trna_charging`, `ppgpp_regulation`, `superhelical_density`, etc.). For the platelet model:
- Override defaults in `PlateletSimulation.__init__()` or
- Define a new `DEFAULT_PLATELET_KWARGS` dict
- At minimum, the engine only requires: `seed`, `lengthSec`, `initialTime`, `timeStepSafetyFraction`, `maxTimeStep`, `outputDir`, `simData`

---

## 14. Data Flow Diagram

```
                    sim_data (knowledge base)
                         |
           +-------------+-------------+
           |             |             |
           v             v             v
     BulkMolecules  UniqueMols   LocalEnvironment
     (counts)       (objects)    (concentrations)
           |             |             |
           |    Views    |    Views    |    Views
           |  +------+   |  +------+  |  +------+
           +->| Proc |<--+->| Proc |<-+->| Proc |
              +------+      +------+     +------+
                 |              |            |
                 v              v            v
              evolveState(): modify allocated counts/objects
                 |              |            |
                 v              v            v
              Merge back into state containers
                 |
                 v
              Listeners observe state, compute derived quantities
                 |
                 v
              Loggers write to shell/disk
```

---

## 15. Critical Path for Phase 2 (Minimal Platelet Stub)

To get a runnable platelet simulation, you need exactly:

1. **`SimulationDataPlatelet`** with:
   - `submass_name_to_index` (e.g. `{'protein': 0, 'lipid': 1, 'small_molecule': 2}`)
   - `compartment_abbrev_to_index` (e.g. `{'c': 0, 'e': 1}`)
   - `constants.n_avogadro`
   - `internal_state.bulk_molecules.bulk_data` (id + mass arrays)
   - `internal_state.unique_molecule.unique_molecule_definitions` (can be empty dict)
   - `internal_state.unique_molecule.unique_molecule_masses` (can be empty)
   - `molecule_groups.*_division` (can be empty lists)
   - `external_state.*` (minimal media definition)

2. **One process** (e.g. `RestingDecay`) that:
   - Creates a `bulkMoleculesView` for a few molecules
   - In `calculateRequest()`: requests some molecules
   - In `evolveState()`: decrements counts (exponential decay)

3. **One listener** (e.g. `Mass`) that:
   - Reads `BulkMolecules` mass and logs it

4. **`PlateletSimulation`** class wiring it all together

5. **`calcInitialConditions()`** that sets starting molecule counts

This is achievable in ~500 lines of new code, proving the framework works for non-E. coli organisms.
