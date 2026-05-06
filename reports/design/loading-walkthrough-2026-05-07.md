# Loading & runtime walkthrough — platelet-wcm — 2026-05-07

A complete trace from `python runscripts/manual/runPlateletSim.py` to
the first row of output values on disk, with file/line references.

Companion to `code-walkthrough-2026-05-07.md` (file-by-file repo tour),
`demo-2026-05-07.md` (feature demo flow), and
`biology-overview-2026-05-07.md` (biology 1-pager). This doc is the
*temporal* view: in what order do things happen, what's executing
when, and where does the very first output value come from.

---

## 60-second summary

```
runPlateletSim.main()                              # runPlateletSim.py:169
  ├── argparse → length_sec, seed, ca_ex_mM, ip3_forced
  └── run_platelet_sim()                           # line 66
      ├── cs_mod·CA_EX_UM = ca_ex_mM × 1000        # module override
      ├── CalciumDynamics._ip3_forced = …        # class override
      │
      ├── SimulationDataPlatelet()                 # build the knowledge base
      ├── persist sim_data under kb/               # for analysis
      │
      ├── PlateletSimulation(simData=…, lengthSec=…, …)
      │     └── Simulation.__init__                # simulation.py:116
      │          ├── kwarg parsing, output-dir prep, RNG seed
      │          └── self._initialize(sim_data)    # line 164
      │               ├── instantiate states, processes, listeners
      │               ├── state.initialize() for each state
      │               ├── state.allocate() for each state
      │               ├── _initialConditionsFunction(sim_data)  # line 197
      │               ├── process.initialize() for each process
      │               ├── listener.initialize() + .allocate()
      │               ├── state.calculateMass()    # initial mass
      │               ├── listener.initialUpdate() # initial attr values
      │               └── logger.initialize(self)  # FIRST ROW (line 234)
      │                    ├── open TableWriter per state/listener
      │                    ├── tableCreate() on each (attrs JSON)
      │                    └── copyData → tableAppend (t=0 row)
      │
      └── sim.run()                                # line 252
            └── run_incremental(lengthSec + initialTime)  # line 267
                  while time < end:
                      _simulationStep += 1
                      _timeTotal += _timeStepSec
                      _pre_evolve_state()          # line 319
                      for processes in _processClasses:
                          _evolveState(processes)  # line 334
                            ├── state.updateQueries()    # current totals
                            ├── process.calculateRequest()  # requests
                            ├── state.partition(processes)  # allocate
                            ├── process.evolveState()    # biology step
                            ├── process.wasTimeStepShortEnough()
                            └── state.merge(processes)
                      _post_evolve_state()         # line 377
                            ├── state.calculateMass()
                            ├── listener.update()
                            └── logger.append → tableAppend  # OUTPUT
```

For the platelet model with default settings the loop runs ~200 times
in ~30 s of wall-clock to produce a 200 s simulation.

---

## Stage by stage with file pointers

### 1. CLI parse — `runscripts/manual/runPlateletSim.py`

`main()` (line 169) builds the `argparse` parser, resolves the output
directory under `out/`, calls `write_metadata()` (saves
`metadata/metadata.json`), then `run_platelet_sim()`.

`run_platelet_sim()` (line 66) does **two attribute-level overrides
before the simulation is constructed**:

```python
cs_mod·CA_EX_UM = float(ca_ex_mM) * 1000.0          # extracellular Ca²⁺ in uM
CalciumDynamics._ip3_forced = bool(ip3_forced)      # Dolan IP3 forcing on/off
```

This is the *only* mutation point for these conditions; the ODE and
listener read them on every step from the same module-level globals.

### 2. Build the knowledge base — `SimulationDataPlatelet()`

`reconstruction/platelet/simulation_data.py` constructs (in order):

```python
self.constants       = Constants()
self.molecule_groups = MoleculeGroups()
self.internal_state  = InternalState()    # builds bulk_data from _MOLECULES
self.external_state  = ExternalState()
self.process         = Process(self)      # RestingDecay + CalciumSignalling
```

`Process(sim_data)` instantiates `CalciumSignalling`, which wires up
the 27-species ODE state vector + every rate-constant table from
`calcium_signalling.py`. **No simulation has started yet** — this is
just data construction. Takes a few ms.

`runPlateletSim.write_sim_data()` (line 56) persists the resulting
object to `kb/` so downstream analysis can load it without rebuilding.

### 3. Construct the simulation — `PlateletSimulation(simData=…, lengthSec=…, …)`

`models/platelet/sim/simulation.py` is just declarations:

```python
class PlateletSimulation(Simulation):
    _internalStateClasses = (BulkMolecules, UniqueMolecules)
    _externalStateClasses = (LocalEnvironment,)
    _processClasses = ((RestingDecay, CalciumDynamics),)  # one group
    _listenerClasses = (Mass, CalciumTrace)
    _initialConditionsFunction = calcInitialConditions
```

The base class's `__init__` does the actual work
(`wholecell/sim/simulation.py:116`):

- **lines 117–127** — validate the subclass declared the required
  attributes.
- **lines 129–139** — merge `kwargs` into instance attrs (each
  `_<name>`); warn on unknown kwargs. This is where `_lengthSec`,
  `_seed`, `_outputDir`, `_simData` etc. get set.
- **line 142** — `self._timeStepSec = min(MAX_TIME_STEP, self._maxTimeStep)`
  → **1.0 s**.
- **line 146** — RNG seeded from `_seed`.
- **lines 153–155** — output dir wiped clean and recreated. *Worth
  knowing during the walkthrough:* every run starts with an empty
  `simOut/`. There's no append-to-existing.
- **line 160** — `self._initialize(sim_data)` does the heavy lifting.

### 4. Wire everything up — `Simulation._initialize(sim_data)` (line 164)

This is where the main execution happens. In execution
order:

```
166-168  collect _processClasses into a flat set
170-172  instantiate each state, process, listener
            (just _orderedAbstractionReference — calls cls() on each)
174      append default listeners (EvaluationTime)
176      _initLoggers()  → creates Shell logger and (if logToDisk) Disk logger

181-186  for each internal_state:
            assign seeded RandomState
            internal_state.initialize(self, sim_data)
            (BulkMolecules reads bulk_data; UniqueMolecules reads
             unique_molecule definitions)
188-189  for each external_state: initialize(sim, sim_data, timeline)
194-195  for each internal_state: allocate()
            (creates the underlying numpy arrays at the right size)

197      self._initialConditionsFunction(sim_data)
            ← FIRST POINT WHERE MOLECULE COUNTS BECOME REAL

199-204  for each process: initialize(sim, sim_data)
            (creates BulkMoleculesView / UniqueMoleculesView /
             EnvironmentView)
206-207  for each listener: initialize(sim, sim_data)
209-210  for each listener: allocate()

212      self._timeTotal = self.initialTime()
            ← simulated clock = 0 (or whatever initialTime kwarg was)

220-222  state.calculateMass()                # initial mass numbers
224-226  external_state.update()              # media timeline applied
228-230  listener.initialUpdate()    # first attribute values

232-234  for each logger: initialize(self)
            ← OPENS TableWriters AND WRITES THE t=0 ROWS
```

### 5. Initial conditions — `models/platelet/sim/initial_conditions.py`

Called at line 197 above with the freshly-constructed `sim_data`. For
the platelet model it copies
`sim_data.internal_state.bulk_molecules.initial_counts` into the
`BulkMolecules` container. After this call:

- `sim.internal_states['BulkMolecules'].container.counts()` matches the
  initial-counts array.
- Cytosolic Ca²⁺ count = 361 (= 100 nM × 6 fL × N_A × 10⁻⁶), DTS Ca²⁺ =
  38 842, etc.

This is *where the model becomes a model* — before this call,
`_internalStateClasses = (BulkMolecules, …)` is just a tuple of class
references with no data behind them.

### 6. Process initialize — `CalciumDynamics.initialize(sim, sim_data)`

(`models/platelet/processes/calcium_dynamics.py`)

```python
def initialize(self, sim, sim_data):
    super().initialize(sim, sim_data)
    self._solver = sim_data.process.calcium_signalling     # ODE module
    # 27-species view; see calcium_signalling·MOLECULE_NAMES
    self._molecules = self.bulkMoleculesView(self._solver.molecule_names)
    self._atp = self.bulkMoleculesView(np.array(['ATP[c]'], …))  # ATP debit
```

The views are **not the data** — they're handles into the
`BulkMolecules` container that the framework uses to route requests.
Same pattern for `RestingDecay`.

### 7. Where the first output values come from

`Simulation._initialize` line 234: `logger.initialize(self)` is the
moment the t=0 values get to disk.

`Disk.initialize(sim)` (`wholecell/loggers/disk.py:43`):

```python
self.mainFile = TableWriter(os.path.join(self.outDir, "Main"))
self.mainFile.writeAttributes(initialTime=…, startTime=…)
self.createTables(sim)     # TableWriter per state/listener; attrs JSON
self.copyData(sim)            # ← writes the FIRST ROW of every column file
```

`copyData()` (line 95) calls `sim.tableAppend(mainFile)` (writes `time`
+ `timeStepSec`) and then for every state and listener:

```python
obj.tableAppend(saveFile)
```

For `CalciumTrace`, this writes the 14 columns at their initial values:
`ca_cyt_nM ≈ 99.9`, `ca_dts_uM ≈ 250.0`, `ip3_nM ≈ 50.1`, etc. — the
values that came out of `listener.initialUpdate()` two lines earlier.

**That's the t=0 row in `simOut/CalciumTrace/ca_cyt_nM` (and every
other column file).**

### 8. The main loop — `Simulation.run() → run_incremental()` (lines 252, 267)

```python
while self.time() < run_until and not self._isDead:
    self._simulationStep += 1
    self._timeTotal += self._timeStepSec      # advance simulated clock
    self._pre_evolve_state()                  # adjust dt, reset mass diffs
    for processes in self._processClasses:  # one group
        self._evolveState(processes)          # v
    self._post_evolve_state()        # mass, listener.update, log
```

Note `_processClasses` is a *tuple of tuples* — the outer tuple is
process *groups* (executed sequentially within one timestep); the
inner tuple is processes within a group (sharing one partition step).
The platelet model has one group containing two processes.

### 9. One step in detail — `_evolveState(processes)` (line 334)

The five-stage per-timestep dance, in order:

1. **`state.updateQueries()`** (line 339) — every view re-reads its
   current state totals so the process sees up-to-date counts when it
   asks "how much do I have?".
2. **`process.calculateRequest()`** (line 346) — for `CalciumDynamics`:
   reads current 27-species counts, calls
   `solver.molecules_to_next_time_step(counts, dt, t_sim, ip3_forced)`
   which integrates the ODE for one second using
   `scipy.integrate.solve_ivp(method='BDF')`, returns `(delta,
   atp_cost)`. The view's `requestIs()` is called for each species the
   process intends to *withdraw*.
3. **`state.partition(processes)`** (line 352) — `BulkMolecules`
   distributes molecules to processes proportional to their requests,
   weighted by priority. Returns each process's actual allocation.
4. **`process.evolveState()`** (line 359) — `CalciumDynamics.evolveState`:
   applies the deltas via `countsInc(delta)`, debits ATP via
   `countsDec(atp_cost)`. *The biology has now happened for this step.*
5. **`process.wasTimeStepShortEnough()`** (line 364) +
   **`state.merge(processes)`** (line 370) — validation + reconcile
   partitioned counts back into the master container.

For the calcium model the heavy lifting is in step 2 (the ODE
integration); steps 3, 4, 5 are bookkeeping.

### 10. Output for this step — `_post_evolve_state()` (line 377)

```python
state.calculateMass()         # update Mass listener inputs
listener.update()             # CalciumTrace re-reads, computes columns
hook.postEvolveState(self)
logger.append(sim) → copyData(sim) → tableAppend on every column file
```

The new row of `ca_cyt_nM`, `ca_dts_uM`, etc. is now on disk. The loop
continues with `time = 2 s`, then `3 s`, …

---

## Why the first row in `CalciumTrace` is `t=0, ca_cyt_nM ≈ 99.9`

That value is the count `361` from the initial-conditions table
converted to nM:

```
361 / (6.022 × 10²³ × 6 × 10⁻¹⁵ × 10⁻⁶) × 1000  ≈  99.9 nM
```

It's *not* the result of any ODE step — it's the initial-conditions
table written before any process has run. The first ODE step delivers
row `t=1` to disk, which (for the +Ca²⁺ stimulated condition) jumps to
~298 nM because of the IP3R-drainage / SERCA-redistribution dynamics
during that single 1 s integration. That single-timestep jump is the
"resting spike" question — see the demo Q&A and design doc §6.8 D7.

---

## Key points

1. **The entry point is `run_platelet_sim()` in
   `runscripts/manual/runPlateletSim.py` (line 66).** Two
   attribute-level overrides (`cs_mod.CA_EX_UM`,
   `CalciumDynamics._ip3_forced`) happen here, before the simulation
   is constructed — they're how the EDTA / resting conditions are
   flipped.
2. **The wiring is in `wholecell/sim/simulation.py:_initialize`
   (lines 164–234).** It's ~70 lines and reads top-to-bottom in
   execution order. If asked "where do molecules come into being?",
   it's line 197 (the initial-conditions function).
3. **The first output value lands at line 234** —
   `logger.initialize(self)` opens the `TableWriter`s and writes the
   t=0 row before any process has ever run.
4. **The main loop is `run_incremental` at line 267**. Per-step body
   is `_pre_evolve_state` → `_evolveState(processes)` × however many
   groups → `_post_evolve_state`.
5. **The biology happens in `process.evolveState()` (line 359)** — for
   the calcium model that's a `scipy.integrate.solve_ivp` BDF step
   over 27 species inside
   `CalciumSignalling.molecules_to_next_time_step()`.
6. **The output writing happens in `logger.append(sim)` at line 397**
   — same `copyData()` pattern as the initial write, just at every
   step instead of once.

