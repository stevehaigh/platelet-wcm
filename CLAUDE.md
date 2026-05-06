# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Whole-cell model of a **human platelet**, built as a dissertation project.
Repo: `stevehaigh/platelet-wcm` (forked and pruned from `CovertLab/wcEcoli`).
The simulation engine in `wholecell/` is reused; all E. coli biology has been removed
and replaced with platelet-specific processes, listeners, and reconstruction code.

Current focus: intracellular Ca²⁺ dynamics (IP3R / SERCA / PMCA / SOCE / CaM),
validated against Dolan & Diamond 2014. Single-cell, no division.

## Build & Run

Python 3.11.5 (pinned via `.python-version`; pyenv recommended). All commands assume
`PYTHONPATH="$PWD"` from the repo root. **No Cython compile step is required** —
the performance-critical `.pyx` modules from wcEcoli have been removed.

```bash
# Install dependencies
pip install -r requirements.txt

# Run a simulation (no ParCa — sim_data is constructed directly)
PYTHONPATH=$PWD python runscripts/manual/runPlateletSim.py [sim_outdir] --length 200

# Generate analysis plots
PYTHONPATH=$PWD python runscripts/manual/analysisPlatelet.py [sim_outdir]

# Pick a specific plot
PYTHONPATH=$PWD python runscripts/manual/analysisPlatelet.py [sim_outdir] --plot calcium_trace

# Run the Phase 3 two-condition validation (with vs without extracellular Ca²⁺)
PYTHONPATH=$PWD python runscripts/manual/runPhase3.py [sim_outdir] --length 200

# Run the Dash webapp locally (http://localhost:8050)
make run     # foreground with hot reload
make stop    # kill it
```

All runscripts support `-h` for full options.

### Run-time conditions

`runPlateletSim.py` exposes three CLI flags that change the biology being
simulated, beyond simulation length and seed:

| Flag | Default | Effect |
|------|---------|--------|
| `--length N` | 60 | Simulation length in seconds |
| `--seed N` | 0 | RNG seed (currently no stochastic processes use it) |
| `--ca-ex-mM X` | 1.2 | Extracellular Ca²⁺ in mM. Set `0` for the Dolan Fig. 4 EDTA / no-extracellular-Ca²⁺ condition (SOCE inactive, PM leak inactive). |
| `--no-ip3-forcing` | (off) | Disable the Dolan Fig. S2 IP3 time curve, so IP3 stays at its 50 nM baseline. Use for an at-rest / un-stimulated sim. |

Where each behaviour is defined in code:

- **Extracellular Ca²⁺** is the module-level constant `CA_EX_UM` in
  `reconstruction/platelet/dataclasses/process/calcium_signalling.py`. The
  runscript overrides it via `cs_mod.CA_EX_UM = ca_ex_mM * 1000.0` before
  the sim is constructed; `_ode_rhs` reads it on every ODE step. Both the
  SOCE current and the basal PM Ca²⁺ leak are gated on `CA_EX_UM > 0`.
- **IP3 forcing** is the class attribute `CalciumDynamics._ip3_forced` in
  `models/platelet/processes/calcium_dynamics.py`. The runscript overrides
  it before constructing `PlateletSimulation`. When `True`, the
  `CalciumDynamics` process passes `ip3_forced=True` into the ODE solver,
  which applies `ip3_forcing_uM(t)` (the Dolan Fig. S2 fit, also in
  `calcium_signalling.py`) to the IP3 state variable each step.

The same conditions are exposed on the **webapp** Configure tab as form
fields (Extracellular Ca²⁺ mM, IP3 forcing checkbox). The three webapp
presets — IP3 transient (+Ca²⁺), EDTA transient, Resting — are defined
in `wholecell/webapp/tabs/configure.py:PRESETS` and differ in exactly
those three knobs (length, `ca_ex_mM`, `ip3_forced`). The Phase 3 driver
`runPhase3.py` runs the +Ca²⁺ and EDTA conditions back-to-back.

## Tests & Type Checking

```bash
# All platelet tests
PYTHONPATH=$PWD python3 -m pytest models/platelet/tests/ -v

# Single test file or method
PYTHONPATH=$PWD python3 -m pytest models/platelet/tests/sim/test_simulation.py
PYTHONPATH=$PWD python3 -m pytest wholecell/tests/utils/test_units.py::TestUnits::test_some_method

# Type checking (platelet paths only)
python -m mypy models/platelet/ reconstruction/platelet/ \
    runscripts/manual/runPlateletSim.py runscripts/manual/analysisPlatelet.py
```

CI (`.github/workflows/ci.yml`) runs pytest + mypy on every push and PR to `main`.
Set `OPENBLAS_NUM_THREADS=1` for reproducible numerics.

## Architecture

The simulation runs in **1-second discrete timesteps**:

```
SimulationDataPlatelet  (constructed in code; no ParCa step)
  -> PlateletSimulation loop -> Listeners record to disk
    -> Analysis plots read from disk
```

### Framework vs. Model

- `wholecell/` — **Framework**: model-agnostic simulation engine, base classes, I/O,
  utilities. Inherited from CovertLab/wcEcoli; division-related code removed/unused.
  Includes `wholecell/webapp/` (Dash app for browsing runs).
- `models/platelet/` — **Platelet model**: processes, listeners, analysis, sim wiring.
- `reconstruction/platelet/` — **Parameters**: `SimulationDataPlatelet` + dataclasses.
  Replaces wcEcoli's heavyweight ParCa with a directly-constructed parameter object.

### Simulation Loop

`PlateletSimulation` (in `models/platelet/sim/simulation.py`) subclasses
`Simulation` (in `wholecell/sim/simulation.py`).

Each timestep:
1. **Processes** declare resource requests via `calculateRequest()`
2. **States** partition molecules to processes
3. **Processes** execute via `evolveState()`
4. **States** merge results back
5. **Listeners** observe and record (read-only)

Processes run in **dependency-ordered groups**. Within a group, processes share a
partition step. The current platelet model has a single group:
`(RestingDecay, CalciumDynamics)`.

### Core Abstractions

| Concept | Base class | Platelet impl | Purpose |
|---------|-----------|---------------|---------|
| **Process** | `wholecell/processes/process.py` | `models/platelet/processes/` (RestingDecay, CalciumDynamics) | Biological submodels that modify state |
| **State** | `wholecell/states/internal_state.py` | `BulkMolecules`, `UniqueMolecules`, `LocalEnvironment` | Cellular state containers |
| **Listener** | `wholecell/listeners/listener.py` | `models/platelet/listeners/` (Mass, CalciumTrace) | Observe and record data each timestep |
| **Analysis** | `models/platelet/analysis/analysisPlot.py` | `single/` (calcium_trace, scaffold_summary) | Post-simulation plots |

There is no platelet equivalent of E. coli variants yet; the simulation runs a single
"platelet_stub" condition, written into the output path as `platelet_stub_{seed:06d}/`.

### Process Lifecycle

Processes implement:
- `initialize(sim, sim_data)` — setup, create views
- `calculateRequest()` — declare molecule needs
- `evolveState()` — update cellular state
- `isTimeStepShortEnough()` / `wasTimeStepShortEnough()` — optional timestep validation

Processes access state through **views**:
- `bulkMoleculesView(moleculeIDs)` — bulk molecule counts (integers)
- `uniqueMoleculesView(moleculeName)` — individual molecule instances
- `environmentView(moleculeIDs)` — external environment concentrations

### Calcium signalling — current focus

The biologically rich process is `CalciumDynamics`. It is a thin wrapper around the
ODE solver and rate constants in
`reconstruction/platelet/dataclasses/process/calcium_signalling.py`, which encodes:

- **IP3R** — 6-state Markov model (Sneyd & Dufour 2002; Purvis & Bhatt 2008 constants)
- **SERCA** — E1/E2 enzymatic cycle (Dode 2002)
- **PMCA** — 5-state CaM-coupled (Caride 2007 Table 3): basal + Ca₄·CaM-activated paths
- **SOCE** — STIM1 dimerisation + Orai1 MWC flux (Dolan & Diamond 2014; Hoover)
- **CaM** — three sub-states (free, Ca₂·CaM, Ca₄·CaM) acting as cytosolic Ca²⁺ buffer
- **IP3** — currently a forced time curve (Dolan 2014 Fig. S2); `_ip3_forced=True`.
  In v0.3 a P2Y1 upstream process will produce IP3 endogenously.

`CalciumTrace` listener records 14 columns (Ca²⁺ pools, CaM/PMCA sub-states, IP3, SOCE flux).
The 5-panel `single/calcium_trace.py` plot is the headline validation figure.

Validation target: Dolan & Diamond 2014 Fig. 4 (Ca²⁺ transients with/without
extracellular Ca²⁺).

### State Partitioning

`BulkMolecules` uses integer counts. Processes request molecules and get partitioned
allocations; priority levels control allocation order.

`UniqueMolecules` tracks individual instances with attributes — currently unused by
the platelet model but available for future granule modelling.

### Data I/O

Simulation output uses a custom binary columnar format (zlib-compressed):
- **Write**: `TableWriter` (in `wholecell/io/tablewriter.py`)
- **Read**: `TableReader` (in `wholecell/io/tablereader.py`)

```python
from wholecell.io.tablereader import TableReader
reader = TableReader(os.path.join(simOutDir, 'CalciumTrace'))
ca_cyt = reader.readColumn('ca_cyt_uM')
```

Each listener/state writes to `simOut/{Name}/` with column files + `attributes.json`.

### Units System

Physical quantities use the `Unum` library (via `wholecell/utils/units.py`):

```python
from wholecell.utils import units
mass = 5 * units.fg
concentration = 100 * units.umol / units.L
value = mass.asNumber(units.g)  # dimensionless float
```

### SimulationDataPlatelet

Built directly in `reconstruction/platelet/simulation_data.py`. Provides the minimal
interface the engine expects: `submass_name_to_index`, `compartment_abbrev_to_index`,
`constants`, `molecule_groups`, `internal_state`, `external_state`, `process`.

Compartments specific to platelets:

| Abbrev | Compartment |
|--------|-------------|
| `c`    | cytoplasm |
| `dts`  | dense tubular system (ER equivalent; Ca²⁺ store) |
| `dg`   | dense granule lumen |
| `ag`   | alpha-granule lumen |
| `m`    | mitochondrial matrix |
| `pl`   | plasmalemma (membrane) |
| `e`    | extracellular / open canalicular system |

Note `pl` (not `pm`) to avoid clashing with `m` (mitochondria).

Rate constants and ODE state for calcium signalling live in
`reconstruction/platelet/dataclasses/process/calcium_signalling.py`. Initial counts
live in `reconstruction/platelet/dataclasses/internal_state.py` and the
`raw_data/molecules.tsv` table.

## Key Conventions

### Style

- **Indentation**: TABs (4-space tab stops) — project-wide convention, not spaces
- **Naming**: `ClassName`, `function_name`, `GLOBAL_CONSTANT`; some camelCase from the
  upstream framework remains, but snake_case is preferred for new code
- **Imports**: absolute from repo root (`from wholecell.sim.simulation import Simulation`)
- **Line length**: soft 79, hard 99
- **Docstrings**: triple-quoted, imperative first sentence; raw strings for LaTeX (`r"""..."""`)
- See `docs/style-guide.md`

### Molecule IDs

Molecules carry compartment tags: `ATP[c]`, `Ca[c]`, `Ca[dts]`, `STIM1_dimer[pl]`, etc.

### Output Directory Structure

```
out/{sim_dir}/
  metadata/                                  # git hash, seed, timestamp
  kb/simData.cPickle                         # serialized SimulationDataPlatelet
  platelet_stub_{seed:06d}/
    {seed:06d}/
      generation_000000/
        000000/
          simOut/                            # Listener data (TableWriter format)
          plotOut/                           # Analysis figures
```

(The variant/seed/generation/cell nesting is preserved from wcEcoli so the webapp's
"Inspect Data" tab can browse platelet runs without modification.)

### Environment

Set `OPENBLAS_NUM_THREADS=1` to avoid threading artifacts in numerical results.

### Webapp

Dash app at `wholecell/webapp/`. `make run` starts it locally. Pushing to the `webapp`
branch (`make deploy`) triggers `.github/workflows/deploy-azure.yml`, which builds
docker images (`docker/runtime/`, `docker/webapp/`) and deploys to Azure Container
Instances at `platelet-wcm.uksouth.azurecontainer.io`.

### Reports & docs

- `platelet-plan/` — high-level project plan (`overview.md`, `architecture.md`)
- `reports/design/` — design docs (calcium pathway, runtime scaffold, reconstruction)
- `reports/lab-books/` — dated session notes (most recent is current state of the work)
- `reports/data/`, `reports/figures/` — calibration data and rendered figures
- `make pdfs` builds PDFs of any `reports/*.md` into `reports/pdf/`
- `docs/style-guide.md` — full style guide
