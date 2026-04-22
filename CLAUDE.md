# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Whole-cell model of *Escherichia coli* from the Covert Lab at Stanford.
Fork: `stevehaigh/wcEcoli` (upstream: `CovertLab/wcEcoli`).

## Build & Run

Python 3.11.3. All commands assume `PYTHONPATH="$PWD"` from the repo root.

```bash
# Compile Cython extensions (required before running)
make clean compile

# Run parameter calculator (ParCa) — ~18 min
python runscripts/manual/runParca.py [sim_outdir]

# Run simulation — ~10 min
python runscripts/manual/runSim.py [sim_dir]

# Run with a variant (e.g., nutrient timeline index 2)
python runscripts/manual/runSim.py -v timelines 2 2

# Analysis plots
python runscripts/manual/analysisSingle.py [sim_dir]
python runscripts/manual/analysisCohort.py [sim_dir]
python runscripts/manual/analysisVariant.py [sim_dir]
```

All runscripts support `-h` for full option docs. Options can be abbreviated (`--cpus 8` or `-c8`).

## Tests & Type Checking

The repo is pinned to Python 3.11.5 via `.python-version`. If pyenv is initialised
in your shell (`eval "$(pyenv init -)"` in `.zshrc`/`.bash_profile`), plain
`python3` and `pytest` will auto-select the right interpreter.

```bash
# Run all tests
python3 -m pytest

# Run a single file or test method
python3 -m pytest models/platelet/tests/
python3 -m pytest wholecell/tests/utils/test_units.py::TestUnits::test_some_method

# Type checking
mypy
```

Tests use pytest (some legacy unittest). No CI linter beyond mypy.

## Architecture

The simulation runs in **1-second discrete timesteps** through this pipeline:

```
Raw Data (72+ TSV files)
  -> ParCa (parameter fitting) -> sim_data
    -> Simulation loop -> Listeners record to disk
      -> Analysis plots read from disk
```

### Framework vs. Model

- `wholecell/` — **Framework**: model-agnostic simulation engine, base classes, I/O, utilities
- `models/ecoli/` — **E. coli model**: processes, listeners, analysis plots, variants
- `reconstruction/ecoli/` — **ParCa**: fits raw experimental data into `sim_data` parameters

The framework is designed to be reusable for other organisms. E. coli-specific code subclasses the framework.

### Simulation Loop

`EcoliSimulation` (in `models/ecoli/sim/simulation.py`) subclasses `Simulation` (in `wholecell/sim/simulation.py`).

Each timestep:
1. **Processes** declare resource requests via `calculateRequest()`
2. **States** partition molecules to processes
3. **Processes** execute via `evolveState()`
4. **States** merge results back
5. **Listeners** observe and record (read-only)

Processes run in **dependency-ordered groups** (e.g., TfUnbinding -> Equilibrium -> TfBinding -> Transcription/Translation -> ...). Within a group, processes run in parallel on shared state.

### Core Abstractions

| Concept | Base class | E. coli impl | Purpose |
|---------|-----------|--------------|---------|
| **Process** | `wholecell/processes/process.py` | `models/ecoli/processes/` (17 processes) | Biological submodels that modify state |
| **State** | `wholecell/states/internal_state.py` | `BulkMolecules`, `UniqueMolecules`, `LocalEnvironment` | Cellular state containers |
| **Listener** | `wholecell/listeners/listener.py` | `models/ecoli/listeners/` (18 listeners) | Observe and record data each timestep |
| **Analysis** | `models/ecoli/analysis/analysisPlot.py` | `single/`, `cohort/`, `multigen/`, `variant/`, `parca/` | Post-simulation plots |
| **Variant** | Functions in `models/ecoli/sim/variants/` | ~32 variant types | Modify `sim_data` for experiments |

### Process Lifecycle

Processes implement:
- `initialize(sim, sim_data)` — setup references
- `calculateRequest()` — declare molecule needs
- `evolveState()` — update cellular state
- `isTimeStepShortEnough()` / `wasTimeStepShortEnough()` — timestep validation

Processes access state through **views**:
- `bulkMoleculesView(moleculeIDs)` — bulk molecule counts (integers)
- `uniqueMoleculesView(moleculeName)` — individual molecule instances (e.g., specific ribosomes)
- `environmentView(moleculeIDs)` — external environment concentrations

### State Partitioning

`BulkMolecules` uses integer counts. Processes request molecules and get partitioned allocations. Request priority levels control allocation order (e.g., metabolism=-10 gets priority, degradation=10 goes last).

`UniqueMolecules` tracks individual instances with attributes (e.g., chromosome fork positions, ribosome elongation progress).

### Data I/O

Simulation output uses a custom binary columnar format:
- **Write**: `TableWriter` (in `wholecell/io/tablewriter.py`) — zlib-compressed column files
- **Read**: `TableReader` (in `wholecell/io/tablereader.py`)

```python
from wholecell.io.tablereader import TableReader
reader = TableReader(os.path.join(simOutDir, 'BulkMolecules'))
counts = reader.readColumn('counts')
ids = reader.readAttribute('moleculeIds')
```

Each listener/state writes to `simOut/{Name}/` with column files + `attributes.json`.

### Analysis Plots

Analysis classes override `do_plot()` and are organized by scope:
- **single/** — one cell, one generation
- **cohort/** — multiple seeds, one variant
- **multigen/** — one lineage across generations
- **variant/** — compare variants
- **parca/** — parameter calculator validation

Each directory's `__init__.py` defines group TAGS (e.g., `CORE`, `ACTIVE`, `METABOLISM`). Run specific plots with `--plot aaCounts` or groups with `--plot CORE`.

### Variants

Variant functions in `models/ecoli/sim/variants/` modify `sim_data` and return metadata:

```python
def variant_name(sim_data, index):
    # Modify sim_data
    return {"shortName": "label", "desc": "..."}, sim_data
```

Common variants: `wildtype`, `gene_knockout`, `condition`, `timelines`, `ppgpp_conc`, `time_step`.

### Units System

Physical quantities use the `Unum` library (via `wholecell/utils/units.py`). All simulation data maintains units until computation, then converts via `.asNumber(units)`:

```python
from wholecell.utils import units
mass = 5 * units.fg
concentration = 100 * units.mmol / units.L
value = mass.asNumber(units.g)  # dimensionless float
```

### ParCa (Parameter Calculator)

Main entry: `reconstruction/ecoli/fit_sim_data_1.py` (~150KB). Orchestrates fitting steps:
raw data -> basal expression -> TF conditions -> condition fitting (CVXPY optimization) -> promoter binding -> final adjustments -> `out/.../kb/simData.cPickle`.

Auto-generates ODE dataclasses in `reconstruction/ecoli/dataclasses/process/`.

## Key Conventions

### Style

- **Indentation**: TABs (4-space tab stops) — this is a project-wide convention, not spaces
- **Naming**: `ClassName`, `function_name`, `GLOBAL_CONSTANT`; camelCase exists in older code and is acceptable but snake_case preferred for new code
- **Imports**: absolute imports from repo root (`from wholecell.sim.simulation import Simulation`)
- **Line length**: soft target 79, harder target 99, no hard limit
- **Docstrings**: triple-quoted, imperative first sentence. Use raw strings for LaTeX: `r"""..."""`
- See `docs/style-guide.md` for full guide

### Molecule IDs

Molecules include compartment tags: `atp[c]` (cytoplasm), `atp[p]` (periplasm), `atp[e]` (extracellular).

### Output Directory Structure

```
out/{sim_dir}/
  kb/simData.cPickle          # Fitted parameters
  {variant}_{index:06d}/      # e.g., wildtype_000000/
    {seed:06d}/               # e.g., 000000/
      generation_{gen:06d}/
        {cell:06d}/
          simOut/             # Listener data (TableWriter format)
          plotOut/            # Analysis figures
```

### Cython Extensions

Three performance-critical modules in `wholecell/utils/`:
- `_build_sequences.pyx` — polymerization sequence building
- `_fastsums.pyx` — fast monomer summation
- `mc_complexation.pyx` — Monte Carlo complexation

Must compile before running: `make compile`.

### Environment

Set `OPENBLAS_NUM_THREADS=1` to avoid threading artifacts in numerical results. Docker requires >=4GB RAM allocation.

### Docker

Two-layer build: `wcm-runtime` (Python + deps) -> `wcm-code` (source + compiled extensions). Dockerfiles in `cloud/docker/`.
