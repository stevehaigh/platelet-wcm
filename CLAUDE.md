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

# 2-D dose-response sweep over ADP × thrombin (one row per cell → heatmaps + 3-D surface)
PYTHONPATH=$PWD python runscripts/manual/runDoseSweep.py [sim_outdir] --grid 9 --length 200

# Run the Dash webapp locally (http://localhost:8050)
make run     # foreground with hot reload
make stop    # kill it
```

All runscripts support `-h` for full options.

### Run-time conditions

`runPlateletSim.py` exposes these CLI flags for the biology being
simulated, beyond simulation length and seed:

| Flag | Default | Effect |
|------|---------|--------|
| `--length N` | 60 | Simulation length in seconds |
| `--seed N` | 0 | RNG seed (currently no stochastic processes use it) |
| `--ca-ex-mM X` | 1.2 | Extracellular Ca²⁺ in mM. Set `0` for the Dolan Fig. 4 EDTA / no-extracellular-Ca²⁺ condition (SOCE inactive, PM leak inactive). |
| `--at-rest` | (off) | Shorthand for `--thrombin-peak-nM 0 --adp-peak-uM 0 --atp-ex-peak-uM 0`. All agonists stay at REST level → cell sits at its endogenous fixed point. |
| `--thrombin-peak-nM X` | (module default: 1.0) | Peak thrombin (nM) during the activation transient. Drives PAR1/PAR4. |
| `--adp-peak-uM X` | (module default: 10.0) | Peak ADP (µM) during the transient. Drives P2Y1. |
| `--atp-ex-peak-uM X` | (module default: 10.0) | Peak extracellular ATP (µM) during the transient. Drives P2X1. |
| `--agonist-delay S` | 0.0 | Seconds the model settles at its fixed point before the agonist time courses start. Useful for ignoring the start-up transient. |

Where each behaviour is defined in code:

- **Extracellular Ca²⁺** is the module-level constant `CA_EX_UM` in
  `reconstruction/platelet/dataclasses/process/calcium_signalling.py`. The
  runscript overrides it via `cs_mod.CA_EX_UM = ca_ex_mM * 1000.0` before
  the sim is constructed; `_ode_rhs` reads it on every ODE step. Both the
  SOCE current and the basal PM Ca²⁺ leak are gated on `CA_EX_UM > 0`.
- **Agonist stimulation** is three optional class attributes on
  `CalciumDynamics` in `models/platelet/processes/calcium_dynamics.py`:
  `_thrombin_peak_nM`, `_adp_peak_uM`, `_atp_ex_peak_uM`. `None` (default)
  → the module-level peak constants are read live at call time
  (`THROMBIN_PEAK_NM`, `ADP_PEAK_UM`, `ATP_EX_PEAK_UM`). `0` → that
  receptor sees only its REST level (a "resting" sim has all three set to
  0). The agonist functions (`thrombin_nM`, `adp_uM`, `atp_ex_forcing_uM`)
  use the `peak_X=None` sentinel pattern so module-constant reassignment
  takes effect immediately — required for dose-sweep work (issue #45).

The same conditions are exposed on the **webapp** Configure tab as form
fields (Extracellular Ca²⁺ mM, "Run at rest" checkbox). The four webapp
presets — Agonist transient (+Ca²⁺), Agonist transient (60 s settle),
EDTA transient, Resting — are defined in
`wholecell/webapp/tabs/configure.py:PRESETS` and differ across four
biology-affecting knobs: `ca_ex_mM`, `at_rest`, `agonist_delay_s`, and
`length_sec` (the last sets how much of the response is observed; the
first three set what biology runs). The Phase 3 driver `runPhase3.py`
runs the +Ca²⁺ and EDTA conditions back-to-back.

## Tests & Type Checking

```bash
# All platelet tests
PYTHONPATH=$PWD python3 -m pytest models/platelet/tests/ -v

# Fast iteration — skip simulation-running tests (~3 s vs ~24 s)
PYTHONPATH=$PWD python3 -m pytest models/platelet/tests/ -m "not slow"

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
`(RestingDecay, CalciumDynamics, GranuleSecretion, ThromboxaneSynthesis)`.

### Core Abstractions

| Concept | Base class | Platelet impl | Purpose |
|---------|-----------|---------------|---------|
| **Process** | `wholecell/processes/process.py` | `models/platelet/processes/` (RestingDecay, CalciumDynamics, GranuleSecretion, ThromboxaneSynthesis) | Biological submodels that modify state |
| **State** | `wholecell/states/internal_state.py` | `BulkMolecules`, `UniqueMolecules`, `LocalEnvironment` | Cellular state containers |
| **Listener** | `wholecell/listeners/listener.py` | `models/platelet/listeners/` (Mass, CalciumTrace, SecretionTrace, ThromboxaneTrace) | Observe and record data each timestep |
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

- **GPCR cascade** — P2Y1 (ADP, reversible), PAR1 (thrombin, fast irreversible cleavage),
  PAR4 (thrombin, ~10× slower cleavage), P2X1 (ATP-gated ionotropic). Receptors drive a
  Gαq exchange/GTPase cycle which activates PLCβ.
- **PI cycle** — PLCβ-catalysed PIP2 → IP3 + DAG; IP3 5-phosphatase / 3-kinase degradation;
  PIP2 resynthesis (lumped). IP3 is now an endogenous state variable, not a forced curve.
- **PKC feedback (v0.6)** — DAG + Ca²⁺ activate PKC (lumped conventional PKCα/β + novel
  PKCδ), closing the previously dead-end DAG branch, via two activity-dependent brakes:
  (1) **P2Y1 desensitisation** (`P2Y1_active[pl]` → `P2Y1_desensitised[pl]`; ADP-arm
  specific, Mundell 2006 / Nicholas 2023) and (2) **PLCβ phosphorylation** (`PLCb_inactive`
  → `PLCb_phosphorylated[c]`, out of the Gq-activatable pool; shared-node brake that lowers
  IP3 toward baseline, Purvis 2008). Both engage after the early Ca²⁺ peak (~10–15 s PKC
  delay), so Dolan 5/5 is preserved; the Ca²⁺ response is store-limited, so the receptor
  desensitisation fraction and IP3 are the clear readouts. See `[pkc.*]` in the calcium TOML
  and the `pkc` / `plcb` experiments in `runPerturbation.py`.
- **IP3R** — Li-Rinzel 1994 reduction of de Young–Keizer 1992 (quasi-steady m∞, one slow ODE for h)
- **SERCA** — E1/E2 enzymatic cycle (Dode 2002)
- **PMCA** — 5-state CaM-coupled (Caride 2007 Table 3): basal + Ca₄·CaM-activated paths
- **SOCE** — STIM1 dimerisation + Orai1 MWC flux (Dolan & Diamond 2014; Hoover)
- **CaM** — three sub-states (free, Ca₂·CaM, Ca₄·CaM) acting as cytosolic Ca²⁺ buffer
- **Stimulus input** — agonist time courses (`thrombin_nM`, `adp_uM`, `atp_ex_forcing_uM`)
  whose peaks are controlled per-call by the `peak_*` kwargs (with `None` →
  read the module default constant live). Passing 0 for a given peak gives
  REST level for that receptor; all three zero is a resting / un-stimulated sim.

`CalciumTrace` listener records 18 columns (Ca²⁺ pools, CaM/PMCA sub-states, IP3, SOCE flux, SERCA+PMCA ATP cost, PKC active, P2Y1 desensitised fraction + PLCβ phosphorylated fraction).
The 5-panel `single/calcium_trace.py` plot is the headline validation figure.

Validation target: Dolan & Diamond 2014 Fig. 4 (Ca²⁺ transients with/without
extracellular Ca²⁺).

### Downstream PKC effects — granule secretion + autocrine ADP (v0.61)

`GranuleSecretion` (in `models/platelet/processes/granule_secretion.py`) is the
first PKC *output* (v0.61), wiring PKC out of its v0.6 brake-only role.

**Slice 1 — secretion.** Each timestep it relocates pre-existing granule cargo —
`ADP[dg]`, `5HT[dg]` (dense), `FGA[ag]` (α) → the extracellular space `[e]`, and
`SELP[ag]` → a surface state `SELP_surface[pl]` (the P-selectin activation
marker). Release is first-order in the remaining pool, scaled by a
`PKC_active × Ca²⁺` coincidence gate that keys off PKC activation *above* a
resting-tone floor, so resting secretion is exactly zero.

**Slice 2 — autocrine ADP loop.** Secreted `ADP[e]` is fed back onto the P2Y1
drive inside the calcium ODE: `_ode_rhs` adds its pericellular concentration
(`secreted_adp_count × _UM_PER_COUNT_EX`, threaded via CalciumDynamics) to the
exogenous ADP forcing, closing PKC → secretion → ADP → P2Y1. `V_EX_L` (effective
pericellular volume, ~66 fL) is a **calibration choice** set so full dense-granule
release ≈ 10 µM (the standard dose). The loop self-limits via ecto-NTPDase
clearance (`ADP[e] → AMP[e]`, first-order `k_ntpdase`, in GranuleSecretion) plus
the v0.6 P2Y1 desensitisation brake and finite cargo. Effect is sub-integer on the
30 s Dolan goldens (P2Y1 is minor vs thrombin/PARs and the response is
store-limited) → **goldens stay byte-identical, Dolan 5/5 preserved, no regen**;
it shows clearly in a thrombin-only sim (zero exogenous ADP) where secreted ADP is
the sole P2Y1 driver.

Rate constants / volumes live in
`reconstruction/platelet/dataclasses/process/granule_secretion.py` and the
volume block of `calcium_signalling.py` (Python, not TOML — the kinetics-as-data
scaffold is still calcium-only). The `SecretionTrace` listener records
secreted-cargo counts, released / surface-exposed fractions, the gate, and the
autocrine `adp_e_uM`.

**Thromboxane Slice A — TXA₂ synthesis (§2, production only).**
`ThromboxaneSynthesis` (in `models/platelet/processes/thromboxane_synthesis.py`)
lumps cPLA₂ → COX-1 → TXA₂-synthase into one Ca²⁺ × PKC-gated production term
(same resting-floor gate → zero at rest), scaled by the **aspirin knob**
`COX1_FACTOR` (module-level in the dataclass, read live; `0` = aspirin knockout,
abolishes TXA₂). De-novo `TXA2[e]` decays first-order (t½ ≈ 30 s) to the stable
ELISA metabolite `TXB2[e]`. The `ThromboxaneTrace` listener records `txa2_uM`,
`txb2`, and the gate. **Additive — no Gq feedback yet, calcium ODE untouched,
goldens byte-identical.** The autocrine **TXA₂ → TP → Gq** loop (Slice B —
`+ tp_a` into the `total_active_R` sum, regen goldens + re-verify Dolan) and
integrin (§3) remain unimplemented. Design:
`reports/design/pkc-downstream-effects-2026-06-12.qmd` §1–2.

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

Rate constants for calcium signalling are externalised to
`reports/params/calcium-v0.6.toml` (issue #32 Phase 2; v0.6 adds `[pkc.*]`)
and loaded at import time by
`reconstruction/platelet/dataclasses/process/_params_loader.py:load_calcium_kinetics()`.
`calcium_signalling.py` consumes the loaded dict as `_KINETICS = load_calcium_kinetics()`
and assigns the remaining ODE state / per-channel scalars; physical constants
(R, T, F, NA), structural integers, and compartment volumes stay in Python.

The molecule inventory (id, mass, initial count, class for all 74 species)
lives in `reports/params/species-v0.6.tsv` and is loaded by
`reconstruction/platelet/dataclasses/_species_loader.py:load_species()`,
exposed in `internal_state.py` as `_MOLECULES`. There is no `raw_data/`
directory in this repo.

**Scope note.** The kinetics-as-data scaffold is currently calcium-only:
the loader, the TSV, and `runscripts/manual/buildKineticsReview.py` are all
hardcoded to the calcium pathway. Extending to a second pathway (e.g.,
mitochondrial metabolism, integrin signalling, cytoskeleton) would require
a new `<pathway>-v0.N.toml`, a parallel `_<pathway>_loader.py`, registration
in whichever process consumes it, and updates to `CHAPTER_TITLES` in the
review renderer so its sections render. None of that scaffolding exists yet
(see `reports/design/kinetics-as-data-2026-05-22.qmd` "Level 2" sketch).

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
- `reports/data/`, `reports/figures/` — calibration data and rendered figures (figures are snapshotted per model release under a version subdir, e.g. `reports/figures/v0.5/`)
- `make pdfs` builds PDFs of any `reports/*.md` into `reports/pdf/` (pandoc + xelatex)
- `make quarto-pdfs` builds PDFs of any `reports/*.qmd` into `reports/pdf-quarto/` (Quarto + xelatex)
- `runscripts/manual/buildDocsSite.py` renders every `.md` / `.qmd` under `reports/{design,lab-books,data,decks,external}` to HTML and writes an auto-generated index at `reports/site/index.html` (lab-book and design-doc listings are produced automatically from the file headers)
- Published site: `https://stevehaigh.github.io/platelet-wcm/reports/site/` — manual publish for now; no GitHub Actions workflow drives it yet
- `reports/design/README.md` — navigator for the design-doc directory (what's current canonical vs historical)
- `docs/style-guide.md` — full style guide

**.md vs .qmd — which to use:**

- **`.md` (pandoc)** — default for prose-only design docs and lab books. No new tool dependencies; build pipeline is the established `make pdfs` target.
- **`.qmd` (Quarto)** — use for **diagram-heavy** docs and anything that benefits from `quarto preview` live-reload during editing. Native mermaid → PDF rendering (no `mmdc` step). Same xelatex pipeline under the hood, same fonts and margins as the pandoc rule, so existing reports/design/ visual style is preserved.
- Quarto frontmatter for this repo: see `reports/design/model-status-2026-05-13.qmd` for the canonical block (format → pdf with `include-in-header: ../pandoc-header.tex` to inherit the project's LaTeX header).
- During iteration: `quarto preview reports/design/<doc>.qmd` watches the file and reloads HTML on save. For PDF iteration, `quarto render <doc>.qmd --to pdf` is faster than going through make.
