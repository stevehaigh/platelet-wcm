# Platelet Whole-Cell Model

A computational whole-cell model of a human platelet, built as a dissertation project.
Models intracellular calcium dynamics, granule release, and receptor signalling using a
discrete-timestep simulation framework.

> **Framework attribution.** The simulation engine (`wholecell/`) is derived from
> [CovertLab/wcEcoli](https://github.com/CovertLab/wcEcoli), the whole-cell model of
> *E. coli* developed by the [Covert Lab](https://www.covert.stanford.edu/) at Stanford.
> Please cite their work if you use anything derived from the framework:
>
> - Karr JR et al. (2012) *Cell* 150(2):389–401. https://doi.org/10.1016/j.cell.2012.05.044
> - Ahn-Horst TA et al. (2022) *npj Syst Biol Appl* 8:30. https://doi.org/10.1038/s41540-022-00242-9


## Quick-start

Python 3.11.5 via pyenv. Run all commands from the repo root.

```bash
# Install dependencies
pip install -r requirements.txt

# Run a 200-second agonist-stimulated simulation (activation transient)
PYTHONPATH=$PWD python runscripts/manual/runPlateletSim.py out/my_run --length 200

# Generate analysis plots
PYTHONPATH=$PWD python runscripts/manual/analysisPlatelet.py out/my_run

# Run tests
PYTHONPATH=$PWD python -m pytest models/platelet/tests/
```

The simulation supports several biologically-distinct conditions, controlled by
the extracellular Ca²⁺ flag and the per-agonist peak flags:

```bash
# Resting / un-stimulated (zero all agonist peaks)
python runscripts/manual/runPlateletSim.py out/resting --length 300 --at-rest

# Phase 1 transient with extracellular Ca²⁺ (Dolan 2014 nominal — default)
python runscripts/manual/runPlateletSim.py out/transient --length 200

# Phase 3 EDTA condition (no extracellular Ca²⁺ — SOCE inactive)
python runscripts/manual/runPlateletSim.py out/edta --length 200 --ca-ex-mM 0

# Dose-sweep variant — override per-receptor peak agonist concentrations
python runscripts/manual/runPlateletSim.py out/sweep --thrombin-peak-nM 3.0 --adp-peak-uM 1.0

# Phase 3 driver: runs both ±extracellular-Ca²⁺ conditions and produces
# the Dolan 2014 Fig. 4 comparison figure
python runscripts/manual/runPhase3.py out/phase3 --length 200
```

The same conditions are available in the webapp Configure tab as form fields
(Extracellular Ca²⁺ mM, "Run at rest" checkbox) and as four presets:
Agonist transient (+Ca²⁺), Agonist transient (60 s settle, +Ca²⁺),
EDTA transient, Resting.

Set `OPENBLAS_NUM_THREADS=1` in your shell profile for consistent numerical results.
See [`docs/create-pyenv.md`](docs/create-pyenv.md) for full environment setup.


## Code layout

```
models/platelet/           Platelet model
  processes/               Biological submodels (CalciumDynamics, GranuleSecretion, ...)
  listeners/               Data recorders (CalciumTrace, ...)
  analysis/single/         Post-simulation plots (calcium_trace.py, ...)
  sim/                     PlateletSimulation — wires processes and listeners
  tests/                   Regression and unit tests

reconstruction/platelet/   Parameter container
  simulation_data.py       SimulationDataPlatelet — constructed directly (no ParCa)
  initialization.py        Helpers that seed BulkMolecules at sim start
  dataclasses/             Per-process / per-state dataclasses + TOML/TSV loaders

wholecell/                 Generic simulation framework (from CovertLab/wcEcoli)
  sim/                     Simulation loop
  states/                  BulkMolecules, UniqueMolecules, LocalEnvironment
  processes/               Process base class
  listeners/               Listener base class + TableWriter/TableReader I/O
  utils/                   Units (Unum), math utilities

wholecell/webapp/          Browser-based web interface (Dash)
runscripts/manual/         Entry points: runPlateletSim.py, analysisPlatelet.py, webapp.py
reports/                   Lab books, design docs, figures, calibration data
docker/                    Dockerfiles for staging deployment
```


## Editing parameters

Calcium-pathway rate constants and the molecule inventory are externalised to
TSV/TOML files under `reports/params/` and loaded at import time. You do **not**
need to edit Python to change a value — edit the data file, re-run the sim.

| File | Purpose | Loader |
|------|---------|--------|
| `reports/params/calcium-v0.6.toml` | Rate constants, calibration scalars, agonist forcing peaks, and `[references.*]` bibliography for the calcium pathway | `reconstruction/platelet/dataclasses/process/_params_loader.py` |
| `reports/params/species-v0.6.tsv` | Molecule inventory: `id`, `mass_fg`, `initial_count`, `molecule_class` for all 83 species | `reconstruction/platelet/dataclasses/_species_loader.py` |

**Change a rate constant.** Open `calcium-v0.6.toml`, find the section
(e.g. `[serca.cycle]`), edit the value, re-run. The inline `# ...` comment
on each row is the per-parameter provenance / citation; update it too if the
new value comes from a different source.

**Add a new rate constant to an existing section.** Add a `key = value`
row in the TOML section, then add the corresponding `K_FOO['new_key']`
reference in `calcium_signalling.py` where the dict is consumed.

**Add a new receptor / sub-pathway within calcium.** (1) Add a new
`[section.subsection]` block to `calcium-v0.6.toml` with its rate
constants and inline citations. (2) Add a `K_FOO = dict(_KINETICS['section']['subsection'])`
line in `calcium_signalling.py` near the existing `K_*` block. (3) Wire it
into `_ode_rhs()`. (4) If the receptor adds a species, append a row to
`species-v0.6.tsv` (id with compartment tag, mass_fg, initial_count,
class). (5) Add a `[references.<key>]` block for any new citations and a
`match = [...]` list so the kinetics review auto-links them.

**Regenerate the clickable review PDF** (renders the TOML to
`reports/design/kinetics-v0.6-review.pdf` with auto-linked citations + a
BibTeX side-output at `reports/params/calcium-v0.6-references.bib`):

```bash
make kinetics-review        # needs quarto + xelatex on PATH
```

CI builds the same artifact on every PR — download it from the Actions run
under "Artifacts → kinetics-review".

> **Scope.** The kinetics-as-data scaffold is currently calcium-only.
> Extending to a non-calcium pathway (e.g. mitochondrial metabolism,
> cytoskeleton) needs a new `<pathway>-v0.N.toml`, a parallel
> `_<pathway>_loader.py`, and updates to `CHAPTER_TITLES` in
> `runscripts/manual/buildKineticsReview.py`. See
> `reports/design/kinetics-as-data-2026-05-22.qmd` for the
> post-dissertation generalisation sketch.


## Web interface

A [Dash](https://dash.plotly.com/) app for launching simulations and exploring results
without the command line.

```bash
# Run locally
make run        # starts webapp at http://localhost:8050/
make stop       # stop the process
```

Deployed to Azure Container Instances on push to the `webapp` branch.


## Infrastructure

| Concern | Tool |
|---------|------|
| CI | GitHub Actions — pytest + mypy on every push |
| Staging deploy | GitHub Actions → Azure Container Instances (`platelet-wcm.uksouth.azurecontainer.io`) |
| Python env | pyenv 3.11.5 + virtualenv |
| Output format | Custom binary columnar format (TableWriter/TableReader, zlib-compressed) |


## Documentation

- [`docs/README.md`](docs/README.md) — setup and run guide
- [`docs/create-pyenv.md`](docs/create-pyenv.md) — Python environment setup
- [`docs/dev-tools.md`](docs/dev-tools.md) — recommended dev tools
- [`docs/style-guide.md`](docs/style-guide.md) — code style conventions
- [`reports/design/`](reports/design/) — model design documents
- [`reports/lab-books/`](reports/lab-books/) — development journal
