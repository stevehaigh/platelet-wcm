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

# Run a 60-second simulation
PYTHONPATH=$PWD python runscripts/manual/runPlateletSim.py out/my_run --length 60

# Generate analysis plots
PYTHONPATH=$PWD python runscripts/manual/analysisPlatelet.py out/my_run

# Run tests
PYTHONPATH=$PWD python -m pytest models/platelet/tests/
```

Set `OPENBLAS_NUM_THREADS=1` in your shell profile for consistent numerical results.
See [`docs/create-pyenv.md`](docs/create-pyenv.md) for full environment setup.


## Code layout

```
models/platelet/           Platelet model
  processes/               Biological submodels (CalciumDynamics, Scaffold, ...)
  listeners/               Data recorders (CalciumListener, ...)
  analysis/single/         Post-simulation plots (calcium_trace.py, ...)
  sim/                     PlateletSimulation — wires processes and listeners
  tests/                   Regression and unit tests

reconstruction/platelet/   Parameter calculator
  sim_data.py              SimulationDataPlatelet — fitted parameters
  fit_sim_data.py          Parameter fitting pipeline

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
| Staging deploy | GitHub Actions → Azure Container Instances (`wcecoli-webapp.uksouth.azurecontainer.io`) |
| Python env | pyenv 3.11.5 + virtualenv |
| Output format | Custom binary columnar format (TableWriter/TableReader, zlib-compressed) |


## Documentation

- [`docs/README.md`](docs/README.md) — setup and run guide
- [`docs/create-pyenv.md`](docs/create-pyenv.md) — Python environment setup
- [`docs/dev-tools.md`](docs/dev-tools.md) — recommended dev tools
- [`docs/style-guide.md`](docs/style-guide.md) — code style conventions
- [`reports/design/`](reports/design/) — model design documents
- [`reports/lab-books/`](reports/lab-books/) — development journal
