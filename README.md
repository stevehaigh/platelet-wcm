# wcEcoli — Whole-Cell Model of *Escherichia coli*

> **Attribution.** This repository is a fork of
> [CovertLab/wcEcoli](https://github.com/CovertLab/wcEcoli), the whole-cell model of
> *Escherichia coli* developed by the
> [Covert Lab](https://www.covert.stanford.edu/) at Stanford University.
> The original authors retain all credit for the simulation engine, model biology,
> parameter fitting pipeline, and analysis infrastructure. Please cite their work if
> you use anything derived from this codebase:
>
> - Karr JR et al. (2012) *Cell* 150(2):389–401. https://doi.org/10.1016/j.cell.2012.05.044
> - Ahn-Horst TA et al. (2022) *npj Syst Biol Appl* 8:30. https://doi.org/10.1038/s41540-022-00242-9
>
> This fork diverges from the upstream project. It is maintained independently and is
> not intended as a contribution back to the original repository.


## What this fork adds

The primary addition is a **browser-based web interface** (the "webapp") that wraps the
existing wcEcoli simulation pipeline and makes it accessible without command-line expertise.
It is also the starting point for future whole-cell modelling work beyond *E. coli*.

All original simulation code, the parameter calculator (ParCa), the analysis scripts, and
the underlying framework are unchanged from upstream.


## The Web Interface

The webapp is a [Dash](https://dash.plotly.com/) application (Plotly's reactive web
framework, built on Flask). It runs as a single Python process serving a browser UI on
`localhost`, and interacts with the simulation code through the same Python API used by the
command-line runscripts.

### Four tabs

| Tab | Purpose |
|---|---|
| **Configure** | Set up and launch simulations. Quick-start presets (wildtype, nutrient shifts, ppGpp sweep, acetate), variant type, number of generations, random seeds, and regulation toggles. |
| **Run Status** | Live job queue. Shows each submitted job, its current phase (parameter fitting → simulating → analysing → done), and elapsed time. |
| **Inspect Data** | Interactive chart. Browse any listener and any data column from any completed run. Overlay multiple runs on the same graph, add/remove traces, and apply transforms (normalise, log scale). |
| **Explore Plots** | Side-by-side gallery of the 80+ pre-generated analysis plots. Compare two runs at once. |

### Running locally

```bash
# Fastest — run directly with hot-reload (no Docker needed)
make run

# Or inside Docker
make docker-run
```

Open [http://localhost:8050/](http://localhost:8050/).

### Docker build targets

```bash
make build        # Build wcm-webapp image (fast — no recompile)
make build-code   # Rebuild wcm-code when simulation code changes
make build-all    # Full rebuild from scratch (~30 min)
make stop         # Kill a directly-running webapp process
```

Simulations launched from the webapp run as Docker subprocesses using the `wcm-code` image,
so the webapp process itself stays responsive. Completed run data is written to `./out/` and
mounted into any Docker container.


## Code layout

```
wholecell/webapp/          Web interface (this fork's primary addition)
  app.py                   Dash application factory; wires all tabs and registers callbacks
  results.py               Data-access layer; scans out/ tree, wraps TableReader
  jobs.py                  SQLite-backed job queue; background thread runs simulations
  tabs/
    configure.py           Simulation setup form with presets
    runs.py                Live job-queue monitor
    inspect_data.py        Interactive listener data browser
    explore.py             Analysis plot gallery
  assets/style.css         Custom CSS

runscripts/manual/webapp.py   Entry point: parses args, calls create_app(), starts server
```

The original codebase layout is unchanged:

```
wholecell/           Generic simulation engine (organism-agnostic)
  processes/         Base classes for biological processes
  states/            BulkMolecules, UniqueMolecules, LocalEnvironment
  listeners/         Base listener class and I/O (TableWriter/TableReader)
  sim/               Simulation loop
  utils/             Units, Cython extensions, math utilities

models/ecoli/        E. coli-specific model
  processes/         17 process implementations
  listeners/         18 listener implementations
  analysis/          80+ analysis plot scripts (single/cohort/multigen/variant/parca)
  sim/variants/      ~32 variant functions

reconstruction/ecoli/  Parameter Calculator (ParCa)
  fit_sim_data_1.py  Main orchestrator; fits raw data into sim_data
  dataclasses/       Structured parameter objects

cloud/docker/        Dockerfiles (runtime → code → webapp, three-layer build)
.github/workflows/   CI/CD: builds and deploys webapp to Azure Container Instances
```


## Quick-start (command line)

Python 3.11. Run all commands from the repo root with `PYTHONPATH="$PWD"`.

```bash
# Compile Cython extensions (required before running)
make clean compile

# Fit parameters — run once per condition (~18 min)
python runscripts/manual/runParca.py out/my_run

# Run a simulation (~10 min)
python runscripts/manual/runSim.py out/my_run

# Run analysis plots
python runscripts/manual/analysisSingle.py out/my_run
```

All runscripts support `-h` for full option documentation.


## Infrastructure

- **CI/CD:** GitHub Actions builds Docker images and deploys to Azure Container Instances
  on push to the `webapp` branch.
- **Job runner:** Simulations are dispatched as Docker subprocesses (`steve-wcm-code`
  image), keeping the webapp process responsive.
- **Data format:** Simulation output uses a custom binary columnar format
  (`TableWriter`/`TableReader`) with zlib compression, unchanged from upstream.


## Future direction

This repository forms the starting point for whole-cell modelling work beyond *E. coli*.
The generic simulation engine in `wholecell/` is designed to be reused for other organisms.
That work is developed on a separate branch.


## Original documentation

The original Covert Lab documentation is preserved in [`docs/`](docs/README.md) and covers:
- Full setup instructions (Docker and pyenv)
- FireWorks workflow configuration
- Google Cloud Platform deployment
- Model architecture and conventions
