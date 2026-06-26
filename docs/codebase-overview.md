# Codebase overview

> A map of the repository for someone (human or AI) who needs to find where
> something lives and how the pieces connect. Pair this with
> [`architecture.md`](architecture.md) (the *why*) and the top-level
> [`README.md`](../README.md) / [`CLAUDE.md`](../CLAUDE.md) (commands + run
> conditions).

## Top-level layout

```
models/platelet/         The platelet model: processes, listeners, analysis, sim wiring
reconstruction/platelet/ Parameters: SimulationDataPlatelet + dataclasses + loaders
wholecell/               Generic engine inherited from CovertLab/wcEcoli (+ tui)
runscripts/manual/       Entry points (run sims, analysis, drivers, TUI)
reports/                 Lab books, design docs, params (TOML/TSV), figures, data, thesis
docs/                    Developer docs (this orientation set, style guide, env setup)
platelet-plan/           High-level plan + the wcEcoli engine architecture map
source-info/             Source PDFs (calcium-papers/) — summarised in docs/papers/
out/                     Simulation output (gitignored runs)
.github/                 CI workflows
```

## `models/platelet/` — the model

```
sim/
  simulation.py          PlateletSimulation(Simulation): declares state classes,
                         the process group, listeners, IC function
  initial_conditions.py  seeds BulkMolecules at sim start (reads sim.run_config)
processes/               biology that MUTATES state (run each timestep)
  resting_decay.py       slow turnover; the only stochastic process (numpy global RNG)
  calcium_dynamics.py    the rich one — thin wrapper around the Ca²⁺ ODE
  granule_secretion.py   PKC×Ca²⁺-gated cargo release + autocrine ADP loop + NTPDase clearance
  thromboxane_synthesis.py  Ca²⁺×PKC-gated TXA2 production (aspirin knob) → TXB2 decay
  integrin_activation.py    αIIbβ3 inside-out switch via Akt/Rap1b; PAC-1 readout
listeners/               OBSERVE state (read-only) and write columns to disk
  mass.py                dry/total mass
  calcium_trace.py       ~18 Ca²⁺/cascade columns — the headline data source
  secretion_trace.py     secreted cargo, gate, autocrine adp_e_uM
  thromboxane_trace.py   txa2/txb2, gate, tp_active_frac
  integrin_trace.py      active/resting αIIbβ3, active_frac, akt_active, rap1b_gtp
analysis/                POST-SIM plots (never run inside the sim)
  analysisPlot.py, singleAnalysisPlot.py   base classes
  phase3_dolan_fig4.py   the Dolan Fig.4 ±Ca_ex comparison figure
  single/
    calcium_trace.py     5-panel headline validation figure
    demo_*.py            per-theme TUI demo figures (calcium/integrin/thromboxane/
                         secretion/reversibility) sharing _demo_common.py
    granule_secretion.py, scaffold_summary.py
tests/                   see docs/development-workflow.md for the testing strategy
```

**Mental model:** *Processes write, Listeners read, Analysis plots.* A process
never plots; a listener never mutates; analysis only reads `simOut/`.

## `reconstruction/platelet/` — parameters & knowledge base

```
simulation_data.py       SimulationDataPlatelet — built directly (no ParCa)
run_config.py            RunConfig frozen dataclass — all per-sim conditions/knobs
initialization.py        initialize_bulk_molecules (applies count_overrides for KOs)
knockouts.py             logical-entity → species map for expression knockouts
dataclasses/
  constants.py           R, T, F, NA, structural integers
  internal_state.py      _MOLECULES (from species TSV), compartments, submasses
  external_state.py      extracellular environment
  molecule_groups.py
  _species_loader.py     loads reports/params/species-v0.6.tsv  → 83 species
  process/
    calcium_signalling.py   THE big one: _ode_rhs (pure fn) + all K_* rate dicts +
                            ODE state assembly. Consumes _KINETICS from the loader.
    _params_loader.py       load_calcium_kinetics() ← reports/params/calcium-v0.6.toml
    granule_secretion.py, thromboxane_synthesis.py, integrin_activation.py,
    resting_decay.py, process.py
```

### Parameters live as data

| File | Purpose | Loader |
|------|---------|--------|
| `reports/params/calcium-v0.6.toml` | calcium rate constants, calibration scalars, agonist peaks, `[references.*]` | `process/_params_loader.py` |
| `reports/params/species-v0.6.tsv` | molecule inventory (id, mass_fg, initial_count, class) for 83 species | `dataclasses/_species_loader.py` |

Edit the data file → re-run; no Python edit needed for a value change. **Scope
caveat:** this kinetics-as-data scaffold is **calcium-only**. The other processes
(secretion, thromboxane, integrin) carry their constants in their Python
dataclasses, not TOML.

## `wholecell/` — the generic engine (inherited)

Used essentially unchanged. Key entry points an assistant will touch:

```
sim/simulation.py        the Simulation base class (the loop)
processes/process.py     Process base class (calculateRequest/evolveState/views)
listeners/listener.py    Listener base class
states/                  bulk_molecules.py, unique_molecules.py, local_environment.py
views/view.py            the request/allocate/merge interface
io/tablewriter.py, tablereader.py   the columnar binary format
utils/units.py           Unum-based physical units
tui/                     Textual "experiment bench" TUI — wholecell/tui/
```

Don't put platelet biology in `wholecell/`. The only platelet-aware code here is
the TUI presentation layer.

## `runscripts/manual/` — entry points

The important ones (all support `-h`; all need `PYTHONPATH=$PWD`):

| Script | What it does |
|--------|--------------|
| `runPlateletSim.py` | run one sim. `run_platelet_sim(...)` is the importable entry used by tests |
| `analysisPlatelet.py` | generate analysis figures from a run (`--plot NAME`, `--out-name`) |
| `runPhase3.py` | the ±extracellular-Ca²⁺ Dolan Fig.4 two-condition driver |
| `runFromConfig.py` | run from a `{length_sec, seed, run_config}` JSON spec (the TUI backend) |
| `runDoseSweep.py` / `runDoseResponse.py` | ADP×thrombin 2-D sweeps → heatmaps/surface |
| `runPerturbation.py` | single-mechanism knockouts (PMCA/MCU/PKC P2Y1/PLCβ) |
| `runSecondWave.py` | autocrine amplification at a weak agonist (v0.6 vs aspirin vs full) |
| `runDemoExperiments.py` | the four write-up experiments (baseline/aspirin/MCU/clopidogrel) |
| `replayTui.py` / `runTui.py` | animated terminal replay / the experiment-bench TUI |
| `plotInhibitoryAxis.py`, `plotIntegrin.py`, `plot*.py` | bespoke figures |
| `buildKineticsReview.py` | render the clickable kinetics-review PDF from the TOML |
| `buildDocsSite.py` | render reports/*.{md,qmd} to the HTML docs site |

## Where the data and prose live (`reports/`)

- `reports/params/` — the TOML/TSV parameter files (+ `.bib` side-outputs).
- `reports/data/` — provenance & calibration data. **Start here:**
  `calcium-data-provenance.md` (verified value-by-value provenance) and
  `platelet-proteome-data-sources.md`.
- `reports/design/` — design docs (per pathway/feature). `reports/design/README.md`
  says which are current-canonical vs historical. Notable:
  `validation-map-2026-06-19.qmd` (the 5/5 rethink),
  `network-reference-2026-06-24.qmd` (the pathway map).
- `reports/lab-books/` — dated session notes; **the most recent is the current
  state of the work.**
- `reports/experiments/` — the demo-experiment write-ups.
- `reports/thesis/` — the dissertation draft + supplementaries.
- `reports/archive/` — superseded docs (kept for history).
- `platelet-plan/` — `overview.md` (the clone→prune→rebuild plan) and
  `architecture.md` (the detailed wcEcoli engine map with per-component
  reusability assessment).

## Conventions you must match

- **Indentation: TABs** (4-space stops) — project-wide, including Python.
- **Naming:** `ClassName`, `function_name`, `GLOBAL_CONSTANT`; some upstream
  camelCase remains, snake_case preferred for new code.
- **Imports:** absolute from repo root.
- **Line length:** soft 79, hard 99.
- **Figures:** matplotlib mathtext for chemistry (`$\mathrm{Ca}^{2+}$`), not raw
  unicode; detailed standalone captions/legends.
- Full guide: [`style-guide.md`](style-guide.md).

## Gotchas

- **`PYTHONPATH=$PWD`** and **`OPENBLAS_NUM_THREADS=1`** are required for correct
  imports and reproducible numerics.
- **Runs are not currently seed-reproducible** — `RestingDecay` draws from
  numpy's *global* RNG, not the sim seed.
- **Use the uv-managed venv, not bare `python3`** — system `python3` may lack the
  deps or be the wrong version. Run via `uv run python …` (or activate `.venv`);
  `make tui` already does. The interpreter is pinned to **3.11.5** via
  `.python-version`, which uv reads automatically.
- The replay TUI's deps (`rich`, `textual`, `textual-plotext`, `plotext`) are an
  **optional extra** (the `viz` extra); model smoke tests `importorskip`.
- "Loop / knockout / perturbation effects are usually invisible under the default
  saturating agonist" (the response is store-limited) — isolate one agonist and
  read IP₃, or use the baseline overlay, to see them.
