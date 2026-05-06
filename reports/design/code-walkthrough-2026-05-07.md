# Code walkthrough cheatsheet — platelet-wcm — 2026-05-07

Companion to `demo-2026-05-07.md`. That doc is for "watch the model
run." This one is for "see how it's built." Use it to navigate the
codebase live during the walkthrough.

---

## The 30-second mental model

Two layers, intentionally separate:

| Layer | Where | What |
|---|---|---|
| **Framework** | `wholecell/` | Generic simulation engine inherited from CovertLab/wcEcoli. Process / state / listener base classes, the time-stepping loop, columnar I/O, the Dash webapp. Untouched by the platelet biology. |
| **Platelet model** | `models/platelet/` + `reconstruction/platelet/` | All the platelet biology — calcium signalling, the molecule inventory, parameters, analysis plots, tests. New code, ~3 500 LOC. |

If a question is about *how a sim runs* (the loop, partitioning, I/O),
look in `wholecell/`. If it's about *what the sim does* (the biology,
parameters, plots), look in `models/platelet/` and
`reconstruction/platelet/`.

---

## Top-level directory map

```
platelet-wcm/
├── wholecell/                      # framework (inherited, mostly untouched)
│   ├── sim/simulation.py           # the time-stepping loop
│   ├── processes/process.py        # Process base class
│   ├── listeners/listener.py       # Listener base class
│   ├── states/                  # BulkMolecules, UniqueMolecules, LocalEnv
│   ├── views/view.py               # how processes ask states for molecules
│   ├── io/tablereader.py           # columnar I/O (zlib + JSON attrs)
│   ├── io/tablewriter.py
│   └── webapp/                     # Dash UI (4 tabs)
│       ├── app.py                  # Dash app factory
│       ├── jobs.py                 # SQLite job queue + subprocess runner
│       ├── results.py              # filesystem scan of out/
│       └── tabs/
│           ├── configure.py        # 3 biologically-distinct presets + form
│           ├── runs.py             # job-status polling table
│           ├── inspect_data.py     # listener data browser
│           └── explore.py          # analysis plot browser
│
├── models/platelet/                # platelet model
│   ├── sim/
│   │   ├── simulation.py        # PlateletSimulation: wires processes/listeners
│   │   └── initial_conditions.py   # sets initial counts on bulk views
│   ├── processes/
│   │   ├── resting_decay.py # protein decay (Burkhart 2012; inert at 200 s)
│   │   └── calcium_dynamics.py     # * the headline process (thin wrapper)
│   ├── listeners/
│   │   ├── mass.py                 # dryMass / proteinMass etc.
│   │   └── calcium_trace.py        # * 14-column Ca²⁺ recording
│   ├── analysis/
│   │   ├── single/                 # one-sim plots
│   │   │   ├── calcium_trace.py    # * 5-panel headline figure
│   │   │   └── scaffold_summary.py
│   │   └── phase3_dolan_fig4.py    # * Phase 3 two-condition comparison
│   └── tests/
│       ├── sim/
│       │   ├── test_simulation.py  # boots the engine
│       │   └── test_regression.py  # 60 s golden-run tolerance suite
│       └── analysis/
│           ├── test_analysis.py    # plot-pipeline smoke
│           └── test_phase3.py      # Phase 3 baseline lock (4 tests)
│
├── reconstruction/platelet/        # platelet "knowledge base"
│   ├── simulation_data.py       # SimulationDataPlatelet (parameter object)
│   ├── raw_data/molecules.tsv      # canonical molecule list
│   └── dataclasses/
│       ├── constants.py
│       ├── molecule_groups.py
│       ├── internal_state.py    # 27-species inventory + initial conditions
│       ├── external_state.py
│       └── process/
│           └── calcium_signalling.py    # *** THE biology file (ODE+rates)
│
├── runscripts/manual/
│   ├── runPlateletSim.py           # * single-condition CLI entry point
│   ├── runPhase3.py                # * Phase 3 two-condition driver
│   ├── analysisPlatelet.py         # invokes the analysis plots
│   └── webapp.py                   # webapp entry point
│
├── reports/                        # docs + artefacts (per-session)
│   ├── design/   # design docs (calcium-dynamics-design.md = as-built ref)
│   ├── data/                       # provenance + saved JSON snapshots
│   ├── figures/                    # rendered figures (incl. Phase 3 result)
│   └── lab-books/                  # dated session journals
│
├── docs/                        # contributor docs (style, env, dev tools)
└── docker/                         # runtime + webapp images
```

Stars (★) mark the files most likely to come up in walkthrough Q&A.

---

## End-to-end trace: one sim, file by file

This is the 8-stop tour to use during the walkthrough — each stop is
one file, opened in order.

### 1. `runscripts/manual/runPlateletSim.py` (entry point)

The user's CLI starts here. Notable bits:
- Builds an `argparse` parser: `--length`, `--seed`, `--ca-ex-mM`,
  `--no-ip3-forcing`. The last two are the run-time biology toggles.
- Creates a `SimulationDataPlatelet()` instance and persists it to
  `kb/` for downstream analysis to load.
- Constructs `PlateletSimulation(simData=…, lengthSec=…, …)` and
  calls `sim.run()`.
- Two module-attribute overrides happen at the top of
  `run_platelet_sim` *before* the simulation is constructed:
  - `cs_mod.CA_EX_UM = float(ca_ex_mM) * 1000.0` — flips
    extracellular Ca²⁺ for the EDTA condition.
  - `CalciumDynamics._ip3_forced = bool(ip3_forced)` — flips the IP3
    forcing curve on/off for resting vs stimulated.
  Module-attribute mutation is intentional: it's the simplest way to
  flip a one-shot run without restructuring the process / sim_data
  contract. Both flags are also recorded in `metadata.json`.

### 2. `reconstruction/platelet/simulation_data.py` (the knowledge base)

`SimulationDataPlatelet` is the equivalent of wcEcoli's
`SimulationDataEcoli`, but tiny — it just wires together a few
dataclasses and exposes:
- `submass_name_to_index` (`'protein'` / `'metabolite'`)
- `compartment_abbrev_to_index` (`c, dts, dg, ag, m, e, pl`)
- `.constants`, `.molecule_groups`, `.internal_state`, `.external_state`,
  `.process`

No ParCa: the sim_data is built directly in code, not fitted from raw
data.

### 3. `reconstruction/platelet/dataclasses/internal_state.py`

The species table — one row per molecule with `(id, mass_fg,
initial_count, class)`. 27 calcium-pathway species + a handful of
metabolites + a few placeholder proteins (granule cargo, integrins).

Two interesting details to point out:
- The SERCA `E1 / E1·Ca` counts are 2963 / 2963 (Phase 2a (iii)
  pre-equilibration), *not* Dolan's verbatim 5920 / 6 — comment
  in-line explains why.
- The CaM ladder initial conditions are derived from detailed-balance at 100 nM cyt,
  not Dolan Table S1.

### 4. `reconstruction/platelet/dataclasses/process/calcium_signalling.py` ★★★

**This is the biology file.** Almost all the dissertation-relevant
maths lives here:

- Module-level constants up top: volumes, `N_A`, `V_IM` / `V_PM`,
  `CA_EX_UM`, `J_PM_LEAK_IONS_S`, `IP3` forcing parameters.
- `K_IP3R`, `K_SERCA`, `K_PMCA`, `K_CAM`, `K_CAM_PMCA`, `K_STIM`,
  `K_MWC`, `PUNCTA` — every rate-constant table, with primary-source
  attribution in comments.
- `_phi_*` helper functions implement the Sneyd & Dufour 2002
  φ-function rate laws for IP3R.
- `_mwc_open_fraction` — Hoover & Lewis 2011 MWC solver.
- `ip3_forcing_uM(t)` — the Dolan Fig. S2 IP3 time curve.
- **`_ode_rhs(t, y, t_sim_start, ip3_forced)`** — the heart of the
  model. ~200-line function that builds `dy` term by term: IP3R flux,
  SERCA cycle, PMCA basal + CaM-activated, CaM ladder, STIM cycle,
  MWC SOCE, PM leak, IP3 forcing.
- `class CalciumSignalling` — wraps `_ode_rhs` with the per-timestep
  entry point `molecules_to_next_time_step(counts, dt, t_sim,
  ip3_forced)`.

If asked "where is the Po⁴ tetramer cooperativity?", point at lines
~493–509. If asked "where is the SOCE current?", lines ~647–662. If
asked "where is the EDTA fix?", line ~652 (the `CA_EX_UM > 0` guard).

### 5. `models/platelet/sim/simulation.py` (the wiring)

40 lines. `PlateletSimulation` subclasses
`wholecell.sim.simulation.Simulation` and declares:
- `_internalStateClasses = (BulkMolecules, UniqueMolecules)`
- `_externalStateClasses = (LocalEnvironment,)`
- `_processClasses = ((RestingDecay, CalciumDynamics),)`
- `_listenerClasses = (Mass, CalciumTrace)`
- `_initialConditionsFunction = calcInitialConditions`

That's it — the framework does the rest.

### 6. `models/platelet/processes/calcium_dynamics.py`

Thin wrapper, ~70 lines. The `CalciumDynamics(Process)` class:
- `initialize`: gets `BulkMoleculesView`s over the 27-species ODE state
  + `ATP[c]`.
- `calculateRequest`: calls
  `sim_data.process.calcium_signalling.molecules_to_next_time_step(
  counts, dt, t_sim, ip3_forced=True)`. Returns the count deltas and
  ATP cost. Requests from the partitioner.
- `evolveState`: applies the deltas, debits ATP.

The biology is **not here** — this file is purely the framework
contract. Useful to point out: this is what makes the biology testable
in isolation (`molecules_to_next_time_step` can be called from a
script, no PlateletSimulation needed — used by the lab-book sweeps).

### 7. `models/platelet/listeners/calcium_trace.py` (the observer)

Records 14 columns per timestep. Notable:
- Imports `cs_mod` (the calcium_signalling module) rather than
  `from … import CA_EX_UM` so the runtime EDTA override propagates
  here too.
- Re-derives SOCE flux for logging (the ODE is the source of truth;
  the listener just needs the same number for the analysis plot).

### 8. `models/platelet/analysis/single/calcium_trace.py` (the plot)

Reads `simOut/CalciumTrace/*` columns + `simOut/BulkMolecules/counts`,
produces the 5-panel headline figure. The Dolan reference band in
panel 1 is the analytical schematic baked into the file (functions
`_dolan_reference_nM` / `_dolan_reference_band_nM`).

### Phase 3 fork: `runPhase3.py` + `analysis/phase3_dolan_fig4.py`

For the Phase 3 walkthrough, replace stops 1, 8 with:

- **`runscripts/manual/runPhase3.py`** — runs `run_platelet_sim` twice
  with `ca_ex_mM=1.2` and `ca_ex_mM=0.0`, then calls the comparison
  plot and writes `phase3_summary.json`.
- **`models/platelet/analysis/phase3_dolan_fig4.py`** — loads
  `reports/data/dolan-2014-fig4-reference.json`, builds the 3-panel
  comparison figure, evaluates the 5 acceptance criteria, returns a
  JSON summary.

---

## "Where do I look for…?" — the index

| Question | File |
|---|---|
| The 27-species ODE state | `reconstruction/platelet/dataclasses/internal_state.py` |
| The IP3R rate constants | `…/calcium_signalling.py` `K_IP3R` |
| The SERCA cycle equations | `…/calcium_signalling.py` `_ode_rhs`, search "SERCA" |
| The MWC SOCE model | `…/calcium_signalling.py` `_mwc_open_fraction` + `K_MWC` |
| The IP3 forcing curve | `…/calcium_signalling.py` `ip3_forcing_uM` |
| The PM leak | `…/calcium_signalling.py` `J_PM_LEAK_IONS_S` (constant) and `_ode_rhs` (term applied with EDTA guard) |
| Why our SERCA E1/E1·Ca differ from Dolan | `internal_state.py` SERCA block; design doc §6.8 D5 |
| Why we depart from `γ_SOC=24 fS` | `calcium_signalling.py` near `GAMMA_SOC_S`; design doc §6.8 D3 |
| The 9 deviations from primary sources | design doc §6.8 |
| The 5-panel plot code | `models/platelet/analysis/single/calcium_trace.py` |
| The Phase 3 comparison plot | `models/platelet/analysis/phase3_dolan_fig4.py` |
| The Phase 3 reference values | `reports/data/dolan-2014-fig4-reference.json` |
| The webapp Configure form | `wholecell/webapp/tabs/configure.py` |
| The webapp jobs queue | `wholecell/webapp/jobs.py` |
| The Phase 1 lab-book diagnosis | `reports/lab-books/lab-book-2026-05-01-phase1-complete.md` |
| The Phase 2a γ-sweep + rate balance | `reports/lab-books/lab-book-2026-05-05-phase2a-investigation.md` |
| The Phase 3 result write-up | `reports/lab-books/lab-book-2026-05-06-phase3-results.md` |

---

## Run-time conditions: what gets flipped, where

Three knobs change the *biology* the model simulates (length and seed
don't). Each follows the same pattern: a CLI flag on
`runPlateletSim.py` overrides a module-level constant or class
attribute *before* the simulation is constructed; the ODE solver and
listener read it on every timestep.

| Knob | CLI flag | Default | Defined at | Read by |
|---|---|---|---|---|
| Extracellular Ca²⁺ | `--ca-ex-mM X` | 1.2 mM (Dolan nominal) | `cs_mod.CA_EX_UM` (module constant) in `reconstruction/platelet/dataclasses/process/calcium_signalling.py` | `_ode_rhs` SOCE current + PM leak block (both gated on `CA_EX_UM > 0`); the `CalciumTrace` listener for the recorded SOCE flux column |
| IP3 forcing | `--no-ip3-forcing` | ON (Dolan Fig. S2) | `CalciumDynamics._ip3_forced` (class attr) in `models/platelet/processes/calcium_dynamics.py` | passed each step into `cs.molecules_to_next_time_step(... ip3_forced=...)`; ODE applies `ip3_forcing_uM(t)` to the IP3 state when True |
| Webapp presets | (form) | preset-dependent | `wholecell/webapp/tabs/configure.py:PRESETS` — list of dicts, each with `length_sec`, `ca_ex_mM`, `ip3_forced`, `description` | the Configure-tab callback fills the form fields; `jobs.py` translates them to the runscript flags above |

End-to-end on a webapp click: **Configure tab preset** → form field
values → **submit** → `jobs.py` builds CLI args including
`--ca-ex-mM <X>` and (if off) `--no-ip3-forcing` → subprocess runs
`runPlateletSim.py` → `cs_mod.CA_EX_UM` and
`CalciumDynamics._ip3_forced` overridden in `run_platelet_sim` →
`PlateletSimulation` constructed with those module-level values
already in place → ODE reads them throughout the run → metadata.json
records both for traceability.

If asked **"why module-level attribute mutation?"**: the alternative —
threading these values through `SimulationDataPlatelet` and the process
constructor — would touch four files for what is, today, a one-line
runscript override. The Phase 3 work demonstrated the pattern; if the
v0.3 receptor cascade introduces dependent state these will probably
move into the sim_data object.

---

## Inherited framework: what to point at if asked

The wcEcoli engine pieces we reuse, in order of "you'll definitely
get asked about this":

1. **`wholecell/sim/simulation.py`** — the simulation loop.
   Per-timestep flow: `_pre_evolve_state` → for each process group:
   `updateQueries` → `calculateRequest` → `partition` → `evolveState`
   → `wasTimeStepShortEnough` → `merge` → `external_state.update` →
   `_post_evolve_state`. Listeners run last.

2. **`wholecell/processes/process.py`** — base class with
   `calculateRequest` / `evolveState` / `bulkMoleculesView` / etc.

3. **`wholecell/states/bulk_molecules.py`** — the integer-count
   state with priority-based partitioning. ATP[c] currently has no
   competing demanders in the platelet model, but the partitioning
   *would* kick in if it did.

4. **`wholecell/io/{tablewriter,tablereader}.py`** — columnar binary
   format with zlib chunks + JSON attrs. Each listener writes one
   directory under `simOut/`.

If asked "what's E. coli left?" — short answer: nothing biological;
the framework is generic by design.

---

## Test-suite tour (if it comes up)

- `models/platelet/tests/sim/test_simulation.py` — boots PlateletSimulation
- `models/platelet/tests/sim/test_regression.py` — locks 60 s golden-run
  numbers (peak Ca²⁺ band, dryMass, etc.)
- `models/platelet/tests/analysis/test_analysis.py` — plot-pipeline smoke
- `models/platelet/tests/analysis/test_phase3.py` — Phase 3 baseline lock
- `wholecell/tests/{containers,io,states,utils}/...` — framework-side tests

Run all: `pytest models/platelet/tests/ -q` (current count: 21 platelet,
89 framework, 110 total).

---

## "Why did you make X choice?" — three answers worth memorising

1. **"Why deterministic ODE, not Gillespie?"** — Cytosolic Ca²⁺ is at
   the continuum-limit boundary (~361 ions at rest), but the SERCA
   cycle has 1 000 s⁻¹ transitions across 11 892 enzymes — Gillespie
   would fire millions of events per simulated second. ODE gives the
   correct mean; quantisation noise on round-trip is bounded ±1 per
   species per step. Design doc §3.1 has the full argument.

2. **"Why integer counts at all, not concentrations?"** — Inherited
   from wcEcoli. Important for the upstream cascade in v0.3 where
   PLC-Gq is ~1 molecule; you need integer-aware partitioning to do
   that honestly. Today's calcium core could equally be done in
   concentrations, but the framework is set up to keep them discrete.

3. **"Why one process for all of calcium signalling?"** — The IP3R /
   SERCA / PMCA / SOCE / CaM subsystem is tightly coupled through
   shared cytosolic Ca²⁺. Splitting it across processes would require
   artificial partitioning at every interface. The ODE captures the
   coupling naturally. Mirrors the wcEcoli `TwoComponentSystem`
   design choice.

---

## Live commands during the walkthrough

```bash
# Find the rate-constant tables in the biology file
grep -n "^K_" reconstruction/platelet/dataclasses/process/calcium_signalling.py

# Show the 27-species inventory
sed -n '38,97p' reconstruction/platelet/dataclasses/internal_state.py

# Show the wiring
cat models/platelet/sim/simulation.py

# Run all three biology conditions while talking — each ~30 s
PYTHONPATH=$PWD OPENBLAS_NUM_THREADS=1 \
    pyenv exec python runscripts/manual/runPlateletSim.py walk_t --length 200
PYTHONPATH=$PWD OPENBLAS_NUM_THREADS=1 \
    pyenv exec python runscripts/manual/runPlateletSim.py \
        walk_edta --length 200 --ca-ex-mM 0
PYTHONPATH=$PWD OPENBLAS_NUM_THREADS=1 \
    pyenv exec python runscripts/manual/runPlateletSim.py \
        walk_resting --length 300 --no-ip3-forcing

# Phase 3 driver — runs both +/-Ca conditions, writes comparison plot
PYTHONPATH=$PWD OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg \
    pyenv exec python runscripts/manual/runPhase3.py walk_phase3 --length 200

# Show the deviations table
grep -A 200 "^### 6.8" reports/design/calcium-dynamics-design.md │ head -200

# Show the override sites for the run-time conditions
grep -n 'CA_EX_UM\|_ip3_forced' \
    runscripts/manual/runPlateletSim.py \
    reconstruction/platelet/dataclasses/process/calcium_signalling.py \
    models/platelet/processes/calcium_dynamics.py
```
