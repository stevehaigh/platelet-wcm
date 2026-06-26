# Architecture and rationale

> What the simulation framework is, how a platelet timestep runs, and *why* it is
> shaped this way. For the file/module layout see
> [`codebase-overview.md`](codebase-overview.md). For the original, detailed
> wcEcoli engine map (with reusability assessment per component) see
> [`../platelet-plan/architecture.md`](../platelet-plan/architecture.md).

## The core decision: reuse the engine, replace the biology

The project forked [CovertLab/wcEcoli](https://github.com/CovertLab/wcEcoli) and
asked one question: *which parts of a whole-cell model are about E. coli, and
which are about "being a whole-cell simulator"?* The answer drove the whole
design:

- **Kept (generic):** the simulation loop, the Process/State/Listener
  abstractions, the View system (how processes request and receive molecules),
  resource partitioning, the columnar binary I/O (TableWriter/TableReader), the
  units system (Unum), and the analysis-plot plumbing.
- **Removed (E. coli-specific):** every E. coli process (transcription,
  translation, replication, metabolism, complexation…), every E. coli listener,
  the entire `reconstruction/ecoli` ParCa pipeline, and **cell division**
  (platelets don't divide — the engine already supports never completing a cell
  cycle).
- **Removed (performance):** the Cython `.pyx` modules. The platelet model is
  small enough to run in pure Python/NumPy/SciPy, so there is **no compile step**.

The payoff: a tested, modular engine with a clean extension point, and a model
codebase small enough to reason about end-to-end.

## The four layers

```
  reconstruction/platelet/         (1) Parameters / knowledge base
     SimulationDataPlatelet  ──────────►  built directly in code, NO ParCa step
            │ sim_data
            ▼
  wholecell/sim/simulation.py      (2) Generic engine (inherited, unchanged)
     PlateletSimulation(Simulation)
            │ wires
            ▼
  models/platelet/                 (3) The platelet model
     processes/  states  listeners      biology written as Processes + Listeners
            │ writes columns
            ▼
  models/platelet/analysis/        (4) Analysis — reads simOut, makes figures
```

1. **Reconstruction** — `SimulationDataPlatelet` is the knowledge base the engine
   reads (`submass_name_to_index`, `compartment_abbrev_to_index`, `constants`,
   `molecule_groups`, `internal_state`, `external_state`, `process`). In wcEcoli
   this is produced by a heavyweight "ParCa" (parameter calculation) step from
   raw data. **We construct it directly in Python** from the TSV/TOML parameter
   files — there is no `raw_data/` directory and no ParCa. This is the single
   biggest simplification versus wcEcoli, and it is appropriate because the
   platelet parameter set is small and hand-curated.

2. **Simulation engine** — `wholecell/sim/simulation.py` is used **as-is**.
   `PlateletSimulation` just declares its state classes, process groups,
   listeners, and an initial-conditions function. The engine is fully generic;
   nothing platelet-specific lives in `wholecell/` except the TUI.

3. **The model** — `models/platelet/` holds the biology as **Processes** (modify
   state) and **Listeners** (observe & record). This is where the science lives.

4. **Analysis** — post-hoc plots read the columnar output. Analysis never runs
   inside the sim.

## The simulation loop (1-second discrete timesteps)

Each timestep, in order:

1. **`calculateRequest()`** — every process declares which molecules (and how
   many) it needs.
2. **partition** — states allocate molecules to processes (proportional +
   priority when demand exceeds supply). Only `BulkMolecules` is partitioned.
3. **`evolveState()`** — each process executes its biology for this `dt`, mutating
   its allocated counts.
4. **merge** — states reconcile the per-process changes back into the global
   counts.
5. **Listeners** observe the merged state (read-only) and write a row to disk.

Processes run in **dependency-ordered groups**; within a group they share one
partition step. The platelet model currently has a **single group**:
`(RestingDecay, CalciumDynamics, GranuleSecretion, ThromboxaneSynthesis, IntegrinActivation)`.

### Why a 1-second discrete step over a stiff continuum ODE?

The interesting Ca²⁺ kinetics are sub-second and stiff. The model resolves this
by keeping the **engine** at 1 s (the bookkeeping/partition granularity) while
`CalciumDynamics` integrates its ODE system **internally with SciPy** across each
1 s engine step. So the engine provides structure, resource accounting, and I/O;
the actual Ca²⁺ math is a proper ODE solve. This is the key adaptation that lets
a discrete copy-number engine carry a continuous biochemical model.

## How the biology is organised: the pure ODE + thin process pattern

`CalciumDynamics` (the biologically rich process) is deliberately **thin**. The
science is a pure function:

```
_ode_rhs(y, t, config, step_inputs)   # in calcium_signalling.py
```

- `y` — the ODE state vector (Ca²⁺ pools, IP₃R/SERCA/PMCA sub-states, CaM
  sub-states, IP₃, GPCR cascade, TP/P2Y12 etc.).
- `config` — the per-sim `RunConfig` (see below), read as a plain argument.
- `step_inputs` — autocrine species fed back by name (secreted ADP[e], TXA₂[e]).

Keeping the RHS a **pure function of (state, config, inputs)** is what makes the
model testable (unit-test the rate law directly), reproducible, and free of
hidden global state.

## Run-time configuration: `RunConfig` (and why it replaced globals)

A simulation's conditions — extracellular Ca²⁺, agonist peaks, feedback gains,
aspirin/knockout knobs, count overrides — are a **frozen dataclass**
(`reconstruction/platelet/run_config.py`) attached as `sim.run_config`.
Processes and listeners read it in `initialize()`; the ODE receives it as an
argument.

This **replaced (v0.62) an earlier pattern of mutating module globals in place**
(`cs_mod.CA_EX_UM`, `CalciumDynamics._adp_peak_uM`, `tx_mod.COX1_FACTOR`, …),
which was process-global and only correct by save/restore discipline. The
rationale: a frozen, per-sim config object is thread-safe, makes every run
reproducible from a single serialisable spec, and removes a whole class of
"previous run leaked into this one" bugs. The TUI exploits this — it runs each
sim as a subprocess of `runFromConfig.py` taking a `{length_sec, seed,
run_config}` JSON spec, so a config doubles as a shareable preset.

## Compartments and molecule IDs

Molecules carry compartment tags, e.g. `Ca[c]`, `Ca[dts]`, `STIM1_dimer[pl]`.
Platelet compartments:

| Abbrev | Compartment |
|--------|-------------|
| `c`   | cytoplasm |
| `dts` | dense tubular system (ER-equivalent Ca²⁺ store) |
| `dg`  | dense granule lumen |
| `ag`  | alpha-granule lumen |
| `m`   | mitochondrial matrix |
| `pl`  | plasmalemma (membrane) — note `pl` not `pm`, to avoid clashing with `m` |
| `e`   | extracellular / open canalicular system |

State is held in three containers: **BulkMolecules** (integer copy numbers,
partitioned), **UniqueMolecules** (per-instance objects — present but currently
unused by the platelet model), and **LocalEnvironment** (extracellular
concentrations).

## Integer counts at the edge of the continuum

At ~100 nM cytosolic Ca²⁺ in 6 fL there are only ~360 ions, and the Gq–PLC
complex can be <1 molecule/cell. The engine's integer-count representation sits
right at the edge of continuum validity. The model is **deterministic** (the ODE
runs in concentration space; counts are the bookkeeping currency), which is an
accepted simplification: it cannot reproduce the cell-to-cell variability that
stochastic models (Purvis 2008, Sveshnikova 2015) show. This is a stated
limitation, not an oversight.

## Output directory structure

The variant/seed/generation/cell nesting is preserved from wcEcoli so the
analysis tooling (e.g. `analysisPlatelet.py`) can locate a run's output
unmodified:

```
out/{sim_dir}/
  metadata/                                # git hash, seed, timestamp
  kb/simData.cPickle                       # serialized SimulationDataPlatelet
  platelet_stub_{seed:06d}/{seed:06d}/generation_000000/000000/
    simOut/                                # listener columns (TableWriter format)
    plotOut/                               # analysis figures
```

There is only one "variant" (`platelet_stub`); there is no platelet equivalent of
E. coli variants yet.

## Design principles, in short

- **Engine generic, biology specific.** Nothing platelet-specific in `wholecell/`
  (bar the TUI).
- **Parameters are versioned data with citations**, not magic numbers in code.
- **The rate law is a pure function** — config in, derivatives out.
- **Each new layer is normalised at rest** so the resting fixed point and the
  Dolan transient are preserved (this is also why "5/5 passes by construction" —
  see [`validation-and-regressions.md`](validation-and-regressions.md)).
- **Reproducibility via `RunConfig`**, not global mutation.
