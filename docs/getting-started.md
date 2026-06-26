# Getting started

A hands-on tour of the repository for someone seeing it for the first time. It
assumes you can read code and follow a bit of cell biology, but **not** that
you're an expert in either — terms are explained as they come up. By the end you
will have run a simulation, found your way around the code, and run a real
*in-silico* knockout experiment (removing a mitochondrial calcium channel) and
read the result.

If you just want the reference map of where everything lives, skip to
[Codebase overview](codebase-overview.md). This page is the gentler, narrative
front door.

---

## 1. What is this, in 30 seconds

This is a **whole-cell-style computational model of a human platelet** — the
small cell fragment in your blood that triggers clotting. The model currently
focuses on **calcium signalling**: when a platelet meets a chemical "activate!"
signal, the concentration of calcium ions (Ca²⁺) inside it spikes, and that
spike is the master switch that drives everything downstream (shape change,
granule release, clumping together). We simulate one platelet over time and
track how its molecules change.

A few terms you'll meet repeatedly:

| Term | Plain meaning |
|------|---------------|
| **Agonist** | A signalling molecule that *activates* the platelet (e.g. thrombin, ADP). |
| **Calcium store / DTS** | An internal compartment (the *dense tubular system*, the platelet's version of the ER) that holds Ca²⁺ and releases it on demand. |
| **Cytosol** | The main body of the cell. "Cytosolic Ca²⁺" is the headline number we track. |
| **Knockout** | Removing or disabling one protein, *in silico* here, to see what breaks — the classic way to test what a component does. |
| **Validation** | Checking the model reproduces a real wet-lab measurement (ours: Dolan & Diamond 2014). |

The simulation engine is reused from [CovertLab's *wcEcoli*](https://github.com/CovertLab/wcEcoli)
(a whole-cell model of *E. coli*); all the bacterial biology was removed and
replaced with platelet biology. For the *why* behind the project, read the
[Domain overview](domain-overview.md); for the design, the
[Architecture](architecture.md).

---

## 2. Set up (about five minutes)

The project uses **[uv](https://docs.astral.sh/uv/)** to manage Python (pinned
to 3.11.5) and dependencies. From the repo root:

```bash
uv sync --extra viz      # core deps + the terminal UI extras
```

That creates a `.venv/` with everything you need. Full detail (and the
`--all-extras` / `dev` variants) is in [environment.md](environment.md).

Two conventions used in every command below:

- **`uv run python …`** runs Python inside the managed venv without you having to
  activate it. (You *can* `source .venv/bin/activate` once and then drop the
  `uv run` prefix.)
- **`PYTHONPATH=$PWD`** — always run from the repo root with this set, so imports
  like `from wholecell.sim.simulation import Simulation` resolve.
- **`OPENBLAS_NUM_THREADS=1`** keeps the numerics reproducible.

---

## 3. Your first run (a two-minute win)

Run a baseline simulation — a platelet hit with a standard agonist transient,
for 200 simulated seconds:

```bash
PYTHONPATH=$PWD OPENBLAS_NUM_THREADS=1 \
  uv run python runscripts/manual/runPlateletSim.py out/first-run --length 200
```

This writes recorded data ("listeners") under
`out/first-run/…/simOut/`. (The deep `platelet_stub_…/generation_…/` nesting is
inherited from wcEcoli — you rarely touch it directly.)

Now turn that into plots:

```bash
PYTHONPATH=$PWD OPENBLAS_NUM_THREADS=1 \
  uv run python runscripts/manual/analysisPlatelet.py out/first-run
```

The headline figure is `calcium_trace` — the cytosolic Ca²⁺ transient. To pull a
number out yourself, every listener column is readable with `TableReader`; the
helper below finds the (nested) `simOut` for you:

```python
import os
from wholecell.io.run_results import find_variants, find_cells
from wholecell.io.tablereader import TableReader

run = "out/first-run"
variant = find_variants(run)[0]
simout = find_cells(run, variant)[0]["simout_path"]
ca = TableReader(os.path.join(simout, "CalciumTrace")).readColumn("ca_cyt_nM")
print(f"resting {ca[0]:.0f} nM  →  peak {ca.max():.0f} nM")
```

`CalciumTrace` has ~25 columns (Ca²⁺ pools, IP₃, the receptor states, ATP costs,
the feedback fractions). `ca_cyt_nM`, `ca_dts_uM`, and `ip3_nM` are the ones
you'll use most.

---

## 4. A tour of the code

The repository has **four layers**. The key idea (see
[Architecture](architecture.md)) is *reuse the generic engine, replace the
biology*:

```
wholecell/                  The generic simulation engine (from wcEcoli).
                            Timestep loop, state containers, binary I/O, units,
                            and the terminal UI (wholecell/tui/). Biology-agnostic.

models/platelet/            THE MODEL. The platelet-specific biology:
  processes/                  submodels that change the cell's state each step
  listeners/                  recorders that observe state and write it to disk
  sim/                        wires the processes + listeners into a simulation
  analysis/single/            the plots (calcium_trace, the demo figures, …)
  tests/                      unit + behavioural tests

reconstruction/platelet/    THE PARAMETERS. The fixed "knowledge base":
  simulation_data.py          the parameter object handed to the engine
  run_config.py               RunConfig — how you configure ONE run (see §5)
  dataclasses/process/        rate constants + the calcium ODE itself

runscripts/manual/          ENTRY POINTS. The scripts you actually run
                            (runPlateletSim.py, analysisPlatelet.py, the
                            experiment drivers).
```

For the full map, see [Codebase overview](codebase-overview.md).

**The simulation loop**, each 1-second timestep, in one breath: the *processes*
declare which molecules they need → the *states* hand out those molecules →
the processes run their biology and update the state → the *listeners* record
everything. Repeat.

**Where the biology actually lives.** Two files do most of the scientific work:

- `reconstruction/platelet/dataclasses/process/calcium_signalling.py` — the
  **calcium ODE** (the system of differential equations for IP3R / SERCA / PMCA /
  SOCE / the GPCR receptor cascade) plus its rate constants. This is the heart of
  the model. The process `models/platelet/processes/calcium_dynamics.py` is just
  a thin wrapper that calls the solver each timestep.
- `reports/params/calcium-v0.6.toml` — those rate constants externalised as
  data, and `reports/params/species-v0.6.tsv` — the inventory of every molecule
  (id, mass, starting count).

So "the model" is: a list of molecules (TSV) + a set of rate constants (TOML) +
the equations that use them (`calcium_signalling.py`), driven by the engine.

---

## 5. Worked example: knocking out the MCU

This is the part worth doing slowly, because it teaches the general pattern for
*every* experiment in the repo.

### The biology

Mitochondria can soak up calcium through a channel called the **MCU**
(mitochondrial calcium uniporter). The interesting question: if you remove the
MCU, what happens to the cytosolic Ca²⁺ spike?

There are two competing intuitions:

- MCU is a **buffer** — it pulls Ca²⁺ *out* of the cytosol. Remove it and you'd
  expect cytosolic Ca²⁺ to go **up**.
- But MCU also sits at mitochondria–store contact sites and **helps the store
  release** Ca²⁺ (it relieves a brake on the IP3R release channel). Remove it and
  release weakens, so cytosolic Ca²⁺ would go **down**.

Two real platelet studies (Ghatge 2026; Ajanel 2025) knocked out MCU and measured
*reduced* cytosolic Ca²⁺ — the second effect wins. Our model includes that
release coupling, so it should reproduce the **reduction**. Let's check.

### Run it (one command)

```bash
PYTHONPATH=$PWD OPENBLAS_NUM_THREADS=1 \
  uv run python runscripts/manual/runPerturbation.py getting-started-mcu --experiment mcu
```

You'll see (numbers are deterministic, seed 0):

```
[mcu] MCU loss reduces the cytosolic transient (v0.7 #76 coupling, +Ca)
  MCU V_max ×0     → peak cyt  435.9 nM, recovery-AUC     93649 nM·s, DTS min   20.1 µM
  MCU V_max ×1     → peak cyt  530.6 nM, recovery-AUC     64746 nM·s, DTS min    0.8 µM
```

Read that off: the **knockout (`×0`) peaks at 436 nM vs the baseline (`×1`) 531
nM — an ~18% reduction**, matching the wet-lab direction. Notice too that the
knockout's store stays fuller (`DTS min` 20.1 vs 0.8 µM) — the loss of release is
exactly why the cytosol sees less. The run also writes a two-panel figure
`out/getting-started-mcu/mcu_traces.png` (knockout in red, baseline in black) and
a `mcu.npz` of the raw traces.

### How it's configured — the transferable bit

Open `runscripts/manual/runPerturbation.py` and find the `EXPERIMENTS` registry.
The `mcu` entry says, in effect: *take the rate constant `K_MITO['V_max_MCU']`
and multiply it by each factor in `[0.0, 1.0]`*. The `0.0` is the knockout.

Under the hood, that factor is just a field on a **`RunConfig`** — the single
object that describes one run. Everything you can vary about an experiment is a
field on it (see `reconstruction/platelet/run_config.py` for the fully-documented
list). The MCU knockout is one line:

```python
from reconstruction.platelet.run_config import RunConfig
from runscripts.manual.runPlateletSim import run_platelet_sim, resolve_sim_path

ko = RunConfig(ca_ex_mM=1.2, mcu_vmax_scale=0.0)   # 0.0 = knockout, 1.0 = normal
run_platelet_sim(resolve_sim_path("out/my-mcu-ko"),
                 length_sec=400, seed=0, run_config=ko)
```

`RunConfig` is **frozen** (immutable) and gets recorded into the run's
`metadata/` for provenance — so a run is fully reproducible from its config.
The same `*_scale` / knob pattern gives you the other experiments:

| Field | What it does |
|-------|--------------|
| `mcu_vmax_scale=0.0` | MCU knockout (this example) |
| `cox1_factor=0.0` | aspirin (abolishes thromboxane) |
| `p2y12_block=1.0` | clopidogrel / ticagrelor (P2Y12 antagonist) |
| `integrin_act_scale=0.0` | Glanzmann thrombasthenia (no integrin activation) |
| `ca_ex_mM=0.0` | the EDTA / no-external-Ca²⁺ validation condition |

> **Try this:** the MCU result *depends on a modelling choice*. Add
> `mito_coupling_gain=0.0` to the `RunConfig` and the release coupling switches
> off — now you get the old "buffer-only" behaviour where the knockout *raises*
> cytosolic Ca²⁺ instead. Flipping that one knob is a good way to see how the
> model encodes a biological hypothesis.

### Read the result yourself

The perturbation run saved the traces into an `.npz`, so you don't need to dig
through `simOut/`:

```python
import numpy as np
d = np.load("out/getting-started-mcu/mcu.npz")
for factor, trace in zip(d["factors"], d["cyt"]):   # cyt[i] is the trace for factors[i]
    print(f"mcu_vmax ×{factor:g}: peak {trace.max():.1f} nM")
```

---

## 6. More experiments you can run

- **`runPerturbation.py --experiment pmca|pkc|plcb`** — the other single-knob
  scans (calcium pump, and two PKC feedback brakes).
- **`runDemoExperiments.py`** — runs the four "thesis demo" experiments
  (baseline, aspirin, MCU knockout, clopidogrel) in one pass with a uniform
  300 s-settle → 300 s-stimulate protocol, and renders their figures. Add
  `--tag 2026-06-26` to archive a dated set without overwriting the committed
  figures. The write-ups are in [`reports/experiments/`](../reports/experiments/).
- **`make tui`** — the **interactive** way: a terminal app where you tick
  knockouts, set conditions, hit run, and watch the Ca²⁺ trace live. Good for
  poking around without writing any code.

A note you'll meet often: under a strong saturating agonist the Ca²⁺ response is
*store-limited*, so many knockouts barely move the peak. The clearest readout is
usually a more specific column (a receptor-desensitised fraction, IP₃, the
integrin PAC-1 fraction) — that's why each experiment harvests its own.

---

## 7. Where to go next

- [Domain overview](domain-overview.md) — the biology and the project's purpose.
- [Architecture](architecture.md) — the four layers, the timestep, the pure-ODE
  pattern, and why `RunConfig` replaced global state.
- [Codebase overview](codebase-overview.md) — the complete directory map.
- [Development workflow](development-workflow.md) — branching, PRs, the test
  layers, and CI.
- [Validation and regressions](validation-and-regressions.md) — an honest
  account of what the model is and isn't validated against.
- [`reports/lab-books/`](../reports/lab-books/) — dated session notes; the most
  recent is the current state of the work.

The root [`CLAUDE.md`](../CLAUDE.md) is the densest single reference for the
whole repo, if you want the firehose.
