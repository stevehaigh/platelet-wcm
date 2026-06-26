# Platelet WCM — Developer Docs

## orientation (start here)

A set of orientation docs aimed at AI assistants and new readers — what the
project is, how it's built, and how we know it works:

1. [Domain overview](domain-overview.md) — the biological problem and how we tackle it
2. [Architecture and rationale](architecture.md) — the engine, the timestep, the *why*
3. [Codebase overview](codebase-overview.md) — where everything lives
4. [Development workflow](development-workflow.md) — branch, PR, the test layers, and the behavioural regression suite
5. [Validation and regressions](validation-and-regressions.md) — an honest account of "Dolan 5/5", the retired byte-identical goldens, and the behavioural regression suite that replaced them
6. [Source papers](papers/) — one summary per primary source (+ index)

Deeper references: [`../platelet-plan/`](../platelet-plan/) (high-level plan +
wcEcoli engine map), [`../reports/data/calcium-data-provenance.md`](../reports/data/calcium-data-provenance.md)
(value-by-value provenance), [`../reports/design/`](../reports/design/) (design
docs), [`../reports/lab-books/`](../reports/lab-books/) (the latest = current state).

## setup

1. [Python environment setup](environment.md) — uv + Python 3.11.5 + packages
2. [Recommended dev tools](dev-tools.md) — git, make, etc. (also has the legacy pyenv setup)

## running

```bash
# Three biologically-distinct run conditions (use --help for the full flag list):

# 1. Agonist transient with extracellular Ca²⁺ (Phase 1, default)
PYTHONPATH="$PWD" OPENBLAS_NUM_THREADS=1 python runscripts/manual/runPlateletSim.py \
    out/transient --length 200

# 2. EDTA transient (Phase 3 no-extracellular-Ca²⁺)
PYTHONPATH="$PWD" OPENBLAS_NUM_THREADS=1 python runscripts/manual/runPlateletSim.py \
    out/edta --length 200 --ca-ex-mM 0

# 3. Resting (zero all agonist peaks — no extracellular stimulus)
PYTHONPATH="$PWD" OPENBLAS_NUM_THREADS=1 python runscripts/manual/runPlateletSim.py \
    out/resting --length 300 --at-rest

# Run analysis plots on any output
OPENBLAS_NUM_THREADS=1 PYTHONPATH="$PWD" python runscripts/manual/analysisPlatelet.py out/transient

# Phase 3 driver — runs (1) + (2) and produces the Dolan 2014 Fig. 4 comparison
PYTHONPATH="$PWD" OPENBLAS_NUM_THREADS=1 python runscripts/manual/runPhase3.py \
    out/phase3 --length 200
```

The condition flags override module-/class-level state before the
simulation is constructed: `--ca-ex-mM` overrides `cs_mod.CA_EX_UM`, and
`--at-rest` (or the per-receptor `--thrombin-peak-nM` / `--adp-peak-uM`
/ `--atp-ex-peak-uM` flags) overrides `CalciumDynamics._thrombin_peak_nM`,
`_adp_peak_uM`, and `_atp_ex_peak_uM`. The same knobs are exposed as
editable fields (with built-in and user-saved presets) in the TUI
experiment bench (`make tui`). See the
top-level [README](../README.md) and each script's `--help` for all
options.

## development

* [Coding style guide](style-guide.md)

## tests

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH="$PWD" python -m pytest models/platelet/tests/ -v
```

CI runs automatically on push/PR via `.github/workflows/ci.yml`.
