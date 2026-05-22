# Platelet WCM — Developer Docs

## setup

1. [Required development tools](dev-tools.md) — pyenv, gcc, make, git
2. [Creating the pyenv runtime environment](create-pyenv.md) — Python 3.11.5 + packages

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
form fields and four presets in the webapp Configure tab. See the
top-level [README](../README.md) and each script's `--help` for all
options.

## development

* [Coding style guide](style-guide.md)

## tests

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH="$PWD" python -m pytest models/platelet/tests/ -v
```

CI runs automatically on push/PR via `.github/workflows/ci.yml`.
