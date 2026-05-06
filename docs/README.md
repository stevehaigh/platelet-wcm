# Platelet WCM — Developer Docs

## setup

1. [Required development tools](dev-tools.md) — pyenv, gcc, make, git
2. [Creating the pyenv runtime environment](create-pyenv.md) — Python 3.11.5 + packages

## running

```bash
# Three biologically-distinct run conditions (use --help for the full flag list):

# 1. IP3 transient with extracellular Ca²⁺ (Phase 1, default)
PYTHONPATH="$PWD" OPENBLAS_NUM_THREADS=1 python runscripts/manual/runPlateletSim.py \
    out/transient --length 200

# 2. EDTA transient (Phase 3 no-extracellular-Ca²⁺)
PYTHONPATH="$PWD" OPENBLAS_NUM_THREADS=1 python runscripts/manual/runPlateletSim.py \
    out/edta --length 200 --ca-ex-mM 0

# 3. Resting (no IP3 stimulus)
PYTHONPATH="$PWD" OPENBLAS_NUM_THREADS=1 python runscripts/manual/runPlateletSim.py \
    out/resting --length 300 --no-ip3-forcing

# Run analysis plots on any output
OPENBLAS_NUM_THREADS=1 PYTHONPATH="$PWD" python runscripts/manual/analysisPlatelet.py out/transient

# Phase 3 driver — runs (1) + (2) and produces the Dolan 2014 Fig. 4 comparison
PYTHONPATH="$PWD" OPENBLAS_NUM_THREADS=1 python runscripts/manual/runPhase3.py \
    out/phase3 --length 200
```

The two condition flags (`--ca-ex-mM`, `--no-ip3-forcing`) override
`cs_mod.CA_EX_UM` and `CalciumDynamics._ip3_forced` respectively before
the simulation is constructed; the same conditions are exposed as form
fields and three presets in the webapp Configure tab. See the top-level
[README](../README.md) and each script's `--help` for all options.

## development

* [Coding style guide](style-guide.md)

## tests

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH="$PWD" python -m pytest models/platelet/tests/ -v
```

CI runs automatically on push/PR via `.github/workflows/ci.yml`.
