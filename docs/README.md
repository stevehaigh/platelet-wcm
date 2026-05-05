# Platelet WCM — Developer Docs

## setup

1. [Required development tools](dev-tools.md) — pyenv, gcc, make, git
2. [Creating the pyenv runtime environment](create-pyenv.md) — Python 3.11.5 + packages

## running

```bash
# Run a simulation (60 s, seed 0, output to out/platelet_manual/)
OPENBLAS_NUM_THREADS=1 PYTHONPATH="$PWD" python runscripts/manual/runPlateletSim.py

# Run analysis plots on the output
OPENBLAS_NUM_THREADS=1 PYTHONPATH="$PWD" python runscripts/manual/analysisPlatelet.py out/platelet_manual
```

See the top-level [README](../README.md) and each script's `--help` for all options.

## development

* [Coding style guide](style-guide.md)

## tests

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH="$PWD" python -m pytest models/platelet/tests/ -v
```

CI runs automatically on push/PR via `.github/workflows/ci.yml`.
