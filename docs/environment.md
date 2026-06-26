# Python environment setup (uv)

The project uses **[uv](https://docs.astral.sh/uv/)** to manage the Python
interpreter and dependencies. Python is pinned to **3.11.5** via the repo's
`.python-version` file, which uv reads automatically — so you don't pick a version
by hand.

> Migrated from pyenv (2026-06-25). pyenv is no longer required; because
> `.python-version` is shared by both tools, an existing pyenv setup still works
> if you have one. The old pyenv/PyCharm instructions are kept for reference in
> [`create-pyenv.md`](create-pyenv.md) and [`dev-tools.md`](dev-tools.md).

## Install uv

macOS (Homebrew): `brew install uv` — or the standalone installer:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

## One-time project setup

From the repo root:

```bash
uv python install 3.11.5      # installs the pinned interpreter (no-op if present)
uv venv                       # creates .venv using 3.11.5 (from .python-version)
uv pip install -r requirements.txt
uv pip install -r requirements-viz.txt   # optional: replay TUI / experiment bench
```

`uv venv` creates a `.venv/` in the repo (gitignored). uv resolves the
interpreter from `.python-version`, so the venv is always 3.11.5.

## Running commands

Two equivalent styles — pick one:

```bash
# A) prefix with `uv run` (no activation needed)
PYTHONPATH=$PWD OPENBLAS_NUM_THREADS=1 uv run python runscripts/manual/runPlateletSim.py out/my_run --length 200

# B) activate the venv once, then use plain `python`
source .venv/bin/activate
PYTHONPATH=$PWD OPENBLAS_NUM_THREADS=1 python runscripts/manual/runPlateletSim.py out/my_run --length 200
```

Always run from the repo root with `PYTHONPATH=$PWD`, and set
`OPENBLAS_NUM_THREADS=1` for reproducible numerics. `make run` and `make tui`
already invoke `uv run` for you.

## Tests & type-checking

```bash
PYTHONPATH=$PWD OPENBLAS_NUM_THREADS=1 uv run python -m pytest models/platelet/tests/ -v
uv run python -m mypy models/platelet/ reconstruction/platelet/ \
    runscripts/manual/runPlateletSim.py runscripts/manual/analysisPlatelet.py
```

## Notes

- **No compile step** — this is pure Python/NumPy/SciPy; the wcEcoli `.pyx`
  modules were removed.
- **CI** does not use uv yet — it installs with `pip` from
  `requirements-platelet-ci.txt` (see `.github/workflows/ci.yml`). Switching CI to
  `astral-sh/setup-uv` is an optional follow-up.
- **Dependencies remain in `requirements*.txt`** (not `pyproject.toml`). Moving
  them into `pyproject.toml` + a committed `uv.lock` (full uv-native, `uv sync`)
  is a deliberate future step, not done here.
- For IDE config, point the interpreter at `.venv/bin/python`.
