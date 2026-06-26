# Development workflow — build, branch, test

> How code gets written, checked, and merged here. For the testing *philosophy*
> and the honest discussion of the "golden / 5/5" regression approach, see its
> companion [`validation-and-regressions.md`](validation-and-regressions.md).

## Environment

- **Python 3.11.5**, pinned via `.python-version`. Managed with
  [**uv**](https://docs.astral.sh/uv/) — it reads `.python-version`, so the right
  interpreter is selected automatically (no pyenv shims / precmd hook).
- Always run from the repo root with **`PYTHONPATH=$PWD`**.
- Set **`OPENBLAS_NUM_THREADS=1`** for reproducible numerics (threading changes
  floating-point reduction order).
- **No compile step** — the wcEcoli `.pyx` Cython modules were removed; this is
  pure Python/NumPy/SciPy.
- The replay-TUI extras are optional (`requirements-viz.txt`). CI uses the
  slimmer `requirements-platelet-ci.txt`.

```bash
# one-time setup (uv installs 3.11.5 if needed and creates .venv)
uv python install 3.11.5
uv venv
uv pip install -r requirements.txt          # + requirements-viz.txt for the TUI

# run things through the venv — no activation needed
PYTHONPATH=$PWD uv run python runscripts/manual/runPlateletSim.py out/my_run --length 200
PYTHONPATH=$PWD uv run python runscripts/manual/analysisPlatelet.py out/my_run
```

> `uv run <cmd>` runs `<cmd>` in the project `.venv`. If you prefer, activate it
> once with `source .venv/bin/activate` and drop the `uv run` prefix. **pyenv is
> no longer required**; `.python-version` is shared by both, so an existing pyenv
> setup keeps working if you have one.

## Branch & PR flow

- **`main`** is the integration branch and the CI target. Work happens on feature
  branches; **don't commit straight to `main`**.
- **`webapp`** is a deploy branch — pushing to it triggers the Azure Container
  Instances deploy (`.github/workflows/deploy-azure.yml`). Don't use it for
  normal development.
- Open a **PR into `main`**; CI must be green before merge. Recent history is
  small, focused PRs (e.g. "#88 Finish the network rename", "#90 Gitignore
  rendered experiment PDFs"), often paired with a lab-book entry or design doc.
- Commits/PRs are made **only when the user asks**.
- `gh pr edit` / `gh issue close` can fail here on a projects-classic GraphQL
  error; use the REST API (`gh api -X PATCH/POST …`) as a workaround.

### Issue-driven, documented work

Substantive changes are tracked by GitHub issues and accompanied by a
**design doc** (`reports/design/*.qmd`) before, and a **lab-book entry**
(`reports/lab-books/lab-book-YYYY-MM-DD-*.md`) after. The most recent lab book is
the source of truth for "where the work is now." Design docs carry an as-built
status note when the implementation diverges from the plan.

## CI (`.github/workflows/ci.yml`)

Runs on every push and PR to `main`. Three jobs:

1. **pytest** — `models/platelet/tests/` then `wholecell/tests/` (Python 3.11.5,
   `OPENBLAS_NUM_THREADS=1`, `MPLBACKEND=Agg`).
2. **kinetics-review PDF** — renders `reports/design/kinetics-v0.6-review.pdf`
   from the calcium TOML (Quarto + TinyTeX) and uploads it as an artifact. Has
   documented retries for transient CTAN/TinyTeX flakes.
3. **mypy** — type-checks platelet paths only (`models/platelet/`,
   `reconstruction/platelet/`, and the two main runscripts).

## Running the tests

```bash
# Everything (~24 s)
PYTHONPATH=$PWD uv run python -m pytest models/platelet/tests/ -v

# Fast iteration — skip the sim-running tests (~3 s)
PYTHONPATH=$PWD uv run python -m pytest models/platelet/tests/ -m "not slow"

# One file / method
PYTHONPATH=$PWD uv run python -m pytest models/platelet/tests/sim/test_simulation.py

# Type checking
uv run python -m mypy models/platelet/ reconstruction/platelet/ \
    runscripts/manual/runPlateletSim.py runscripts/manual/analysisPlatelet.py
```

Tests requiring a full simulation are marked `@pytest.mark.slow`. Use
`-m "not slow"` while iterating; run the full suite before pushing.

## The test layers (what each kind protects)

The suite (`models/platelet/tests/`) is layered by intent:

| Layer | Examples | Protects against |
|-------|----------|------------------|
| **Unit** | `test_calcium_signalling.py`, `test_resting_decay.py` | a rate law / pure helper computing the wrong thing |
| **Integration** | `test_simulation.py`, `test_dose_sweep.py`, `test_perturbation.py`, `test_knockouts.py`, `test_count_overrides.py` | a sim that won't run, or a knob with no effect / the wrong effect |
| **Subsystem validation** | `test_validation_targets.py`, `test_thromboxane.py`, `test_integrin.py`, `test_inhibitory_axis.py`, `test_second_wave.py`, `test_secretion.py` | a module losing its *expected biological direction/magnitude* |
| **Acceptance — behavioural** | `test_acceptance.py` (`@slow`) | a headline biological result drifting out of its band — resting equilibrium, the Dolan ±Ca²⁺ transient, the MCU/P2Y12/COX-1 knockouts, the resting-quiescence invariant |
| **Regression — structural band** | `test_regression.py` (`@slow`) | metrics drifting out of a physiological tolerance band (peak Ca²⁺ within ±30% of baseline; resting IP₃ 40–60 nM; dry mass within 1%; DTS never below cytosol; SOCE flux ≥ 0) |
| **Analysis** | `analysis/test_analysis.py`, `analysis/test_phase3.py`, `test_figures.py` | plot code crashing on real output |

### The acceptance + regression mechanisms

- **`test_acceptance.py`** is the headline behavioural contract: one readable,
  biologically-anchored band per result (resting equilibrium; the Dolan ±Ca²⁺
  transient and SOCE differential; MCU / P2Y12 / COX-1 knockouts; the
  resting-quiescence invariant that every gated output is exactly 0 at rest). A
  failure *names the biology that moved*, not just a column. It replaced the
  former byte-identical goldens (`test_byte_identical.py` + `golden/*.npz`) and
  the all-or-nothing "Dolan 5/5" gate — see
  [`validation-and-regressions.md`](validation-and-regressions.md).

- **`test_regression.py`** asserts lower-level *tolerance bands* and structural
  invariants (initial Ca²⁺ ~100 nM, peak 200–800 nM, DTS never drains below
  cytosol, SOCE flux ≥ 0, dry mass within 1% of a stated baseline). Each band has
  a comment recording the baseline value, when it last changed, and **why** — so a
  legitimate biology change updates the band with a paper-trail, and an accidental
  one fails loudly.

## Changing parameters vs changing behaviour (and the regression contract)

A recurring distinction in this repo — almost every commit message states which
case it is:

- **Behaviour-preserving change** (refactor, parameter externalisation, adding a
  *terminal output* that doesn't touch the Ca²⁺ ODE): the acceptance and
  regression bands stay green *unchanged* — the result does not move. Examples:
  granule secretion Slice 1, thromboxane Slice A, integrin αIIbβ3.
- **Intentional biology change** (adding ODE states, changing a rate law): a band
  may legitimately shift. Re-run `test_acceptance.py` / `test_regression.py`, and
  where a result has moved for a known reason, **update the band and its comment**
  with the new value and the reason — so the change carries a paper-trail and an
  *un*intended move still fails loudly.

Treat a band change as a reviewable decision, not a chore: the band comment is the
record of why the number is what it is.

## Reproducibility caveats

- **Runs are not currently seed-reproducible**: `RestingDecay` draws from numpy's
  *global* RNG, not the sim seed. The deterministic golden scenarios avoid this
  by exercising paths that don't depend on that draw.
- Use `RunConfig` (or the `runFromConfig.py` JSON spec) to capture a run's exact
  conditions — that spec is the reproducible/shareable unit, not a pile of CLI
  flags.

## Docs & figures pipeline

- `make pdfs` → PDFs of `reports/*.md` (pandoc + xelatex).
- `make quarto-pdfs` → PDFs of `reports/*.qmd` (Quarto + xelatex).
- `make kinetics-review` → the clickable kinetics-review PDF from the TOML.
- `runscripts/manual/buildDocsSite.py` → the HTML docs site under `reports/site/`.
- **Prefer `.qmd` (Quarto)** for new design docs/reports (diagram-heavy, live
  preview); `.md` (pandoc) for prose-only. See the project conventions in
  `CLAUDE.md`.
- Figures: matplotlib **mathtext** for chemistry, **detailed standalone
  captions** (figures are thesis artefacts).

## A sensible loop for a change

1. State the goal as a verifiable criterion (a test, a band, a figure).
2. Branch off `main`.
3. Make the surgical change; match TAB indentation and existing style.
4. Run `pytest -m "not slow"` + `mypy` while iterating; full suite before pushing.
5. Decide the regression contract: behaviour-preserving (bands stay green
   unchanged) or intentional biology change (update the moved band + its comment).
6. Write/refresh the lab-book entry; update `CLAUDE.md`/README if a flag, path, or
   symbol changed (grep the whole repo for the old name first).
7. Open a PR into `main`; merge on green CI.
