---
title: "Lab book — 2026-06-16: TUI experiment bench (P0–P2) + count-override knockouts"
---

# Lab book — 2026-06-16: TUI experiment bench (P0–P2)

## Context

Built a **terminal UI** (Textual) that turns the model into an interactive
*experiment bench*: edit run conditions and pathway knobs, knock out
receptors/pathways, run, and watch the Ca²⁺ trace stream live. Scoped first as a
design doc (`reports/design/tui-tinkering-dashboard-2026-06-15.qmd`), then built
P0 → P2. All work is on the `new_UI` worktree (branch `worktree-new_UI`); not yet
pushed / PR'd.

The design frames every tunable value as **three tiers**: Tier 0 = `RunConfig`
(live and clean since v0.62), Tier 1 = ~200 TOML rate constants, Tier 2 = 78 TSV
initial counts. Locked decisions (doc §12): subprocess per run; `textual-plotext`
line charts; the TUI is intended to eventually replace the Dash webapp.

## What landed

**P0 — skeleton.** `wholecell/tui/app.py` (`PlateletBenchApp`), entry
`runscripts/manual/runTui.py`, `make tui`. Edit the core run conditions → run as a
subprocess → poll `simOut/live.csv` → stream cytosolic + DTS Ca²⁺ as plotext line
charts.

**P1 — the real bench.**
- Full `RunConfig` surfaced as a compact, collapsible grouped form (Stimulus /
  Feedback loops / Pumps & brakes), schema + assembly in `wholecell/tui/runspec.py`.
- New config-file runner `runscripts/manual/runFromConfig.py` (`{length_sec, seed,
  run_config}` JSON → `RunConfig`) — needed because the `runPlateletSim.py` CLI only
  exposes a subset of `RunConfig`; the spec doubles as a preset.
- Inline **KO** checkboxes on the 8 Tier-0 scale/gain knobs (force → 0).
- **Presets** (`wholecell/tui/presets.py`): built-ins (Agonist / 60 s settle / EDTA /
  Resting / Aspirin / Glanzmann) + user save/load.
- **Compare-to-baseline** overlay (`b`), **modified-from-defaults** indicator +
  **stale-preset** marker, and **on-demand 5-panel figure** (`f` → `analysisPlatelet.py`
  → opens the PDF).

**P2 — expression knockouts (the model refactor).**
- `RunConfig.count_overrides: Dict[str,int]` applied at state-seeding
  (`reconstruction/platelet/initialization.py:initialize_bulk_molecules`, wired from
  `models/platelet/sim/initial_conditions.py` via `sim.run_config`). Keeps `RunConfig`
  separate from `sim_data` — they meet at seeding. Empty → byte-identical; unknown id
  raises.
- Logical-entity map `reconstruction/platelet/knockouts.py` (6 receptors/integrin +
  SERCA → all conformational sub-states, so a knockout zeroes the whole entity and
  preserves conservation). Surfaced in the TUI as a "Knockouts (remove protein)"
  checkbox group → `count_overrides`.

## Verification

- **Tests:** 28 in `wholecell/tests/tui/` (incl. a slow end-to-end that drives a real
  5 s sim through the TUI, ticks aspirin + a PAR1 expression-KO, and checks both reach
  the written spec) + 4 each for `count_overrides` and the entity map in
  `models/platelet/tests/sim/`. Full fast suite **155 passed, no regressions**.
- **End-to-end biology:** under a thrombin-only stimulus, a PAR1+PAR4 expression
  knockout collapses the response to its resting fixed point — peak cyt 316 → 107 nM
  (−66 %), peak IP₃ 164 → 50 nM (−69 %). Confirms count→0 genuinely removes function
  (the ODE reads receptor amounts from the state vector, not a TOML constant).

## Gotchas (so the next session doesn't relive them)

1. **Textual `App` attribute collisions.** `App` owns `_running` *and* `_ready` (it does
   `await self._ready()`). Assigning a bool to either silently breaks the app (first as
   "run never starts", then as `'bool' object is not callable` on mount). Namespace every
   custom App instance attribute (we use `_sim_*`, `_ui_ready`, `_applied_*`, …) and
   `hasattr(App(), name)`-check before adding new ones.
2. **Knockouts are invisible under the default saturating agonist.** Peak cytosolic Ca²⁺
   is store-limited and the other agonists compensate, so a receptor KO barely moves it.
   Isolate the agonist (e.g. thrombin-only, ADP/ATP = 0) and read IP₃, or use the
   baseline overlay. Same lesson as the v0.6 PKC brake being invisible on cyt Ca²⁺.
3. **`python3` ≠ the sim interpreter.** On the dev machine `python3` resolves to a system
   Python with `textual` but not `textual-plotext` / the sim deps. `make tui` therefore
   uses `pyenv exec python` (pinned 3.11.5 has everything); `make`'s `/bin/sh` can reach
   pyenv at `/opt/homebrew/bin/pyenv`.

## Outcome

Opened **PR #62** (`worktree-new_UI` → `main`); CI green (mypy + pytest 3.11.5 +
kinetics-review PDF). Requested a GitHub Copilot review — positive summary plus two
trivial nits (a stale `app.py` docstring; an `open().read()` handle leak in a test),
both fixed and pushed.

## P3 (rate overrides) — dropped

Decided **not** to build Tier-1 rate overrides in the TUI. Two reasons:

1. **Not a biological intervention.** The ~200 TOML constants are mostly *elementary*
   rates of multi-state schemes (SERCA's E1↔E2 cycle, the PMCA/CaM states, the SOCE MWC
   parameters, buffer on/off rates). No real perturbation edits one of those in
   isolation; real interventions act on *expression* (copy number — P2) or *whole-protein
   activity* (the Tier-0 scale knobs), both already exposed.
2. **It breaks the calibrated resting state.** The resting fixed point and the
   pre-equilibrated buffer / sub-state counts are pinned to the original rate ratios (see
   `reports/design/calibration-coupling-2026-05-25.qmd`). Editing a rate without
   re-deriving the dependent counts/rates gives a cell that isn't at rest at t=0 — a
   silently non-physiological run that's easy to over-interpret.

**Positioning:** the TUI is a **demo / explanation** tool, not an analysis workbench.
Sensitivity analysis and parameter sweeps stay **scripted** (`runPerturbation.py`,
`runDoseSweep.py` → reproducible figures); future TUI polish should favour
clarity/storytelling over analysis depth. **P0–P2 is the intended endpoint.**
