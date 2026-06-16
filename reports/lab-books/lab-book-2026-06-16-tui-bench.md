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

## Next

- **P3 — `rate_overrides`** (Tier 1): edit any of the ~200 TOML rate constants. The most
  invasive piece — needs the build-time refactor that moves `_KINETICS` off import-time
  and threads a per-sim kinetics bundle into `_ode_rhs`. Deferred.
- Otherwise: open a PR for the P0–P2 work.
