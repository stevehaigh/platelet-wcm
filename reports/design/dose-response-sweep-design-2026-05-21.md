---
title: "Dose-Response Sweep — Design Plan"
author: "Steve Haigh"
date: "2026-05-21"
---

# Dose-Response Sweep — Design Plan

**Status:** design plan, no code yet. Captures what needs to change in
the model to support a 2-D agonist dose sweep cleanly (no global
mutation, no monkey-patching, no flag double-duty).

**Update 2026-05-22:** the original draft proposed *splitting*
`ip3_forced` into two flags with a legacy shim. Decision since
revised: **expunge** the Dolan Fig. S2 IP3-overlay path entirely —
the v0.4 GPCR cascade produces IP3 endogenously, so the forced
overlay is redundant. This shrinks the refactor and removes a
deprecated stimulus channel. Sections 3, 4, 5, 6, 8 rewritten below;
section 2 framing kept so future readers see why the two channels
were coupled to begin with.

**Audience:** future-me picking this up next week; the dissertation
supervisor wanting to know what blocks the feature.

**Trigger:** platelet-researcher feedback — *"can we do a 2-way dose
response sweep that compares cyt Ca spikes based on different values of
ADP and Thrombin antagonists? Output would be a 3-D plot ADP / thrombin
/ cyt Ca max peak."*

---

## 1. What we want to deliver

A runscript that sweeps `(ADP_peak, thrombin_peak)` on an N × N grid
and reports peak cytosolic Ca²⁺ per cell, plus:

- `sweep.csv` — one row per grid cell (`adp_peak_uM,
  thrombin_peak_nM, peak_ca_cyt_nM`)
- `sweep_surface.png` — matplotlib 3-D surface (headline figure)
- `sweep_heatmap.png` — 2-D heatmap with cell-value annotations
- A `Dose Response` tab in the Dash webapp rendering an interactive
  Plotly surface from the same CSV

The semantic mapping is *"lower agonist peak = stronger competitive
antagonism"* (Option 1 from the scoping discussion — no new antagonist
species, no Kᵢ values needed). Default grid: 5 × 5, log-spaced, 200 s
sims, 25 runs total (~5–10 min).

## 2. Why the current code can't do this cleanly

Two unrelated knobs are bolted onto one flag in
`reconstruction/platelet/dataclasses/process/calcium_signalling.py`,
and the agonist peaks are baked into function default args at import.
Both have to be addressed before the sweep can be expressed without
hacks.

### 2.1 `ip3_forced` does double duty

The same boolean gates two independent stimulus channels in
`_ode_rhs`:

```python
# Channel A — Dolan Fig. S2 IP3 overlay (the original meaning)
if ip3_forced:
    dy[_IDX['IP3[c]']] += ip3_forcing_uM(t)  # forced curve onto IP3[c]

# Channel B — agonist stimulus protocol (added in v0.4 / v0.4.1)
if ip3_forced:
    adp_um_now = adp_uM(t_sim_start + t, delay=ip3_delay)
    thr_nm_now = thrombin_nM(t_sim_start + t, delay=ip3_delay)
else:
    adp_um_now = ADP_REST_UM        # 0
    thr_nm_now = THROMBIN_REST_NM   # 0

# (parallel gate for atp_ex_forcing_uM around line 1238)
```

This was a pragmatic shortcut when the GPCR cascade was added — the
existing flag was the obvious lever — but it permanently couples them.
The three states that exist today are:

| `ip3_forced` | Dolan IP3 overlay | GPCR agonists |
|---|---|---|
| `True`  | applied | applied   |
| `False` | skipped | snap to rest (no stimulation at all) |

The state we *need* for a dose sweep — *"GPCR agonists yes, Dolan
overlay no"* — is unreachable. Turning the flag off zeroes the
agonists; turning it on overlays the Dolan curve, which dominates the
peak and washes out the dose-dependence.

### 2.2 Peak constants are captured at import

```python
ADP_PEAK_UM         = 10.0
THROMBIN_PEAK_NM    = 1.0

def adp_uM(t, delay=0.0, peak_uM=ADP_PEAK_UM):   # captured here
    ...

def thrombin_nM(t, delay=0.0, peak_nM=THROMBIN_PEAK_NM):  # captured
    ...
```

Default-arg expressions evaluate once at function definition, so
`cs_mod.ADP_PEAK_UM = X` from a runscript has no effect on what
`adp_uM(t)` returns — the closure already holds the original value.
The existing `CA_EX_UM` override pattern works because that constant
is read live inside `_ode_rhs`, not captured.

The only way to vary the peak *today* is to monkey-patch the function
itself. That's what the prototype sweep driver did. It works
mechanically but it's an anti-pattern — invisible to readers, fragile
under refactoring, hostile to tests.

## 3. The refactor — three independent fixes

Each fix is small, mechanically obvious, and individually testable.

### 3.1 Live-read peak constants

In `calcium_signalling.py`, change the three peak-bearing functions
from baked defaults to live module reads:

```python
def adp_uM(t, delay=0.0, peak_uM=None):
    if peak_uM is None:
        peak_uM = ADP_PEAK_UM        # resolved at call time
    ...

def thrombin_nM(t, delay=0.0, peak_nM=None):
    if peak_nM is None:
        peak_nM = THROMBIN_PEAK_NM
    ...

def atp_ex_forcing_uM(t, delay=0.0, peak_uM=None):
    if peak_uM is None:
        peak_uM = ATP_EX_PEAK_UM
    ...
```

After this change, `cs_mod.ADP_PEAK_UM = X` actually takes effect at
the next call — matching the existing `cs_mod.CA_EX_UM` override
pattern exactly. No new globals introduced.

**Verification:** unit test that sets `cs_mod.ADP_PEAK_UM = 0.1`, runs
a 10 s sim, and asserts the trajectory differs measurably from the
default-peak run.

### 3.2 Expunge the Dolan Fig. S2 IP3-overlay path

The v0.4 GPCR cascade (P2Y1 / PAR1 / PAR4 → Gαq → PLCβ → IP3) is now
the canonical IP3 source. The forced overlay channel is a v0.3 legacy
that no current scientific question needs. Remove it outright:

- Delete `ip3_forcing_uM(t)` and any helpers it pulls in
- Delete the IP3-overlay `if ip3_forced: dy[IP3] += ...` block in
  `_ode_rhs`
- Rename the remaining gate (which controls ADP / thrombin / ATP\_ex
  application) from `ip3_forced` → `apply_agonist_protocol` to match
  what it actually does

The three call sites collapse to one flag:

| Today | After |
|---|---|
| IP3 forcing overlay block | **deleted** |
| ATP\_ex gate for P2X1 (~line 1238) | `apply_agonist_protocol` |
| ADP / thrombin gate for GPCRs (~line 1283) | `apply_agonist_protocol` |

Single flag propagates up through:

```
CalciumSignalling.molecules_to_next_time_step(
    counts, dt, t_sim, apply_agonist_protocol, ip3_delay,
)
    └── _ode_rhs(t, y, t_sim_start, apply_agonist_protocol, ip3_delay)
```

The dose sweep then runs with `apply_agonist_protocol=True` — the
GPCR cascade alone is the stimulus, exactly as the biology demands.
"At rest" sims use `apply_agonist_protocol=False`.

### 3.3 Thread agonist peaks through the call chain

Mirror the existing `ca_ex_mM` pattern in `run_platelet_sim`:

```python
def run_platelet_sim(
        ...,
        ca_ex_mM=DEFAULT_CA_EX_MM,
        adp_peak_uM=None,
        thrombin_peak_nM=None,
        apply_agonist_protocol=True,
        ip3_delay=DEFAULT_IP3_DELAY,
        ):
    cs_mod.CA_EX_UM = float(ca_ex_mM) * 1000.0
    if adp_peak_uM is not None:
        cs_mod.ADP_PEAK_UM = float(adp_peak_uM)
    if thrombin_peak_nM is not None:
        cs_mod.THROMBIN_PEAK_NM = float(thrombin_peak_nM)
    CalciumDynamics._apply_agonist_protocol = apply_agonist_protocol
    CalciumDynamics._ip3_delay = float(ip3_delay)
    ...
```

The sweep driver then becomes ~3 lines per cell, no monkey-patching:

```python
run_platelet_sim(
    cell_dir,
    adp_peak_uM=adp,
    thrombin_peak_nM=thr,
    apply_agonist_protocol=True,
    length_sec=200,
)
```

## 4. Migration — direct caller updates

No legacy shim. Existing callers update in place:

| Caller | Before | After |
|---|---|---|
| `runPlateletSim.py` CLI | `--no-ip3-forcing` | `--at-rest` (renamed; sets `apply_agonist_protocol=False`) |
| `runPlateletSim.py` defaults | `ip3_forced=True` | `apply_agonist_protocol=True` |
| `runPhase3.py` | `ip3_forced=True` (Dolan IP3 overlay drives Ca²⁺) | `apply_agonist_protocol=True` with GPCR agonist profile that produces a comparable Ca²⁺ transient (see risk §8.1) |
| `wholecell/webapp/tabs/configure.py:PRESETS` | `ip3_forced` in three presets | `apply_agonist_protocol`; "IP3 transient" preset renamed to "Agonist stimulation" |
| `wholecell/webapp/jobs.py` | `ip3_forced` config field | `apply_agonist_protocol` |
| `CLAUDE.md` run-time conditions table | row for `--no-ip3-forcing` | row for `--at-rest` |
| Phase 3 acceptance JSON / lab books | `ip3_forced` key | `apply_agonist_protocol` |

The metadata schema in `runscripts/manual/runPlateletSim.py:write_metadata`
also loses the `ip3_forced` field. Old runs on disk keep their metadata
as-is; the analysis tooling shouldn't care (it reads the trace, not the
metadata).

## 5. Files touched

Split across two issues — see §9 below.

### 5a. Prep refactor (issue A)

| File | Change | Est. LoC |
|---|---|---|
| `reconstruction/platelet/dataclasses/process/calcium_signalling.py` | Live-read peaks (3 fns); delete `ip3_forcing_uM`; delete IP3-overlay block in `_ode_rhs`; rename flag `ip3_forced` → `apply_agonist_protocol` | ~40 |
| `models/platelet/processes/calcium_dynamics.py` | Rename `_ip3_forced` → `_apply_agonist_protocol`; pass new name into solver | ~10 |
| `runscripts/manual/runPlateletSim.py` | New `adp_peak_uM` / `thrombin_peak_nM` kwargs; rename CLI flag; drop `ip3_forced` metadata field | ~25 |
| `runscripts/manual/runPhase3.py` | Switch from IP3-overlay stimulation to GPCR-driven stimulation; possibly retune acceptance criteria (see §8.1) | ~15 |
| `wholecell/webapp/tabs/configure.py` | Rename preset field + label; drop "IP3 forcing" checkbox | ~10 |
| `wholecell/webapp/tabs/runs.py` | Update condition-summary rendering (it currently prints "IP3 ON" / "rest") | ~5 |
| `wholecell/webapp/jobs.py` | Rename config field | ~3 |
| `CLAUDE.md` | Update run-time conditions table; remove `--no-ip3-forcing` row | ~15 |
| `models/platelet/tests/` | Update any tests that pass `ip3_forced`; add coverage that peak constants are live-read | ~30 |

Total: ~150 changed, single PR.

### 5b. Dose sweep (issue B, blocked on A)

| File | Change | Est. LoC |
|---|---|---|
| `runscripts/manual/runDoseSweep.py` | **New** — sweep driver, calls `run_platelet_sim` per cell, writes CSV / NPZ / static plots | ~180 |
| `wholecell/webapp/tabs/dose_response.py` | **New** — Plotly surface + heatmap tab reading `sweep_summary.json` | ~120 |
| `wholecell/webapp/app.py` | Register new tab | ~5 |
| `models/platelet/tests/` | Smoke test the sweep driver on a tiny grid | ~20 |

Total: ~325 new, single PR.

## 6. Order of work

**Prep refactor (issue A) — merge first:**

1. **Live-read peaks** + a 4-line test that `cs_mod.ADP_PEAK_UM = 0.1`
   actually changes the trajectory. Smallest, lowest risk, verifies
   the mechanism end-to-end before any deletion.
2. **Delete IP3 overlay path** — `ip3_forcing_uM` and its call in
   `_ode_rhs`. Run existing tests; they should still pass except for
   ones that explicitly probe the overlay (those get updated to drive
   via agonists or get deleted).
3. **Rename + plumb** the remaining flag through
   `molecules_to_next_time_step` → `CalciumDynamics` →
   `run_platelet_sim` → CLI. Update all callers (`runPhase3.py`,
   webapp presets, jobs config, metadata writers, lab-book notes).
4. **Re-validate Phase 3** — confirm the SOCE differential and dts-min
   acceptance criteria still pass under GPCR-driven stimulation. If
   not, retune the ADP / thrombin profile (not the SOCE biology) so
   the Ca²⁺ transient shape matches the Dolan target. Capture the
   choice in `phase3_summary.json`.
5. **Docs** — `CLAUDE.md` run-time conditions table, this design doc,
   any lab-book references to `ip3_forced`.

**Dose sweep (issue B) — depends on A landing:**

6. **Write the sweep driver** against the new kwargs. No
   monkey-patch. CSV + NPZ + static PNGs.
7. **Webapp tab** + register in `app.py`.
8. **Smoke test** + a 5×5 200 s reference run; commit the resulting
   PNGs to `reports/figures/` so the dissertation has the headline
   figure under version control.

## 7. What this design explicitly does not do

- **No `StimulusConfig` dataclass.** Would be cleaner — config-as-data,
  easy to log / diff / persist — but it's a much bigger touch on every
  runscript and `SimulationDataPlatelet` wiring. Worth doing later if
  the number of stimulus knobs keeps growing; not justified by today's
  four flags. See `kinetics-as-data-sketch.md` for the broader
  config-as-data direction; the same logic applies to stimulus config.
- **No competitive-antagonist species.** Sweeping the agonist peak
  directly is mathematically equivalent to a competitive antagonist at
  fixed Kᵢ, and it requires no new reactions. Adding real antagonists
  (`cangrelor[e]`, `vorapaxar[e]`) with binding kinetics is a separate
  v0.5 question and needs Kᵢ values from literature.
- **No scrape of PlateletWeb.** Discussed and rejected earlier — see
  GH issue #42.
- **No model-level expansion** (P2Y12 / desensitisation / GPVI). Out
  of scope here; tracked under their own issues.

## 8. Risks

1. **Phase 3 validation may need retuning.** Dolan 2014 Fig. 4 was
   reproduced by *forcing* IP3 directly with the Fig. S2 curve. After
   the refactor the only IP3 source is the GPCR cascade, so the Ca²⁺
   transient must come from an ADP / thrombin profile that drives the
   cascade to produce a comparable IP3 spike. The acceptance criteria
   (peak ≥ X nM, DTS min ≤ Y µM, SOCE differential ≥ 100 nM) are
   biology-level and should still hold, but the *measured* peak may
   shift slightly. Plan: re-run Phase 3 as the last step of issue A;
   if criteria fail, tune the agonist time-course parameters
   (`THROMBIN_PEAK_NM`, `ADP_PEAK_UM`, their τ values), **not** the
   downstream rate constants. Document the chosen profile in
   `phase3_summary.json` and a lab-book entry.
2. **Other legacy forcing paths.** `gq_signal_uM` and any other v0.3
   forcing helpers should be reviewed for the same expunge-or-keep
   decision. If they're already dead code post-v0.4, delete them in
   the same PR.
3. **No legacy shim — old runscripts in branches will break** when
   they merge. Acceptable cost for a dissertation-stage repo, but
   worth flagging in the issue so any open work-in-progress branches
   are rebased rather than merged after issue A lands.
4. **Webapp jobs DB** persists the `ip3_forced` field across restarts
   if SQLite. Old job rows will have `ip3_forced` in their config
   JSON; new code should tolerate that key being absent. One-line
   `.get('apply_agonist_protocol', cfg.get('ip3_forced', True))` at
   read time covers it without a migration.

## 9. Tracking issues

- **#44 — Expunge `ip3_forced` + Dolan Fig. S2 IP3 overlay path.**
  Prep refactor; must merge before #45. Covers §3 (all three fixes),
  §4 (caller migration), §5a (files), §6 steps 1–5, §8 (risks).
- **#45 — 2-D agonist dose sweep.** Blocked on #44. Covers §5b
  (files), §6 steps 6–8, the user-facing deliverable in §1.

## 10. Pointers

- GH issue #42 — annotation ingestion (referenced PlateletWeb
  discussion; this dose-sweep work is the trigger for "expanding
  beyond Ca²⁺" mentioned there)
- `reconstruction/platelet/dataclasses/process/calcium_signalling.py`
  lines ~640–730 (GPCR cascade definitions), ~1238 (P2X1 gate),
  ~1283–1288 (GPCR agonist gate), ~922 (`_ode_rhs` signature)
- `runscripts/manual/runPhase3.py` — the existing two-condition
  driver pattern the sweep driver mirrors
- `wholecell/webapp/tabs/explore.py` — existing tab pattern the new
  Dose Response tab mirrors
