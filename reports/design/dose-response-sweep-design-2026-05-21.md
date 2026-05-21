---
title: "Dose-Response Sweep — Design Plan"
author: "Steve Haigh"
date: "2026-05-21"
---

# Dose-Response Sweep — Design Plan

**Status:** design plan, no code yet. Captures what needs to change in
the model to support a 2-D agonist dose sweep cleanly (no global
mutation, no monkey-patching, no flag double-duty).

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

### 3.2 Split `ip3_forced` into two flags

Introduce two orthogonal booleans:

| Flag | Controls |
|---|---|
| `apply_dolan_ip3_curve` | The Fig. S2 IP3 overlay onto `IP3[c]`. Original meaning of `ip3_forced`. |
| `apply_agonist_protocol` | Whether ADP / thrombin / ATP\_ex are applied at all (otherwise they snap to rest). |

Three call sites in `_ode_rhs` migrate:

| Today | New flag |
|---|---|
| IP3 forcing block (Fig. S2 overlay) | `apply_dolan_ip3_curve` |
| ATP\_ex gate for P2X1 (~line 1238) | `apply_agonist_protocol` |
| ADP / thrombin gate for GPCRs (~line 1283) | `apply_agonist_protocol` |

Both flags propagate up through:

```
CalciumSignalling.molecules_to_next_time_step(
    counts, dt, t_sim,
    apply_dolan_ip3_curve, apply_agonist_protocol, ip3_delay,
)
    └── _ode_rhs(t, y, t_sim_start,
                 apply_dolan_ip3_curve, apply_agonist_protocol, ip3_delay)
```

### 3.3 Thread agonist peaks through the call chain

Mirror the existing `ca_ex_mM` pattern in `run_platelet_sim`:

```python
def run_platelet_sim(
        ...,
        ca_ex_mM=DEFAULT_CA_EX_MM,
        adp_peak_uM=None,
        thrombin_peak_nM=None,
        apply_dolan_ip3_curve=True,
        apply_agonist_protocol=True,
        ip3_delay=DEFAULT_IP3_DELAY,
        ):
    cs_mod.CA_EX_UM = float(ca_ex_mM) * 1000.0
    if adp_peak_uM is not None:
        cs_mod.ADP_PEAK_UM = float(adp_peak_uM)
    if thrombin_peak_nM is not None:
        cs_mod.THROMBIN_PEAK_NM = float(thrombin_peak_nM)
    CalciumDynamics._apply_dolan_ip3_curve = apply_dolan_ip3_curve
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
    apply_dolan_ip3_curve=False,
    apply_agonist_protocol=True,
    length_sec=200,
)
```

## 4. Backward compatibility

Existing callers — `runPlateletSim.py`, `runPhase3.py`, webapp presets
in `wholecell/webapp/tabs/configure.py:PRESETS` — all use `ip3_forced`
with both meanings coupled. Preserve their behaviour with a legacy
shim in `run_platelet_sim`:

```python
def run_platelet_sim(..., ip3_forced=None, ...):
    if ip3_forced is not None:
        # Legacy kwarg — both channels follow the single flag.
        apply_dolan_ip3_curve = ip3_forced
        apply_agonist_protocol = ip3_forced
```

So:

- The `--no-ip3-forcing` CLI flag still produces an at-rest sim
- The Phase 3 driver (Dolan Fig. 4) still runs with both channels on
- The Resting webapp preset still maps to both channels off

Add a deprecation note in the docstring. The legacy kwarg can be
removed later — it's a one-line edit when the time comes.

## 5. Files touched

| File | Change | Est. LoC |
|---|---|---|
| `reconstruction/platelet/dataclasses/process/calcium_signalling.py` | Live-read peaks (3 fns); split flag in `_ode_rhs` + `molecules_to_next_time_step` | ~30 |
| `models/platelet/processes/calcium_dynamics.py` | Add `_apply_dolan_ip3_curve` / `_apply_agonist_protocol` class attrs; pass both into solver | ~10 |
| `runscripts/manual/runPlateletSim.py` | New kwargs + CLI flags + legacy `ip3_forced` shim | ~25 |
| `runscripts/manual/runPhase3.py` | None (legacy kwarg still works) | 0 |
| `runscripts/manual/runDoseSweep.py` | **New** — sweep driver, calls `run_platelet_sim` per cell, writes CSV / NPZ / static plots | ~180 |
| `wholecell/webapp/tabs/dose_response.py` | **New** — Plotly surface + heatmap tab reading `sweep_summary.json` | ~120 |
| `wholecell/webapp/app.py` | Register new tab | ~5 |
| `wholecell/webapp/tabs/configure.py` | *Optional* — separate checkbox for the new flag, or keep coupled in presets | ~5 |
| `models/platelet/tests/` | Cover the four flag combinations × at least one agonist peak | ~30 |
| `CLAUDE.md` | Update run-time conditions table to describe both flags | ~15 |

Total: ~420 new + ~85 changed, single PR.

## 6. Order of work

1. **Live-read peaks** + a 4-line test that `cs_mod.ADP_PEAK_UM = 0.1`
   actually changes the trajectory. Smallest, lowest risk, verifies
   the mechanism end-to-end before any plumbing.
2. **Split the flag** in `_ode_rhs` only (no upstream plumbing yet) +
   a unit test enumerating the four flag combinations.
3. **Plumb through** `molecules_to_next_time_step` → `CalciumDynamics`
   → `run_platelet_sim`. Legacy `ip3_forced` shim included. Confirm
   `runPhase3.py` and existing tests still pass byte-for-byte.
4. **Write the sweep driver** against the new kwargs. No monkey-patch.
5. **Webapp tab** + register in `app.py`.
6. **Docs** — update `CLAUDE.md` run-time conditions table.

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

1. **The `gq_signal_uM` legacy forcing path** (deprecated in v0.4 / #9
   but possibly still referenced) may also need migrating. Check before
   step 2.
2. **Webapp presets** currently express stimulus combinations as a
   single `ip3_forced` value in `PRESETS`. After the split, the presets
   should map to explicit `(apply_dolan_ip3_curve,
   apply_agonist_protocol)` tuples to make the new combinations
   reachable from the UI. If we don't update them, the new flag is
   shell-only.
3. **Existing CLAUDE.md docs and lab books** reference
   `--no-ip3-forcing` semantics. The legacy shim preserves behaviour
   so the docs remain accurate, but the new conditions table needs
   adding so future-me knows both knobs exist.

## 9. Pointers

- GH issue #42 — annotation ingestion (referenced PlateletWeb
  discussion; this dose-sweep work is the trigger for "expanding
  beyond Ca²⁺" mentioned there)
- `reconstruction/platelet/dataclasses/process/calcium_signalling.py`
  lines ~640–730 (GPCR cascade definitions), ~1238 (P2X1 gate),
  ~1283–1288 (GPCR agonist gate), ~922 (`_ode_rhs` signature)
- `runscripts/manual/runPhase3.py` — the existing two-condition
  driver pattern this sweep driver mirrors
- `wholecell/webapp/tabs/explore.py` — existing tab pattern the new
  Dose Response tab mirrors
