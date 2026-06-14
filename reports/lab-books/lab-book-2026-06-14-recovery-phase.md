---
title: "Lab book — 2026-06-14: recovery-phase validation — the sustained plateau is the autocrine second wave"
---

# Lab book — 2026-06-14: recovery-phase validation

## Context

The 30-s "Dolan 5/5" criteria score only the **peak** of the Ca²⁺ transient.
The **recovery phase** (t > 30 s) was pinned by a strict `xfail` in
`test_validation_targets.py::TestDolanRecoveryPhase` that read: *"the sustained
+Ca plateau (~430 nM) sits above Dolan's 200–275 nM band … the recovery phase
is not yet calibrated."* This session investigated **why**, choosing the
biology-first route (understand the mechanism before tuning constants) rather
than a multi-knob constant fit.

This follows the integrin §3 work (same session); branch `integrin-v0.63`.

## Diagnosis (read-only sweeps)

The peak is already correct in **both** conditions (+Ca 435 ∈ Dolan 400–500;
EDTA 317 ∈ 225–325). Two *distinct* recovery-phase gaps:

1. **+Ca plateau too high** (model ~432 vs Dolan 200–275 nM). The model's DTS
   store empties to ~0 µM, whereas Dolan's only partially depletes
   (120–180 µM). An empty store pins SOCE wide open (~19 nM/s), holding cyt at
   430.
2. **EDTA doesn't return to baseline** (model ~235 vs Dolan ~75 nM by ~80 s).

Strengthening the PLCβ brake lowered IP₃ dramatically (223 → 89 nM) but barely
moved the plateau (432 → 392) — so the plateau is **SOCE / empty-store
limited, not IP₃ limited** once the store is drained. PMCA could not be raised
to clear EDTA without breaking the resting fixed point (×8 PMCA → resting cyt
collapses 100 → 34). No single knob closed both gaps.

## Mechanism (reading the ODE)

- **The IP₃R is a Nernst-driven channel** (`calcium_signalling.py` ~L941):
  whenever it is open it equilibrates DTS ⇌ cyt, so an open IP₃R *cannot hold a
  store gradient*. The store empties because IP₃ stays high; the `m∞⁴`
  dependence makes refill a sharp threshold (~IP₃ < 75 nM → store refills).
- **The sustained IP₃ is driven by the v0.61 autocrine loops, not the
  agonist.** Even under a transient, reversible ADP stimulus (ADP clears with
  τ ≈ 30 s), IP₃ stayed pinned at ~250 nM. Isolating the loops showed the
  culprit: with **TXA₂ → TP → Gq** (and secreted ADP → P2Y1) disabled, IP₃
  finally decays to its 50 nM baseline. The autocrine amplifiers sustain
  Gq → IP₃ and hold the cell activated — the platelet **"second wave."**
- **EDTA clearance is buffer-limited.** With the IP₃R closed and the store
  isolated, cytosolic free Ca²⁺ falls only as the large CaM / gelsolin
  bound-Ca pool releases while extrusion is capacity-limited — a structurally
  long time constant (cyt/store recovery lags IP₃ by ~100 s).

## The conclusion

**Dolan (2014) drove a *transient* IP₃ dose and modelled neither thromboxane
nor granule secretion.** Our model adds those autocrine positive-feedback loops
— biology Dolan omitted — and they sustain the response. So the high sustained
plateau is a model **prediction (extra biology)**, not a recovery-phase
calibration defect. The decisive test: with the loops off and a transient
reversible-ADP stimulus (the **Dolan-equivalent** configuration), the model
returns to baseline (cyt → 100 nM, IP₃ → 50 nM, the DTS store refills),
reproducing Dolan's recovery (Fig. 4C). The Ca-handling machinery
(IP₃R / SERCA / SOCE / PMCA / NCX) is sound; the divergence is the amplifiers.

| condition (400 s) | cyt_end | IP₃_end | DTS_end |
|---|---|---|---|
| full v0.61, +Ca (loops on)        | 353 nM | 224 nM | 0 µM  |
| full v0.61, EDTA (loops on)       | 178 nM | 264 nM | 0 µM  |
| Dolan-equiv., +Ca (loops off)     | 100 nM | 50 nM  | 68 µM |
| Dolan-equiv., EDTA (loops off)    | 100 nM | 51 nM  | 43 µM |

## Changes landed

1. **NCX bug fix** (`calcium_signalling.py`). The Na⁺/Ca²⁺ exchanger — a Ca²⁺
   *extruder* to the (infinite) extracellular sink, driven by the Na⁺ gradient
   — was gated *inside* the `if ca_ex_uM > 0` block, i.e. switched off under
   EDTA. It had been lumped with the genuine influx pathways (SOCE, PM leak,
   P2X1, which *do* need external Ca²⁺). Moved out of the conditional so it runs
   in both conditions; its allosteric gate keeps it silent at rest (5.8 ions/s
   at 100 nM) and strong at high cyt (236 ions/s at 300 nM). **Both +Ca goldens
   stay byte-identical** (NCX already ran there); EDTA now extrudes correctly.
2. **Re-scoped the validation** (`test_validation_targets.py`). Removed the
   misleading `xfail` ("not yet calibrated") and added a passing
   characterisation test `test_autocrine_loops_sustain_else_ip3_recovers`: the
   full model sustains IP₃ (second wave) while the Dolan-equivalent recovers it
   toward baseline. The class docstring now states the framing.
3. **Figure** (`runscripts/manual/plotRecoveryPhase.py`). 3-panel recovery
   contrast (cyt Ca²⁺, IP₃, DTS store) × {full, Dolan-equivalent} × {+Ca, EDTA}
   with a detailed caption — the headline recovery-phase artefact.

## What was *not* done, and why

**No PLCβ-brake strengthening.** `k_plcb_phos` is explicitly **model-choice**
(the TOML header notes Purvis 2008 Table 1 has the known 100× transcription
erratum, so its rows are not transferred verbatim) — there is no specific
literature value we are under-shooting, so strengthening it has no biological
basis. It was also rendered moot by the framing: the sustained plateau is the
autocrine second wave (correct biology), so suppressing it with a stronger brake
would hide a genuine feature rather than fix a defect.

## Status

`test_validation_targets.py` 7/7 pass (no xfail). Byte-identical goldens
preserved → Dolan 5/5 intact. Resting drift (cyt 100 → 55 over 200 s at rest)
is pre-existing and noted for a later session.
