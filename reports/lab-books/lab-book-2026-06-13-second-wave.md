---
title: "Lab book — 2026-06-13: v0.61 downstream PKC effects complete + the autocrine second-wave experiment"
---

# Lab book — 2026-06-13: v0.61 downstream effects + second wave

## Context

v0.6 wired PKC's two *braking* feedbacks (P2Y1 desensitisation, PLCβ
phosphorylation). This session built the *output*/amplifying side — the v0.61
downstream effects scoped in
`reports/design/pkc-downstream-effects-2026-06-12.qmd` — and then answered the
obvious question: **do any of these changes actually do anything to the
response?**

All work is on branch `PKC-P2Y1-desensitisation` (not yet pushed; PR #54 still
covers v0.6 only).

## What landed (4 + 1 commits)

1. **Granule secretion (§1 Slice 1).** `GranuleSecretion` process relocates
   pre-existing cargo (`ADP[dg]`, `5HT[dg]`, `FGA[ag]` → `[e]`; `SELP[ag]` →
   `SELP_surface[pl]`) on a `PKC_active × Ca²⁺` coincidence gate with a
   resting-tone floor (resting secretion exactly zero). Additive — Ca²⁺ ODE
   untouched, goldens byte-identical.
2. **Autocrine ADP loop (§1 Slice 2).** Secreted `ADP[e]` fed back onto the
   P2Y1 drive in the ODE (`_ode_rhs` adds `secreted_adp_count × _UM_PER_COUNT_EX`).
   Pericellular volume `V_EX_L` (~66 fL) calibrated so full dense-granule release
   ≈ 10 µM (a **modelling choice** — local concentration is geometry/density
   dependent). Self-limited by ecto-NTPDase clearance (`ADP[e] → AMP[e]`).
3. **Thromboxane synthesis (§2 Slice A).** `ThromboxaneSynthesis` lumps
   cPLA₂→COX-1→synthase into one Ca²⁺×PKC-gated production term; `COX1_FACTOR`
   is the **aspirin knob** (0 → no TXA₂). `TXA2[e]` decays (t½≈30 s) to the
   stable ELISA metabolite `TXB2[e]`. Additive.
4. **Autocrine TXA₂ → TP → Gq (§2 Slice B).** TP receptor
   (`TP_inactive/active[pl]`, `[gpcr.tp]`) added to the ODE; active TP joins
   `total_active_R` (`+ tp_a`). Goldens regenerated (the 2 new ODE states
   perturb `at_rest` ~0.003 % — benign solver jitter; `default_activation` was
   byte-identical). **Dolan 5/5 preserved.**

**PKC is now the hub of four feedback loops** (2 brakes + 2 amplifiers) — the
push–pull substrate of the platelet second wave the design argued for.

## Key finding — store-limit hides the loops under a strong agonist

Under the **saturating default agonist** (thrombin + 10 µM ADP + ATP, 300 s),
the full v0.61 model is essentially indistinguishable from v0.6 on the headline
calcium outputs:

| output | v0.61 vs v0.6 (max diff) |
|---|---|
| cytosolic Ca²⁺ | 0.06 % (≈ 0.3 nM) — none |
| DTS store Ca²⁺ | 0.00 % |
| SOCE flux, pump ATP | < 0.3 % |
| IP₃ | ~4 % (+8.5 nM) |
| **P2Y1 desensitised fraction** | **~20 %** |

The store empties either way, so the cytosolic transient is pinned by the
SOCE/pump balance regardless of the extra Gq drive. The amplifiers' effect lives
in the receptor / messenger readouts (IP₃, P2Y1 desensitisation) and in the
entirely new outputs (P-selectin, secreted ADP/5-HT, TXA₂/TXB₂, TP activation).

## The second wave — where the loops *do* change Ca²⁺

The amplifiers bite when the primary stimulus is **weak and transient**, so
there is a recovery phase for them to act on. Scanned thrombin (always empties
the store — irreversible PAR cleavage is cumulative → no submaximal regime) and
ADP (reversible → graded). At **ADP 0.5 µM, ADP-only, 300 s**:

- All conditions reach the same store-limited **peak** (~315 nM).
- In the recovery phase they split: the open-loop model (v0.6) decays toward
  baseline (**117 nM** by 300 s); the closed-loop model stays activated
  (**210 nM**) — a **+93 nM (~80 %) sustained second wave** that does not exist
  in v0.6.
- **Decomposition:** "aspirin" (TXA₂ off, ADP loop on) sits on top of "full" in
  the Ca²⁺ panel → the sustained Ca²⁺ is carried by the **autocrine ADP loop**
  (secreted dense-granule ADP sustains P2Y1 → store stays depleted → SOCE holds
  Ca²⁺ up). TXA₂ adds Gq drive visible in **IP₃** (full > aspirin) but not in
  Ca²⁺ (store/SOCE-limited).

This is the classic biology: dense-granule/autocrine amplification and the
aspirin-sensitive second wave matter for *weak* agonists, not saturating ones.

## Committed as an experiment

`runscripts/manual/runSecondWave.py` — open vs closed-loop comparison at a weak
transient ADP pulse (3 conditions: `v06` / `aspirin` / `full`), writing
`second_wave_traces.png`, `.npz`, `_summary.json`. Smoke test in
`models/platelet/tests/sim/test_second_wave.py`.

Clean loop toggle added: **`cs_mod.AUTOCRINE_ADP_GAIN`** (module-level, read
live in `_ode_rhs`; default 1.0, 0 = open loop) — mirrors the `CA_EX_UM` /
`COX1_FACTOR` override pattern, so the experiment needs no monkeypatching.
`AUTOCRINE_ADP_GAIN × 1.0` is exact, so goldens stay byte-identical.

## Notes / next

- The TXA₂ loop gain (TP count, `[gpcr.tp]` affinity, `k_prod`) and the autocrine
  volume `V_EX_L` are the tunable knobs; current values are conservative,
  literature-anchored model choices.
- **Deferred:** a `runPerturbation.py` aspirin experiment (needs the runner
  generalised to override the scalar `COX1_FACTOR`, not a `cs_mod` dict);
  integrin §3 (per-cell PAC-1 affinity state only — aggregation is inter-cellular
  and out of single-cell reach).
- Branch has 5 unpushed commits — consider push / a v0.61 PR.
- Tests 68 → all pass, mypy clean. Goldens regenerated once (Slice B).
