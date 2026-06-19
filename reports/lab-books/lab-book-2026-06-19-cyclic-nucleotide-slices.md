---
title: "Lab book — 2026-06-19 (cont.): inhibitory-axis Slices 1 + 4 — the cAMP-raising drugs and the secretion brake"
---

# Lab book — 2026-06-19 (cont.): cyclic-nucleotide pharmacology (Slices 1 + 4)

## Context

Continuation of the same day's inhibitory-axis work (after #10 / PRs #66, #68
merged). Pulled forward the two highest value-per-effort deferred slices from
issue #70: **Slice 4** (PKA brake on granule secretion) and **Slice 1** (the
cAMP-*raising* Gs arm — PGI₂/Gs/AC + explicit PDE3A — i.e. the prostacyclin and
PDE3-inhibitor drugs). Together with #10's P2Y12 (cAMP-lowering) arm, the model
now holds the **bidirectional** cAMP/PKA axis and the headline antiplatelet /
vasodilator pharmacology. The deeper Ca²⁺-core recalibration (Slice-1-mechanistic
and the clamp rework) stays parked on #70.

---

## 1. The PKA-headroom recalibration (the enabling change)

The blocker for cAMP-*raising* drugs: resting cAMP (1 µM) sat well above the PKA
EC50 (0.3 µM), so resting PKA was ~92 % active and could barely rise — PGI₂ could
only move it a few %, giving negligible suppression. Fix: **raise the PKA EC50 to
≈ resting cAMP** (`k_pka_uM` 0.3 → 1.0), so resting PKA is ~50 % with headroom in
**both** directions. This is also more biologically honest (resting PKA isn't
saturated). Because the brake is normalised to 1.0 at resting cAMP regardless of
`k_pka_uM` (`PKA_REST_FRAC = pka_active_frac(CAMP_REST_COUNT)`), **the resting
fixed point and Dolan 5/5 are still preserved by construction** — only the
drug-effect magnitudes change. Re-verified: #10's clopidogrel/VASP tests still
pass; VASP initial counts re-seeded to the new resting fraction (0.45).

## 2. Slice 1 — the cAMP-raising Gs arm (lean rate terms)

Per the design doc, AC/PDE entered as **rate terms, not species** (no PTGIR/Gs
state). The cAMP node became:

```
d[cAMP]/dt = V_AC·(1 + ac_gs_max·gs(PGI2) + forskolin)·(1 − i_gi_max·p2y12_frac)
             − k_pde·(1 − pde3_block)·cAMP
```

- `gs(PGI2) = pgi2_nM/(pgi2_nM + K_pgi2_nM)` — saturating Gs/AC stimulation.
- `RunConfig.pgi2_nM` (PGI₂/iloprost), `forskolin` (direct-AC fold), `pde3_block`
  (cilostazol/dipyridamole ∈ [0,1]). All default 0 → rest unchanged.
- Registered in the TUI RunSpec schema (the `test_schema_keys` invariant) under a
  new "Inhibitory / drugs (cAMP-raising)" group.

## 3. Slice 4 — PKA brake on secretion

Mirrored the #10 integrin brake: `GranuleSecretion` reads `cAMP[c]`, computes
`pka_brake_factor`, and multiplies it into the SNARE release rate (normalised to
1.0 at rest, so resting secretion stays exactly zero). `SecretionTrace` gains a
`pka_brake` column. Now the cAMP-raising drugs visibly brake dense-granule
release, and clopidogrel slightly dis-inhibits it — extending the drug readouts
from integrin onto secretion (lumi-aggregometry / P-selectin targets).

## 4. Result — bidirectional pharmacology (weak-agonist regime, ADP 0.5 µM)

| condition | cAMP | PAC-1 (t150) | VASP-P | ADP released |
|---|---|---|---|---|
| rest | 1.00 | 0 % | 0.45 | 0 % |
| control (ADP lowers cAMP) | ↓0.33 | 65.8 % | 0.29 | 99.6 % |
| clopidogrel (P2Y12 block) | 1.00 | 50.2 % | 0.45 | 97.8 % |
| **PGI₂ 50 nM** | ↑4.9 | **13.8 %** | 0.61 | **53 %** |
| **forskolin ×5** | ↑4.6 | 13.8 % | 0.61 | 53 % |
| **cilostazol (PDE3 0.8)** | ↑3.8 | 20.8 % | 0.59 | 64 % |

So the axis now runs both ways around the resting tone: autocrine ADP lowers
cAMP (more activation), clopidogrel blocks that (back toward rest), and the
cAMP-elevating drugs push above rest → strong suppression of **both** integrin
and secretion + a rise in VASP-P. New figure
`reports/figures/v0.7/cyclic_nucleotide_drugs.png`
(`plotInhibitoryAxis.py --figure drugs`).

## 5. Tests / checks

- New `TestCyclicNucleotideDrugs` (5 tests): PGI₂/PDE3-block raise cAMP;
  cAMP-raising suppresses PAC-1 and secretion; VASP-P rises. The 6 existing #10
  tests still pass (11 total in the file).
- Byte-identical goldens regenerated (the `k_pka` change alters the activated
  IP3R brake; rest unchanged). Dry mass unchanged (VASP total still 50 000).
- Full suite: platelet **122** + wholecell **117** pass; mypy clean.

## 6. Experiments & reproduction (the inhibitory-axis experiment set)

All experiments are driven by per-run `RunConfig` knobs (no monkeypatching);
each figure script runs its own sims. Drug knobs introduced across the axis:
`p2y12_block` (clopidogrel, #10), `pgi2_nM` / `forskolin` / `pde3_block`
(Slice 1), plus `cox1_factor` (aspirin, v0.61).

| Experiment | Command | Artifact / readout |
|---|---|---|
| **Mechanism** — ADP→P2Y12→cAMP↓→PKA↓→VASP-P↓ (standard agonist) | `python runscripts/manual/plotInhibitoryAxis.py --figure mechanism` | `reports/figures/v0.7/inhibitory_axis_mechanism.png` |
| **Antiplatelet treatments** — control vs aspirin vs clopidogrel vs both (weak agonist; PAC-1, Gαq, cAMP, VASP/PRI, TXA₂, Ca²⁺) | `… --figure treatments` | `reports/figures/v0.7/antiplatelet_treatments.png` |
| **Cyclic-nucleotide drugs** — control vs clopidogrel vs PGI₂ / forskolin / cilostazol (PAC-1, cAMP, VASP/PRI, secretion) | `… --figure drugs` | `reports/figures/v0.7/cyclic_nucleotide_drugs.png` |
| **Granule-secretion kinetics** (#8) | `python runscripts/manual/analysisPlatelet.py --plot granule_secretion <run>` | `analysis/single/granule_secretion.py` |
| **Behavioural assertions** (acceptance criteria, both directions) | `pytest models/platelet/tests/sim/test_inhibitory_axis.py` | `TestInhibitoryAxis` (#10) + `TestCyclicNucleotideDrugs` (Slices 1+4) |
| **Calibration probe** (rest invariant + drug dose panel) | inline `RunConfig(pgi2_nM=…, forskolin=…, pde3_block=…, p2y12_block=…)` → read `CalciumTrace`/`IntegrinTrace`/`SecretionTrace` | §4 table above |

Readout provenance (single source of truth): cAMP / PKA / VASP-P / P2Y12 →
`CalciumTrace`; PAC-1 → `IntegrinTrace.active_frac`; ADP release →
`SecretionTrace.adp_released_frac`; the PKA brake values → `IntegrinTrace` /
`SecretionTrace` `pka_brake`. Earlier same-day work (the P2Y12 axis itself, the
cytosolic-Ca²⁺ clamp probes that motivated braking the *outputs*, and the
"5/5" → validation-portfolio reframe) is in
`lab-book-2026-06-19-p2y12-inhibitory-axis.md`; the validation portfolio is
`reports/design/validation-map-2026-06-19.qmd`.

## State / next

- v0.7 inhibitory axis now spans **Slices 1, 2, 4** — the full cAMP/PKA "off"
  module (Gs raise + Gi lower) plus brakes on IP₃R, integrin, and secretion.
- Remaining on #70: **Slice 3** (NO/cGMP/PKG) and the deeper **mechanistic
  resting cAMP tone + Ca²⁺-clamp recalibration** (the only route to making the
  inhibitory axis visible on cytosolic Ca²⁺ itself). Post-dissertation scope.
