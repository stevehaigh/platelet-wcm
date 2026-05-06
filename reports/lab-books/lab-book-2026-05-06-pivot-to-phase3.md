---
title: "Lab book — 2026-05-06: Path B falsified, pivot to path A for Phase 3"
---

# Lab book — 2026-05-06

## Session summary

Tried path B (tune SERCA `k_release_r` to recover Dolan resting DTS at
250 µM). A sweep across 8 orders of magnitude (4 × 10⁻³ → 0) showed
zero effect on the ODE's resting fixed point: every value gives the
same end state (cyt ≈ 3 nM, DTS = 0 at γ_IP3R = 10 pS). A direct
rate-balance dump at the IC pinned down why — `k_release_r` controls
only the reverse leak from DTS into E2P·Ca when the store is full;
forward delivery is hard-capped by E2P·Ca occupancy upstream.

Pivoting to **path A**: accept the current resting state as the v0.2
working point, frame Phase 3 as a *transient-shape* comparison (with
vs without extracellular Ca²⁺), and document the deviation rather
than calibrating it away.

---

## Path B finding

### Sweep at γ = 10 pS, no IP3 forcing, 300 s integrations

| `k_release_r` (µM⁻²·s⁻¹) | cyt t=300 (nM) | DTS t=300 (µM) | E2P·Ca | E2P |
|---|---|---|---|---|
| 4 × 10⁻³ (Purvis primary) | 3.32 | 0.0 | 0 | 12 |
| 1 × 10⁻³ | 3.32 | 0.0 | 0 | 12 |
| 1 × 10⁻⁴ | 3.32 | 0.0 | 0 | 12 |
| 1 × 10⁻⁵ | 3.32 | 0.0 | 0 | 12 |
| 1 × 10⁻⁷ | 3.32 | 0.0 | 0 | 12 |
| 1 × 10⁻¹⁰ | 3.32 | 0.0 | 0 | 12 |

(All values run to the same fixed point; the sweep was effectively a
no-op on the long-time behaviour.)

### Rate balance at the Dolan IC, varying `k_release_r`

| `k_release_r` | d[cyt]/dt (ions/s) | d[DTS]/dt (ions/s) | net SERCA into DTS |
|---|---|---|---|
| 4 × 10⁻³ | +111 913 | −118 288 | **−6 000** (reverse) |
| 1 × 10⁻³ | +111 913 | −107 788 | +4 500 |
| 1 × 10⁻⁴ | +111 913 | −104 638 | +7 650 |
| 1 × 10⁻⁵ | +111 913 | −104 323 | +7 965 |
| 1 × 10⁻⁷ | +111 913 | −104 288 | **+8 000 (ceiling)** |
| 1 × 10⁻¹⁰ | +111 913 | −104 288 | +8 000 |
| 0 | +111 913 | −104 288 | +8 000 |

(Net SERCA = d[DTS]/dt + 112 288 since IP3R drains DTS at ~112 k ions/s
at the IC — see lab-book 2026-05-05 §"Rate balance at the Dolan IC".)

**The cycle's maximum forward delivery at full DTS is ~8 000 ions/s**,
set by `k_release_f × E2P·Ca = 1 000 × 4 = 4 000 events/s = 8 000
ions/s`. E2P·Ca itself is fixed by upstream supply (the binding step
at the IC delivers tracer levels of E1·Ca through the cycle). Setting
`k_release_r` to zero buys exactly the 6 000 ions/s reverse leak that
the original value introduced — nothing more.

IP3R inflow at γ = 10 pS at the IC = 112 000 ions/s. SERCA's ceiling
of 8 000 ions/s is **14× short** — `k_release_r` cannot bridge this
gap regardless of magnitude.

### What would actually work, and why we're not doing it

Path B as a calibration would require *at least two* simultaneous
moves:

1. **Reduce γ_IP3R** so IP3R inflow ≤ 8 000 ions/s at full DTS
   (lab-book 2026-05-01 originally proposed γ ≈ 0.6 fS). Phase 2a
   (lab-book 2026-05-05) showed any γ low enough to balance SERCA
   gives DTS *overshoot* (up to ~451 µM at γ = 6 × 10⁻¹⁶) and cyt
   *undershoot* (~3 nM). Same magnitude problem, opposite imbalance.
2. **Increase SERCA cycle throughput** by tuning `k_release_f`,
   `k_phos_f`, `k_conf_f` upward — multiple Purvis primary-source
   rate constants moved off-published. Multi-parameter calibration
   without a coherent biological rationale.
3. **Add a missing biology compartment** — mitochondrial Ca²⁺ buffer
   (issue #22), absorbing the cyt-DTS imbalance the published rates
   leave unaccounted for. Bigger scope, ~1–2 working days.

(1) is a multi-parameter scan with a thin defence. (2) departs from
primary sources without a published reason. (3) is the most
biologically defensible but blocks Phase 3 by the scope expansion it
implies. None is a single defensible lever, especially against the
mid-July deadline.

---

## Path A plan (re-affirmed)

Per the v0.2 design doc §6.8 D7, this deviation is now documented as
*the* known v0.2 limitation: the model's natural fixed point sits at
lower cyt and higher DTS than Dolan, motivating future work on
mitochondrial buffering (#22). The dissertation framing becomes
"validated transient peak amplitude and qualitative SOCE-dependent
shape against Dolan 2014 Fig. 4; resting absolute baseline differs by
factor X for reasons traced to SERCA cycle throughput vs IP3R
conductance, addressed in the v0.2.5 MCU extension."

### Phase 3 acceptance criteria (reframed for path A)

| Criterion | Target | Source |
|---|---|---|
| Peak cyt Ca²⁺ with `CA_EX = 1.2 mM` | 200–800 nM (Dolan ±30%) | already passing at 299.5 nM (lab-book 2026-05-01) |
| Peak cyt Ca²⁺ with `CA_EX = 0` (EDTA) | similar magnitude, slightly smaller | Dolan 2014 Fig. 4 qualitative |
| SOCE plateau present with `CA_EX = 1.2 mM` | sustained elevation 30–60 s | Dolan 2014 Fig. 4 |
| SOCE plateau absent with `CA_EX = 0` | faster decay, no plateau | Dolan 2014 Fig. 4 |

The peak passes one criterion already; the new validation is that the
*shape* responds correctly to extracellular Ca²⁺ removal.

### Phase 3 workflow

1. CLI flag `--ca-ex-mM <FLOAT>` (default 1.2) on `runPlateletSim.py`,
   threading through to override `CA_EX_UM` in the calcium signalling
   ODE.
2. Digitised Dolan 2014 Fig. 4 reference traces, committed under
   `reports/data/`.
3. New analysis plot `models/platelet/analysis/single/phase3_dolan_fig4.py`
   — two-panel comparison of cyt + DTS with Dolan reference overlay.
4. Phase 3 driver script `runscripts/manual/runPhase3.py` that runs
   both conditions into one output directory and triggers the analysis.
5. Pass/fail evaluation against the four criteria above; numerical
   comparison committed alongside the plot.
6. Regression test in `models/platelet/tests/` that locks the Phase 3
   numbers against future drift.
7. Lab-book entry for the Phase 3 results; design-doc §7.2 / §7.3
   updates from "OPEN / FUTURE" to "PASS / DEVIATION DOCUMENTED";
   close issue #20.

Effort estimate: ~2 working days end-to-end. Mostly mechanical — no
new biology debug.

---

## Files touched / inspected this session

| File | Change |
|---|---|
| `reconstruction/platelet/dataclasses/process/calcium_signalling.py` | Read-only — confirmed `K_SERCA['k_release_r'] = 4.0e-3` |
| `reports/lab-books/lab-book-2026-05-05-phase2a-investigation.md` | Reference for the Phase 2a diagnosis |

No source changes from this session. Path A implementation begins in a
follow-up commit.

---

*Branch:* `main` · *Last commit:* `a8239011` (Cleanup, 2026-05-06) ·
*Working tree:* `reports/design/calcium-dynamics-design.md` design-doc
rewrite still uncommitted.
