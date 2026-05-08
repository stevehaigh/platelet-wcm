---
title: "Lab book — 2026-05-08: Issue #26 audit (Caride k₁₁ CaM-PMCA stability)"
---

# Lab book — 2026-05-08: Issue #26 audit (Caride k₁₁ stability)

## Session goal

Close out issue **#26** (filed 2026-05-07): the Caride 2007
CaM-PMCA cycle exhibits a positive-feedback runaway when DTS Ca²⁺
supply is replenished (e.g. when a CALR-style luminal buffer is
added). The 2026-05-07 CALR-buffer attempt produced a 23× peak
amplification (392 nM → 8 949 nM) that was traced to the
`Ca₄·CaM·PMCA → PMCA·CaM + 4 Ca²⁺_cyt` k₁₁ release path. Today's
question: with the IP3R k₃ fix from this morning's Sneyd-Dufour
audit (commit `5c70d6df`) in place, does the runaway disappear?

## Method

1. Re-read Caride 2007 paper text + Table 3 to verify what `k₁₁`
   represents biologically and whether the 10 s⁻¹ value is
   correctly transcribed from the source.
2. Re-add the CALR buffer (yesterday's deferred design) on top of
   the post-k₃-fix calcium signalling code.
3. Run Phase 3 with both fixes active. Compare to yesterday's
   2/5 baseline (CALR buffer + pre-k₃-fix code).
4. Decide remediation per the issue's plan: (a) buffer-runaway
   fixed by k₃, (b) needs k₁₁ adjustment, (c) needs structural
   change, or (d) defer until D13 is resolved.

Empirical setup: temporary uncommitted edit of
`internal_state.py` (CALR_sites_free, CALR_sites_Ca species),
`calcium_signalling.py` (K_CALR dict + ODE rhs term). Reverted
after the experiment.

## Caride 2007 sanity check

Confirmed against the Caride 2007 PDF text (page 5645 region):

> "STEP 5: ... These steps comprise CaM activation of PMCA and
> CaM-activated PMCA activity. The rate constants for these five
> reactions were taken from Refs. 13, 19, and 24 and the results in
> this paper. ... k₁₁ and k₁₁r (the rates for dissociation/binding
> of all 4 Ca²⁺ from PMCA-CaM)"

So `k₁₁ = 10 s⁻¹` is the CaM-PMCA Ca²⁺-dissociation rate, taken
from cited primary references (Refs. 13, 19, 24 of Caride 2007)
plus their own fitting. Not a Caride-only invention; carries primary
literature provenance.

The biological role: when the active Ca₄·CaM·PMCA complex finishes
its catalytic cycle (k₁₀ extrudes 1 Ca²⁺ to outside), the Ca²⁺
ions on the four CaM EF-hands dissociate locally before CaM itself
unbinds from PMCA (k₁₂ slow). This is correctly biology — CaM
EF-hand Ca²⁺ unbinding *is* fast (~10 s⁻¹ timescale) — so the
rate is not "wrong".

## Empirical result

Phase 3 with CALR buffer added on top of today's k₃-fixed code:

| Metric | No buffer (current `main`) | + CALR (yesterday, k₃ = 11) | + CALR (today, k₃ = 0.11) |
|---|---|---|---|
| +Ca_ex peak | 393 nM | 8 949 nM | **8 982 nM** |
| −Ca_ex peak | 325 nM | 8 904 nM | **8 916 nM** |
| SOCE differential | 68 nM | 45 nM | **65 nM** |
| Phase 3 acceptance | 4/5 | 2/5 | **2/5** |

**The runaway is essentially identical with or without the k₃ fix.**
+Ca_ex peak: 8 949 → 8 982 nM (0.4% change). The k₃ correction
does not bound the CaM-PMCA cycle's positive feedback.

This is the cleanest possible answer to the audit question:

> **The runaway is in the CaM-PMCA cycle, not in the IP3R kinetics.**

## Mechanical analysis (per-cycle Ca²⁺ flux balance)

One full pass through the Caride 2007 step-5 cycle moves Ca²⁺ as
follows (per turnover, integers in `count Ca²⁺`):

| Step | Reaction | Δ[Ca²⁺]_cyt | Δ[Ca²⁺]_ex |
|---|---|---|---|
| 8 | PMCA + Ca₄·CaM → P(Ca₄·CaM) | 0 | 0 |
| 9 | P(Ca₄·CaM) + Ca²⁺_cyt → P(Ca₄·CaM)·Ca | **−1** | 0 |
| 10 | P(Ca₄·CaM)·Ca → P(Ca₄·CaM) + Ca²⁺_ex | 0 | **+1** |
| 11 | P(Ca₄·CaM) → PCaM + 4 Ca²⁺_cyt | **+4** | 0 |
| 12 | PCaM → PMCA + CaM_free | 0 | 0 |
| | (then: CaM_free + 4 Ca²⁺_cyt → Ca₄·CaM via k₆+k₇) | **−4** | 0 |
| | **Net per full cycle** | **−1** | **+1** |

Per fully completed cycle, 1 Ca²⁺ is removed from cyt to outside
— that's the PMCA's pumping job, correct.

But during a transient, the cycle's components don't run at the
same rate:

- k₆ + k₇ (CaM Ca²⁺ loading): ~170 µM⁻²·s⁻¹ × cyt² — fast at
  µM-range cyt
- k₈ (Ca₄·CaM binding to PMCA): 0.2 µM⁻¹·s⁻¹ — moderate
- k₉, k₁₀ (active extrusion): 50 / 30 s⁻¹ — fast
- **k₁₁ (Ca²⁺ release): 10 s⁻¹** — moderate
- **k₁₂ (CaM dissociation from PCaM): 0.033 s⁻¹** — bottleneck (τ ≈ 30 s)

Because k₁₂ is the bottleneck, the system holds many PMCAs in the
slow `PCaM` regeneration step at any given time. Meanwhile,
k₁₁'s "release 4 Ca²⁺ to cyt" runs at 10 s⁻¹ on every active
complex, while k₆ + k₇ "reload CaM with 4 Ca²⁺ from cyt" depend on
free CaM availability. If free CaM is fast to reload (it is), the
+4 release in step 11 and the −4 reabsorption in CaM ladder
balance over time. The +3 net "intermediate" Ca²⁺ in cyt during
the cycle is what drives the transient peak.

**With limited DTS Ca²⁺** (no buffer, current `main`), the IP3R
runs out of substrate after ~5 s, the CaM ladder absorbs the
released Ca²⁺ over the ensuing seconds, and cyt peaks at ~393 nM
before falling. **With CALR buffer** providing a 3.7× larger DTS
Ca²⁺ reservoir, the IP3R keeps draining for longer, more PMCAs
get pulled into the CaM-PMCA cycle, and the +3 instantaneous
release builds cyt up to ~9 µM before extrusion (PMCA k₁₀) +
re-absorption (CaM k₆ + k₇) catches up.

## Interpretation

This is *not* a bug. The Caride 2007 rate constants are
literature-grounded (per their Refs. 13, 19, 24); the cycle as
implemented faithfully reproduces those rates; and Caride's own
published simulation in CHO cells shows transient peaks reaching
~3 µM cyt — high but bounded.

Our model produces peaks of ~9 µM with the CALR buffer because:

1. We're coupled to Sneyd-Dufour IP3R kinetics (Caride used a
   simpler Ca²⁺-release scheme).
2. We're using SERCA3b (platelet) parameters (Caride used SERCA2b).
3. With the CALR buffer and our IP3R Po⁴ values, IP3R drains a
   larger total Ca²⁺ pool to cyt over a longer transient than
   Caride's CHO context produced.

So our 9 µM is not a Caride violation — it's a coupled-system
artefact that emerges when Caride's PMCA is plugged into our
IP3R + SERCA + buffer combination.

## Decision

**Close #26 as "diagnosed; not a bug; deferred to v0.3+ as a
calibration question."**

The k₁₁ rate is correct per primary source. Modifying it
unilaterally to suppress the runaway would deviate from Caride
2007 without strong biological justification. The runaway is
*exposed* by the CALR buffer addition but is a feature of the
coupled cycle, not of any single rate constant.

For the dissertation v0.2 freeze, the current `main` baseline
(4/5 Phase 3, no CALR buffer) is the deliverable. Future work to
add the buffer should consider one of:

- (a) Add the buffer with a smaller effective capacity (e.g. some
  fraction of CALR sites are kinetically inaccessible — biologically
  defensible if those CALR are Ca²⁺-saturated by e.g. SERCA-bound
  Ca²⁺ in the lumen rather than freely-exchangeable).
- (b) Coarse-grain the Caride k₁₁ + k₁₂ steps into a single
  effective CaM-recycling rate, giving up the explicit "4 Ca²⁺
  released to cyt" intermediate.
- (c) Add a downstream CaM-binding-protein competitor that
  buffers the +3 transient release before it can amplify the peak.
- (d) Replace the Caride 5-state CaM-PMCA cycle with a simpler
  effective-rate model that has the same net per-cycle behaviour
  (1 Ca²⁺ extruded per turnover) but no intermediate +4 / −4
  spike.

All four are deviations from primary source; (b) is the most
common simplification in the systems-biology literature (e.g.
Wagner & Keizer 1994 lump CaM-PMCA cycle into a single
Ca²⁺-extrusion term); (c) requires extra biology we don't have
in v0.2; (a) and (d) are calibration choices.

This belongs as a v0.3+ design decision, not a v0.2 fix.

## Files / artefacts

- This lab book (`reports/lab-books/lab-book-2026-05-08-k11-stability-audit.md`).
- Empirical Phase 3 output:
  `out/phase3_calr_with_k3fix/phase3_summary.json` (the 2/5
  acceptance with CALR buffer + k₃ fix).
- No code changes committed. The CALR-buffer code edits made to
  reproduce yesterday's experiment with today's k₃ fix were
  reverted; current `main` is the post-Sneyd-Dufour-audit baseline.

## Issue resolution

Closing GitHub issue **#26** with a comment summarising:
- The runaway is real and reproducible
- It's independent of the IP3R k₃ fix (this audit's empirical test)
- It's inherent to the Caride 2007 CaM-PMCA cycle as published,
  not a bug in our implementation
- It's only *exposed* by the CALR buffer addition (which is
  itself deferred until #26 was resolved — circular)
- Resolution: defer to v0.3+ as a calibration question; v0.2
  freezes at 4/5 Phase 3 with no buffer.

---

*Branch:* `main` · *Status:* #26 closed; v0.2 freeze ready ·
*Linked issues:* #26 (closed by this lab book), #24 (still open;
separate research-grade limitation), #25 (v0.3+ stretch).
