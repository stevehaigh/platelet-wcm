# Validation and regressions — an honest account

> This doc is deliberately candid. It explains what "Dolan 5/5" and the "golden"
> tests actually buy us, where they fall short, and why **preventing regressions
> as the model grows is an open, unsolved problem** in this project. It expands
> on `reports/design/validation-map-2026-06-19.qmd` (the canonical write-up) for
> a code-oriented reader. For the mechanics of running the tests see
> [`development-workflow.md`](development-workflow.md).

## TL;DR

- **"Dolan 5/5"** = the model reproduces Dolan & Diamond 2014 Fig. 4 (the
  cytosolic Ca²⁺ transient, with and without extracellular Ca²⁺). It was the
  right headline when the model *was* the Ca²⁺ transient.
- It is **no longer a sufficient measure of correctness** for a model that now
  spans the GPCR cascade, PKC feedback, secretion, thromboxane, integrin, and the
  inhibitory axis. We keep it as a **regression invariant**, not as evidence the
  new biology is right.
- The **regression suite** — the behavioural acceptance tests
  (`test_acceptance.py`) plus the structural bands (`test_regression.py`) — is
  good at catching *change*. It is **not** evidence of *correctness*, and it gets
  weaker as the model grows. **This is the open challenge.**

## What validation vs. regression mean concretely

Two different things often get conflated:

1. **The Dolan validation** — a *validation* claim: the Ca²⁺ core matches a
   specific published experiment (Fig. 4: resting level, peak, ±Ca_ex difference,
   store depletion/refill). Asserted behaviourally in
   `test_acceptance.py::TestDolanTransient` (paired peaks in band + the SOCE
   differential); `phase3_dolan_fig4.py` still computes the old five criteria to
   *annotate the plot*, but they are no longer a pass/fail gate.

2. **The regression suite** — *regression* guards, not validation:
   - `test_acceptance.py` — the headline **behavioural** contract: one band per
     biological result (resting equilibrium, the Dolan transient, the MCU /
     P2Y12 / COX-1 knockouts, the resting-quiescence invariant).
   - `test_regression.py` — lower-level physiological **tolerance bands** (peak
     Ca²⁺ within ±30% of a stated baseline, resting IP₃ 40–60 nM, dry mass within
     1%, DTS never drains below cytosol, SOCE flux ≥ 0).

A passing regression test says "this run behaves as it did last time," not "this
run is biologically correct."

> The former **byte-identical goldens** — bit-exact NPZ snapshots of two
> scenarios — were retired 2026-06-26: brittle, platform-sensitive (NumPy/SciPy
> version, FP hardware), and blind to whether a drift mattered. The behavioural
> bands above replace them. The all-or-nothing **"Dolan 5/5" count** was retired
> with them, in favour of the explicit per-result assertions.

## Why 5/5 is no longer the headline

Three reasons, each demonstrable from the current model (see the validation-map
doc for the probe data):

1. **It passes by construction.** Every layer added since the Ca²⁺ core is
   deliberately *normalised at rest* so the resting fixed point and the Dolan
   transient don't move (e.g. the PKA brake is exactly 1.0 at resting cAMP;
   terminal outputs like integrin don't touch the Ca²⁺ ODE at all). A criterion
   engineered to keep passing is a guard against accidental breakage — not
   evidence the new biology is correct.

2. **Low discriminating power (the "Ca²⁺ clamp").** Cytosolic Ca²⁺ is
   store-limited and SOCE-clamped, so it is nearly insensitive to any single
   release/clearance lever. Probes moved the peak by single-digit nM across large
   parameter swings (IP₃R open-probability gain 1→8: 293→298 nM; PMCA k_cat
   1.0→0.3×: 408→412 nM). The buffers (CaM/gelsolin/calreticulin) plus the
   SOCE↔clearance balance pin the level, so **many different parameterisations
   pass 5/5** — it can't tell them apart.

3. **Scope mismatch.** The clinically interesting biology has its *own*, stronger
   assays that the Ca²⁺ transient simply cannot speak to (VASP/PRI for P2Y12
   inhibition; PAC-1 flow cytometry for integrin; lumi-aggregometry for
   secretion; antiplatelet drug dose-response).

## The Ca²⁺ clamp as a stated limitation

The clamp is a model property worth stating plainly in the thesis, not hiding:

- **What it means.** The model currently **cannot quantitatively predict changes
  in cytosolic Ca²⁺ amplitude** driven by upstream modulation (PKC, P2Y12/cAMP).
  Those signals are real in the model but are absorbed by buffering and the SOCE
  clamp before they reach free cytosolic Ca²⁺.
- **How we work with it.** Feedback/inhibitory signals are validated on the
  readouts they actually move — the functional outputs (integrin, secretion,
  thromboxane) and signalling-node readouts (IP₃, cAMP, VASP-P) — not on free
  cytosolic Ca²⁺. This is why, e.g., the PKA brake acts on the αIIbβ3 rate
  (visible: clopidogrel lowers PAC-1) as well as on IP₃R (real but ~invisible on
  cyt Ca²⁺).
- **Practical consequence for developers.** "Knockout/loop/perturbation effects
  are invisible under the default saturating agonist." To *see* a change, isolate
  one agonist and read IP₃, or use the baseline overlay — don't expect cyt Ca²⁺
  to move.
- **Future work.** Revisiting the cytosolic buffer load and the SOCE-clamp
  calibration so upstream modulation reaches free Ca²⁺ is a candidate Ca²⁺-core
  project — and it would need its own re-validation against Dolan, which is
  exactly the regression-invariant role 5/5 is well suited to.

## The validation portfolio (what we judge the model on now)

Each subsystem is judged against the assay a wet lab would use.
Status: ✓ working comparison; ◐ partial/qualitative; ○ aspirational.

| Subsystem | Validation target (assay) | Status | Where |
|---|---|---|---|
| Ca²⁺ core | Dolan & Diamond 2014 Fig. 4 (±Ca_ex) | ✓ now a **regression invariant** | `test_acceptance.py`, `phase3_dolan_fig4.py` |
| PI cycle / IP₃ | resting IP₃ ≈ 50 nM; transient shape (Purvis 2008; Dolan Fig. S2) | ✓ | `test_regression.py` |
| PKC feedback | rapid P2Y1 desensitisation; IP₃ ↓ (Mundell 2006; Purvis 2008) | ◐ direction + band | `test_validation_targets.py` |
| Granule secretion | lumi-aggregometry (ATP/5-HT kinetics); P-selectin flow | ◐ shape | `granule_secretion.py` plot |
| Thromboxane | aspirin abolishes TXA₂; autocrine second wave | ✓ direction | `test_thromboxane.py` |
| Integrin (αIIbβ3) | PAC-1 flow cytometry; graded dose-response; drug shifts | ✓ | `plotInhibitoryAxis.py`, `test_integrin.py` |
| Inhibitory axis | VASP/PRI (P2Y12); PGI₂/forskolin Ca²⁺ suppression | ✓ VASP/PRI; ○ PGI₂ | `test_inhibitory_axis.py` |

This portfolio is the right way to *present* the model. It is **not yet a
solution to the regression problem** — see below.

## The open challenge: regression-proofing a growing model

This is the honest, unresolved part. As the model gains subsystems, "how do we
know a change didn't quietly break something?" gets harder, and our current tools
each have a real weakness:

- **Byte-identical goldens were brittle and uninformative — now retired.** They
  caught *any* change, including legitimate ones, and a failure told you *that*
  output moved, not *whether it mattered*; they were platform-sensitive
  (NumPy/SciPy version, FP hardware) and covered only **two** scenarios and
  **four** columns. `test_acceptance.py`'s per-result behavioural bands replaced
  them (candidate directions 2 and 6 below). The deeper weaknesses below persist.

- **Bands encode a baseline, not a truth.** A band around "the value we happened
  to get" drifts every time we intentionally change biology and re-baseline; over
  enough re-baselines it can ratchet away from physiology unnoticed. This still
  applies to `test_acceptance.py`'s *absolute-value* bands (resting level, Dolan
  peak). Its *relative* checks (knockout vs wild type) and *exact-zero* checks
  (resting quiescence) dodge it — they assert a structural relationship, not a
  remembered number — but literature-anchored magnitudes (candidate 3) remain the
  missing piece.

- **"Passes by construction" hides real coupling.** Because each new layer is
  normalised at rest to keep 5/5 green, 5/5 can't *detect* a new layer that's
  subtly wrong — it was built specifically not to perturb the thing 5/5 measures.

- **The clamp eats discriminating power.** The most-tested observable (cyt Ca²⁺)
  is the least sensitive one, so the strongest regression test covers the part of
  the model least able to reveal a problem.

- **Subsystem tests assert direction, not magnitude.** Many portfolio tests check
  "aspirin lowers TXA₂" or "clopidogrel lowers PAC-1" — directional sanity, which
  a badly-miscalibrated-but-still-monotonic model would still pass.

### What would actually help (candidate directions, not yet built)

These are options to discuss, not decisions:

1. **Per-subsystem golden traces** (not just the Ca²⁺ columns) — freeze
   representative outputs for secretion/thromboxane/integrin/inhibitory runs, so a
   change anywhere is caught where it shows.
2. **Tolerance-based, not byte-exact, regression** on those traces (e.g.
   `assert_allclose` with a documented rtol) to survive benign FP/library drift
   while still catching real moves.
3. **Validation tests with magnitude bands tied to literature**, not to our own
   last run — e.g. PAC-1 reaches a literature-anchored fraction, secretion
   half-time in a measured range — so the target is external.
4. **A "characterisation test" snapshot suite** that records many summary metrics
   per scenario into a committed table; diffs become a reviewable artifact rather
   than a pass/fail.
5. **Sensitivity/identifiability as a first-class check** — given the clamp, know
   which parameters each observable is actually sensitive to, and test there.
6. **Retire byte-identical — done (2026-06-26).** It was always a scaffold for a
   specific migration (#32), not a permanent fixture; `test_acceptance.py`'s
   behavioural bands replace it.

Most of this is still unimplemented. **The working contract now is:**
`test_acceptance.py` is the headline behavioural suite (the Dolan transient inside
it is the cheap regression invariant on the Ca²⁺ core); judge new biology on its
own subsystem assay; and treat every band re-baseline as a reviewed decision with
a recorded reason (the baseline-comment discipline in `test_regression.py`). The
byte-identical goldens and the all-or-nothing 5/5 gate are retired.

## Recommendation (current stance)

1. Keep the **Dolan transient** (in `test_acceptance.py`) as the cheap
   **regression invariant** on the Ca²⁺ core — explicit behavioural bands, not a
   5/5 count.
2. Present validation as the **portfolio table**, each subsystem vs its own assay.
3. State the **Ca²⁺ clamp** as a known limitation with a future-work hook.
4. Treat **regression-proofing the whole model as an explicit open work item** —
   `test_acceptance.py` resolves the brittleness of byte-identical goldens, but
   literature-anchored magnitudes and per-subsystem golden traces (candidates 1,
   3, 4 above) remain open.
