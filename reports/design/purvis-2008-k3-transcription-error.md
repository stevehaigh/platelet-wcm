---
title: "A 100× transcription error in Purvis 2008 Table 1 (IP3R k₃) propagated through the platelet-Ca²⁺ modelling literature"
date: 2026-05-08
author: platelet-wcm dissertation project (Steve Haigh)
---

# A 100× transcription error in Purvis 2008 Table 1 (IP3R k₃)

## TL;DR

The widely-cited Purvis *et al.* 2008 platelet-signalling model
(*Blood* 112:4069–4079) Table 1 lists the IP3R closing-rate
constant `k₃ = 11 s⁻¹`. The original Sneyd & Dufour 2002 *PNAS*
paper that Purvis cites as the source for these rate constants
gives **`k₃ = 0.11 s⁻¹`** in its Fig 4 caption, with body-text
confirmation that `φ₃ ≈ 0.1 s⁻¹` at low cytosolic Ca²⁺. Purvis's
value is **100× too large**. This error has propagated into at
least one downstream model (Dolan & Diamond 2014, *Biophys J*
106:2049–2060, the canonical platelet-Ca²⁺ Monte Carlo model) and
into the platelet-wcm v0.2 implementation reported in this
dissertation, before being identified by direct audit of the
Sneyd-Dufour primary source on 2026-05-08.

The mechanical effect of the error is small in transient
simulations (it primarily affects the steady-state IP3R `s`
sub-state population), but the methodological lesson is large: a
quantitative parameter widely re-used across two decades of
platelet-Ca²⁺ modelling was wrong by two orders of magnitude.

## The primary source

Sneyd J., Dufour J.-F. 2002. "A dynamic model of the type-2
inositol trisphosphate receptor." *PNAS* 99(4):2398–2403.

The model defines a 6-state Markov ladder for one IP3R subunit.
The transition `O → S` (open subunit → "shut" inactivated state)
has the rate function

    φ₃(c) = k₃ · L₅ / (L₅ + c)

where `c = [Ca²⁺]_cyt` and `L₅ = 54.7 µM`.

Two pieces of evidence pin `k₃` to **0.11 s⁻¹**:

1. **Fig 4 caption** (the parameter listing): "k₃ = 0.11 s⁻¹·µM⁻¹"
   (the `µM⁻¹` in the unit annotation is a documentation typo —
   the formula `k₃·L₅/(L₅+c)` is dimensionally consistent only
   when `k₃` has units of `s⁻¹`).
2. **Body text** (page 2401, `Results and Discussion`): "Our model
   predicts that Ca²⁺-independent inactivation of the IPR by IP3
   has a time constant of around 10 s (φ₃ ≈ 0.1 s⁻¹)".

Numerically, with `k₃ = 0.11`, `φ₃(c → 0) = k₃·L₅/L₅ = k₃ = 0.11
s⁻¹` ✓ — matches the body text exactly. With `k₃ = 11`,
`φ₃(c → 0) = 11 s⁻¹` — 100× too fast, time-constant ~91 ms instead
of the paper's 10 s.

### Dimensional analysis (why the caption's `µM⁻¹` annotation must be the typo)

The dimensional argument is worth making explicit, because the
caption's unit annotation looks plausible at first glance.

**Forward derivation — what units must `k₃` have?**

In the Sneyd-Dufour ODE, the `O → S` transition contributes:

    do/dt ⊃ −φ₃(c) · o
    ds/dt ⊃ +φ₃(c) · o

The state variable `o` is a fraction (or count) — dimensionless,
or `count`. For `do/dt` to have the correct units (`fraction/s`
or `count/s`), `φ₃` must have units of `s⁻¹`.

Then, walking the formula `φ₃ = k₃·L₅ / (L₅ + c)` outward:

| Quantity | Units |
|---|---|
| `(L₅ + c)` (denominator) | `µM` |
| `φ₃` (required) | `s⁻¹` |
| Therefore numerator `k₃·L₅` | `s⁻¹ · µM` |
| `L₅` | `µM` |
| Therefore **`k₃`** | **`s⁻¹`** |

So `k₃` must be in `s⁻¹`. The Fig 4 caption's `s⁻¹·µM⁻¹`
annotation is dimensionally inconsistent with the formula.

**Cross-check — are the other `k_i` annotations also typos?**

Worth verifying, because if `k₃`'s `µM⁻¹` is a typo, you might
suspect the same of `k₁`, `k₂`, `k₄`. They aren't.

| Constant | Caption units | Formula it appears in | Numerator structure | Verdict |
|---|---|---|---|---|
| `k₁` | s⁻¹·µM⁻¹ | `(k₁L₁ + l₂)·c / (L₁ + c(1 + L₁/L₃))` | `(s⁻¹·µM⁻¹)·µM = s⁻¹`, then `× c (µM)` = `s⁻¹·µM` | ✓ correct |
| `k₂` | s⁻¹·µM⁻¹ | `(k₂L₃ + l₄c) / (...)` | same as k₁ | ✓ correct |
| `k₃` | s⁻¹·µM⁻¹ | `k₃·L₅ / (L₅ + c)` | `(s⁻¹·µM⁻¹)·µM = s⁻¹`, but **no `c` factor** | ✗ wrong; needs k₃ in s⁻¹ |
| `k₄` | s⁻¹·µM⁻¹ | `(k₄L₅ + l₆)·c / (L₅ + c)` | same as k₁ | ✓ correct |

The structural difference: `φ₁`, `φ₂`, `φ₄` all have a `c`
factor in the numerator (because they describe transitions that
involve binding a Ca²⁺ ion or IP3 ion — the rate is proportional
to the binding-partner concentration). `φ₃` does not have a `c`
factor — the `O → S` transition is a unimolecular conformational
change, no second molecule binding. So `k₃` cannot legitimately
carry the same `s⁻¹·µM⁻¹` units as its neighbours; the caption
appears to have copy-pasted the same unit annotation across all
four `k_i` without checking each formula individually.

This is a well-defined, single-line typo in the published Sneyd-Dufour
2002 caption. The body text "`φ₃ ≈ 0.1 s⁻¹`" is the correct
specification.

## The transcription

Purvis J.E., Chatterjee M.S., Brass L.F., Diamond S.L. 2008.
"A molecular signaling model of platelet phosphoinositide and
calcium regulation during homeostasis and P2Y1 activation."
*Blood* 112(10):4069–4079.

Table 1 of Purvis 2008 (page 4072) lists the IP3R rate constants
for the closing transition as:

    k₃ = 11 s⁻¹     L₅ = 54.7 µM     k₋₃ = 29.8 s⁻¹

The 17 other K_IP3R rate constants (k₁, k₋₁, k₂, k₋₂, k₄, k₋₄,
l₂, l₋₂, l₄, l₋₄, l₆, l₋₆, L₁, L₃, L₅) match Sneyd-Dufour 2002
Fig 4 exactly. Only k₃ is mis-transcribed by 100×.

The most likely explanation is a numerical-typo cascade: the
Sneyd-Dufour Fig 4 caption has the value `0.11` in close adjacency
to the `L₅ = 54.7 µM` constant, and the unit annotation
`s⁻¹·µM⁻¹` is itself a typo (per the body text). Whoever
transcribed the value into Purvis Table 1 may have either dropped
the leading `0.` or implicitly multiplied by a characteristic
[Ca²⁺] (×100 µM gives 11). Without an erratum from Purvis *et al.*
or correspondence with the authors we can't be certain of the
exact mechanism, only the magnitude (factor 100×).

### Could the ×100 have been intentional?

It's worth considering charitably: a careful transcriber, noticing
the caption's `s⁻¹·µM⁻¹` annotation but reading the formula as
written, might have multiplied by a "characteristic" `[Ca²⁺]` to
convert into `s⁻¹` and reasoned the result was the rate they
should record in their own Table 1. With `[Ca²⁺] = 100 µM`,
`0.11 µM⁻¹·s⁻¹ × 100 µM = 11 s⁻¹` — exactly what Purvis records.

Three reasons this is unlikely to be the actual mechanism:

1. **`100 µM` is not a meaningful characteristic [Ca²⁺] for the
   o→s transition.** Cytosolic Ca²⁺ ranges over five orders of
   magnitude (~10 nM resting, ~1 µM transient peak, ~10 µM patch-
   pipette saturating). 100 µM is far above any cytosolic regime
   and would only correspond to luminal (DTS / SR) concentrations
   — but the o→s transition is gated by *cytosolic* Ca²⁺ via the
   `(L₅ + c)` denominator. Picking 100 µM here would be a category
   error, not a deliberate calibration. The reader would naturally
   reach for `[Ca²⁺] ≈ 0.1–1 µM` (resting/transient cyt range);
   that gives `0.011–0.11 s⁻¹`, not `11 s⁻¹`.
2. **Sneyd-Dufour's body text is unambiguous.** The "`φ₃ ≈ 0.1 s⁻¹`,
   time constant ~10 s" statement on page 2401 directly specifies
   the φ₃ rate. A deliberate reinterpretation by Purvis would have
   to also override that — i.e., the transcriber would have to
   read the caption value, multiply it by 100, *and* believe the
   body-text estimate was wrong. That's a much stronger claim
   than just "I assumed the caption units."
3. **No footnote or explicit annotation in Purvis Table 1.** A
   deliberate calibration adjustment from a primary source would
   conventionally be flagged with a "rate constants converted from
   Sneyd-Dufour 2002 by multiplication by [...]" footnote.
   Purvis Table 1 carries no such annotation; the row is
   presented as a verbatim transcription. The 17 other K_IP3R rate
   constants match Sneyd-Dufour exactly with no conversion factor,
   and Purvis's own `k₋₄` row correctly transcribes Sneyd-Dufour's
   value (0.54 s⁻¹) despite Sneyd-Dufour's caption labelling
   `k₋₄` units inconsistently with the formula too — so Purvis
   was reading values, not "fixing" units.

The simplest hypothesis remains a digit-shift transcription error
(`0.11 → 11`) — exactly the kind of typo that survives proof-reading
because the surrounding rate constants (4, 11.4, 4707, etc.) are
all in the same order of magnitude and `11 s⁻¹` doesn't look
unreasonable on its face.

## The propagation

| Year | Reference | k₃ value | Source |
|---|---|---|---|
| 2002 | Sneyd & Dufour, *PNAS* 99:2398 | **0.11 s⁻¹** | primary (correct) |
| 2008 | Purvis *et al.*, *Blood* 112:4069 Table 1 | **11 s⁻¹** | transcription error |
| 2014 | Dolan & Diamond, *Biophys J* 106:2049 Table 1 | **same as Purvis** | inherited via Purvis |
| 2026 | platelet-wcm v0.2 | **same as Purvis** | inherited via Dolan |

Dolan & Diamond 2014's Table 1 lists `IP3R closing` with rate law
`[IP3R_o] · k₃·L₅/(L₅ + [Ca²⁺]_cyt)` and the same numerical
values as Purvis 2008. The platelet-wcm v0.2 model under audit
here uses the Dolan numbers via the
`reports/data/calcium-data-provenance.md` 2026-04-23 cross-check.

So the error appears to have been "frozen in" at the Purvis 2008
step and propagated through the literature without being
re-checked against the Sneyd-Dufour primary source for at least
18 years.

## Mechanical impact in our model

Quantitative effect of correcting `k₃: 11 → 0.11`:

### Markov-chain equilibrium at IP3 = 50 nM, ca_cyt = 100 nM (resting)

| Sub-state | Dolan Table S1 (filtered IC) | Pre-fix (k₃ = 11) | Post-fix (k₃ = 0.11) |
|---|---|---|---|
| n  | 60.96% | 44.65% | 46.39% |
| o  | 19.67% | 10.27% | 10.67% |
| a  |  4.90% | 14.18% | 14.73% |
| i1 | 12.58% | 13.49% | 14.02% |
| i2 |  1.88% | 13.63% | 14.16% |
| **s**  | **0.00%** |  **3.78%** |  **0.04%** ← matches Dolan |

The headline mechanical effect: the `s` sub-state population
collapses from ~3.78 % to ~0.04 %, finally matching Dolan Table
S1's reported `s ≈ 0`. Other sub-states barely move because the
ladder's overall flux balance is dominated by the
o ↔ a equilibrium (φ₄ = 4707 s⁻¹ at saturating Ca²⁺), not by
the much-slower o ↔ s shutting transition.

### Po⁴ at the same conditions

| | Po | Po⁴ |
|---|---|---|
| Dolan Table S1 IC | 0.0638 | 1.65×10⁻⁵ |
| Pre-fix (k₃ = 11) | 0.1378 | 3.61×10⁻⁴ (22× Dolan) |
| Post-fix (k₃ = 0.11) | 0.1432 | 4.21×10⁻⁴ (25× Dolan) |

Po⁴ actually rises slightly with the fix (because the shutting
flux that previously diverted material into `s` no longer fires,
leaving more population in `o` and `a`). The 22-25× ratio to
Dolan's filtered-IC Po⁴ is unchanged in order of magnitude — that
gap is *not* attributable to this transcription error and is now
firmly attributed to **Dolan's Monte Carlo filtering procedure**
(which required only macro-concentration stationarity and did not
constrain IP3R sub-state populations to the equilibrium of the
Sneyd-Dufour rate laws), per the 2026-05-08 audit.

### Phase 3 acceptance

| | Pre-fix | Post-fix |
|---|---|---|
| +Ca_ex peak | 392 nM ✓ | 393 nM ✓ |
| −Ca_ex peak | 325 nM ✓ | 325 nM ✓ |
| SOCE differential | 67 nM ✗ | 68 nM ✗ |
| Total | **4/5** | **4/5** |

No regression and no improvement at the macroscopic acceptance
level. The fix is mathematically necessary (matches primary
source) but practically minor in transient amplitudes.

## Why the error persisted

Three contributing factors made this 100× error survive 18 years
of citation:

1. **The primary Sneyd-Dufour 2002 Fig 4 caption itself contains
   a unit-annotation typo** (`s⁻¹·µM⁻¹` instead of `s⁻¹` for
   `k₃`). This makes the listed value `0.11` look like it might
   need to be multiplied by some characteristic [Ca²⁺] to get a
   meaningful rate, plausibly motivating an "effective `k₃`" of
   `0.11 × 100 µM = 11`.
2. **The body text correction (`φ₃ ≈ 0.1 s⁻¹`) is buried** in the
   `Results and Discussion` section as part of a comparison of
   model predictions to other published rate measurements. Easy
   to miss when transcribing rate constants from the figure
   caption.
3. **The macroscopic effect is small** — the IP3R `o ↔ s`
   transition is a minor side-branch of the 6-state ladder; the
   dominant `o ↔ a` rates are 4-5 orders of magnitude faster and
   carry essentially all of the gating signal. A 100× error in
   `k₃` shifts the `s` sub-state population by a few percent but
   doesn't materially change Po, peak Ca²⁺, or any of the
   dose-response observables that get plotted in Sneyd-Dufour
   Figs 4–8 or in subsequent platelet-Ca²⁺ models. So the error
   *can survive* without being noticed in any of the downstream
   validation work.

## Implications for the broader platelet-Ca²⁺ modelling community

Anyone who uses the Sneyd-Dufour 2002 type-2 IP3R kinetics via
Purvis 2008 / Dolan 2014 / platelet-wcm v0.2 should re-check
their `k₃` value. The fix is trivial — change one constant from
`11` to `0.11` — and the impact on macroscopic predictions is
small. The methodological value is that this is the kind of
parameter error that's worth catching.

Models that may be affected (not exhaustive — derived from the
direct citation chain in `reports/data/calcium-data-provenance.md`
and forward-citation searches not yet performed):

- Purvis *et al.* 2008 — original transcription site.
- Dolan & Diamond 2014 — inherited via Purvis Table 1.
- Sveshnikova *et al.* 2015 (*Nature Comm. Biol.*), 2025 (review),
  if they use Purvis or Dolan as their kinetics source.
- platelet-wcm v0.2 (this dissertation) — corrected 2026-05-08
  in commit `5c70d6df`.

A targeted forward-citation search (Google Scholar "cited by"
Purvis 2008 / Dolan 2014, filtered for IP3R kinetics
implementations) would identify any other affected models.

## Bottom line

> *Purvis et al. 2008 Table 1 has `k₃ = 11 s⁻¹` for the IP3R
> `o → s` transition. The Sneyd-Dufour 2002 primary source has
> `k₃ = 0.11 s⁻¹`. Purvis's value is 100× too large.
> Mechanically the error has small effects on transient
> dose-response predictions (which is why it was missed) but
> shifts the resting `s` sub-state population by ~100×.
> Identified by direct audit of the Sneyd-Dufour 2002 PDF on
> 2026-05-08; corrected in platelet-wcm v0.2 in commit
> `5c70d6df`.*

## References

- Sneyd J., Dufour J.-F. 2002. A dynamic model of the type-2
  inositol trisphosphate receptor. *PNAS* 99(4):2398–2403.
- Purvis J.E., Chatterjee M.S., Brass L.F., Diamond S.L. 2008.
  A molecular signaling model of platelet phosphoinositide and
  calcium regulation during homeostasis and P2Y1 activation.
  *Blood* 112(10):4069–4079.
- Dolan A.T., Diamond S.L. 2014. Systems modeling of Ca²⁺
  homeostasis and mobilization in platelets mediated by IP3 and
  store-operated Ca²⁺ entry. *Biophys J* 106(9):2049–2060.

## Audit trail in this repository

- Discovery: lab-book entry at
  `reports/lab-books/lab-book-2026-05-08-sneyd-dufour-audit.md`.
- Corrected source code:
  `reconstruction/platelet/dataclasses/process/calcium_signalling.py`
  (`K_IP3R['k3'] = 0.11`).
- Updated provenance documentation:
  `reports/data/calcium-data-provenance.md` (IP3R section).
- Numerical evidence: equilibrium-distribution comparison at
  `reports/data/ip3r-equilibrium-2026-05-07.json` (pre-fix) and
  by re-running `runscripts/manual/checkIP3REquilibrium.py`
  (post-fix; reproducible).
- Commit: `5c70d6df`, 2026-05-08.
