---
title: "Possible Typographical Error in 2008 Paper by Purvis et al."
date: 2026-05-08
author: "Steve Haigh"
keywords: IP3 receptor, Sneyd-Dufour, Purvis, Dolan, platelet calcium signalling, type-2 IP3R, rate constant, transcription error
---

# Possible Typographical Error in 2008 Paper by Purvis et al.

## Abstract
The widely-cited platelet calcium-signalling model of Purvis et al.
(2) records the IP3 receptor (type-2) closing rate constant as
`k₃ = 11 s⁻¹` in its Table 1.

The original Sneyd & Dufour (1) source paper that (2) cites for
these kinetics gives `k₃ = 0.11 s⁻¹` (Fig 4 caption, with body-text
confirmation that `φ₃ ≈ 0.1 s⁻¹` at low cytosolic Ca²⁺).

The Purvis value is 100× too large.

The error has propagated into at least one downstream
modelling paper (3) and into independent re-implementations of
the same kinetic scheme. The mechanical effect on transient
dose-response simulations is small, which is why the error has
been able to persist; the IP3R `s` sub-state population at
resting cytosolic Ca²⁺ is the most affected quantity, shifting
from ~3.8 % (with `k₃ = 11`) to ~0.04 % (with `k₃ = 0.11`).

A short dimensional argument fixes the value unambiguously to
`k₃ = 0.11 s⁻¹`. This note describes the error, its propagation,
and the dimensional reasoning that resolves it, to assist anyone
using the Sneyd-Dufour 2002 type-2 IP3R kinetic scheme via (2)
or downstream papers.

## 1. The rate constant in question

Sneyd & Dufour (1) introduce a 6-state Markov-chain model for
one subunit of the type-2 IP3 receptor, with states `R` (rest),
`O` (open), `A` (activated), `S` (shut), `I₁`, `I₂`
(inactivated). The constant `k₃` is the rate-law parameter for
the unimolecular `O → S` transition (Ca²⁺-independent
inactivation by IP3). It appears in exactly one place in the
rate-law formalism: the φ₃ function

$$\varphi_3(c) = \frac{k_3 \cdot L_5}{L_5 + c}$$

where `c = [Ca²⁺]_cyt` (in µM) and `L₅ = 54.7 µM` is a
Sneyd-Dufour-style equilibrium-constant parameter. The reverse
rate `k₋₃ = 29.8 s⁻¹` is a constant, not a φ-function.

## 2. Dimensional analysis: why `k₃` must have units `s⁻¹`

The state `o` is a count (or fraction); for the ODE term
`do/dt ⊃ −φ₃(c) · o` to give a dimensionally correct rate
(`count/s` or `1/s`), `φ₃` must have units of `s⁻¹`.

Working outward from the formula:

| Quantity | Units |
|---|---|
| `(L₅ + c)` (denominator) | µM |
| `φ₃` (required) | s⁻¹ |
| Numerator `k₃ · L₅` (forced) | s⁻¹ · µM |
| `L₅` | µM |
| Therefore `k₃` (forced) | **s⁻¹** |

So `k₃` must be in `s⁻¹`. There is no consistent reading of the
formula in which `k₃` carries a `µM⁻¹` factor.

### 2.1 Cross-check: are the other `k_i` annotations also typos?

Worth verifying, because in the Sneyd-Dufour (1) Fig 4 caption
the constants `k₁, k₂, k₃, k₄` are all annotated with the same
`s⁻¹·µM⁻¹` units. They are not all typos; only `k₃` is.

| Constant | Caption units | Where it appears | Numerator structure | Verdict |
|---|---|---|---|---|
| `k₁` | s⁻¹·µM⁻¹ | `(k₁L₁ + l₂) · c / [...]` | `(s⁻¹·µM⁻¹)·µM = s⁻¹`, then `× c (µM)` = `s⁻¹·µM`; ÷ µM denominator → `s⁻¹` | ✓ correct |
| `k₂` | s⁻¹·µM⁻¹ | `(k₂L₃ + l₄c) / [...]` | analogous to `k₁` | ✓ correct |
| `k₃` | s⁻¹·µM⁻¹ | `k₃ · L₅ / (L₅ + c)` | `(s⁻¹·µM⁻¹)·µM = s⁻¹`, but **no `c` factor** in numerator; ÷ µM denominator → `µM⁻¹·s⁻¹` | ✗ wrong dimensions for a rate |
| `k₄` | s⁻¹·µM⁻¹ | `(k₄L₅ + l₆) · c / (L₅ + c)` | analogous to `k₁` | ✓ correct |

The structural difference is that `φ₁`, `φ₂`, `φ₄` describe
transitions involving binding of a Ca²⁺ or IP3 ion, so the rate is
proportional to the binding-partner concentration and a `c`
factor appears in the numerator. `φ₃` describes a unimolecular
conformational change with no second molecule binding, so its
formula has no `c` factor and `k₃` cannot legitimately carry the
same units as its neighbours. The Fig 4 caption appears to have
copy-pasted a single unit annotation across all four `k_i` without
checking each formula individually.

## 3. The primary source

Two pieces of evidence in the original Sneyd & Dufour (1) paper
pin `k₃` to **0.11 s⁻¹** numerically:

1. **Fig 4 caption** (the parameter listing): "k₃ = 0.11 s⁻¹·µM⁻¹"
   (the `µM⁻¹` annotation is the typo identified in §2; the
   numerical value is correct).
2. **Body text** (Results and Discussion): "Our model predicts that
   Ca²⁺-independent inactivation of the IPR by IP3 has a time
   constant of around 10 s (φ₃ ≈ 0.1 s⁻¹)."

The body text is unambiguous about both the magnitude and the
units of `φ₃`. Plugging `φ₃ ≈ 0.1 s⁻¹` back into the formula at
low `c`:

$$\varphi_3(c \to 0) = \frac{k_3 \cdot L_5}{L_5} = k_3$$

so `k₃ ≈ 0.11 s⁻¹` — consistent with the value listed in the
caption when interpreted with the corrected units.

## 4. The transcription in Purvis 2008

Purvis et al. (2) Table 1 (page 4072) records the IP3R rate
constants for the closing transition as

> `k₃ = 11 s⁻¹·µM⁻¹    L₅ = 54.7 µM    k₋₃ = 29.8 s⁻¹`

i.e. with the **same `s⁻¹·µM⁻¹` unit annotation** as the
Sneyd & Dufour (1) Fig 4 caption (carrying forward the
documentary typo identified in §2), but with the value listed as
**11** rather than the original **0.11**.

The 17 other Sneyd-Dufour rate constants in (2) Table 1 (`k₁`,
`k₋₁`, `k₂`, `k₋₂`, `k₄`, `k₋₄`, `l₂`, `l₋₂`, `l₄`, `l₋₄`,
`l₆`, `l₋₆`, `L₁`, `L₃`, `L₅`) match the (1) Fig 4 values
exactly. **Only `k₃` is mis-transcribed**, and only in the
numerical value (by a clean factor of 100×).

## 5. Could the ×100 have been intentional?

It is worth considering whether Purvis et al. (2)
deliberately multiplied by 100. The most natural such hypothesis:
a careful transcriber reads the (1) caption's `s⁻¹·µM⁻¹`
annotation literally, recognises that `0.11 µM⁻¹·s⁻¹` cannot
plausibly be a unimolecular rate constant (the `O → S`
transition has no second molecule to bind), multiplies by a
"characteristic" Ca²⁺ concentration to convert into `s⁻¹`, and
records the result. With `[Ca²⁺] = 100 µM`,
`0.11 µM⁻¹·s⁻¹ × 100 µM = 11 s⁻¹` — exactly the (2) numerical
value.

This hypothesis is unlikely as the resulting rate would be `s⁻¹`
(not `s⁻¹·µM⁻¹`), and a transcriber would have updated the
unit annotation accordingly. The fact that the units are kept
unchanged rules out a deliberate unit-conversion mechanism.

A deliberate calibration adjustment from a primary source
would conventionally be flagged with an annotation. The (2)
table carries no such annotation.

The simplest hypothesis is therefore a digit-shift transcription
error (`0.11 → 11`).

## 6. Why the error has persisted

Four factors plausibly contribute to this 100× error surviving 18
years of citation:

1. **The primary (1) Fig 4 caption itself contains a
   unit-annotation typo** (`s⁻¹·µM⁻¹` instead of `s⁻¹` for `k₃`),
   making the value `0.11` look as though it might need to be
   multiplied by some characteristic [Ca²⁺] to give a rate.
2. **The (1) body-text correction is buried** in the Results
   and Discussion section. The unambiguous statement
   "`φ₃ ≈ 0.1 s⁻¹`" is part of a larger paragraph comparing
   model predictions to other published rate measurements; a
   reader pulling rate constants from the figure caption can
   easily miss the qualitative confirmation that follows.
3. **The (2) Table 1 is rotated 90° on the printed page.**
   A reader wanting to verify `k₃` has to either rotate the
   printed page, rotate the screen view, or bend to an awkward
   angle, then read across cramped formula cells in a small
   typeface. The cognitive friction makes both the original
   transcription error more likely and subsequent catches less
   likely.
4. **The macroscopic effect of the error is small.**
   A 100× error in `k₃` shifts the `S` sub-state population
   by a few percent at typical resting cytosolic Ca²⁺ levels
   and changes the channel open probability `Po⁴` by less than
   ~20 %. None of the dose-response observables plotted in (1)
   Figs 4–8, or in subsequent platelet-Ca²⁺ models that use the
   same kinetics, are sensitive enough to the `O → S` rate to
   expose the error without a deliberate sub-state audit.

## 7. Mechanical impact: what changes when `k₃` is corrected

The headline mechanical effect of correcting `k₃: 11 → 0.11` is
on the IP3R `S` sub-state population at the Markov-chain
equilibrium, evaluated at typical resting conditions
(`[IP3] = 50 nM`, `[Ca²⁺]_cyt = 100 nM`):

| Sub-state | Pre-fix (`k₃ = 11`) | Post-fix (`k₃ = 0.11`) |
|---|---|---|
| `n` | 44.7 % | 46.4 % |
| `o` | 10.3 % | 10.7 % |
| `a` | 14.2 % | 14.7 % |
| `i₁` | 13.5 % | 14.0 % |
| `i₂` | 13.6 % | 14.2 % |
| `s` | **3.78 %** | **0.04 %** |

The `s` population collapses by approximately two orders of
magnitude, finally matching the Dolan & Diamond (3) Table S1
representative initial-condition value of `S ≈ 0`. Other
sub-state populations move by less than two percentage points,
because the ladder's overall flux balance is dominated by the
`O ↔ A` equilibrium rather than the much slower `O ↔ S`
shutting transition.

Channel open probability `Po = (0.9·a + 0.1·o)⁴` changes from
`Po = 0.138`, `Po⁴ = 3.6×10⁻⁴` to `Po = 0.143`, `Po⁴ = 4.2×10⁻⁴`
— a ~17 % rise, because the shutting flux that previously
diverted material into `S` no longer fires, leaving slightly more
population in `O` and `A`. This is a small enough macroscopic
change that transient dose-response simulations are essentially
unaffected, which is consistent with the error having survived
without notice.

## 8. Verification

A reader wishing to verify the analysis above can do so without
running any simulation:

1. **Check the (1) Fig 4 caption value**: it lists
   `k₃ = 0.11` (with the `s⁻¹·µM⁻¹` annotation discussed in §2).
2. **Check the (1) body text** (Results and Discussion,
   page 2401): "`φ₃ ≈ 0.1 s⁻¹`".
3. **Check (2) Table 1** (page 4072, in the IP3R closing
   row): `k₃ = 11 s⁻¹`.
4. **Confirm the dimensional argument** in §2: the formula
   `φ₃ = k₃·L₅/(L₅+c)` is dimensionally consistent with `φ₃`
   in `s⁻¹` only when `k₃` has units `s⁻¹`.

## 9. Implications for downstream models

Anyone using the (1) type-2 IP3R kinetics by way of (2) Table 1
should re-check their `k₃` value against the primary source.
Models known to be in the citation chain include at least
(3) (which inherits the (2) Table 1 `k₃ = 11 s⁻¹`); a forward-
citation search ("cited by (2)" filtered for IP3R kinetics
implementations) would likely identify others.

## 10. Suggested correction

For models built on the (2) Table 1 or (3) Table 1 IP3R
kinetics:

```
- k₃ = 11 s⁻¹      (Purvis et al. (2) Table 1)
+ k₃ = 0.11 s⁻¹    (Sneyd & Dufour (1) Fig 4 + body text;
                    100× correction)
```

No other rate constants are affected.

## 11. AI assistance
An AI assistant (Anthropic Claude) was used during the broader
platelet calcium-signalling modelling work in which this
discrepancy was identified, including to cross-reference
rate-constant values from primary-source PDFs against
intermediate-source tables.

The same assistant was also used to help draft this note.
All factual claims, the dimensional analysis, and the
bibliographic citations have been verified by the author
against the cited primary sources.

## References

(1) Sneyd, J. & Dufour, J.-F. (2002). A dynamic model of the
type-2 inositol trisphosphate receptor. *Proceedings of the
National Academy of Sciences*, 99(4):2398–2403.
https://doi.org/10.1073/pnas.032281999

(2) Purvis, J. E., Chatterjee, M. S., Brass, L. F. & Diamond, S.
L. (2008). A molecular signaling model of platelet
phosphoinositide and calcium regulation during homeostasis and
P2Y1 activation. *Blood*, 112(10):4069–4079.
https://doi.org/10.1182/blood-2008-05-157883

(3) Dolan, A. T. & Diamond, S. L. (2014). Systems modeling of
Ca²⁺ homeostasis and mobilization in platelets mediated by IP3
and store-operated Ca²⁺ entry. *Biophysical Journal*,
106(9):2049–2060.
https://doi.org/10.1016/j.bpj.2014.03.028
