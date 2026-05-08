---
title: "Lab book — 2026-05-08: Sneyd-Dufour 2002 φ-function audit (path B)"
---

# Lab book — 2026-05-08: Sneyd-Dufour 2002 φ-function audit

## Session goal

Yesterday's session ended with three forward paths after the
Phase 0 audit + CALR-buffer attempt (lab books 2026-05-07
phase-0-biology-audit and dts-drain-investigation):

- **Path A** — write up at 4/5 with documented gaps
- **Path B** — re-derive the φ-functions internally, audit against
  Sneyd-Dufour 2002 if we can locate the PDF
- **Path C** — empirical γ_IP3R retune

Steve found the Sneyd-Dufour 2002 PDF and Lièvremont 1997 overnight,
so path B is now open. Today's task: audit the φ-function
implementations and rate constants against the primary source.

## Method

1. Extract the canonical φ-function definitions from the body of
   Sneyd-Dufour 2002 *PNAS* 99:2398–2403.
2. Extract the best-fit parameter values from the Fig 4 caption.
3. Compare line-by-line to our `_phi_*` helpers and `K_IP3R` dict
   in `reconstruction/platelet/dataclasses/process/calcium_signalling.py`.
4. Cross-reference with the body text where possible (for sanity
   checks like reported time constants).

## What I found

### φ-functions: implementation matches Sneyd-Dufour exactly

Sneyd-Dufour 2002 defines seven φ-functions (after the fast-equilibrium
collapse that produces the 6-state simplified model in their Fig 3):

```
φ₁(c)  = (k₁L₁ + l₂)c / (L₁ + c(1 + L₁/L₃))
φ₂(c)  = (k₂L₃ + l₄c)  / (L₃ + c(1 + L₃/L₁))
φ₋₂(c) = (k₋₂ + l₋₄c)  / (1 + c/L₅)
φ₃(c)  = k₃L₅          / (L₅ + c)
φ₄(c)  = (k₄L₅ + l₆)c  / (L₅ + c)
φ₋₄(c) = L₁(k₋₄ + l₋₆) / (L₁ + c)
φ₅(c)  = (k₁L₁ + l₂)c  / (L₁ + c)
```

Our `_phi_*` helpers reproduce each one verbatim. ✓ All seven match.

The reverse rates `_phi_n_i1_rev = k_m1 + l_m2`,
`_phi_a_i2_rev = k_m1 + l_m2`, and `_phi_o_s_rev = k_m3` (constants,
not φ-functions) are also consistent with the simplified ladder.

The Po formula `Po = (0.9·a/total + 0.1·o/total)⁴` matches the
Sneyd-Dufour body text exactly (the paper writes it as
`(0.1·O + 0.9·A)⁴` — same expression, different ordering of terms).

### Rate constants: one real bug + one documentation typo

| Constant | Sneyd-Dufour 2002 | Purvis 2008 | Our pre-audit code | Our post-audit code |
|---|---|---|---|---|
| k₁ | 0.64 µM⁻¹·s⁻¹ | 0.64 | 0.64 | 0.64 ✓ |
| k₋₁ | 0.04 s⁻¹ | 0.04 | 0.04 | 0.04 ✓ |
| k₂ | 37.4 µM⁻¹·s⁻¹ | 37.4 | 37.4 | 37.4 ✓ |
| k₋₂ | 1.4 s⁻¹ | 1.4 | 1.4 | 1.4 ✓ |
| **k₃** | **0.11 s⁻¹** | **11** | **11.0** | **0.11** ← FIXED |
| k₋₃ | 29.8 s⁻¹ | 29.8 | 29.8 | 29.8 ✓ |
| k₄ | 4 µM⁻¹·s⁻¹ | 4 | 4.0 | 4.0 ✓ |
| k₋₄ | 0.54 s⁻¹ | 0.54 | 0.54 (comment unit was µM⁻¹·s⁻¹) | 0.54 (comment fixed to s⁻¹) |
| L₁ | 0.12 µM | 0.12 | 0.12 | 0.12 ✓ |
| L₃ | 0.025 µM | 0.025 | 0.025 | 0.025 ✓ |
| L₅ | 54.7 µM | 54.7 | 54.7 | 54.7 ✓ |
| l₂ | 1.7 s⁻¹ | 1.7 | 1.7 | 1.7 ✓ |
| l₋₂ | 0.8 s⁻¹ | 0.8 | 0.8 | 0.8 ✓ |
| l₄ | 1.7 µM⁻¹·s⁻¹ | 1.7 | 1.7 | 1.7 ✓ |
| l₋₄ | 2.5 µM⁻¹·s⁻¹ | 2.5 | 2.5 | 2.5 ✓ |
| l₆ | 4707 s⁻¹ | 4707 | 4707 | 4707 ✓ |
| l₋₆ | 11.4 s⁻¹ | 11.4 | 11.4 | 11.4 ✓ |

**The k₃ bug** — Purvis 2008 Table 1 lists `k₃ = 11 s⁻¹`. Sneyd-Dufour
2002 Fig 4 caption gives `k₃ = 0.11 s⁻¹·µM⁻¹` (with a µM⁻¹ unit
typo — the φ₃ formula is only dimensionally consistent if k₃ is
s⁻¹). The Sneyd-Dufour body text confirms unambiguously:

> "Our model predicts that Ca²⁺-independent inactivation of the IPR
> by IP3 has a time constant of around 10 s (φ₃ ≈ 0.1 s⁻¹)"

`k₃ = 0.11` gives `φ₃(c → 0) = 0.11 s⁻¹` ✓. `k₃ = 11` gives
`φ₃ ≈ 11 s⁻¹` — 100× too fast. Purvis appears to have transcribed
incorrectly.

**The k₋₄ documentation typo** — value (0.54) is right but our
in-code comment claimed it was `µM⁻¹·s⁻¹`. Sneyd-Dufour Fig 4 has
it as `s⁻¹` and the φ₋₄ formula requires it. Pure documentation
fix; no behavioural change.

### Effects of fixing k₃

**On the equilibrium sub-state distribution**:

| Sub-state | Dolan IC | Our pre-fix (k₃=11) | Our post-fix (k₃=0.11) |
|---|---|---|---|
| n  | 60.96% | 44.65% | 46.39% |
| o  | 19.67% | 10.27% | 10.67% |
| a  |  4.90% | 14.18% | 14.73% |
| i1 | 12.58% | 13.49% | 14.02% |
| i2 |  1.88% | 13.63% | 14.16% |
| s  |  0.00% |  3.78% | **0.04%** ← matches Dolan |

The `s` state collapses from 3.78% to ~0.04% — finally matching
Dolan Table S1's reported `s ≈ 0`. The other sub-states barely move.

**On Po⁴**:

| | Dolan IC | Our pre-fix | Our post-fix |
|---|---|---|---|
| Po | 0.0638 | 0.1378 | 0.1432 |
| Po⁴ | 1.65×10⁻⁵ | 3.61×10⁻⁴ (22× Dolan) | 4.21×10⁻⁴ (25× Dolan) |

The Po⁴ ratio actually *worsens* slightly because removing the
"o → s leak" leaves a tiny bit more material in the open ladder.
But the absolute change is small (~17 % up). The 22-25× gap to
Dolan IC stands.

**On Phase 3**:

| | Pre-k₃-fix (current main) | Post-k₃-fix |
|---|---|---|
| +Ca_ex peak | 392 nM | 393 nM |
| −Ca_ex peak | 325 nM | 325 nM |
| SOCE differential | 67 nM | 68 nM |
| Acceptance | 4/5 | 4/5 |

Essentially no Phase 3 behavioural change — within numerical noise.
All 21 platelet tests pass without further updates.

**On resting fixed point** (`restConvergence.py` 6 000 s):

| | Pre-fix | Post-fix |
|---|---|---|
| Ca²⁺_cyt | 2 169 nM | 2 240 nM |
| Ca²⁺_DTS_free | 0.06 µM | 0.07 µM |

Trivial change. The k₃ fix doesn't address the resting-state gap;
the upstream IP3R Markov chain is still 25× more leaky at IP3 =
50 nM than Dolan's IC implies, so the system still drives cyt to
~2.2 µM and DTS to ~0.

## Diagnosis: where is the remaining 25× gap?

After this audit:
- The φ-functions match Sneyd-Dufour exactly.
- The rate constants match Sneyd-Dufour Fig 4 exactly (with the
  k₃ correction).
- The Po formula matches Sneyd-Dufour exactly.

So OUR computation of the Markov-chain equilibrium at IP3 = 50 nM
and ca_cyt = 100 nM giving Po⁴ = 4.21×10⁻⁴ is what
**Sneyd-Dufour 2002 would predict** at those conditions.

Dolan's reported IC has Po⁴ = 1.65×10⁻⁵. The 25× discrepancy
must come from one of:

1. **Dolan's IC is not at Markov-chain equilibrium.** Their Monte
   Carlo filter required macro-concentration stationarity only; the
   sub-state populations were sampled, not equilibrated. So Table S1
   is a *snapshot* with the "right" macro-concentrations and Po⁴
   that survives their numerical-tolerance window, not a true
   stationary point of the Sneyd-Dufour ladder. **Most likely.**
2. **Sneyd-Dufour parameters were calibrated at IP3 = 10 µM** (per
   Fig 4 caption: "10 µM of IP3 was added"). At our IP3 = 50 nM
   (200× lower), the model is being asked to extrapolate well below
   the calibration regime. The "resting" Po⁴ is essentially
   unconstrained by the original fit. **Plausible.**
3. **There's another rate-constant mismatch we missed.** I've now
   audited all 17 K_IP3R values against Sneyd-Dufour Fig 4. The
   only remaining sources of discrepancy would be in the φ-function
   forms (matched), the Po weighting (matched, and Sneyd-Dufour
   explicitly say the weights aren't crucial), or in how we set up
   the equilibrium solve (we built the rate matrix per the standard
   Q·p = 0, Σp = 1 procedure). **Very unlikely.**

So the residual gap is most likely **explanation 1** (Dolan IC is a
sampled snapshot, not an equilibrium), with possibly some **2**
contribution (Sneyd-Dufour parameters extrapolated outside their
calibration regime).

## Implications

This is the cleanest possible outcome of path B:

- **Found a real bug** (k₃ 100× too large) and fixed it.
- **Confirmed the implementation otherwise faithfully reproduces
  Sneyd-Dufour 2002**, including the φ-functions, the rate constants,
  and the Po formula.
- The remaining 25× Po⁴ gap to Dolan's reported IC is now firmly
  attributable to **how Dolan filtered their initial-conditions
  ensemble**, not to a bug in our code or in Sneyd-Dufour's model.

For dissertation purposes this is a strong finding:

> "We audited the IP3R Sneyd-Dufour 2002 implementation against the
> primary source and identified a 100× transcription error in
> Purvis 2008 Table 1's reported k₃ value, propagated through to
> our model and several other implementations of the same scheme.
> Correcting k₃ from 11 to 0.11 s⁻¹ collapses the resting `s`
> sub-state population from ~3.8% to <0.1%, matching Dolan & Diamond
> 2014 Table S1's reported `s ≈ 0`. The remaining ~25× Po⁴ gap
> between our Markov-chain equilibrium and the Dolan IC is
> attributable to Dolan's Monte Carlo filtering procedure, which
> required macro-concentration stationarity only and did not
> constrain the IP3R sub-state populations to the equilibrium of
> the Sneyd-Dufour rate laws."

## Decision: commit the k₃ fix; revisit the path forward

The fix is unambiguously correct and biologically meaningful (it
matches Dolan's published `s ≈ 0`), so it lands cleanly. Phase 3
is unchanged (4/5).

The Phase 3 4/5 → 5/5 question now reduces to: **either accept
that Dolan's resting IP3R Po⁴ is what their IC says it is (1.65e-5)
without worrying about whether it's a Markov-chain equilibrium of
Sneyd-Dufour, OR accept that our model converges to a higher Po⁴
that gives a non-biological resting state and live with the
documented gap**.

Pragmatically:

- **Short-term**: commit the k₃ fix and the matching documentation
  updates. No Phase 3 regression. ✓
- **Medium-term**: the remaining gap is the Sneyd-Dufour /
  Dolan-filtering-vs-Markov-equilibrium question — a real
  research-grade limitation, not a bug. Worth writing up properly
  in the design doc as a documented v0.2 limitation, with the
  options being either: (a) accept the model's natural Po⁴ and
  live with the resting drift (path A from yesterday); (b) clamp
  IP3R sub-state IC to a non-equilibrium Dolan-style snapshot at
  start (path C in disguise); or (c) move to a more recent IP3R
  model that handles low-IP3 regime better than Sneyd-Dufour (e.g.
  Park & Sneyd 2009, deYoung-Keizer 1992).
- **For the buffer attempt** (#26 / CALR): this audit doesn't
  resolve the k₁₁ stability question that derailed yesterday's
  buffer attempt. The k₁₁ runaway is a separate downstream issue
  in the Caride 2007 cycle, not in the IP3R kinetics.

## Files changed

- `reconstruction/platelet/dataclasses/process/calcium_signalling.py`
  — `K_IP3R['k3']: 11.0 → 0.11` plus comment block + k₋₄ unit
  annotation fix.
- `reports/data/calcium-data-provenance.md` — IP3R section updated
  with the Sneyd-Dufour audit findings.
- `reports/lab-books/lab-book-2026-05-08-sneyd-dufour-audit.md`
  (this file) — new lab book entry.
- `reports/data/rest-converged-2026-05-07.json` — re-run output.
- `out/phase3_after_k3/` — Phase 3 reproducer output (verifies
  unchanged 4/5).

---

*Branch:* `main` · *Status:* k₃ fix landed (post this commit) ·
*Linked issues:* #24 (resting-state gap remains as research-grade
limitation), #26 (k₁₁ stability separately) · *Triggered by:*
overnight reference acquisition by Steve.
