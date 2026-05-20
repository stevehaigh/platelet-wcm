---
title: "Lab book — 2026-05-15: Mike's recent-platelet collection, Sveshnikova/Lenoci/Mazet state comparison, PIP2 transcription fix"
---

# Lab book — 2026-05-15: paper triage + PIP2 fix + follow-up issues

## Context

Pre-Mike-chat (Thursday 2026-05-14) triage of the "Recent Platelet"
Zotero collection. Question: do the recent papers validate, invalidate,
or open up the v0.4.1 model, and what should we read for parameters?

## Papers reviewed

| Paper | Verdict for v0.4.1 |
|---|---|
| **Sveshnikova & Panteleev 2025** *Hämostaseologie* 45:49 (review) | **High value.** Sibling ODE model; same compartment topology (cytosol/DTS/mito/PM). Methodologically the canonical write-up of the Russian-group lineage. |
| **Ajanel et al. 2025** *Circ. Res.* 137:474 (MCU / ITAM) | Low priority for GPCR-only scope. MCU dispensable for PAR/P2Y-driven mitochondrial Ca²⁺; matters for GPVI/CLEC-2 and procoagulant transition. Cite if/when adding GPVI or procoagulant state. |
| **Balabin et al. 2024** *Biomed. Khim.* 70:394 (personalisation) | Abstract-only; key takeaway is empirical SERCA/IP3R proteomic correlation across 15 healthy donors — they compensate. Useful sensitivity-analysis citation. |
| **Starodubtseva et al. 2024** *PLOS ONE* 19:e0308679 (point-source ADP) | Population-level phenomenological model; not directly comparable to a single-cell ODE. Notes log-normal ADP sensitivity threshold (median 10⁻⁷-10⁻⁶ M) for future P2Y1/P2Y12 work. |
| **Cardillo & Barakat 2025** *Biomech. Model. Mechanobiol.* 24:1465 (2D plug) | Thrombus scale — different model class; not relevant to single-cell biochemistry. |
| **Anjum et al. 2025** *Blood* 145:1568 (aging platelets) | Reminds us a single Ca²⁺ trace is implicitly a "young/reticulated" answer. One-line dissertation caveat. |

## Initial-state comparison vs Sveshnikova / Lenoci / Mazet

Converting Sveshnikova/Lenoci concentrations to copies (at their volumes)
revealed that **their protein concentrations don't always correspond to
physically plausible copy numbers** — they're effective ODE
concentrations tuned for kinetics. The repo (whole-cell tradition from
wcEcoli) tracks real copies, which should match proteomics. Headline
findings:

| Species | Repo | Sveshnikova 2025 (3 fL) | Lenoci 2011 | Verdict |
|---|---|---|---|---|
| Ca²⁺ cyt rest | 100 nM | 13 nM | 55 nM | Repo defensible (Dolan canonical) |
| Ca²⁺ DTS | 250 µM | 1000 µM | — | Repo at low end of 300-1000 µM range |
| IP3 rest | 50 nM | 50 nM | 5 nM | Match ✓ |
| **PIP2** | **31 µM (count 1.12e5)** | 200 µM | 200 µM | **Bug: Mazet 2020 says 1.12 × 10⁶, repo had 1.12 × 10⁵** — see fix below |
| PLCβ total | 1 000 | 60 nM ≈ 108 copies | 3 nM ≈ 5 copies | Repo at low end; Burkhart 2012 ~5000. See issue. |
| Gαq total | 5 000 | — | 4.3 nM ≈ 8 copies | Match Mazet 2020 ✓ |
| PAR1 | 2 500 | 6 nM ≈ 11 copies | 6 nM | Repo matches Brass 1992 / Burkhart 2012 proteomics ✓ |
| P2Y1 | 150 | — | — | Match Coller 1995 ✓ |

## PIP2 transcription fix (this session)

`internal_state.py` had `PIP2 = 112_000` with citation "Mazet 2020
Supp Table S6 cited platelet PIP2 ~ 1.12 × 10⁵". Mazet 2020 main text
actually says *"a value of 1.12 × 10⁶ per cell"*. Off by 10×.

**Behaviour-neutral fix:**

- `internal_state.py` PIP2 count: 112_000 → 1_120_000 (now ≈310 µM at
  6 fL, consistent with Sveshnikova/Lenoci 200 µM).
- `calcium_signalling.py` `K_PLCB['k_cat']`: 2.26e-7 → 2.26e-8
  (preserves `k_cat × PIP2` product, so IP3 production rate at rest
  stays at 3.62 PIP2/s and Dolan match is preserved).
- `test_regression.py` dry-mass baseline: 245.97 fg → 247.72 fg (the
  1.008e6 extra PIP2 × 1.74e-6 fg = 1.75 fg extra lipid mass).

All 21 platelet tests pass; Phase 3 Dolan 5/5 still green.

## Follow-up issues opened

- **#36 — Stochastic simulation of low-copy GPCR intermediates.**
  Sveshnikova 2015 argues PLCβ-Gq-GTP / PAR1-active sit at single-digit
  counts → stochastic needed for subpopulation effects. Cheap first cut:
  Poisson draws around deterministic means for active sub-states.
- **#37 — Resting cytosolic Ca²⁺ baseline.** Repo 100 nM vs Sveshnikova
  13 nM. Repo is defensible (Dolan / Heemskerk canonical); citation
  needs tightening in `internal_state.py`.
- **#38 — DTS free Ca²⁺ baseline.** Repo 250 µM vs Sveshnikova 1 mM.
  Check whether new CALR + HSP90B1 + BiP + CREC buffering captures the
  difference (total bound + free) and decide whether to bump.
- **#39 — PLCβ copy number.** Current 1000 is documented as a stub;
  Burkhart 2012 reports PLCβ-3 ~ 5000. Either bump (and rescale `k_act`
  / `k_cat`) or pin down a tighter citation for 1000.

## What to read next

For PI-cycle parameters the most direct source is now in hand:

- **Mazet, Tindall, Gibbins & Fry 2020** *Sci. Rep.* 10:13244 — the
  framework citation for the PI cycle in `calcium_signalling.py`. Worth
  re-reading their Supp. Table S6 + S7 to double-check the other PI
  cycle copies (PIP, PI4P, PI, OCRL1, PIP5K) before adding the full
  PI/PI4P chain (out of v0.4 scope).
- **Lenoci et al. 2011** *Mol. BioSyst.* 7:1129 — full rate-constant
  table for PAR1 → Gq → PLCβ. Most useful for cross-checking
  the `K_PAR1` / `K_GQ` / `K_PLCB` "calibrated" values, all of which
  are currently uncited in the repo.
- **Sveshnikova et al. 2015** *Mol. BioSyst.* 11:1052 — base 27-ODE
  PAR1 model. Useful as a sibling ODE for sanity-checking dynamic
  responses, although parameter tables are in ESI.

## Status snapshot

- Model: v0.4.1 (unchanged); regression baseline bumped to 247.72 fg.
- Tests: 21/21 pass; Phase 3 Dolan 5/5.
- 8 issues closed across this session lineage; 4 new issues opened
  (above).
