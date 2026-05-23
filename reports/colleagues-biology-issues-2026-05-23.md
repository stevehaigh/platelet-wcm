---
title: "Biology issues to share with colleagues"
date: 2026-05-23
---

# Open biology issues — `stevehaigh/platelet-wcm`

This is a curated, grouped view of the open issues that are about
**model biology** rather than engineering. Generated 2026-05-23 from
the `biology` GitHub label. Live filter (always up to date):

<https://github.com/stevehaigh/platelet-wcm/issues?q=is%3Aopen+label%3Abiology>

Engineering / infrastructure issues are tagged `tech` and filtered
out here. If you want everything:

<https://github.com/stevehaigh/platelet-wcm/issues>

## Quick orientation

The model is currently at **v0.5.x**: GPCR cascade (PAR1 / PAR4 /
P2Y1 / P2X1) → Gαq → PLCβ → PI cycle → IP3 → IP3R/SERCA/DTS-release →
cyt Ca²⁺ → PMCA / SOCE / NCX / MCU / NCLX / CaM / multi-buffer DTS.
Validated against Dolan & Diamond 2014 Fig. 4 (Phase 3 acceptance
5/5). Dose-response surface across ADP × thrombin reproduced (issue
#45, closed). Dissertation submission ~August 2026.

The biological items below are roughly grouped by where input from a
collaborator would be most useful. The grouping is mine; the issues
themselves are at the URLs.

## Calibration & baselines — pointed questions where biology input matters most

These are open calibration questions where literature values disagree
and the model's chosen value affects everything downstream. **Highest
value group to discuss.**

| # | Question | Notes |
|---|---|---|
| [#37](https://github.com/stevehaigh/platelet-wcm/issues/37) | Is resting cytosolic Ca²⁺ 100 nM (Purvis / Dolan) or 13 nM (Sveshnikova)? | The model currently sits at ~100 nM. ~10× difference between the two literatures. |
| [#38](https://github.com/stevehaigh/platelet-wcm/issues/38) | Is DTS *free* Ca²⁺ 250 µM (Dolan / Purvis) or 1 mM (Sveshnikova)? | Affects IP3R driving force; the model has 250 µM at rest. |
| [#39](https://github.com/stevehaigh/platelet-wcm/issues/39) | PLCβ copy number — 1000 (placeholder) vs ~5000 (Burkhart 2012 proteomics) | Affects PI cycle gain. |

## Receptors — what to add next

The model has PAR1, PAR4, P2Y1, P2X1. Three obvious next-receptor
candidates; each is a real Gq / IP3R pathway in human platelets.

| # | Receptor | Pathway | Notes |
|---|---|---|---|
| [#10](https://github.com/stevehaigh/platelet-wcm/issues/10) | **P2Y₁₂** | ADP / Gi-coupled / inhibitory amplification | The canonical "missing receptor" alongside P2Y1. Clinically central (clopidogrel target). |
| [#33](https://github.com/stevehaigh/platelet-wcm/issues/33) | **TP** (thromboxane A2) | Gq positive-feedback arm | TXA2 release reinforces activation; aspirin's mechanism. |
| [#34](https://github.com/stevehaigh/platelet-wcm/issues/34) | **GP6 / GPVI** | Collagen → Syk → PLCγ2 → IP3 | Different cascade family (ITAM not Gq); collagen-mediated activation. |

## Granule release — the next platelet-biology arc beyond Ca²⁺

Once the Ca²⁺ scope is solid, granule secretion is the natural next
biology layer. Dense granules and α-granules each carry distinct
cargos with their own clinical readouts. **Currently scoped but not
implemented** in the v0.x line.

| # | Piece | Notes |
|---|---|---|
| [#5](https://github.com/stevehaigh/platelet-wcm/issues/5) | Granules as `UniqueMolecule` instances | Engine-side modelling decision; biology question is "are granules discrete-event-like enough to merit per-granule tracking?" |
| [#6](https://github.com/stevehaigh/platelet-wcm/issues/6) | `GranuleRelease` process — Ca²⁺-triggered exocytosis | The exocytosis trigger threshold + kinetics are biology questions. |
| [#7](https://github.com/stevehaigh/platelet-wcm/issues/7) | `GranuleState` listener | Records cargo inventory + cumulative secretion. |
| [#8](https://github.com/stevehaigh/platelet-wcm/issues/8) | `granule_secretion.py` analysis | Plots comparable to lumiaggregometry data. |

## Methodology / cascade dynamics

| # | Question | Notes |
|---|---|---|
| [#36](https://github.com/stevehaigh/platelet-wcm/issues/36) | **Stochastic** simulation of low-copy GPCR intermediates? | Sveshnikova et al. argue PAR1·Gαq·GTP and active PAR1 sit at single-digit copies per platelet; ODEs can't capture the platelet-to-platelet heterogeneity that produces the procoagulant subpopulation. MVP: stochastic P2Y1 only. |
| [#40](https://github.com/stevehaigh/platelet-wcm/issues/40) | Full hybrid **SSA + ODE** upstream cascade | The natural extension of #36 — entire upstream cascade stochastic, downstream Ca²⁺ machinery deterministic. Reproduces Sveshnikova's procoagulant-subpopulation framework. |

## Future / stretch — out of dissertation scope; want signal on relative interest

| # | Scope | Notes |
|---|---|---|
| [#12](https://github.com/stevehaigh/platelet-wcm/issues/12) | `IntegrinActivation` — αIIbβ3 inside-out signalling | Outside-in / inside-out coupling. Important for thrombus formation but a major scope expansion. |
| [#13](https://github.com/stevehaigh/platelet-wcm/issues/13) | `Metabolism` process — mitochondrial ATP + glycolysis | Currently the model assumes infinite ATP. Real cell has finite reserves; matters for activation duration. |
| [#35](https://github.com/stevehaigh/platelet-wcm/issues/35) | Multi-cell thrombus formation | Couples this single-cell WCM to a multi-cell aggregation framework (Diamond / Brass lineage). |
| [#41](https://github.com/stevehaigh/platelet-wcm/issues/41) | Spatial Ca²⁺ microdomains — submembrane subcompartment | Currently a well-mixed cytosol; sub-membrane microdomains around SOCE / IP3R puncta are biologically real but require a non-trivial PDE/lattice approach. |

## Writing

| # | Scope |
|---|---|
| [#43](https://github.com/stevehaigh/platelet-wcm/issues/43) | Perspective paper: "Different cells, different models" — comparing WCM frameworks across cell types (platelet here vs E. coli wcEcoli) |

## Mixed (UI surface for biology scenarios)

| # | Scope |
|---|---|
| [#11](https://github.com/stevehaigh/platelet-wcm/issues/11) | Webapp: named-scenario library (resting / thrombin_1nM / collagen_5ug_ml / …). Biology question: what scenarios? Tech question: where do they live in the UI / file layout? |

## Notes on what's been done recently (for context when reading the issues)

- **v0.5.x (May 2026)**: 2-D ADP × thrombin dose-response sweep (#45, closed). Headline finding: cascade is **graded upstream (peak IP3), binary at the IP3R→DTS-release step (peak Ca²⁺ saturates), graded again in the integrated response (Ca²⁺ AUC)**. See `reports/figures/dose-sweep-9x9*-panel.png` (one figure per range) and the live interactive surface at `reports/figures/dose-sweep-9x9-transition-interactive.html` (served via Pages).
- **v0.4.x (April–May 2026)**: GPCR cascade added (PAR1 / PAR4 / P2Y1 / P2X1 → Gαq → PLCβ → PI cycle → IP3). Replaces the previous forced-IP3 curve.
- **v0.3.x (March–April 2026)**: PMCA / SOCE / NCX / MCU / NCLX / CaM / multi-buffer DTS calibration. Phase 3 acceptance landed (5/5 against Dolan & Diamond 2014 Fig. 4).
- **Engineering plumbing** (May 2026): byte-identical regression suite, kinetics-as-data design doc, terminal replay TUI for dissertation demo.

## How to comment

GitHub issues directly is best (the discussion is then version-
controlled alongside the model). Comments here, threaded under the
right issue, are easier for me to act on than email. If GitHub
account barriers are an obstacle let me know and we'll work something
else out.
