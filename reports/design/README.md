# `reports/design/` — navigator

Design docs for the platelet WCM. This README signposts what's current.

For *current work state* — what's being touched this week, today's calibration
results — read the most recent file in [`../lab-books/`](../lab-books/) instead.
The lab books are dated session notes and are always more current than anything
in this directory.

> Superseded v0.2–v0.5 design/planning docs have been moved to
> [`../archive/`](../archive/) (see its README). Only current material is listed
> below.

---

## Read this first (current model status)

| File | What it is |
|---|---|
| [`validation-map-2026-06-19.qmd`](validation-map-2026-06-19.qmd) | **The current model-status / validation overview.** What each mechanism is validated against, and the Dolan 5/5 placed in the context of the broader validation portfolio. Start here for "what does the model do now, and how do we know it's right". |
| [`version-comparison-v0.5-v0.6-v0.61-2026-06-14.qmd`](version-comparison-v0.5-v0.6-v0.61-2026-06-14.qmd) | Model evolution v0.5 → v0.6 → v0.61: biology added, model changes, and results per release (pre-v0.7; the inhibitory axis and the #73/#76 arms came after). |

## Pathway designs (current biology)

| File | What it is |
|---|---|
| [`pkc-p2y1-feedback-design-2026-06-11.qmd`](pkc-p2y1-feedback-design-2026-06-11.qmd) | **v0.6** — DAG → PKC → P2Y1 desensitisation negative feedback (Purvis 2008; Mundell 2006 / Nicholas 2023). |
| [`pkc-downstream-effects-2026-06-12.qmd`](pkc-downstream-effects-2026-06-12.qmd) | **v0.61** — design-of-record for the built downstream PKC outputs: granule secretion (autocrine ADP), thromboxane (autocrine TP→Gq, aspirin target), integrin αIIbβ3. PKC as the hub of four feedback loops. |
| [`inhibitory-axis-design-2026-06-15.qmd`](inhibitory-axis-design-2026-06-15.qmd) | **v0.7** — the cAMP/PKA inhibitory axis: P2Y12/Gi (clopidogrel), PGI₂/Gs + forskolin + PDE3 (cilostazol), VASP readout. |
| [`pi3k-akt-rap1b-arm-2026-06-22.qmd`](pi3k-akt-rap1b-arm-2026-06-22.qmd) | **#73** — PI3K/Akt → Rap1b sustain arm; Rap1b-GTP is the integrin's proximal driver (built; carries an as-built status note). |
| [`mcu-coupling-2026-06-23.qmd`](mcu-coupling-2026-06-23.qmd) | **#76** — bounding the mitochondrial Ca²⁺ sink + the MCU → IP3R-relief coupling that flips the MCU-KO direction toward the data (built). |
| [`dts-depletion-literature-2026-06-14.qmd`](dts-depletion-literature-2026-06-14.qmd) | Literature check on whether the model's complete DTS Ca²⁺ depletion is a real discrepancy (incl. §8 V_IM=0, §9 γ_IP3R recalibration). |
| [`model-driven-experiment-design-2026-06-14.qmd`](model-driven-experiment-design-2026-06-14.qmd) | Model-driven experimental design around the DTS-depletion question. |

## Calibration & kinetics

| File | What it is |
|---|---|
| [`calibration-coupling-2026-05-25.qmd`](calibration-coupling-2026-05-25.qmd) | One-page matrix of "if you change X you must re-derive Y" chains in `calcium_signalling.py`. Skim before touching any rate constant. |
| [`kinetics-v0.6-review.qmd`](kinetics-v0.6-review.qmd) | The clickable kinetics review rendered from `calcium-v0.6.toml` (auto-linked citations). Built in CI. |

## Tooling / infrastructure designs

| File | What it is |
|---|---|
| [`tui-tinkering-dashboard-2026-06-15.qmd`](tui-tinkering-dashboard-2026-06-15.qmd) | Terminal-TUI experiment-bench design (edit biology, knock out pathways, run, watch the Ca²⁺ trace live). P0–P2 built in `wholecell/tui/`. |
| [`web-ui-rewrite-2026-06-18.qmd`](web-ui-rewrite-2026-06-18.qmd) | Web-UI rewrite over a shared bench-core (FastAPI + JS). Parked post-dissertation; tracked in issue #64. |
| [`perturbation-figures.qmd`](perturbation-figures.qmd) | Perturbation-figure design (PMCA V_max & MCU). *Note: predates #76; its MCU claim is superseded by the MCU-coupling work.* |
| [`model-status-graphviz.qmd`](model-status-graphviz.qmd) | Graphviz vs Mermaid layout comparison; source of the v0.5 architecture figure (being replaced by the BioRender diagram). |

## Feasibility / external-facing

| File | What it is |
|---|---|
| [`mike-report-2026-05-14.qmd`](mike-report-2026-05-14.qmd) | PNAS-style feasibility write-up on mechanistic whole-cell platelet modelling (v0.4.1-era; some authoring stubs remain). |

## Errata / footnotes

| File | What it is |
|---|---|
| [`purvis-2008-k3-transcription-error.md`](purvis-2008-k3-transcription-error.md) | 100× transcription error in Purvis 2008 Table 1 (IP3R k₃ = 11 s⁻¹ vs the original Sneyd & Dufour 2002 0.11 s⁻¹) and its propagation through downstream models. |

---

## How to add a new doc

- New design docs default to **`.qmd`** (Quarto). See `CLAUDE.md` "Reports & docs" section for the canonical frontmatter, or copy the block from any current `.qmd` here (e.g. [`mcu-coupling-2026-06-23.qmd`](mcu-coupling-2026-06-23.qmd)).
- Date in the filename: `<topic>-YYYY-MM-DD.qmd`. The publish script (`runscripts/manual/buildDocsSite.py`) indexes by this date.
- Add the new file to the relevant section above. Superseded docs go to [`../archive/`](../archive/).
