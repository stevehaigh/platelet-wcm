# `reports/design/` — navigator

Design docs for the platelet WCM. ~20 files accumulated over the v0.2 → v0.4
arc; this README signposts what's current canonical vs historical context.

For *current work state* — what's being touched this week, today's
calibration results — read the most recent file in
[`../lab-books/`](../lab-books/) instead. The lab books are dated session
notes and are always more current than anything in this directory.

---

## Read this first (current canonical)

| File | What it is |
|---|---|
| [`model-status-2026-05-13.qmd`](model-status-2026-05-13.qmd) | v0.4.1 model status. Twelve coupled mechanisms agonist → cytosolic Ca²⁺; Phase 3 Dolan validation 5/5. The single best summary of the model as it stands today. |
| [`kinetics-as-data-2026-05-22.qmd`](kinetics-as-data-2026-05-22.qmd) | Issue #32 audit and implementation plan. Move ~150 numerical constants out of Python into TOML/TSV. Four-level spectrum from "Python-with-TOML-constants" to "fully data-driven biology". |
| [`codebase-review-2026-05-23.qmd`](codebase-review-2026-05-23.qmd) | Dev-readiness audit (max-effort review of `main`). Priority-ordered punch list, calibration-coupling map, test-coverage gap analysis, "not to do" list. |
| [`calibration-coupling-2026-05-25.qmd`](calibration-coupling-2026-05-25.qmd) | One-page matrix of "if you change X you must re-derive Y" chains in `calcium_signalling.py`. Skim before touching any rate constant. |

## Recent work (last ~3 weeks)

| File | What it is |
|---|---|
| [`pkc-p2y1-feedback-design-2026-06-11.qmd`](pkc-p2y1-feedback-design-2026-06-11.qmd) | **v0.6 design.** DAG → PKC → P2Y1 desensitisation negative feedback; closes the dead-end DAG branch. Anchored to Purvis 2008 (PKC module) + Mundell 2006 / Nicholas 2023 (PKC→P2Y1 in platelets). **Implemented** on branch `PKC-P2Y1-desensitisation`: both feedback routes — P2Y1 desensitisation (Mundell/Nicholas) **and** PLCβ phosphorylation (Purvis route, added on Mike's request) — plus `runPerturbation.py` `pkc`/`plcb` figures. calcium-v0.6 TOML/TSV, ODE, listener columns, goldens regenerated, Dolan 5/5 preserved. |
| [`pkc-downstream-effects-2026-06-12.qmd`](pkc-downstream-effects-2026-06-12.qmd) | **v0.61 design.** Scopes the downstream PKC effects deferred from v0.6: granule secretion (autocrine ADP), thromboxane A₂ (autocrine TP→Gq, aspirin target), integrin αIIbβ3 affinity state. Recommends sequencing (secretion → thromboxane → integrin); flags the single-cell limit (aggregation out of reach). Frames PKC as the hub of four feedback loops (2 brakes + 2 amplifiers). Design only — no code. |
| [`pathway-diagram-review-2026-05-19.qmd`](pathway-diagram-review-2026-05-19.qmd) | Annotation list for the BioRender pathway figure used in lab-meeting slide 7. |
| [`lab-meeting-2026-05-14.qmd`](lab-meeting-2026-05-14.qmd) | Lab-meeting presentation: methodology, results, AI-assisted validation. |
| [`mike-report-2026-05-14.qmd`](mike-report-2026-05-14.qmd) | PNAS-style write-up on feasibility of mechanistic whole-cell platelet modelling. |
| [`model-status-graphviz.qmd`](model-status-graphviz.qmd) | Graphviz vs Mermaid layout-engine comparison for pathway visualisations. |
| [`kinetics-as-data-sketch.qmd`](kinetics-as-data-sketch.qmd) | Single-reaction sketch of the post-dissertation refactor — companion to `kinetics-as-data-2026-05-22`. Also exists as `.md`. |
| [`diagrams-v0.3.0.md`](diagrams-v0.3.0.md) | Mermaid versions of the v0.2 ASCII diagrams; updated post Phase 4 / #31 (PI cycle replaces forced IP3). |

## Errata / footnotes

| File | What it is |
|---|---|
| [`purvis-2008-k3-transcription-error.md`](purvis-2008-k3-transcription-error.md) | 100× transcription error in Purvis 2008 Table 1 (IP3R k₃ = 11 s⁻¹ vs original Sneyd & Dufour 2002 0.11 s⁻¹). Documents the discovery and propagation through downstream models. |

## Historical context (v0.2 era — 2026-05-07)

These describe the model *before* the v0.3 / Phase 3+ work (SERCA cycle, DTS
buffers, SOCE, P2X1, CaM ladder). Useful for the demo flow and
walkthroughs but **not** the current biology.

| File | What it is |
|---|---|
| [`biology-overview-2026-05-07.md`](biology-overview-2026-05-07.md) | Single-cell deterministic Ca²⁺ model overview at v0.2/v0.3.0. IP3-mediated transient calibrated against Dolan 2014. |
| [`demo-2026-05-07.md`](demo-2026-05-07.md) | Demo cheatsheet — Dolan 3/5 acceptance criteria at v0.2; SOCE-differential gap is the named root cause for the two failures (now resolved in v0.4). |
| [`code-overview-2026-05-07.md`](code-overview-2026-05-07.md) | Companion to the demo: file/line navigation for the live walkthrough; framework (`wholecell/`) vs model (`models/platelet/`) split. |
| [`code-walkthrough-2026-05-07.md`](code-walkthrough-2026-05-07.md) | Temporal trace from `runPlateletSim.py` to first row of output. |

## Pre-v0.2 design docs (undated)

The original design work that bootstrapped the model. Mostly of historical interest.

| File | What it is |
|---|---|
| [`calcium-dynamics-design.md`](calcium-dynamics-design.md) | v0.2 calcium dynamics design — first real biochemistry targeting Dolan Fig. 4. |
| [`calcium-next-steps-plan.md`](calcium-next-steps-plan.md) | Plan dated 2026-04-29 for reproducing Dolan Fig. 4 — superseded by the model-status doc. |
| [`calcium-signalling-pathway-design.md`](calcium-signalling-pathway-design.md) | Architecture for ADP receptor signalling (P2Y1, P2Y12). Uses the wcEcoli `TwoComponentSystem` pattern with stoichiometry matrix + SymPy-generated ODEs. Predates the current `calcium_signalling.py` dataclass approach. |
| [`platelet-runtime-scaffold-summary.md`](platelet-runtime-scaffold-summary.md) | First runnable `models/platelet/` scaffold — proves the wholecell engine can host a platelet namespace. |
| [`reconstruction-platelet-design.md`](reconstruction-platelet-design.md) | Design for the `reconstruction/platelet/` namespace; keeps useful reconstruction architecture from wcEcoli while replacing bacterium-specific parts. |

---

## How to add a new doc

- New design docs default to **`.qmd`** (Quarto). See `CLAUDE.md` "Reports & docs" section for the canonical frontmatter (or copy from [`model-status-2026-05-13.qmd`](model-status-2026-05-13.qmd)).
- Date in the filename: `<topic>-YYYY-MM-DD.qmd`. The publish script (`runscripts/manual/buildDocsSite.py`) indexes by this date.
- Add the new file to the relevant section above.
