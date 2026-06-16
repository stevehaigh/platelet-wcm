---
title: "Lab book — 2026-06-16: replay-TUI modernisation, model-diagram refresh, and the DTS V_IM correction"
---

# Lab book — 2026-06-16: TUI, diagram, and the DTS V_IM correction

## Context

A mixed session spanning three threads: (1) finishing and landing the terminal
replay-TUI modernisation, (2) reviewing the lab-meeting model diagram and
producing an update plan, and (3) — the substantive science — a
biology-first investigation of how far the DTS Ca²⁺ store depletes, which ended
in a model change (`V_IM = 0`). Notes captured here because they were spread
across PRs, a sweep probe, and a long design discussion.

---

## 1. Replay TUI modernised and merged (PR #60)

Brought `runscripts/manual/replayTui.py` up to the v0.63 biology and landed it.

- **New biology in the schematic** — PKC feedback line (active PKC, P2Y1-desens %,
  PLCβ-phos %) and a downstream-outputs block (secretion %s + autocrine ADP[e],
  TXA₂/TP, integrin PAC-1), reading the optional SecretionTrace/ThromboxaneTrace/
  IntegrinTrace listeners with graceful fallback.
- **Auto-sizing cell box** — the schematic interior now grows from a 96-col floor
  to a 112-col target when the terminal allows and clamps so it never crops
  (`_fit_inner_w`); added vertical breathing rows; DTS/Mito boxes derive their
  widths from one number instead of hand-calibrated literals.
- **`?` help overlay** (`HelpScreen` modal) documenting every on-screen field with
  units + literature source; keyboard-driven (no mouse needed). We deliberately
  did *not* go widget-per-field: the schematic is one ASCII-art `Static`, and
  decomposing it would rewrite the renderer and reintroduce the rich.Live
  blank-region fragility we removed.
- **rich + textual are now an optional extra** (`requirements-viz.txt`), out of the
  core install; CI still covers them; tests `importorskip`.
- Two rounds of Copilot review addressed (dead-code removal; sparkline width/
  alignment via `_visible_len`; docstrings; launcher installs the pinned viz set).
- **Merged to `main`** (squash `c142a44c`); superseded and closed #48.

## 2. Model-diagram refresh (issue #61)

Reviewed the pre-v0.5 BioRender schematic
(`haigh.bio/decks/lab-meeting-2026-05-14/img/image4.png`) against the current
model and produced an actionable update plan:

- `reports/figures/model-diagram-update-checklist.md` — a tick-list (un-grey PKC
  as the hub of 2 brakes + 2 autocrine amplifiers + 3 terminal outputs; granules/
  secretion; thromboxane TXA₂→TP; new integrin αIIbβ3/PAC-1; repaint the greyed
  region as the v0.7 inhibitory axis), an alignment-checked ASCII layout map, and a
  BioRender icon reference (search terms + placeability flags from a live library
  search; PKC/COX-1/cAMP/sGC/PGI₂/VASP have no clean placeable icon).
- Filed **issue #61** (assigned) for the manual BioRender work; checklist on branch
  `model-diagram-checklist`.

## 3. The DTS V_IM correction (the science)

**Question (biology-first).** Does the DTS deplete only to equilibrium with the
cytosol, or empty completely? Complete emptying looked non-physical.

**Finding.** Reading the store trace, the DTS settles at **~0.03 µM free Ca²⁺ —
*below* the cytosolic free Ca²⁺ (~0.42 µM)**. A passive IP₃R cannot pump uphill,
so this was suspicious. The IP₃R is a Nernst flux,
`−γ·N·Po·(V_IM − E_Ca)` with `E_Ca = (RT/zF)·ln([Ca]_dts/[Ca]_cyt)`; zero net flux
sits at `[Ca]_dts/[Ca]_cyt = exp(V_IM·zF/RT)`. The model carried **`V_IM = −60 mV`**
(borrowed from Dolan's −100..−60 mV *plasma-membrane* sampling), giving an
equilibrium ratio of ~0.011 — it actively pulled the store to ~1 % of cytosol.

**Cross-check vs Dolan.** Dolan's own Fig 4A bottoms the store at **120–180 µM
(+Ca) / 80–130 µM (EDTA)** — *partial* depletion. Our ~0 contradicts Dolan's
store trace while still passing the 5/5 criteria, because those score only the
cytosolic peak/SOCE, never store depth.

**V_IM sensitivity sweep** ({−60, −40, −20, 0, +40} mV):

| V_IM | +Ca peak | DTS min (+Ca) | SOCE Δ | Dolan 5/5 |
|------|----------|---------------|--------|-----------|
| −60 mV (old) | 425 nM | 0.03 µM (*below cyt*) | 111 | 5/5 |
| 0 mV (new)   | 417 nM | 1.9 µM | 143 | 5/5 |
| +40 mV       | 398 nM | 30.9 µM | 148 | 5/5 |

Two results: (i) **5/5 is insensitive to V_IM** (peak is buffer/PMCA/SOCE-clamped);
(ii) **V_IM does not control depletion *depth*** — even at 0 mV the store still
depletes ~99 %, because under sustained IP₃ the open IP₃R + SERCA equilibrate it
low. A follow-up probe confirmed IP₃ never returns to its ~50 nM rest
(stays 175–220 nM across configs), so the channel stays open. Dolan's store stops
at 150 µM only because his IP₃ was a *transient* pulse.

**Discussion / decision.** Steve's position: don't chase Dolan's number — under a
*sustained* agonist, deep depletion is the more realistic biology (and drives the
SOCE we validate). Agreed. The one genuinely unphysical part was the
*sub-cytosolic undershoot* — a passive channel can't drive the store below
cytosolic free Ca²⁺, and the ER/DTS membrane holds ≈0 mV (counter-ion permeable),
not −60 mV. **Decision: set `V_IM = 0`.** The store now bottoms at the cytosolic
equilibrium (~1.9 µM +Ca / ~0.4 µM EDTA) — still ~99 % depleted but
thermodynamically honest. `V_PM` stays −60 mV (the plasma membrane *does* hold a
resting potential; only the ER membrane was wrong).

**Implementation.**
- `reconstruction/platelet/dataclasses/process/calcium_signalling.py`: `V_IM_V = 0.0`
  with the thermodynamic rationale in the comment.
- Dolan 5/5 re-verified (SOCE differential *improved* 111→143 nM); byte-identical
  goldens regenerated (`at_rest_30s`, `default_activation_30s`).
- Regression guard added:
  `test_regression.py::test_dts_never_drains_below_cytosol` (DTS free Ca²⁺ ≥
  cytosolic free Ca²⁺ across the run) so a future negative V_IM can't silently
  reintroduce sub-cytosolic draining.
- Design doc updated: `reports/design/dts-depletion-literature-2026-06-14.qmd` §8.
- On branch `dose-response-v0.63` (PR #59 — the DTS-depletion investigation; the
  V_IM fix is its resolution).

**Still open (separable, not forced):** whether to also hold a *typical-stimulus*
residual store (~30–50 %) via stronger IP₃R inactivation (design-doc §6 option b).
V_IM fixed the *floor*; depth remains a defensible-range calibration best framed by
the dose-response result.

---

## Status / next

- PR #60 (TUI) merged to main; #48 closed.
- PR #59 (dose-response + DTS-depletion) now also carries the V_IM=0 fix.
- Issue #61 open for the BioRender diagram work.
- Noted (not fixed): pre-existing mypy `RunConfig(**dynamic_dict)` typing
  complaint in `runDoseResponse.py` (outside the CI mypy path).
