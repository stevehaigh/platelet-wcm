---
title: "Lab book — 2026-06-22: TUI demo figures, thesis Discussion, doc fixes, and a Zou-2022 reversibility follow-up"
---

# Lab book — 2026-06-22: demo figures + reversibility (Zou 2022)

## One-line summary

A consolidation-and-communication day rather than a new-biology day: a critical
appraisal of the model, three factual corrections, a drafted thesis **Discussion**,
a rebuilt **TUI demo-figure system** (focused per-theme figures with baseline
overlays + per-run output dirs), and a **reversibility demo + design doc** prompted
by mapping the model onto Zou et al. 2022. Shipped as **PR #74** (4 commits) and
**issue #73** (a design to implement next). No change to the calcium core →
**Dolan 5/5 preserved** throughout.

## 1. Critical appraisal of the model (what's solid / what's thin)

Started by stress-testing the model against three questions — what gives
confidence it's useful, what does it get right, where is it incomplete/wrong.
Verdict (grounded in a four-front audit of validation, parameters, results,
limitations):

- **Useful as a mechanistic integrator / hypothesis generator, not (yet) a
  quantitative predictor of cytosolic Ca²⁺ amplitude.** The genuinely emergent,
  literature-matched results carry the weight: MCU KO *raising* cytosolic Ca²⁺
  (Ghatge 2026), the autocrine commitment switch, the second wave, graded IP₃.
- **Central limitation:** cytosolic free Ca²⁺ is store-limited / SOCE-clamped, so
  Dolan 5/5 has low discriminating power and now functions as a regression
  invariant; feedback/inhibition is read on functional *outputs*, not free Ca²⁺.
- **Parameter weak spots:** the 6 fL cytosol volume (universal denominator), the
  calibration-coupled γ_IP3R↔SERCA↔γ_SOC↔PM-leak chain, SERCA Vmax possibly
  2–5× high, the downstream module's uncited bare-number rate constants.
- **Structural:** single-cell (no aggregation), `--seed` not actually
  reproducible, no metabolism, well-mixed.

This appraisal framed everything below (it is also the spine of the new thesis
Discussion).

## 2. Three factual fixes (comments/docs only — no behaviour change)

1. `calcium_signalling.py` docstring **γ_IP3R 0.35 → 0.135 pS** (the live,
   recalibrated value from `calcium-v0.6.toml`; 0.35 was never in use).
2. **`N_IP3R = 1328` re-attributed to Dolan 2014 Table S1**, not Burkhart 2012
   (whose total IP₃R is ~4850) — corrected in both the code comment and the TOML.
3. CLAUDE.md `--seed` row corrected: the one stochastic process (`RestingDecay`)
   draws from numpy's **global** RNG, so runs are **not** currently
   seed-reproducible (was "no stochastic processes use it").

## 3. Thesis Discussion drafted

Replaced the planning-note stub in `reports/thesis/draft-thesis.qmd` with a
~1,750-word **Discussion**: feasibility as the principal finding; extensibility
(four downstream modules on the Ca²⁺ core); emergent vs construct-by-design
results; the store-limited Ca²⁺ clamp as the central validation caveat;
parameter and single-cell limitations; the Purvis k₃ / AI-extraction reflection;
future directions. Citations use verified `references.bib` keys; unverified ones
(Ghatge 2026, Sveshnikova 2025, Zou 2022) are `[cite: …]` placeholders to pull
from Zotero. PDF re-renders clean (0 undefined refs; gitignored artifact).

## 4. TUI demo-figure system reworked

**Problem (raised by Steve):** the single 6-panel `demo_overview` figure forced
every demo to show four irrelevant panels.

**Fix:** dropped the monolith for **five focused per-theme figures** —
`demo_calcium` (Ca_cyt, DTS+SOCE, IP₃), `demo_integrin` (PAC-1 + PKA brake),
`demo_thromboxane` (TXA₂/TXB₂), `demo_secretion` (autocrine ADP + cargo release),
`demo_reversibility` (see §5) — sharing `single/_demo_common.py`. Each supports a
grey **baseline overlay** via `PLATELET_BASELINE_SIMOUT`; the "Demo figure" button
renders the set and opens the run's plot folder.

**Supporting changes:**

- **Per-run output dirs:** the TUI now writes each run to
  `out/<save-as-name | timestamp>/` instead of overwriting `tui_run/`, so runs no
  longer clobber each other and a pinned baseline survives to be overlaid.
- `analysisPlatelet.py --out-name` to name a figure file (single `--plot`).
- All `single/demo_*` use mathtext for chemistry (no glyph-drops — the older
  `calcium_trace`/`phase3` figures still warn on raw-unicode `Ca²⁺`).

**Demo set evolution (`tui-demos.md`).** Steve pushed back that the Glanzmann
demo ("knock out the integrin → PAC-1 = 0") is **true-by-construction** — a
multiplication by zero, not biology. Agreed. Replaced/augmented with:

- **Demo 2 — MCU knockout (counterintuitive, emergent):** knocking out
  mitochondrial Ca²⁺ *uptake* **raises** the cytosolic peak (**+15 %, 407 → 468
  nM**, verified) and the store still empties — matches Ghatge 2026 MCU⁻/⁻
  platelets, never tuned to it.
- **Demo 3 — clopidogrel throttles PAC-1 with the integrin intact:** P2Y12 block
  lowers PAC-1 (**≈ 64 % → 53 %**, −11 pts) via the cAMP→PKA brake; the figure's
  PKA-brake trace shows the mechanism (control dis-inhibited ~1.8 vs blocked 1.0).
- **Glanzmann kept as a one-line construct-check foil** inside Demo 3.

Tests updated (per-run dirs, the figure set, a `_run_root` unit test);
**38 TUI+analysis tests pass, mypy clean.**

## 5. Zou 2022 — reversible αIIbβ3 activation (the day's bit of science)

Mapped the model onto **Zou, Swieringa, de Laat et al. 2022** (*IJMS* 23:12512,
"Reversible Platelet Integrin αIIbβ3 Activation and Thrombus Instability"; a
review). Its thesis: αIIbβ3 activation is an intrinsically **reversible**
equilibrium inside-out switch whose persistence needs sustained autocrine ADP and
is opposed by cAMP/PKA; reversibility underlies thrombus instability.

**Empirically tested what the model robustly shows** (long transient-agonist runs):

- **Reversible activation — YES, robust.** PAC-1 rises then falls (ADP-only:
  **0 → 0.67 @ ~150 s → 0.41 @ 600 s**), the relaxation lagging the autocrine
  ADP[e] being cleared by the model's ecto-NTPDase (the paper's apyrase/CD39
  mechanism).
- **PGI₂ / cAMP antagonism — YES, clean.** +PGI₂ 50 nM crushes the peak
  (**0.67 → 0.18**) — the prostacyclin arm.
- **Weak-(ADP)-vs-strong-(thrombin) reversibility — NOT robust.** Initially I
  framed a contrast (ADP reverses more), but a head-to-head test showed
  thrombin-only reverses about as much as ADP-only: the model's **PARs internalise
  rather than latch**, so thrombin's drive also decays. The earlier −17 % was a
  combined high-dose artefact. **Claim dropped** from the demo and thesis para.

**Built:** `single/demo_reversibility.py` (Demo 5; PAC-1 rise-then-fall + its
autocrine-ADP driver) and a thesis Discussion paragraph — both honestly bounded
(reversibility is a *designed* 2-state feature, not a prediction; the agonist-class
distinction isn't reproduced; aggregation is out of single-cell scope).

**The gap, written up as a design + issue.** The paper's *primary* (ir)reversibility
node is **PI3K → Akt → Rasa3 ⊣ Rap1b**, downstream of P2Y12. The model has only
the cAMP/PKA arm of P2Y12 — so clopidogrel **lowers the peak** but does **not
collapse the sustained phase** (the endpoint is unchanged). Wrote a minimal design
to close this (`reports/design/pi3k-akt-rap1b-arm-2026-06-22.qmd`): two lumped
states (`Akt_active`, `Rap1b_GTP`), fast Ca²⁺/CalDAG-GEFI + slow P2Y12→Akt routes
(Stolla 2011), integrin driver swapped from the PKC×Ca gate to Rap1b-GTP, calcium
ODE untouched (Dolan preserved by construction). Tracked as **issue #73**.

## 6. Shipped

Branch `tui-demos-and-doc-fixes`, **PR #74** (4 commits, pushed; CI green
locally — 38 TUI+analysis tests, mypy clean across CI scope):

1. `54172334` — three doc/fact fixes (§2)
2. `ebe3af04` — thesis Discussion draft (§3)
3. `967c90ac` — TUI per-theme demo figures + per-run output dirs (§4)
4. `73111456` — Zou 2022 follow-up: reversibility demo + thesis para + PI3K/Akt
   design (§5)

Plus **issue #73** (implement the PI3K/Akt→Rap1b arm).

## 7. Honest corrections worth remembering

- The integrin **weak-vs-strong agonist reversibility distinction is not robust**
  in the model (PARs internalise) — don't claim it; the PI3K/Akt arm + PAR
  latching would be needed.
- **Glanzmann KO is true-by-construction** — useful only as a foil, not a headline.
- **MCU KO** is the cleanest emergent integrin-adjacent demo; **PGI₂ suppression**
  is the cleanest cAMP/integrin contrast.

## 8. Next

1. Review / merge **PR #74**.
2. Implement **issue #73** (PI3K/Akt→Rap1b arm) — would make the clopidogrel and
   reversibility demos genuinely P2Y12-dependent.
3. Still live: inhibitory-axis Slice 3 (NO/cGMP) + Ca²⁺-clamp recalibration (#70).
