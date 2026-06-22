# Platelet WCM — TUI demo script

Quick walk-through for showing off the terminal experiment bench
(`wholecell/tui/`). Everything below is driven from **inline form widgets** —
no JSON editing. The TUI generates the `run_config.json` for you on each run;
hit **Save** to keep any setup as a named preset.

## Launch

```bash
make tui
```

Uses `pyenv exec python` (pinned 3.11.5) so the viz deps resolve. If it
complains about missing packages:

```bash
pip install -r requirements-viz.txt   # rich, textual, textual-plotext, plotext
```

Keys: **r** run · **b** set baseline · **f** demo figures (per-theme) · **q** quit.
The right pane streams cytosolic + DTS Ca²⁺ live from `live.csv`.

**Output, naming & overlays.** Each run writes to its own directory under
`out/` — named by the **save as…** field if you've typed one, otherwise a
timestamp — so runs no longer overwrite each other and a pinned run survives to
be compared against. **Set baseline** (`b`) pins the last run; the next **Demo
figure** (`f`) then renders **four focused per-theme figures** into the run's
plot folder — `demo_calcium`, `demo_integrin`, `demo_thromboxane`,
`demo_secretion` — each with that baseline overlaid in **grey**, and opens the
folder. Each demo below names the figure it uses, so you only look at the one
that matters (no more hunting for the relevant panel in a wall of plots).

> **Why the demos are framed this way:** under the default *saturating* agonist
> the Ca²⁺ response is **store-limited**, so most knockouts/loops barely move the
> cytosolic trace. The compelling demos either read a *different* output or use a
> *weak/isolated* agonist so the pathway under test is the dominant driver — the
> exception is **MCU** (Demo 2), which competes for cytosolic Ca²⁺ directly during
> the transient and so *does* move it.

---

## Demo 1 — The validation figure: ±extracellular Ca²⁺ (Dolan & Diamond 2014, Fig. 4)

Shows the model's primary validation claim directly on the live trace.

1. Preset dropdown → **Agonist transient** (200 s, Ca²⁺_ex = 1.2 mM). Press **r**.
2. When done, press **b** to freeze it as the grey baseline.
3. Preset dropdown → **EDTA (no Ca_ex)** (Ca²⁺_ex = 0 → SOCE + PM leak off).
   Press **r**.

**What to point at:** with external Ca²⁺ the cytosolic trace holds a sustained
plateau; without it the DTS store depletes and the signal collapses back to
baseline — the canonical Dolan contrast, side-by-side via the baseline overlay.

*(Manual equivalent of the same set: `ca_ex_mM` 1.2 → 0.0.)*

---

## Demo 2 — MCU knockout: mitochondria buffer the cytosol (a counterintuitive prediction)

The model *predicting* something you can't read off the knob — the payoff of
having every Ca²⁺ pathway present and mass-balanced, not an isolated submodel.

1. Preset dropdown → **Agonist transient** (standard +Ca²⁺ stimulus). Press **r**,
   then **b** (wild-type baseline).
2. Tick **MCU V_max → KO** (Pumps / brakes group). Press **r**.

**What to point at:** knocking out mitochondrial Ca²⁺ *uptake* makes the cytosolic
peak go **up**, not down (~+15%; ≈ 407 → 468 nM here) — in **`demo_calcium`** the
**cytosolic Ca²⁺ panel** shows the KO run sitting *above* the grey wild-type
baseline. Meanwhile the **DTS panel** (same figure) shows the store still empties
either way:
mitochondria *redistribute* Ca²⁺ away from the cytosol during the transient
rather than returning it to the store. The naive expectation — remove an uptake
pathway, get less Ca²⁺ — is wrong, and the direction matches the elevated
cytosolic Ca²⁺ in MCU-knockout platelets (Ghatge 2026), a phenotype the model was
never tuned to. It's also the *one* lever that visibly moves the otherwise
store-clamped cytosolic peak.

---

## Demo 3 — Clopidogrel throttles integrin activation *without touching the integrin*

The opposite of a knock-out-the-thing demo: αIIbβ3 is fully present, yet an
antiplatelet drug turns its output down — because the effect arrives through the
cascade (P2Y12 → Gi → cAMP↑ → PKA brake → integrin).

1. **Weak ADP agonist** so the autocrine ADP → P2Y12 axis is in play (it's
   swamped at saturation): set `Thrombin (nM)` = **0**, `ADP (uM)` = **0.5**,
   `ATP_ex (uM)` = **0**, `Length` = **200**. Press **r**, then **b** (control).
2. Add clopidogrel: set `P2Y12 block (clopidogrel)` = **1.0** (Feedback-loops
   group). Press **r**.

**What to point at:** in **`demo_integrin`** the **PAC-1 trace** drops (≈ 64% →
53% here, −11 points) against the grey control — with αIIbβ3 copy number
**unchanged** (the PKA-brake line, right axis, shows the mechanism re-engaging).
Blocking P2Y12 stops autocrine ADP from lowering cAMP, so the PKA brake
re-engages and the terminal output falls. This is the real clopidogrel / VASP
mechanism (the basis of the clinical platelet-reactivity assay), and it shows
several pathways integrating onto one readout.

*Contrast (construct check): the **Glanzmann (aIIbb3 KO)** preset zeroes PAC-1 by
removing the integrin itself — same readout, but trivially, by construction.
Clopidogrel reaching that same readout through the cascade, integrin intact, is
the more informative result.*

---

## Demo 4 — Aspirin and the second wave (autocrine amplification at a weak agonist)

Shows the v0.61 feedback loops in the regime they were built for.

1. **Weaken the stimulus** so the loops matter (they're swamped at saturation):
   set `Thrombin (nM)` ≈ **0.2** and `ADP (uM)` = **0**, `Length` = **300**.
   Press **r**, then **b** (full-loop baseline).
2. Knock out the thromboxane amplifier: tick **COX-1 (aspirin) → KO**
   (*or* pick the **Aspirin (COX-1 KO)** preset, then re-apply the weak
   thrombin/ADP edits). Press **r**.

**What to point at:** with the loops intact, secreted TXA₂ feeds back through TP
and lifts the late ("second-wave") Ca²⁺; aspirin flattens it. With the full-loop
run pinned as baseline, **`demo_thromboxane`** drops to ~0 under aspirin against
the grey baseline, **`demo_secretion`** shows the autocrine ADP[e] / second-wave
difference, and **`demo_calcium`** shows the lifted late cytosolic Ca²⁺.

**Bonus — autocrine ADP in isolation:** thrombin-only run (`ADP (uM)` = 0) makes
*secreted* ADP the sole P2Y1 driver. Toggle **Autocrine ADP gain → KO** to show
the loop's contribution vanish.

---

## Reset

Preset dropdown → **Defaults** clears everything back to schema defaults; the
"Δ from defaults" line under the form always tells you what's currently edited.
