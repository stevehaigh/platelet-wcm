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
figure** (`f`) then renders **five focused per-theme figures** into the run's
plot folder — `demo_calcium`, `demo_integrin`, `demo_thromboxane`,
`demo_secretion`, `demo_reversibility` — each with that baseline overlaid in
**grey**, and opens the folder. Each demo below names the figure it uses, so you only look at the one
that matters (no more hunting for the relevant panel in a wall of plots).

**Match the lengths.** The grey baseline trace stops where *its* run ended, so a
60 s baseline (the default `Length`) against a 200 s run looks truncated at 60 s.
Run the baseline and the comparison at the **same `Length`**; each baseline's
duration is now shown in the legend (e.g. "… (baseline · 60 s)") so a mismatch is
obvious rather than mysterious.

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

## Demo 2 — MCU knockout: the one lever that moves the store-clamped cytosolic peak

A worked example of confronting a model prediction with data — the payoff (and
the honesty test) of having every Ca²⁺ pathway present and mass-balanced, not an
isolated submodel.

1. Preset dropdown → **Agonist transient** (standard +Ca²⁺ stimulus). Press **r**,
   then **b** (wild-type baseline).
2. Tick **MCU V_max → KO** (Pumps / brakes group). Press **r**.

**What to point at:** knocking out mitochondrial Ca²⁺ *uptake* makes the cytosolic
peak go **up**, not down (~+15%; ≈ 407 → 468 nM here) — in **`demo_calcium`** the
**cytosolic Ca²⁺ panel** shows the KO run sitting *above* the grey wild-type
baseline. Meanwhile the **DTS panel** (same figure) shows the store still empties
either way:
mitochondria *redistribute* Ca²⁺ away from the cytosol during the transient
rather than returning it to the store. Removing a buffer and getting *more* of the
buffered species is the **expected** direction (not counterintuitive). **Honest
caveat:** real platelet MCU-knockout data show agonist-evoked cytosolic Ca²⁺
*reduced*, not raised (Ghatge 2026; Ajanel 2025) — the model treats MCU as
buffer-only and so **diverges** here; see `reports/experiments/3-mcu-knockout.qmd`
and issue #76. As a demo it remains the *one* lever that visibly moves the
otherwise store-clamped cytosolic peak.

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
**unchanged**. The **lower panel** shows why: under clopidogrel the PKA brake
stays at ~1.0 (engaged), versus the control's dis-inhibited hump to ~1.8.
Blocking P2Y12 stops autocrine ADP from lowering cAMP, so the PKA brake stays
engaged and the terminal output falls. This is the real clopidogrel / VASP
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

## Demo 5 — Reversible αIIbβ3 activation (Zou et al. 2022)

Maps to Zou et al. 2022 (*Int. J. Mol. Sci.* **23**:12512): αIIbβ3 activation is
an intrinsically **reversible** inside-out switch, not a permanent on-state, and
that reversibility underlies thrombus *instability*. The model's integrin is a
reversible two-state switch, so under a *transient* agonist it shows this directly.

Use a **long run with a transient agonist** so both the rise *and* the fall show:

1. `Thrombin (nM)` = **0**, `ADP (uM)` = **5.0**, `ATP_ex (uM)` = **0**,
   `Length` = **600**. Press **r**, then **f**.

**What to point at:** in **`demo_reversibility`** the **PAC-1 trace rises then
falls** (≈ 0 → 0.67 @ ~150 s → 0.41 @ 600 s) — reversible activation, not a latch.
The right panel shows **autocrine ADP[e]** peaking (~2.9 µM) then cleared to zero
by ~250 s (ecto-NTPDase = the paper's apyrase / CD39 mechanism); the PAC-1
relaxation *lags* the ADP clearance, because the integrin's return to the resting
conformation is slow. This is the single-platelet affinity-state basis of the
paper's thrombus-instability argument.

**Bonus — prostacyclin antagonism:** press **b** to pin that run, then set
`PGI2 / iloprost (nM)` = **50** and re-run. PAC-1 is strongly suppressed (peak
≈ 0.67 → 0.18 vs the grey baseline) — the cAMP/PKA arm antagonising integrin
activation, as prostacyclin / iloprost do in the paper.

*Caveat: reversibility here is a designed feature of the 2-state module, not a
surprising prediction. The paper's primary (ir)reversibility node, PI3K→Akt→Rap1b,
is absent (the model has only the cAMP/PKA arm of P2Y12); the model does **not**
robustly reproduce the paper's weak-(ADP)-vs-strong-(thrombin) reversibility
distinction, because its PARs internalise rather than latch; and disaggregation /
embolization are inter-cellular, out of single-cell reach. The model shows the
proposed cause (affinity reversal), not the observed effect.*

---

## Reset

Preset dropdown → **Defaults** clears everything back to schema defaults; the
"Δ from defaults" line under the form always tells you what's currently edited.
