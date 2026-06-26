# Current model truth table (thesis-sync reference)

**Purpose.** A one-page snapshot of the *as-built* model, used to synchronise the
thesis scaffold and supplements against the code before writing prose. The
cardinal rule: **the thesis must not contradict the actual model.** Compiled
2026-06-26 by reading the code, tests, and `docs/validation-and-regressions.md`
directly — not from memory or older lab books (several of which describe an
earlier model state). Keep this file current; treat it as the source of truth a
claim is checked against.

If a number below feeds a thesis figure legend or a results sentence, it carries
a `file:line` so it can be re-verified.

---

## 1. Version / file names

| Item | Current | Stale label to avoid |
|------|---------|----------------------|
| Calcium kinetics table | `reports/params/calcium-v0.6.toml` | ~~calcium-v0.5.toml~~ |
| Species inventory | `reports/params/species-v0.6.tsv` (83 species) | ~~species-v0.5.tsv~~, ~~63 species~~ |
| Loader default | `_DEFAULT_VERSION = 'v0.6'` (`_params_loader.py`, `_species_loader.py`) | — |
| Biological model "version" | **v0.6 line** (param-file version). The biology has advanced well past the June "v0.4.1 / frozen" state. | ~~v0.4.1~~, ~~"biology frozen"~~ |

The kinetics-as-data scaffold is still **calcium-only**; downstream modules keep
their constants in Python dataclasses, not TOML.

---

## 2. Built modules / processes

| Process | What it models | Touches Ca²⁺ ODE? |
|---------|----------------|-------------------|
| `CalciumDynamics` | GPCR (P2Y1/PAR1/PAR4/P2X1) → Gαq → PLCβ → PI-cycle → IP3R / SERCA / PMCA / SOCE / CaM / **MCU+NCLX** / **NCX** | **Yes — the core ODE** |
| `RestingDecay` | First-order decay of resting protein pool (one stochastic process) | No |
| `GranuleSecretion` | Dense/α-granule cargo release + **autocrine ADP** loop + ecto-NTPDase clearance | Feeds back via `step_inputs` (autocrine ADP → P2Y1) |
| `ThromboxaneSynthesis` | cPLA₂→COX-1→TXA₂ (Slice A) + **autocrine TXA₂→TP→Gq** (Slice B) | Slice B feeds back (TXA₂ → TP) |
| `IntegrinActivation` | αIIbβ3 resting⇌active inside-out switch; **PI3K/Akt→Rap1b arm built (#73)**; PAC-1 readout | No (terminal output) |

Plus the **P2Y12 / Gi / cAMP / PKA inhibitory axis** inside `CalciumDynamics`
(cAMP node, PKA brake on IP3R + integrin, VASP-P readout) and cAMP-raising drugs
(PGI₂ / forskolin / PDE3-block). NCX (`SLC8A1/3`) and NCLX are genuine model
mechanisms (`calcium_signalling.py:574`, `:789`) — the abstract may list them.

PKC is the hub of the feedback loops and three terminal outputs (secretion,
thromboxane, integrin).

---

## 3. MCU — the single most out-of-date thesis claim

**Current behaviour (coupling ON, the shipped default `mito_coupling_gain=1.0`):**
MCU knockout (`mcu_vmax_scale=0`) **REDUCES** agonist-evoked cytosolic Ca²⁺
peak and AUC, **matching the direction** of Ghatge 2026 and Ajanel 2025.

- Mechanism (`ip3r_relief_factor`, `calcium_signalling.py:817`): MCU uptake at
  mitochondria–DTS (MAM) contacts relieves the IP3R's Ca²⁺-dependent
  inactivation; the relief scales with functional MCU capacity and is gated by an
  activation function of cytosolic Ca²⁺ (engages only during the evoked
  transient, **zero at rest** → resting state preserved). KO → relief lost →
  reduced release → reduced cyt Ca²⁺; the fuller store then lowers SOCE
  indirectly (SOCE is *not* gated directly). RunConfig note: `run_config.py:115`.
- Capacity back-pressure (`mito_fill`, `:1469`) bounds matrix accumulation
  (Part 1).
- **Toggling coupling OFF** (`mito_coupling_gain=0`) recovers the old
  **buffer-only** behaviour, where KO *raises* the peak — i.e. the previous
  divergence is reproducible on demand and localises the missing biology.
- Test that pins this: `test_validation_targets.py:~202` asserts
  `ko.max() < wt.max()`, `ko.max() < 0.92*wt.max()`, `ko.sum() < wt.sum()`,
  and `ko_off.max() > wt.max()` (buffer-only sign flip).

**Honest caveats (do not overclaim the fix):** the coupling is a *model choice*
(not a measured constant); it lumps relief as a whole-flux scale (no spatial
microdomain / explicit h-gate); only the endpoints (WT, full KO) are calibrated,
not intermediate MCU levels; the magnitude is modest (~15–20 %) because losing
the MCU *buffer* partly offsets the lost *relief*. The model reproduces the
**direction** of the experimental effect, not a validated magnitude or a unique
mechanism (MCU's support of SOCE is an alternative not separately resolved).

**Numbers for the figure/legend:** to be filled from a fresh measured run (WT vs
KO peak + AUC, coupling on; and KO buffer-only) — the old draft's 436→651 nM
*raise* is the **decoupled** behaviour and must not be presented as the model's
result. The `figures/v0.5/perturbation-mcu-2026-06-10.png` asset shows the old
buffer-only result and **must be regenerated** for the new (coupled) narrative.

**Stale code comment (out of scope, flagged):** `calcium_signalling.py:800-802`
still says the model "is buffer-only and so RAISES the cyt peak on MCU loss" —
predates Part 2; only true with coupling off. Worth a one-line fix in a later
housekeeping pass (not part of this thesis/test rework).

---

## 4. Other drift the thesis must fix

| Thesis claim (old) | Current truth | Fix |
|--------------------|---------------|-----|
| PI3K→Akt→Rap1b arm "is absent" (Discussion) | **Built (#73)** — `akt_rap_step`, `Rap1b_GTP`, `rap1b_scale` knob | State it is built; keep honest scope (P2Y12 block slows rise / lowers stimulated-phase integrin; converges once ADP clears — *not* "reverses at high ADP") |
| `calcium-v0.5.toml` / `species-v0.5.tsv` (Methods) | v0.6 | Update names |
| "model is at v0.4.1; future work = granule release + P2Y12 arm" (Conclusion) | Both **built**; v0.6 line | Rewrite seed: future work is de-clamping + multi-assay validation + NO/cGMP + GPVI + metabolism |
| "21-test unit suite" / "5/5" as the acceptance gate | Being replaced by a behavioural **regression suite** (Task 2) | Describe the regression suite (resting equilibrium · Dolan ±Ca²⁺ transient · drug/KO responses), not a hard-coded count |

---

## 5. Claims — allowed vs not allowed

**Allowed (with caveats):**

- A platelet Ca²⁺ pathway can be modelled mechanistically in a whole-cell
  framework and validated against Dolan & Diamond 2014 (±extracellular Ca²⁺),
  with no hand-fitted forcing downstream of agonist input.
- The architecture is **demonstrably extensible**: secretion, thromboxane,
  integrin, and a P2Y12/cAMP/PKA inhibitory arm compose over the shared
  mass-balanced state without re-fitting the calcium core. *(This is now
  DEMONSTRATED, not merely asserted — an upgrade over the 2026-06-09 review,
  where only one pathway existed.)*
- MCU knockout reduces evoked cytosolic Ca²⁺ (direction matches Ghatge/Ajanel).
- Drug-like knockouts behave correctly in *direction*: aspirin (COX-1) abolishes
  TXA₂; clopidogrel (P2Y12 block) lowers integrin PAC-1 via the cAMP/PKA brake.

**Not allowed (contradicted by code or explicitly disclaimed):**

- ~~"No published whole-cell models of any mammalian cells."~~ → qualify:
  no published **mechanistic / dynamical / mass-balanced** whole-cell model of a
  mammalian cell (Karr→Covert sense); note the Covert lab H1-hESC effort as a
  first step. Don't fly a primacy flag.
- ~~"With further work a whole-cell model is feasible."~~ → tiered framing:
  **major pathway validated (proved) · architecture extensible (now shown) ·
  full WCM not proved/not claimed**; the binding constraint is data, not
  computation.
- ~~"Every reaction carries a primary-source rate constant."~~ → true for the
  **calcium core** (Dolan/Purvis/Caride/Hoover/Burkhart); the downstream modules
  (integrin kinetics, PKC, PKA brake) are **model choices**. Scope the claim.
- Quantitative prediction of how upstream modulation (PKC, P2Y12/cAMP) changes
  *free cytosolic Ca²⁺* — blocked by the **Ca²⁺ clamp**; those effects are shown
  on functional outputs (integrin, secretion, VASP-P), not the Ca²⁺ trace.

---

## 6. Validation status (one line per subsystem)

From `docs/validation-and-regressions.md`. Dolan "5/5" is now a **regression
invariant on the Ca²⁺ core**, not a whole-model correctness claim ("passes by
construction" + low discriminating power under the Ca²⁺ clamp).

| Subsystem | Status |
|-----------|--------|
| Ca²⁺ core (Dolan Fig. 4 ±Ca_ex) | ✓ regression invariant |
| PI cycle / IP₃ | ✓ direction + band |
| PKC feedback | ◐ direction + band |
| Granule secretion | ◐ shape (not kinetics) |
| Thromboxane | ✓ direction (aspirin abolishes) |
| Integrin αIIbβ3 / PAC-1 | ✓ graded dose-response + drug shifts (relative) |
| Inhibitory axis (P2Y12/cAMP/PKA) | ✓ VASP/PAC-1; ○ PGI₂ effect on Ca²⁺ (future) |
| MCU | ◐ direction matches; magnitude modest, mechanism one of several |
