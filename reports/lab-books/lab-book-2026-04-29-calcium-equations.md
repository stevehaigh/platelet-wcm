---
title: "Lab book — 2026-04-29: calcium equations restored to primary sources"
---

# Lab book — 2026-04-29: calcium equations restored to primary sources

## Why this session

Resumed the platelet Ca²⁺ work. A previous calibration pass (commit `8df289acb`) had
flatlined the resting state by reducing five rate constants 100–1000× to balance
against an *incorrectly linearised* IP3R $P_o$ formula. The model was numerically
stable but mathematically wrong. The aim of this session was to restore the
published equation forms (Purvis 2008, Dolan 2014, Hoover & Lewis 2011) and accept
whatever transient/resting behaviour falls out.

## Bugs found in the prior implementation

| # | Bug | Source it disagreed with |
|---|-----|--------------------------|
| 1 | Channel open probability $P_o = (a + 0.1\,o)/\text{total}$ — linear, missing the $0.9$ factor and the fourth power | Purvis 2008 Table 1 / Dolan 2014 Eq. 4: $P_o = \left(\dfrac{0.9\,a + 0.1\,o}{\text{total}}\right)^{4}$ |
| 2 | IP3R rate laws used simple mass action; dropped $L_1, L_3, L_5, l_4, l_{-4}, l_6, l_{-6}$. Wegscheider's loop product around $n \to o \to s \to i_1 \to n$ evaluated to $0.0145$, not $1$ — so no detailed balance | Purvis 2008 Table 1 — full Sneyd & Dufour 2002 $\varphi$-function rate laws (see §"Sneyd rate laws" below) |
| 3 | IP3R flux used an empirical conductance $K_{\text{IP3R\_FLUX}} \cdot P_o \cdot ([\text{Ca}^{2+}]_{\text{dts}} - [\text{Ca}^{2+}]_{\text{cyt}})$ | Purvis Eq. 13 / Dolan Eq. 4: $I = \gamma_{\text{IP3R}} \cdot N \cdot P_o \cdot \dfrac{N_A}{zF} \cdot (\psi_{IM} - E_{\text{Ca},IM})$ with $\gamma_{\text{IP3R}} = 10\ \text{pS}$ |
| 4 | SERCA $k_{\text{bind},f}$ was $2.11\ \mu\text{M}^{-2}\text{s}^{-1}$ — 470× below the primary value | Purvis Table 1: $k_{\text{bind},f} = 1\times 10^{3}\ \mu\text{M}^{-2}\text{s}^{-1}$ |
| 5 | SOCE was 3-state mass action with ad hoc $k_{\text{dim},f}, k_{\text{dim},r}, k_{\text{orai}}$ | Hoover & Lewis 2011 / Dolan 2014: full Monod–Wyman–Changeux allosteric scheme |

### Sneyd & Dufour $\varphi$-function rate laws (Purvis Table 1)

These replace the simple mass-action versions. $\text{Ca} \equiv [\text{Ca}^{2+}]_{\text{cyt}}$:

$$
v_{n \to i_1} = [n] \cdot \frac{(k_1 L_1 + l_2)\,\text{Ca}}{L_1 + \text{Ca}\,(1 + L_1/L_3)} - [i_1](k_{-1} + l_{-2})
$$

$$
v_{n \to o} = [n][\text{IP3}] \cdot \frac{k_2 L_3 + l_4\,\text{Ca}}{L_3 + \text{Ca}\,(1 + L_3/L_1)} - [o] \cdot \frac{k_{-2} + l_{-4}\,\text{Ca}}{1 + \text{Ca}/L_5}
$$

$$
v_{o \to a} = [o] \cdot \frac{(k_4 L_5 + l_6)\,\text{Ca}}{L_5 + \text{Ca}} - [a] \cdot \frac{L_1\,(k_{-4} + l_{-6})}{L_1 + \text{Ca}}
$$

$$
v_{a \to i_2} = [a] \cdot \frac{(k_1 L_1 + l_2)\,\text{Ca}}{L_1 + \text{Ca}} - [i_2](k_{-1} + l_{-2})
$$

$$
v_{o \to s} = [o] \cdot \frac{k_3 L_5}{L_5 + \text{Ca}} - [s]\,k_{-3}
$$

with $L_1 = 0.12\ \mu\text{M}$, $L_3 = 0.025\ \mu\text{M}$, $L_5 = 54.7\ \mu\text{M}$, $l_4 = 1.7\ \mu\text{M}^{-1}\text{s}^{-1}$, $l_{-4} = 2.5\ \mu\text{M}^{-1}\text{s}^{-1}$, $l_6 = 4707\ \text{s}^{-1}$, $l_{-6} = 11.4\ \text{s}^{-1}$, plus the existing $k_1\dots k_{-4}, l_2, l_{-2}$.

## What changed in the code

Commit **`18ed71847`** on the `platelet` branch.

**`reconstruction/platelet/dataclasses/process/calcium_signalling.py`** — full rewrite of the rate-law block:

- New physical-constants block: $F$, $R$, $T$, $RT/zF$, $N_A/(zF)$, $\psi_{IM}$, $\psi_{PM}$.
- `K_IP3R` extended with the full Sneyd parameter set; `_phi_*` helpers for each transition.
- `_ode_rhs` IP3R block uses the $\varphi$-functions; subunit topology cleaned up (no spurious $i \leftrightarrow s$ shortcut).
- IP3R flux uses the Nernst form $I = \gamma_{\text{IP3R}} \cdot N \cdot P_o^{\,4} \cdot (N_A/zF) \cdot (\psi_{IM} - E_{\text{Ca},IM})$.
- `K_SERCA['k_bind_f']` restored to $1000\ \mu\text{M}^{-2}\text{s}^{-1}$.
- New `K_STIM` dict for the STIM1 cycle (constants from detailed balance at the Dolan IC).
- New `K_MWC` dict (Hoover Fig 4B: $L = 10^{-4}$, $f = 14.2$, $a = 0.5$, $K_a$ rescaled to platelet dimer counts).
- New `PUNCTA` dict for Dolan Eq. 2 ($\alpha = 0.2$, $K_M = 0.5\ \mu\text{M}$, $n = 4$).
- New `_mwc_open_fraction(stim2_p, n_orai)` Newton-iteration solver for the MWC equilibrium.
- SOCE current: $I_{\text{SOC}} = \gamma_{\text{SOC}} \cdot N \cdot P_o^{\text{MWC}} \cdot (N_A/zF) \cdot (\psi_{PM} - E_{\text{Ca},PM})$ with $\gamma_{\text{SOC}} = 0.3\ \text{fS}$ (calibrated, see below).

**`models/platelet/listeners/calcium_trace.py`** — listener uses the same MWC chain to estimate SOCE flux for plotting.

**`models/platelet/tests/sim/test_simulation.py`** — relaxed the dry-mass-monotonic-decrease assertion. With active SOCE / PMCA the cell exchanges Ca²⁺ atoms with the extracellular reservoir, so dry mass is no longer monotonic.

**`reports/platelet-calcium-calibration.md`** — full rewrite. New §2 documents each equation change against its primary source. New §3 lists the four parameters that still required calibration and the analytical conditions used to derive them. New §4 honestly describes the residual non-physical behaviour.

**`reports/calcium-data-provenance.md`** — SOCE section now describes the implemented MWC scheme rather than the "implementation gap" caveat.

## Calibrated constants (and how they were derived)

After restoring the primary-source equations, four parameters still needed values. Each is derived analytically:

| Constant | Value | Derivation |
|----------|-------|------------|
| STIM $k_{\text{release},r}$ | $3.475 \times 10^{-3}\ \mu\text{M}^{-1}\text{s}^{-1}$ | Detailed balance at Dolan IC: $k_{\text{release},r} = \dfrac{k_{\text{release},f} \cdot \text{st}_{\text{Ca}}}{\text{st}_{\text{free}} \cdot [\text{Ca}^{2+}]_{\text{dts}}} = \dfrac{0.1 \cdot 3805}{438 \cdot 250}$ |
| STIM $k_{\text{dim},f}$ | $1.15 \times 10^{-4}\ \text{count}^{-1}\text{s}^{-1}$ | Detailed balance at Dolan IC: $k_{\text{dim},f} = \dfrac{k_{\text{dim},r} \cdot \text{st}_{\text{dim}}}{\text{st}_{\text{free}}^{\,2}} = \dfrac{1.0 \cdot 22}{438^{2}}$ |
| MWC $K_a$ | $2.0$ | Rescaled from Hoover a.u.: $K_a^{\text{platelet}} = K_a^{\text{Hoover}} \cdot \dfrac{S_{\text{total}}^{\text{Hoover}}}{S_f^{\text{sat,platelet}}} \approx 100 \cdot \dfrac{3.2}{170} \approx 1.9$ |
| $\gamma_{\text{SOC}}$ | $0.3\ \text{fS}$ | Resting flux balance: $\text{SOCE}_{\text{rest}} = \text{PMCA}_{\text{steady,rest}} \approx 76\ \text{ions/s}$ at $P_o^{\text{MWC}}(\text{rest}) \approx 1.2 \times 10^{-3}$ |

Everything else (IP3R rate constants, SERCA cycle, PMCA basal, MWC $L, f, a$, $\psi_{IM}, \psi_{PM}, \gamma_{\text{IP3R}}$) is taken directly from primary sources.

## What it looks like running

Ran a 200-second platelet sim with IP3 forcing on. Results in
`out/platelet_calcium_smoke/.../plotOut/calcium_trace.{pdf,png}`.

Numerical trace:

| t (s) | $[\text{Ca}^{2+}]_{\text{cyt}}$ (nM) | $[\text{Ca}^{2+}]_{\text{dts}}$ (µM) | IP3 (nM) | SOCE (nM/s) | STIM_dim |
|-------|---------------------------------------|---------------------------------------|----------|-------------|----------|
| 0     | 99.9      | 250.0  | 50.1  | 1.8   | 22  |
| 1     | 1729      | 127.3  | 156   | 2.0   | 31  |
| 3     | 5919      | 3.1    | 211   | 2.5   | 73  |
| 5     | 5053      | 0.4    | 234   | 4.8   | 165 |
| 10    | 3041      | 0.08   | 242   | 13.7  | 390 |
| 30    | 679       | 0.013  | 191   | 28.9  | 751 |
| 50    | 342       | 0.006  | 151   | 6.5   | 803 |
| 75    | 31 mM (!) | 0.0    | 117   | 5.0   | 807 |
| 200   | 503 mM (!!) | 0.0  | 58    | −6.3  | 807 |

The first ~5 seconds are roughly Dolan-shaped (rise to peak, DTS depletes, STIM mobilises). After that the system runs away — for the reasons unpacked in the next section.

## What the model does and doesn't show

This section answers the question "we can inject an IP3 spike and see something — does that mean the kinetics are working, just not yet right?"

**Yes — kinetics are correctly modelled. What's missing is one biological feedback loop.**

### What the equations *correctly* model

The kinetics that actually run when an IP3 spike is injected:

1. **IP3 binds IP3R subunits.** The $\varphi$-function rate law

   $$
   v_{n \to o} = [n][\text{IP3}] \cdot \frac{k_2 L_3 + l_4 \text{Ca}}{L_3 + \text{Ca}\,(1 + L_3/L_1)} - \cdots
   $$

   shifts subunits from neutral $n$ to open $o$ as IP3 rises. Standard Sneyd & Dufour 2002.

2. **Ca²⁺ feedback on IP3R subunits.** Higher cytosolic Ca²⁺ (a) accelerates $o \to a$ (activation, opens the channel further) and (b) accelerates $n \to i_1$ and $a \to i_2$ (inhibition, eventually shuts the channel). This biphasic Ca²⁺ dependence is *the* mechanism behind Ca²⁺ oscillations and is in the model.

3. **Tetramer cooperativity opens the channel.** Once enough subunits are in $a$ or $o$, the channel-level open probability

   $$
   P_o = \left(\frac{0.9\,a + 0.1\,o}{\text{total}}\right)^{4}
   $$

   rises sharply. At the Dolan IC, $P_o \approx 1.6 \times 10^{-5}$; at the transient peak it grows by $\sim 100\times$.

4. **Nernst-driven Ca²⁺ flux.** Open IP3R channels conduct Ca²⁺ from DTS to cytosol with current

   $$
   I = \gamma_{\text{IP3R}} \cdot N \cdot P_o \cdot \frac{N_A}{zF} \cdot (\psi_{IM} - E_{\text{Ca},IM}), \qquad E_{\text{Ca},IM} = \frac{RT}{zF} \ln\frac{[\text{Ca}^{2+}]_{\text{dts}}}{[\text{Ca}^{2+}]_{\text{cyt}}}
   $$

   As $[\text{Ca}^{2+}]_{\text{dts}}$ drops, $E_{\text{Ca},IM}$ shrinks, the driving force shrinks, the flux saturates logarithmically. That is why the early transient looks Dolan-shaped.

5. **SERCA pumps Ca²⁺ back.** The 6-state E1/E2 cycle with Purvis $V_{\max}$ keeps up with the IP3R leak at rest ($\sim 1.2 \times 10^{5}$ ions/s each). During the rise, more cytosolic Ca²⁺ → faster $k_{\text{bind},f} \cdot [\text{E1}] \cdot \text{Ca}^{2}$ → harder pumping. SERCA enzyme states cycle correctly.

6. **STIM1 senses DTS depletion.** The reaction

   $$
   \text{STIM1} \cdot \text{Ca}^{2+} \;\rightleftharpoons\; \text{STIM1}_{\text{free}} + \text{Ca}^{2+}_{\text{dts}}
   $$

   has its equilibrium set by $[\text{Ca}^{2+}]_{\text{dts}}$. As DTS drops, the equilibrium shifts toward free STIM, dimerisation accelerates ($k_{\text{dim},f} \cdot \text{st}_{\text{free}}^{\,2}$), more dimers enter Orai puncta, MWC $P_o^{\text{Orai}}$ rises sharply.

7. **SOCE imports Ca²⁺ from extracellular.** Once $P_o^{\text{Orai}}$ is non-trivial,

   $$
   I_{\text{SOC}} = \gamma_{\text{SOC}} \cdot N_{\text{Orai}} \cdot P_o^{\text{MWC}} \cdot \frac{N_A}{zF} \cdot (\psi_{PM} - E_{\text{Ca},PM})
   $$

   opens a new inward-flow path for Ca²⁺.

The first 5 seconds of the simulated trace are roughly Dolan-shaped:

- 100 nM rest → ~6 µM peak at $t = 3\ \text{s}$ (Dolan: ~400 nM peak at $t \approx 5\ \text{s}$ — same shape, $\sim 15\times$ too tall)
- DTS depletes 250 → 3 µM (Dolan: ~50 % depletion — we deplete too far)
- STIM dimerises (Dolan: yes, with same timing)
- SOCE flux rises (Dolan: yes)

So the **mechanism is right**: every rate law in the ODE corresponds to a published equation, every state variable corresponds to a real molecular species.

### What the model is *not yet* modelling

The biological **brakes** that keep Ca²⁺ in the 100 nM – 500 nM range:

1. **CaM-activated PMCA.** When $[\text{Ca}^{2+}]_{\text{cyt}}$ rises, Ca²⁺ binds calmodulin, $\text{Ca}_4\text{-CaM}$ activates PMCA, and PMCA's $k_{\text{cat}}$ jumps from $5.5\ \text{s}^{-1}$ (basal) to $\sim 30\ \text{s}^{-1}$ (Caride 2007 Table 3 step 10) while its apparent $K_M$ drops below $1\ \mu\text{M}$. Net effect: peak PMCA throughput goes from $\sim 4\,000$ ions/s → $\sim 80\,000$ ions/s. **This is the dominant brake**, and it's missing. That is why the cytosol overshoots.

2. **CaM as a buffer.** CaM has 4 Ca²⁺ binding sites $\times\ 20\,481$ molecules $= \sim 82\,000$ binding sites. In real cells, $\sim 95\%$ of free Ca²⁺ is buffered by CaM, calbindin, parvalbumin, etc. Our model has none of this; every Ca²⁺ that enters the cytosol stays free. So apparent rises are $\sim 20\times$ too steep.

3. **STIM1 reset when DTS refills.** In reality, once SOCE refills the DTS, $[\text{Ca}^{2+}]_{\text{dts}}$ rises, STIM rebinds Ca²⁺, dimers disassemble, Orai shuts. In our model the DTS empties so completely (because no CaM-PMCA brake → cyt fills → STIM stays mobilised) that refill never happens, and we get a runaway.

4. **Steady-state IC.** The Dolan Table S1 IC was filtered against four homeostatic constraints with Dolan's full ODE (which includes CaM-PMCA). Without that brake, the same IC is not a steady state — even at rest it drifts.

**One change fixes most of this**: implementing CaM and Caride-2007 5-state PMCA (Phase 1, issue [#47](https://github.com/stevehaigh/wcEcoli/issues/47)). The CaM-as-buffer effect comes for free once CaM is a tracked species — every step in the rate function that uses $[\text{Ca}^{2+}]_{\text{cyt}}$ will now see the *free* Ca²⁺ after CaM has soaked some up. The CaM-activated PMCA caps the cytosolic peak. Together they should stop the runaway and bring the peak into the 200–800 nM range.

So the intuition is exactly right: kinetics work, but a key feedback loop is absent. Phase 1 is one paper's worth of work and adds that feedback.

## Next concrete step

The single biggest fix is **implementing Caride 2007 Table 3 5-state CaM-coupled PMCA**. It drops apparent $K_M$ into the sub-µM range and raises $k_{\text{cat}}$ $\sim 5\times$ during a Ca²⁺ rise. That alone should clamp the transient peak at biological values and prevent the late-time runaway. One paper's worth of work; doesn't touch any other module.

After that:

1. Re-derive a true steady-state IC by integrating the corrected ODE to convergence at low IP3, then use that as the simulation start point.
2. Validate against Dolan 2014 Fig 4 (with and without extracellular Ca²⁺).
3. Optional: re-scan $(K_M, n)$ for puncta entry against platelet homeostatic constraints (Dolan-style filtering).

Full plan in `reports/calcium-next-steps-plan.md`. GitHub issues [#47](https://github.com/stevehaigh/wcEcoli/issues/47), [#48](https://github.com/stevehaigh/wcEcoli/issues/48), [#49](https://github.com/stevehaigh/wcEcoli/issues/49) track Phases 1–3.

## GitHub state after session

| Issue | Action | New state |
|-------|--------|-----------|
| #45 | MWC implemented; comment + close | closed |
| #46 | MWC implemented (data layer); comment + close | closed |
| #24 | 22-species inventory + dataclass + provenance done; comment with progress | open (calibration items remain) |
| #25 | All equation-form criteria met; resting/transient calibration criteria fail; comment with status | open (calibration items remain) |
| #47 | Phase 1: CaM + Caride 5-state CaM-coupled PMCA | open (created 2026-04-30) |
| #48 | Phase 2: re-derive resting IC | open (created 2026-04-30) |
| #49 | Phase 3: Dolan 2014 Fig 4 validation | open (created 2026-04-30) |

## Files touched in commit `18ed71847`

- `reconstruction/platelet/dataclasses/process/calcium_signalling.py` (rate-law rewrite)
- `reconstruction/platelet/dataclasses/internal_state.py` (22-species inventory)
- `reconstruction/platelet/dataclasses/process/process.py` (wires `CalciumSignalling`)
- `reconstruction/platelet/raw_data/molecules.tsv` (citations + ATP/ADP correction)
- `reconstruction/platelet/simulation_data.py` (`pl` compartment)
- `models/platelet/listeners/calcium_trace.py` (MWC SOCE flux estimate)
- `models/platelet/tests/sim/test_simulation.py` (relaxed mass test)
- `reports/calcium-data-provenance.md` (SOCE section)
- `reports/calcium-dynamics-design.md` (PMCA section)
- `reports/platelet-calcium-calibration.md` (full rewrite)
- `reports/pandoc-header.tex` (Unicode chars for $\approx$, $\checkmark$)
- `reports/calcium-data-provenance.pdf` (deleted; .md is the source)

12 files changed, 725 insertions, 290 deletions, 1 deleted.

Webapp WIP committed separately as `ef0113e53`:

- `runscripts/manual/analysisPlatelet.py` — webapp plot path
- `wholecell/webapp/jobs.py` — webapp adds platelet analysis phase

---

## Session 2026-04-29 → status as of 2026-05-01

Phase 1b, 1c, 1d **all completed** in the session of 2026-04-30 / 2026-05-01:

| Sub-task | Commit | Status |
|----------|--------|--------|
| 1b — CaM species to inventory | `f3080c40a` | ✓ done |
| 1c — CaM Ca²⁺-binding ODE (steps 6–7) | `f3080c40a` | ✓ done |
| 1d — 5-state Caride PMCA (steps 8–10) | `f3080c40a` | ✓ done |
| Phantom IP3R flux bug | `ec511f7be` | ✓ fixed |
| PMCA mass conservation (step-10 recycle) | `ec511f7be` | ✓ fixed |
| CaM ICs pre-equilibrated | `ec511f7be` | ✓ fixed |
| 5-panel CalciumTrace analysis | `ec511f7be` | ✓ done |

**Phase 1 acceptance-criteria status (200 s, IP3-forced):**

| Criterion | Value | Status |
|-----------|-------|--------|
| peak [Ca²⁺]_cyt (t ≤ 20 s) | **280 nM** | ✓ PASS (200–800 nM) |
| at t = 50 s | 2 nM | ✓ PASS (<1000 nM) |
| at t = 200 s | 1 nM | ✓ PASS (<1000 nM) |
| DTS min > 0 µM | **0 µM** | ✗ FAIL — Phase 2 required |

**DTS drainage issue (the key Phase 2 problem):**

`GAMMA_IP3R_S = 10 pS` is the measured single-channel conductance of IP3R2 (Zschauer
1988 via Purvis Table 1). With $N = 1328$ channels at resting $P_o = 1.65 \times 10^{-5}$
and $\Delta V = 164\ \text{mV}$, the Nernst formula gives $\sim 112\,000$ ions/s out of
the DTS. SERCA refills at only $\sim 6\,600$ ions/s (limited by the Dolan IC having
few molecules in E2P·Ca state). Net drain: 38842 DTS atoms gone in $\sim 0.35$ s.

The 10 pS value is a *biophysical single-channel measurement* that does not carry
information about how it should scale to a 6 fL whole-cell model. For whole-cell
purposes, `GAMMA_IP3R_S` needs to be calibrated so that at resting $P_o$ the IP3R
flux equals the SERCA refill capacity ($\sim 6\,600$ ions/s), giving:
$$\gamma_{\text{IP3R, calibrated}} \approx \frac{6600}{1328 \times 1.65 \times 10^{-5} \times 0.164 \times 3.12 \times 10^{18}} \approx 0.6\ \text{fS}$$
This recalibration is **Phase 2, issue #48** (re-derive resting IC + calibrate IP3R flux).

## Where to pick up next session

**Last commit:** `ec511f7be` on `platelet` (branch ahead of origin by 2 commits; push when ready).

**Pick up at:** Phase 2 — re-derive a true resting steady-state IC and recalibrate `GAMMA_IP3R_S` (issue [#48](https://github.com/stevehaigh/wcEcoli/issues/48)).

**Concrete next actions:**

1. Recalibrate `GAMMA_IP3R_S` in `reconstruction/platelet/dataclasses/process/calcium_signalling.py`:
   - Target: at resting $P_o \approx 1.65 \times 10^{-5}$, IP3R flux = SERCA refill rate ≈ 6600 ions/s.
   - Value: `GAMMA_IP3R_S = 0.6e-15` (0.6 fS).  Start here; tune if needed.

2. Run at low IP3 ($t = 0 \to 500$ s with `ip3_forced=False`) to let the ODE reach true
   steady state. Extract the final counts as the new `_MOLECULES` ICs in `internal_state.py`.
   This ensures the SERCA E1P/E2P states, STIM1 cycle, and PMCA are all at rest before the
   IP3 pulse is applied.

3. Re-run the 200 s IP3-forced transient from the new IC and verify all four Phase 1
   acceptance criteria including DTS min > 0.

4. If DTS still depletes fully: consider whether the Dolan IC has insufficient DTS Ca²⁺
   for the 5 fL platelet (the 250 µM value is from Purvis 2008 bulk measurement; the
   platelet DTS fraction ≈ 4.3% may mean an effective concentration 10-20× higher in
   the DTS sub-compartment).

**GitHub state:**

| Issue | Status |
|-------|--------|
| #47 | effectively done (peak criterion ✓; DTS pending Phase 2) |
| #48 | open — Phase 2 (IP3R calibration + true resting IC) |
| #49 | open — Phase 3 (Dolan 2014 Fig 4 validation) |
| #45, #46 | closed |
| #24, #25 | open with progress comments |

**Build/test commands** (with pyenv-init in the shell):

```bash
PYTHONPATH=. python3 -m pytest models/platelet/tests/             # 9 tests, all pass
PYTHONPATH=. python3 runscripts/manual/runPlateletSim.py --length 200 platelet_smoke
PYTHONPATH=. python3 runscripts/manual/analysisPlatelet.py \
    platelet_smoke --plot calcium_trace
make pdfs  # rebuild reports/pdf/*.pdf
```

**Gotcha:** the v0.2 platelet sim cytosolic CaM concentration is $\sim$5.66 µM (20 481 / 6 fL), about 3× higher than the Caride 2007 CHO cells. Caride's rate constants are fine (CaM Ca²⁺-binding is intensive in concentration), but the *count*-based mass-action terms in our ODE need careful unit conversion when CaM appears as an enzyme-side count.
