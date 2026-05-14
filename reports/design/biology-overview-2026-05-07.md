# Platelet WCM — biology overview (v0.3.0)

> *Last updated: 2026-05-12 (Phase 4 — PI cycle replaces forced IP3; first-principles model)*

A single-cell, deterministic model of intracellular Ca²⁺ dynamics in a
resting / activating human platelet. Reproduces the IP3-mediated
calcium transient that initiates platelet activation, calibrated
against Dolan & Diamond 2014 *Biophys. J.* 106:2049–2060 Fig. 4.

## The biological story we tell

A platelet at rest holds cytosolic Ca²⁺ at ~100 nM against an
extracellular concentration of ~1.2 mM and a cytoplasmic-store
(dense tubular system, DTS) free concentration of ~250 µM, by balancing
pumps and channels across two membranes. Most of the cellular Ca²⁺ is
*bound*, not free: calmodulin and gelsolin buffer the cytosol;
calreticulin buffers the DTS lumen. On stimulation (collagen, thrombin,
ADP), receptors generate inositol-1,4,5-trisphosphate (IP3); IP3 opens
IP3 receptors on the DTS and Ca²⁺ floods the cytosol, where it drives
granule secretion, integrin activation and shape change. The transient
resolves over seconds-to-minutes via SERCA reuptake, plasma-membrane
extrusion (PMCA), and store-operated entry (Orai1/STIM1) refilling the
DTS from the extracellular space.

## The pathway, as wired in v0.2.7

```
              +-- P2X1 (PM, ATP-gated) ──────────────── + 1.2 mM Ca_ex
              │
              +-- SOCE: STIM1 (DTS) -→ Orai1 (PM) ────  + 1.2 mM Ca_ex
              │                                            │
              │             +-- CALR / CALR-P (DTS luminal buffer)
              │             │
[IP3 spike] -→ IP3R (DTS) -→ Ca_cyt rises --------→ +-- PMCA (PM) -→ Ca_ex
                                  │                    │
                                  +-- CaM ladder (buffer + activates PMCA)
                                  │
                                  +-- Gelsolin / GSN (cyt buffer; coarse-grained)
                                  │
                                  +-- SERCA (DTS) -→ refill DTS store
```

Eight mechanisms coupled through cytosolic Ca²⁺. ATP is consumed by
both ATPase pumps (SERCA, PMCA). P2X1 is gated by extracellular ATP
(released from dense granules during activation; CD39 clears it).

## Components, with the published kinetic model adopted for each

| Mechanism | Kinetic model | Source |
|---|---|---|
| **IP3R** | deYoung-Keizer 1992 / Li-Rinzel 1994 reduction: one slow inactivation variable `h` + quasi-steady activation `m∞(IP3, Ca)`; Po = m∞⁴ × h tetrameric cooperativity; Nernst flux with γ_IP3R = **0.175 pS** (calibrated, see dissertation-notes §3.1) | deYoung & Keizer 1992 *PNAS* 89:9895; Li & Rinzel 1994 *J Theor Biol* 166:461; Bezprozvanny 1991 / Mak & Foskett 1997 (γ range) |
| **SERCA** | 6-state E1/E2 enzymatic cycle (E2 ⇌ E1 ⇌ E1·Ca → E1P·Ca ⇌ E2P·Ca → E2P → E2); k_bind_f = 500 µM⁻²·s⁻¹ (Purvis 2008 halved in Phase 2 retune); 1 ATP / 2 Ca²⁺ | Dode 2002 (isoform 3b kinetics), Purvis 2008 (rate constants) |
| **PMCA** | 5-state CaM-coupled scheme: basal path (Ca²⁺ binding → extrusion, V_max = 5.5 s⁻¹) **plus** Ca₄·CaM-activated path (V_max = 30 s⁻¹, ~5.5× faster); 1 ATP / Ca²⁺ | Caride et al. 2007 Table 3 |
| **Calmodulin** | Two-lobe Ca²⁺-binding ladder: CaM_free ⇌ Ca₂·CaM ⇌ Ca₄·CaM; ~20 481 molecules; cytosolic Ca²⁺ buffer + PMCA activator | Caride et al. 2007 steps 6–7 |
| **Cytosolic buffer (gelsolin proxy)** | Coarse-grained 1:1 site binding; 800 000 effective Ca²⁺ sites at Kd = 1 µM, k_off = 100 s⁻¹ (fast equilibrium). Represents gelsolin + annexins + Ca·ATP combined | Burkhart 2012 + Yin & Stossel 1979 (copy range); Sage & Rink 1985 (buffer ratio target) |
| **CALR C-domain** | 508 100 low-affinity Ca²⁺ sites at Kd = 1 mM, k_off = 1 000 s⁻¹ (20 324 CALR × 25 sites). Fast equilibrium. The dominant DTS luminal buffer | Burkhart 2012; Vassilakos 1998; Baksh & Michalak 1991 |
| **CALR P-domain** | 20 324 high-affinity Ca²⁺ sites at Kd = 1 µM, k_off = 1 s⁻¹ (1 site per CALR). Slow release; provides ~20 k Ca²⁺ reserve during transient DTS depletion | Vassilakos 1998 |
| **HSP90B1 medium-affinity** | 40 000 sites at Kd = 2 µM, k_off = 1 s⁻¹ (slow). 10 000 HSP90B1 × 4 medium sites. Holds ~40 k Ca²⁺ during the transient (k_off matched to transient timescale). | Argon & Simen 1999; Burkhart 2012 |
| **HSP90B1 low-affinity** | 110 000 sites at Kd = 600 µM, k_off = 600 s⁻¹ (fast). 10 000 HSP90B1 × 11 low sites. ~32 k bound at rest. | Argon & Simen 1999 |
| **BiP / HSPA5** | 50 000 sites at Kd = 500 µM, k_off = 1 000 s⁻¹ (fast). 1 effective site per BiP. ~17 k bound at rest; ~25 % of ER Ca²⁺ store. | Lièvremont 1997; Burkhart 2012 |
| **CREC pool (CALU + RCN1 + RCN2)** | 60 000 aggregated sites at Kd = 1 mM, k_off = 500 s⁻¹ (fast). ~15 000 molecules × 4 sites. ~12 k bound at rest. | Honoré & Vorum 2000; Vorum 1998 |
| **STIM1 sensor cycle** | DTS-bound (Ca-loaded, inactive) ⇌ free monomer ⇌ dimer (active sensor); detailed-balance rate constants | Dolan 2014 + Hoover & Lewis 2011 |
| **Orai1 / SOCE** | Monod–Wyman–Changeux allosteric model: STIM1 dimers translocate into puncta (Hill function on cytosolic Ca²⁺), bind Orai1 tetramers cooperatively; channel opening as fraction of bound STIM1 | Hoover & Lewis 2011, Dolan 2014 eq. 2 + eq. 4 |
| **P2X1 ATP-gated channel** | 3-state coarse kinetics (Closed → Open → Desensitised → Closed); ionotropic trimeric channel; Ca²⁺ flux gated on extracellular Ca²⁺ availability and extracellular ATP forcing. 1 000 channels; γ_Ca = 0.6 fS (calibrated). Distinct from P2Y1 GPCR — see dissertation-notes §7.1 | Mahaut-Smith 2000/2004; Vial & Evans 2002; Hechler 2003 |
| **PI cycle / IP3 production** | **v0.3.0: model now produces IP3 from PI cycle dynamics**, replacing the v0.2.x forced curve. Coarse-grained scheme: Gq activity signal → PLCβ activation → PIP2 hydrolysis → IP3 + DAG. PIP2 resynthesis lumped; IP3 degradation lumped (3-kinase + 5-phosphatase). 5 new species (PIP2, DAG, PLCb_i, PLCb_a) + IP3 retained. Replaces `ip3_forcing_uM` with `gq_signal_uM` (the v0.4 receptor work will replace this with explicit GPCR cascades). | Mazet, Tindall, Gibbins & Fry 2020 *Sci. Rep.* 10:13889 |
| **Resting protein decay** | Exponential decay of all non-calcium-pathway proteins, t½ = 7 days. Operates on platelet-lifespan timescales — *inert on the 200 s transient horizon* — retained for v0.5+ multi-day-scope work | Burkhart et al. 2012 |

## Compartments and copy numbers

| Compartment | Volume | Key species & counts |
|---|---|---|
| Cytosol | 6 fL | 361 Ca²⁺ (100 nM free), 181 IP3 (50 nM), **20 481 CaM** (ladder of 3 states), **800 000 GSN sites** (727 k free / 73 k bound at rest), 10.8 M ATP |
| DTS (cytoplasmic store) | 0.258 fL (4.3 % of cyt) | 38 842 Ca²⁺ (250 µM free), **508 100 CALR C-domain sites** (406 k free / 102 k bound), **20 324 CALR P-domain sites** (~81 free / 20 k bound), 1 328 IP3R (gated via single `h` variable), 11 892 SERCA (6 sub-states), 4 265 STIM1 (3 sub-states) |
| Plasma membrane | (surface) | 769 PMCA (5 sub-states), 1 447 Orai1 monomers (≈ 360 tetrameric channels), **1 000 P2X1 channels** (all closed at rest; 3 sub-states) |
| Extracellular | infinite reservoir | fixed 1.2 mM Ca²⁺ |

Total cellular Ca²⁺ accounting at rest (including bound):

| Pool | Count | µM-equivalent |
|---|---|---|
| Cytosol free | 361 | 100 nM |
| Cytosol bound (CaM ladder) | ~1 280 | ~355 nM |
| Cytosol bound (GSN, coarse) | ~73 000 | ~20 µM |
| DTS free | 38 842 | 250 µM |
| DTS bound (STIM1) | ~3 800 | ~25 µM |
| DTS bound (CALR C-domain) | ~101 600 | ~656 µM |
| DTS bound (CALR P-domain) | ~20 200 | ~131 µM |

Resting bound:free ratio: cyt ≈ 200:1; DTS ≈ 3.2:1 (free + bound;
buffering is dominant) → ~73 % of DTS Ca²⁺ is buffer-bound at rest.

## What we can ask the model

Three biologically-distinct run conditions, exposed in the webapp and
on the CLI:

| Condition | What it tests |
|---|---|
| **IP3 transient (+Ca²⁺_ex)** | Canonical activation. IP3R-driven release, SERCA reuptake, SOCE refilling. Acceptance criterion: peak in 200–800 nM band. |
| **EDTA transient (no Ca²⁺_ex)** | Isolates IP3R contribution by removing extracellular Ca²⁺. SOCE inactive (no source); PM leak inactive. Compares against the +Ca_ex condition to test SOCE dependence (Dolan 2014 Fig. 4). |
| **Resting (no stimulus)** | IP3 stays at 50 nM baseline; no transient driven. Inspects long-time stability. |

Plus a Phase 3 driver that runs the +Ca_ex and EDTA conditions
back-to-back and produces the Dolan Fig. 4 comparison figure, and a
60 s-baseline + 240 s-IP3 transient plot (`plotCaBoundFree.py`) showing
free vs bound Ca²⁺ in both compartments.

## What v0.2.6 reproduces

Validated against Dolan 2014 Fig. 4 + Fig. 3B filtering criteria, after
the Phase 2 buffer biology commit (`7f4a9ffd`, 2026-05-11).

### Resting state (no IP3 stimulus, integrated to 600 s)

| Quantity | Dolan target | Model |
|---|---|---|
| Cytosolic [Ca²⁺] free | 100 nM | **109 nM** ✓ |
| DTS [Ca²⁺] free | 200–300 µM | **264 µM** ✓ |
| IP3 | 50 nM baseline | 50 nM ✓ |

For the first time in this project, the model has a **stable
biologically realistic resting fixed point**. Previous versions
(pre-Phase-2) drifted to cyt > 1 µM and DTS → 0 with no IP3 forcing,
because the DTS had no luminal buffer to absorb SERCA's pumping.

### Phase 3 transient response (200 s, ±Ca²⁺_ex)

![Phase 3 validation figure — 2026-05-11 (Phase 2.5, with P2X1)](/Users/steve/github/platelet-wcm/reports/figures/phase3-dolan-fig4-2026-05-11.png)

| Acceptance criterion | Result |
|---|---|
| Active (+Ca_ex): peak Ca²⁺_cyt > 200 nM | ✓ 493 nM |
| Active (−Ca_ex): peak Ca²⁺_cyt > 200 nM | ✓ 345 nM |
| Peak (+Ca_ex) in Dolan ±30 % band (315–585 nM) | ✓ 493 nM |
| Peak (−Ca_ex) in Dolan ±30 % band (192–358 nM) | ✓ 345 nM |
| SOCE differential: \|peak(+) − peak(−)\| ≥ 100 nM | ✓ **147 nM** |

**All 5 Dolan acceptance criteria pass** (Phase 2.5, 2026-05-11). The
SOCE-differential criterion was closed by adding P2X1: it's a +Ca_ex–
specific fast Ca²⁺ entry pathway that contributes ~150 nM to the peak
in +Ca_ex but is silent in −Ca_ex (no extracellular Ca²⁺ → no flux even
when the channel is open). This is mechanistically distinct from
Orai1/SOCE, which is too slow to contribute to the early peak.

The peak now also occurs **later** in the +Ca_ex run (t ≈ 100–150 s)
rather than the artificial t ≈ 1 s instant-peak of pre-Phase-2.5 — a
qualitatively more realistic transient with a sustained P2X1 + SOCE
plateau, then resolution.

### Free / bound Ca²⁺ during transient

![Free vs bound Ca²⁺ during 60 s baseline + 240 s IP3 transient](/Users/steve/github/platelet-wcm/reports/figures/ca-bound-free-2026-05-11.png)

Generated by `runscripts/manual/plotCaBoundFree.py`. Shows free Ca²⁺,
CaM-bound, GSN-bound, CALR C-domain, CALR P-domain, STIM1-bound, and
IP3 forcing curve all in one figure.

## What v0.2.7 does *not* yet model

The dissertation gap-catalogue is in `reports/dissertation-notes.md`;
this is the short version, ordered by impact.

| Biology | Status | Tracked as |
|---|---|---|
| ~~Other DTS luminal buffers~~ | **Closed Phase 3 / #25 (2026-05-12)** — HSP90B1 (M+L), BiP, CREC pool added. DTS bound:free ratio went from 3.2:1 → 5.8:1 (73% → 85% bound). Free DTS still drops to ~0.02 µM at peak (buffer release rate cannot match IP3R drain), but bound pool retention is now visibly biological. | Issue #25 closed |
| **Mitochondrial Ca²⁺ uniporter (MCU + mNCX)** | Not modelled. Captures Ca²⁺ during the spike, slowly releases over minutes. Three platelet MCU papers in `source-info/calcium-papers/` (Ajanel 2025, Ghatge 2026, Shehwar 2025). | Issue **#22** |
| **Dense granule Ca²⁺ store** | Not modelled. Platelet-specific high-concentration acidic Ca²⁺ store (50–100 mM total); NAADP/TPC-sensitive release. | dissertation-notes §7.2; v0.3+ |
| **Surface-Connected Canalicular System (SCS)** | Not modelled as a separate compartment. Doubles–triples the effective PM surface area for SOCE/PMCA. PM rate constants implicitly absorb this. | dissertation-notes §6.2; v0.3+ |
| **IP3R clustering / microdomain Ca²⁺** | Not modelled. Treated as a well-mixed population of 1 328 independent channels. | dissertation-notes §6.3; v0.3+ |
| **P2Y₁ / Gq / PLCβ upstream cascade** — endogenous IP3 production from ADP signalling | v0.2 uses a forced curve. Canonical replacement: Mazet et al. 2020. | Phase v0.3 (#9, #10) |
| **P2Y₁₂ / Gi / cAMP / PKA modulation** — negative regulation of activation | Not modelled | Phase v0.4 |
| **GPVI signalling cascade** — collagen-mediated activation pathway | Not modelled (Dunster 2015 reference available) | Phase v0.4 |
| **Granule release** (dense and α-granules) | Not modelled (cargo molecules in inventory but no exocytosis kinetics) | Phase v0.3 (#15) |
| **Gelsolin's cytoskeletal dual role** | Modelled only as a passive Ca²⁺ buffer; in real biology, Ca²⁺-bound gelsolin actively severs actin filaments during platelet shape change | dissertation-notes §7.5; v0.3+ |
| **Integrin αIIbβ3 inside-out signalling** | Not modelled | Phase v0.6 (optional) |
| **cAMP / cGMP / NO suppression at rest** | Not modelled (Kleppe 2018 reference available) | Future |

## Known calibration questions

Documented in `reports/dissertation-notes.md`:

- **§3.2** — SERCA cycle flux is still ~2–3× too high at rest vs the
  literature SERCA3b kinetic prediction (the Phase 2 halving brought
  it from 4–5× to 2–3×, but Purvis 2008's k_bind_f is still above
  Dode 2002's measured Vmax / Km would predict). Closing this further
  would let the cyt buffer ratio drop from the current ~200:1 toward
  Sage & Rink 1985's measured ~50:1.
- **§3.1** — γ_IP3R = 0.175 pS is a *calibration anchor* coupled to
  the SERCA rate constants and the cyt buffer load, not an
  independent measurement. Any of those changing requires γ to be
  re-derived.
- **§4.1** — Po = m∞⁴ × h chosen over Li-Rinzel original m∞³ × h
  for consistency with Dolan's tetrameric Po⁴ convention; switching
  would shift γ_IP3R by ~6.5×.

## Design philosophy

Reuse a validated published model where one exists; deviate only when
the published model is genuinely incompatible with the framework or
the biology. Major deviations from primary sources are catalogued in
the design doc + dissertation-notes, each with the deviation's value,
the reason, and a pointer to the lab-book entry that diagnosed it.
Examples:

- IP3R replaced from Sneyd-Dufour to deYoung-Keizer (Phase 1, #27) —
  the Sneyd-Dufour calibration regime (IP3 = 10 µM) didn't extrapolate
  to resting IP3 = 50 nM.
- γ_IP3R reduced from Zschauer 1988's 10 pS to 0.175 pS (Phase 4, #30
  → Phase 2) — bilayer measurements don't transfer to physiological
  conditions; calibrated against analytical SERCA cycle flux instead.
- CALR + cyt buffer added (Phase 2, #28) — Dolan's model was severely
  under-buffered on both sides; our model now reflects luminal +
  cytosolic biology with documented residual gaps.

The dissertation lead is *"reproduce Dolan where possible, deviate
honestly where necessary, document every deviation."*

## Primary references

- **Dolan & Diamond 2014** *Biophys. J.* 106:2049–2060 — the validation target
- **Purvis & Bhatt 2008** *PLoS Comp Biol* 4:e1000050 — kinetic parameterisation foundation
- **deYoung & Keizer 1992** *PNAS* 89:9895; **Li & Rinzel 1994** *J Theor Biol* 166:461 — IP3R model
- **Dode et al. 2002** *J Biol Chem* — SERCA3b kinetics
- **Caride et al. 2007** *J Biol Chem* — PMCA isoform 4b 5-state scheme
- **Vassilakos et al. 1998**; **Baksh & Michalak 1991** — CALR Ca²⁺-binding
- **Burkhart et al. 2012** *Blood* — platelet proteome reference (copy numbers)
- **Sage & Rink 1985** — platelet cytosolic Ca²⁺ buffering ratio measurement
- **Yin & Stossel 1979** — original platelet gelsolin characterisation
- **Hoover & Lewis 2011** *PNAS* — Orai/STIM CRAC channel MWC framework
- **Mazet, Tindall, Gibbins & Fry 2020** *Sci. Rep.* 10:13889 — canonical PI cycle reference for v0.3

All in `source-info/calcium-papers/`; per-value provenance in
`reports/data/calcium-data-provenance.md`. Comprehensive limitations
and assumptions for the dissertation write-up in
`reports/dissertation-notes.md`.
