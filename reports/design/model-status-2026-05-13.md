---
title: "Platelet WCM — Model status (v0.4.1)"
author: "Steve Haigh"
date: "2026-05-13"
---

# Platelet WCM — Model status (v0.4.1)

End-to-end mechanistic cascade from physiological agonist (thrombin
1 nM / ADP 10 µM) to cytosolic Ca²⁺ peak. **No hand-fitted forcing
remains in the calcium pathway.** Twelve coupled mechanisms, validated
against Dolan & Diamond 2014 Fig. 4 at a resting state of cyt = 100 nM,
DTS = 250 µM, IP3 = 50 nM.

---

## Overview diagram

See Platelet WCM v0.4.1 — Ca²⁺ pathway overview

Portrait-oriented as a starting point for a hand-drawn BioRender redraw.
Mermaid source: `reports/figures/v0.5/model-status-2026-05-13.mmd`.

---

## 1. Receptors — P2Y1 + PAR1 + PAR4

**Biology.** Platelet activation is initiated
by three Gq-coupled GPCRs in the plasma membrane. **P2Y1** binds ADP
(K_d ~1 µM) reversibly and dominates the early Ca²⁺ response to ADP.
**PAR1** is cleaved by thrombin at high affinity (EC₅₀ ~0.5 nM); the
N-terminal cleavage exposes a tethered ligand and activation is
essentially irreversible, terminated only by receptor internalisation
over minutes. **PAR4** is the low-affinity thrombin receptor (EC₅₀
~5 nM) and gives the sustained response. All three couple to Gαq —
P2Y₁₂ (Gi-coupled, inhibitory) and GPVI (PLCγ2 spine) are deliberately
out of scope.

**Implementation.** Each receptor is two species — `inactive[pl]`,
`active[pl]` — with first-order activation kinetics in
`calcium_signalling.py` (`K_P2Y1`, `K_PAR1`, `K_PAR4` blocks). PAR1/PAR4
additionally have an `internalized[pl]` sink: as of **v0.4.1** the
cleavage step is *one-way*, with active receptors decaying into the
internalised pool rather than recycling. Initial counts come from
Coller 1995 (P2Y1 = 150, PAR1 = 2 500, PAR4 = 500) and are loaded
through `internal_state.py`. The agonist concentrations are supplied
either as constants in the runscript or via the webapp Configure tab.

## 2. Gαq cycle

**Biology.** The heterotrimeric G-protein cycle converts receptor
activation into a downstream signal. Resting Gαq is GDP-bound and
sequestered with Gβγ; agonist-activated receptors catalyse GDP→GTP
exchange, freeing Gαq-GTP to activate PLCβ. Gαq has intrinsic GTPase
activity (~slow), greatly accelerated by RGS proteins, which terminates
the signal. Total platelet Gαq is ~5 000 molecules (Mazet, Tindall,
Gibbins & Fry 2020); at rest about 2 % is active (Gq_active ≈ 100).

**Implementation.** A single dynamic species `Gq_active[c]` with
inactive Gq implicit as `N_GQ_TOTAL − Gq_active`. The rate law has two
terms: a basal exchange rate and a receptor-catalysed activation
proportional to the summed active receptor count, balanced by
first-order RGS-accelerated GTPase deactivation. Parameters live in
the `K_GQ` block of `calcium_signalling.py` and are sourced from
Mazet 2020 Table 1.

## 3. PLCβ activation

**Biology.** Phospholipase C β is the enzyme that hydrolyses PIP₂ to
IP₃ + DAG. It exists in inactive and Gαq-bound (active) forms. PLCβ
activation is reversible (no proteolytic step), with the active pool
proportional to Gαq-GTP at quasi-steady state. Platelet PLCβ copy
number is small (~857 molecules; Burkhart 2012), but turnover is fast
enough to drive the entire IP₃ transient.

**Implementation.** Two species (`PLCb_inactive[c]`,
`PLCb_active[c]`) coupled by `K_PLCB['k_act']` (Gq-dependent) and
`K_PLCB['k_inact']` (first-order). The catalytic step
(PIP₂ hydrolysis) sits in the same block and is driven by `PLCb_active`
and the substrate count, producing IP₃ and DAG stoichiometrically.

## 4. PI cycle — PIP₂ → IP₃ + DAG

**Biology.** The phosphoinositide cycle is the substrate side of the
IP₃/DAG generating reaction. PIP₂ (PI-4,5-bisphosphate) is the
membrane-tethered phosphoinositide that PLCβ cleaves; the products are
the soluble messenger IP₃ (which diffuses to the DTS) and the
membrane-tethered DAG (which activates PKC, not modelled). Resting
PIP₂ in a platelet is ~112 000 copies (Mazet 2020). This module
**replaces the Dolan 2014 Fig. S2 IP₃ forcing curve** used in v0.2 with
a real upstream mechanism — IP₃ becomes a model *output*, not an
input.

**Implementation.** Three species (`PIP2[c]`, `IP3[c]`, `DAG[c]`) with
the PLCβ catalytic flux from §3 driving the hydrolysis. IP₃ has a
first-order degradation back to baseline (lumped IP₃-kinase /
IP₃-5-phosphatase), tuned so that resting IP₃ = 50 nM. The Mazet 2020
framework (Sci. Rep. 10:13889) supplies the rate constants in the
`K_PI_CYCLE` block. The `_ip3_forced` flag in `CalciumDynamics` is
retained for backwards compatibility with the Dolan validation
condition but defaults to off after v0.4.

## 5. IP3R — deYoung-Keizer with Nernst flux

**Biology.** The IP₃ receptor is a tetrameric Ca²⁺-release channel in
the DTS membrane (the platelet equivalent of the ER). Open
probability depends on both IP₃ (activation) and cytosolic Ca²⁺ (bell-
shaped: activation at low cyt Ca, inhibition at high cyt Ca, giving
calcium-induced calcium release). The deYoung-Keizer 1992 model
captures this with one slow inactivation ODE (`h`) and a quasi-steady
activation gate `m∞(IP3, Ca)`; Li & Rinzel 1994 reduce it to its
modern form. Total platelet IP3R = **1 328 channels** (Burkhart 2012
ITPR2; Dolan 2014 Table S1).

**Implementation.** State variable `IP3R_h[dts]` (slow inactivation
gate) with the Li-Rinzel ODE; `m∞` is computed from cyt Ca and IP₃ on
the fly. Open probability is `m∞⁴ × h` (tetrameric cooperativity).
Ca²⁺ flux uses the **Nernst form** (Purvis 2008 eq. 13 / Dolan 2014
eq. 4): `I = γ·N·Po·(ψ_IM − E_Ca,IM)·NA/(zF)`. The single-channel
conductance `GAMMA_IP3R_S = 0.075 pS` is **back-derived** to balance
SERCA at rest — this is a *calibration-coupled* parameter and any
change to SERCA rate constants requires re-deriving it (see
`reports/dissertation-notes.md §3.1`). Parameters in `K_DYK` and the
preceding constant block.

## 6. SERCA pump — Purvis E1/E2 6-state cycle

**Biology.** SERCA (sarco/endoplasmic reticulum Ca²⁺-ATPase) pumps
cytosolic Ca²⁺ back into the DTS lumen against a ~2 500× concentration
gradient, consuming 1 ATP per 2 Ca²⁺. The platelet isoform is
SERCA3b. Modelled mechanistically as a 6-state enzymatic cycle:
`E2 ⇌ E1 → E1·Ca → E1P·Ca → E2P·Ca → E2P → E2`, capturing
phosphorylation, conformational switch, and Ca²⁺ release into the
lumen.

**Implementation.** Six species (`SERCA_E1`, `SERCA_E2`, `SERCA_E1Ca`,
`SERCA_E1PCa`, `SERCA_E2PCa`, `SERCA_E2P`, all `[dts]`) with the rate
constants of Purvis 2008 Table 1 (`K_SERCA` block). The ATP cost is
counted per cycle in `CalciumDynamics.evolveState`. There is an
**open biology question** about whether the Purvis rate constants
over-estimate the SERCA3b Vmax at resting cyt Ca²⁺ (~4.7 cycles/s/pump
vs. literature 2 % of saturating Vmax ≈ 0.6–1 cycle/s); documented in
the file header and `dissertation-notes.md §3.2`.

## 7. PMCA — Caride 5-state with in-vivo k₁₂ correction

**Biology.** The plasma-membrane Ca²⁺-ATPase extrudes cytosolic Ca²⁺
to the extracellular space, also ATP-driven. It has two parallel
catalytic pathways: a **basal** path (PMCA + Ca ⇌ PMCA·Ca →
PMCA + Ca_ex) and a **Ca₄·CaM-activated** path that runs ~5× faster.
The CaM-activated path makes PMCA the dominant late-phase extruder
once a Ca²⁺ transient has built up cytosolic Ca₄·CaM. Caride 2007
gives the 5-state scheme (steps 4–5 basal, 8–11 CaM-activated).

**Implementation.** Five species (`PMCA`, `PMCA_Ca`, `Ca4_CaM_PMCA`,
`Ca4_CaM_PMCA_Ca`, `PMCA_CaM`, all `[pl]`) with rate constants in the
`K_PMCA` and `K_CAM_PMCA` blocks. v0.3.1 replaced the in-vitro k₁₂
value with an in-vivo correction (Penniston & Enyedi 1998) — without
this, PMCA over-extrudes at low cyt Ca and the resting state drifts
downward. This is one of three places where the Mazet *mosaic data*
critique surfaces.

## 8. Cytosolic buffering — CaM ladder + gelsolin proxy

**Biology.** ~85 % of cytosolic Ca²⁺ is bound to proteins at any
moment; only the free fraction participates in signalling. The two
dominant buffer classes in a platelet are **calmodulin** (two-lobed
Ca²⁺-binding ladder: free → Ca₂·CaM → Ca₄·CaM; Caride 2007 steps 6–7)
and **gelsolin** (a multi-site Ca²⁺-binding actin-severing protein,
~100 000 copies × ~5 sites = ~500 000 sites). Smaller contributions
from annexins, Ca·ATP, and Ca²⁺-binding kinases are lumped together.

**Implementation.** CaM is three species (`CaM_free`, `Ca2_CaM`,
`Ca4_CaM`, all `[c]`) with cooperative two-Ca-at-a-time binding from
the `K_CAM` block — biologically motivated and dynamically important
because **Ca₄·CaM activates PMCA**. Gelsolin is collapsed to a single
coarse-grained 1:1 buffer (`GSN_free[c]`, `GSN_Ca[c]`; `K_GSN` block;
Kd = 1 µM). The `N_GSN = 1 400 000` effective sites is *calibrated*
against the Phase 3 peak band — splitting gelsolin / annexin / Ca·ATP
into separate buffers is on the wishlist.

## 9. DTS buffering — CALR, HSP90B1, BiP, CREC

**Biology.** ~95–99 % of DTS Ca²⁺ is bound to luminal proteins (real
biology). At rest the platelet DTS holds ~250 µM *free* Ca²⁺ in a
~0.26 fL volume, but the **total** Ca²⁺ content is an order of
magnitude higher. The dominant luminal buffer is **calreticulin
(CALR)**, with two domains: the C-domain has ~25 low-affinity sites
(Kd ~1 mM — the dynamic buffer), and the P-domain has one high-affinity
site per molecule (Kd ~1 µM — always saturated, slow exchange).
Additional buffers (HSP90B1 medium- and low-affinity sites, BiP/HSPA5,
and the lumped CREC pool of CALU + RCN1 + RCN2) make up the remainder.

**Implementation.** Eight species across four buffer classes — `CALR_free/Ca`,
`CALR_P_free/Ca`, `HSP90B1_M_free/Ca`, `HSP90B1_L_free/Ca`, `BiP_free/Ca`,
`CREC_free/Ca`, all `[dts]`. Each is a 1:1 reversible binding with on/off
rates in its own K-block (`K_CALR`, `K_CALR_P`, `K_HSP90B1_M/L`, `K_BIP`,
`K_CREC`). Total site count is ~800 000. This was the **Phase 2 / #28**
work landed earlier this month and is what enabled the DTS to behave like a
real Ca²⁺ store rather than overfilling under sustained SERCA pumping.

## 10. SOCE — STIM1 + Orai1, Hoover-Lewis MWC

**Biology.** Store-operated Ca²⁺ entry is the mechanism that refills
the DTS after a Ca²⁺ release event. **STIM1** is a single-pass DTS-
membrane protein whose luminal EF-hand senses DTS Ca²⁺; when [Ca²⁺]_DTS
drops, STIM1 unbinds Ca²⁺, dimerises, and translocates to ER-PM
junctions where it gates **Orai1** in the plasma membrane. Orai1
allows external Ca²⁺ entry into the cytosol, which SERCA then pumps
back into the DTS. Open-probability of Orai1 is a Monod-Wyman-
Changeux (MWC) function of STIM2-puncta occupancy (Hoover & Lewis 2011
framework, calibrated by Dolan 2014 Fig. 4 EDTA-vs-Ca²⁺ paired
condition).

**Implementation.** STIM1 has three states (`STIM1_free`, `STIM1_Ca`,
`STIM1_dim`, all `[dts]`) with Ca²⁺-binding and dimerisation in the
`K_STIM` block; Orai1 open probability uses the MWC formula in
`K_MWC` and `GAMMA_SOC_S`. The SOCE current is gated on
`CA_EX_UM > 0` so the EDTA validation condition disables it. This
module replaced an earlier ad hoc 3-state mass-action model (closed
issues #45/#46).

## 11. P2X1 + NCX — plasma-membrane Ca²⁺ entry/exit

**Biology.** Two additional plasma-membrane pathways. **P2X1** is an
ATP-gated cation channel that opens fast on the millisecond scale when
ATP is released from dense granules; it provides a brief early Ca²⁺
spike that primes the cell before IP3R release. The channel has three
states (closed → open → desensitised → recovered). **NCX** (the
Na⁺/Ca²⁺ exchanger) is a low-affinity Ca²⁺ extruder driven by the Na⁺
gradient — it complements PMCA for high-Ca²⁺ extrusion and dominates
when cyt Ca²⁺ peaks above a few hundred nM.

**Implementation.** P2X1 is three species (`P2X1`, `P2X1_O`, `P2X1_D`,
all `[pl]`) with rate constants and a Ca²⁺-specific conductance
(`GAMMA_P2X1_S = 1.3 fS`) calibrated against Phase 3 timing.
ATP_ex forcing supplies the gating ligand. NCX is parameterised
allosterically (cytosolic Ca²⁺-dependent activation via the `K_NCX`
Hill formulation) and contributes a separate flux term in the ODE
RHS. Both landed in **v0.3.4**.

## 12. Mitochondrial Ca²⁺ — MCU + NCLX

**Biology.** Mitochondria provide a fast cytosolic Ca²⁺ sink during
transients (via the **mitochondrial Ca²⁺ uniporter, MCU**) and a slow
release pathway (via **NCLX**, the mitochondrial Na⁺/Ca²⁺ exchanger).
MCU has a sharp activation curve — half-saturated at ~1 µM with high
Hill cooperativity — meaning it is effectively off at resting cyt
(100 nM) and switches on only during a transient peak. This buffering
bypasses the PMCA-rate-limited extrusion bottleneck.

**Implementation.** Single species `CA2_MITO[m]` with MCU uptake
(Hill function of cyt Ca²⁺, `K_MITO['K_MCU'] = 1 µM`, `n_MCU` high) and
NCLX release (first-order in mito Ca²⁺). Block: `K_MITO`. Landed in
**v0.4.0** (issue #22). One useful negative finding: MCU buffering
*alone* does not accelerate DTS recovery — testable model prediction
that pharmacological MCU activation may actually *prolong* the DTS
overshoot.

---

## Validation summary

| Metric | Target (Dolan 2014) | Model (v0.4.1) |
|---|---|---|
| Resting cyt Ca²⁺ | 100 ± 10 nM | 104 nM |
| Resting DTS Ca²⁺ | 250 µM | 235 µM |
| Resting IP₃ | 50 nM | 50 nM |
| Resting Gαq-active | ~100 | 100 |

Two documented honest findings remain open:

1. **Long-recovery cyt collapse** (>2000 s sustained stim): cyt drops
   to ~6 nM. Caused by deYoung-Keizer m∞ → 0 below d₅ = 82 nM
   (bistability artefact). Phase 3 (200 s window) unaffected.
2. **DTS overshoot to ~1 mM** during 3000 s recovery: PMCA + NCX
   combined extrusion vs. SOCE + P2X1 entry. May be real biology
   (matches long Ca²⁺ tails reported in real activated platelets).

Both are documented in `reports/lab-books/` and are candidates to be
addressed.
