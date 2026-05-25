"""
Calcium signalling dataclass for the platelet whole-cell model.

Holds species ordering, rate constants, compartment volumes, and the
ODE right-hand-side used by the CalciumDynamics process. The ODE covers:

  * IP3R deYoung-Keizer 1992 model (Li-Rinzel 1994 simplification): one
    slow inactivation ODE for h (fraction of non-inhibited channels) plus
    quasi-steady activation m∞(IP3, Ca). Replaces the Sneyd-Dufour 2002
    6-state Markov model, which was correctly implemented but calibrated
    at IP3 = 10 µM and extrapolated poorly to resting IP3 = 50 nM.
  * IP3R Ca²⁺ flux via the Nernst form (Purvis 2008 eq. 13 / Dolan 2014
    eq. 4): I = γ·N·Po·(ψ_IM − E_Ca,IM)·NA/(zF), with γ_IP3R = 0.35 pS
    (calibrated to Dolan 2014 platelet resting state; Phase 4 #30),
    Po = m∞(IP3, Ca)⁴ × h.
  * SERCA E1/E2 cycle (Purvis 2008 Table 1, Dode 2002 kinetics) with the
    primary-source rate constants — including k_bind_f = 1×10¹⁵ M⁻²s⁻¹.
  * PMCA 5-state CaM-coupled scheme (Caride 2007 Table 3): basal path
    (steps 4–5) plus CaM-activated path (steps 8–11).
  * Calmodulin Ca²⁺ binding (Caride 2007 Table 3 steps 6–7): two-lobe
    cooperative ladder CaM → Ca₂·CaM → Ca₄·CaM. Ca₄·CaM activates PMCA
    (5× higher k_cat) and acts as a cytosolic Ca²⁺ buffer.
  * SOCE: Dolan 2014 MWC allosteric scheme (Hoover & Lewis 2011 framework)
    parameterising channel open probability as a function of STIM2 in the
    Orai puncta. Replaces the prior ad hoc 3-state mass-action model
    (issues #45/#46).

State layout: integer counts per species, indexed by `MOLECULE_NAMES`.
ODE works in count units; rate laws convert to concentration internally
where needed (volumes for cytosol vs DTS).

Numerical regime: 1-second outer timestep, scipy.integrate.solve_ivp with
the BDF method (the system is stiff — SERCA cycle rates run up to 1000 s⁻¹,
IP3R φ-function rate laws produce stiff loops).
"""

import math

import numpy as np
from scipy.integrate import solve_ivp

from reconstruction.platelet.dataclasses.process._params_loader import (
	load_calcium_kinetics)


# Externalised rate constants (issue #32 Phase 2). Sections of this dict
# are spliced into the module-level K_* names below as the kinetics-as-data
# refactor lands. Today only `[ip3r.k_dyk]` is sourced from the TOML; all
# other K_* dicts remain inline literals until their slice arrives.
_KINETICS = load_calcium_kinetics()


# ── Compartment volumes ───────────────────────────────────────────────────
# Source: Purvis 2008 (direct measurement; 6 fL cytosol, 4.3% DTS).
V_CYT_L = 6.0e-15      # cytosol, litres
V_DTS_L = 0.258e-15    # DTS,    litres (4.3% of cytosol)
N_A = 6.022e23         # Avogadro

# Conversion factors: concentration_uM = count / (N_A × volume × 1e-6)
_UM_PER_COUNT_CYT = 1.0 / (N_A * V_CYT_L * 1e-6)
_UM_PER_COUNT_DTS = 1.0 / (N_A * V_DTS_L * 1e-6)


# ── Species ordering ──────────────────────────────────────────────────────
# This must match reconstruction/platelet/dataclasses/internal_state.py
# for the calcium-pathway subset. Listed in the order used inside the ODE
# state vector y[].
MOLECULE_NAMES = (
	'CA2_CYT[c]',
	'CA2_DTS[dts]',
	'IP3[c]',
	'IP3R_h[dts]',
	'SERCA_E1[dts]',
	'SERCA_E2[dts]',
	'SERCA_E1Ca[dts]',
	'SERCA_E1PCa[dts]',
	'SERCA_E2PCa[dts]',
	'SERCA_E2P[dts]',
	'PMCA[pl]',
	'PMCA_Ca[pl]',
	'STIM1_free[dts]',
	'STIM1_Ca[dts]',
	'STIM1_dim[dts]',
	'ORAI1[pl]',
	# CaM Ca²⁺-binding ladder (Caride 2007 steps 6–7; Phase 1c)
	'CaM_free[c]',
	'Ca2_CaM[c]',
	'Ca4_CaM[c]',
	# PMCA–CaM complex sub-states (Caride 2007 steps 8–11; Phase 1d)
	'Ca4_CaM_PMCA[pl]',
	'Ca4_CaM_PMCA_Ca[pl]',
	'PMCA_CaM[pl]',
	# Coarse-grained cytosolic Ca²⁺ buffer (gelsolin proxy; scaffold-only,
	# see K_GSN block below for the biology / scope disclosure).
	'GSN_free[c]',
	'GSN_Ca[c]',
	# Calreticulin DTS Ca²⁺-binding sites (Phase 2 / #28; see K_CALR block).
	'CALR_free[dts]',
	'CALR_Ca[dts]',
	# CALR high-affinity P-domain (1 site per CALR; slow release).
	'CALR_P_free[dts]',
	'CALR_P_Ca[dts]',
	# Additional DTS luminal buffers (Phase 3 / #25; see K_HSP90B1_*,
	# K_BIP, K_CREC blocks below for biology / provenance).
	'HSP90B1_M_free[dts]',  # medium-affinity sites (Kd ~ 2 µM)
	'HSP90B1_M_Ca[dts]',
	'HSP90B1_L_free[dts]',  # low-affinity sites (Kd ~ 600 µM)
	'HSP90B1_L_Ca[dts]',
	'BiP_free[dts]',        # HSPA5/BiP, low-affinity
	'BiP_Ca[dts]',
	'CREC_free[dts]',       # aggregated CALU + RCN1 + RCN2 pool
	'CREC_Ca[dts]',
	# P2X1 ATP-gated cation channel (Phase 2.5; see K_P2X1 block).
	'P2X1[pl]',     # closed (resting)
	'P2X1_O[pl]',   # open (conducting)
	'P2X1_D[pl]',   # desensitised
	# PI cycle — replaces forced-IP3 (Phase 4 / #31; Mazet 2020 framework).
	'PIP2[c]',           # PI(4,5)P2 — PLCβ substrate
	'DAG[c]',            # diacylglycerol — PLCβ co-product (PKC substrate)
	'PLCb_inactive[c]',  # PLCβ inactive pool
	'PLCb_active[c]',    # PLCβ active (Gq-bound); catalyses PIP2 hydrolysis
	# Mitochondrial Ca²⁺ (issue #22, 2026-05-12) — provides fast
	# Ca²⁺ uptake during transients (via MCU) and slow release (via
	# NCLX), bypassing the PMCA-rate-limited extrusion bottleneck.
	'CA2_MITO[m]',       # free mito matrix Ca²⁺ count
	# GPCR receptor cascade (issue #9) — P2Y1 + PAR1 + PAR4 → Gαq drives
	# PLCβ activation, which produces IP3 endogenously through the PI cycle.
	'P2Y1_inactive[pl]', # P2Y1 ADP receptor (Gq-coupled), inactive
	'P2Y1_active[pl]',   # P2Y1 active (ADP-bound)
	'PAR1_inactive[pl]', # PAR1 thrombin receptor (high-affinity), inactive
	'PAR1_active[pl]',   # PAR1 cleaved/active (proteolytic, ~irreversible)
	'PAR1_internalized[pl]', # v0.4.1: one-way sink (was 2-state recycling)
	'PAR4_inactive[pl]', # PAR4 thrombin receptor (low-affinity)
	'PAR4_active[pl]',
	'PAR4_internalized[pl]', # v0.4.1: one-way sink
	'Gq_active[c]',      # Gαq-GTP (active); inactive Gq implicit (total - active)
)
# Index lookups for readability inside the rate function.
_IDX = {name: i for i, name in enumerate(MOLECULE_NAMES)}
N_SPECIES = len(MOLECULE_NAMES)


# ── Resting concentrations ────────────────────────────────────────────────
CA_EX_UM = 1200.0       # extracellular Ca²⁺, fixed reservoir (Dolan 2014)
IP3_REST_UM = 0.05      # cytosolic IP3 baseline (50 nM; Purvis 2008)


# ── Physical constants ────────────────────────────────────────────────────
# Used by the Nernst-based IP3R and SOCE flux equations (Purvis Table 1
# row "Ca²⁺ release from DTS" / Dolan eq. 4).
F_FARADAY        = 96485.0          # C/mol
R_GAS            = 8.314            # J/(mol·K)
T_KELVIN         = 310.0            # 37 °C; Purvis/Dolan
RT_OVER_zF_V     = R_GAS * T_KELVIN / (2.0 * F_FARADAY)   # ≈ 0.01334 V (z=2 for Ca²⁺)
NA_OVER_zF       = N_A / (2.0 * F_FARADAY)                # ions per ampere-second (z=2)

# Membrane potentials (Dolan 2014 Methods §"Membrane potentials"):
#   V_IM responsive cluster sits at the upper end of the −100..−60 mV
#   sampling range (V_IM > −70 mV); use the Dolan upper bound.
#   V_PM measured 60..70 mV inside-negative; Dolan uses −60 mV.
V_IM_V = -0.060          # DTS-membrane potential (V)
V_PM_V = -0.060          # plasma-membrane potential (V)


# ── IP3R: deYoung-Keizer 1992 / Li-Rinzel 1994 ───────────────────────────
# Replaces the Sneyd-Dufour 2002 6-state Markov model (see
# lab-book-2026-05-08-sneyd-dufour-audit.md for the audit confirming the
# prior implementation was correct but poorly conditioned at resting IP3).
#
# Li & Rinzel 1994 reduce the 8-state deYoung-Keizer model to one slow ODE:
#
#   m∞(IP3, Ca) = [IP3 / (IP3 + d₁)] × [Ca_cyt / (Ca_cyt + d₅)]
#   dh/dt       = a₂ × [d₂ − (Ca_cyt + d₂) × h]
#   Po          = m∞⁴ × h
#
# h → h∞ = d₂/(Ca_cyt + d₂) in quasi-steady state; τ_h = 1/(a₂(Ca_cyt+d₂)).
# The ⁴ exponent encodes tetrameric cooperativity (all four subunits must be
# non-inhibited and in the m∞ activation state).
#
# Parameters from deYoung & Keizer 1992 PNAS 89:9895-9899 Table 1, via the
# derived dissociation constants in Li & Rinzel 1994 J Theor Biol 166:461-473.
# Values live in `reports/params/calcium-v0.5.toml [ip3r.k_dyk]` (issue #32
# Phase 2). Edit there to change. `dict(...)` is a defensive copy so a
# caller can't mutate the loader's underlying dict.
K_DYK = dict(_KINETICS['ip3r']['k_dyk'])

# IP3R channel ensemble (count + conductance). Values + full
# calibration-coupling commentary live in `reports/params/calcium-v0.5.toml
# [ip3r.channel]` (issue #32 Phase 2 slice 8).
# ⚠ CALIBRATION-COUPLED with K_SERCA below — see TOML and
# calibration-coupling-2026-05-25.qmd chains 1 + 2.
N_IP3R       = _KINETICS['ip3r']['channel']['n_total']
GAMMA_IP3R_S = _KINETICS['ip3r']['channel']['gamma_s']


# SERCA cycle (Purvis 2008 Table 1, Dode 2002 isoform 3b kinetics).
# Values + ⚠ CALIBRATION-COUPLED commentary + open-biology-question
# notes live in `reports/params/calcium-v0.5.toml [serca.cycle]`
# (issue #32 Phase 2 slice 8).
K_SERCA = dict(_KINETICS['serca']['cycle'])


# ── PMCA4b basal path (Caride 2007 Table 3 steps 4–5) ────────────────────
# Steps 4–5 are unchanged; the CaM-activated path (steps 8–11) is below.
# Values live in `reports/params/calcium-v0.5.toml [pmca.basal]`
# (issue #32 Phase 2 slice 3).
K_PMCA = dict(_KINETICS['pmca']['basal'])

# ── CaM Ca²⁺ binding (Caride 2007 Table 3 steps 6–7) ─────────────────────
# Two-lobe cooperative scheme: slow N-lobe (step 6) then fast C-lobe (step 7).
# Ca²⁺ concentrations in µM; rates in µM⁻²·s⁻¹ (forward) or s⁻¹ (reverse).
# Values live in `reports/params/calcium-v0.5.toml [cam.binding]`
# (issue #32 Phase 2 slice 3).
K_CAM = dict(_KINETICS['cam']['binding'])

# ── Coarse-grained cytosolic Ca²⁺ buffer (gelsolin proxy) ────────────────
# Real platelet cytosolic Ca²⁺ buffering ratio is ~50:1 (bound:free; Sage
# & Rink 1985), dominated by gelsolin (~50 000–280 000 copies × multi-site
# Ca²⁺ binding; Burkhart 2012; Yin & Stossel 1979), with smaller
# contributions from annexins, Ca·ATP, and Ca²⁺-binding kinases. Our
# model represents this as a single coarse-grained 1:1 buffer:
#
#   GSN_free + Ca²⁺  ⇌  GSN_Ca       (k_on, k_off; Kd = 1 µM)
#
# `N_GSN` counts effective Ca²⁺-binding sites, not gelsolin molecules.
# Calibrated jointly with CALR (Phase 2 / #28): with the DTS buffered,
# IP3R transients deliver substantially more Ca²⁺ to the cytosol, and
# the cytosolic buffer is needed to keep peaks in the Dolan Fig 4 band.
# The Phase 2 calibration intentionally over-shoots the Sage & Rink
# ~50:1 ratio (current rest ~200:1) — this single tunable buffer
# absorbs slack from all other cytosolic buffers not yet broken out.
#
# See `reports/dissertation-notes.md §1.1` for the literature gap and
# v0.3+ plan to split this into explicit gelsolin / annexin / Ca-ATP.
# Values live in `reports/params/calcium-v0.5.toml [buffers.gsn]` and
# [buffers.gsn_pool] (issue #32 Phase 2 slice 7).
K_GSN = dict(_KINETICS['buffers']['gsn'])

N_GSN = _KINETICS['buffers']['gsn_pool']['n_total']


# ── Calreticulin (CALR) DTS Ca²⁺ buffer — Phase 2 / issue #28 ────────────
# Calreticulin is the dominant luminal Ca²⁺-binding protein in the ER/SR
# (and the DTS, which is the platelet equivalent). It has three domains:
#   N (lectin)
#   P (proline-rich): 1 high-affinity Ca²⁺ site, Kd ~ 1 µM (always saturated
#     at physiological DTS [Ca²⁺] — folded into the existing CA2_DTS pool)
#   C (acidic):      ~25 low-affinity Ca²⁺ sites, Kd ~ 1 mM — this is the
#     dominant *dynamic* DTS Ca²⁺ buffer and what this block models.
#
# Platelet copy number: 20 324 CALR molecules (Burkhart 2012; Dolan 2014
# Table S1). With 25 low-affinity sites per molecule the total binding
# capacity is N_CALR = 508 100 sites. At [Ca²⁺]_DTS = 250 µM (Dolan
# resting): fractional occupancy = 250/(250+1000) = 0.20 →
#   CALR_Ca   = 101 620 sites occupied
#   CALR_free = 406 480 sites empty
#
# Effect: real biology has 95–99 % of DTS Ca²⁺ buffered; previously the
# model had ~9 % (STIM1 EF-hand only) and the resting DTS overfilled to
# > 1 mM in 6 000-s integrations because SERCA pumped into a near-
# unbuffered lumen. With CALR the buffering ratio rises to ~72 % at rest
# and the DTS gains a large "spring" that absorbs SERCA excess and
# releases during IP3R drainage.
#
# Each binding site is treated as an independent 1:1 Ca²⁺ binder. This is
# the standard coarse-graining for calreticulin and matches the Sneyd
# 2014 / Hofer 1998 buffer formulation. v0.3+ may refine to include the
# high-affinity P-domain site or extend to HSP90B1 / CALU (issue #25).
#
# See `reports/dissertation-notes.md §2.1` for the full biological context.
# DTS luminal buffer cluster — Phase 2/3 issues #25, #28.
# ⚠ CALIBRATION-COUPLED with resting DTS Ca²⁺ level (Dolan 250 µM
# target). Values + full literature commentary live in
# `reports/params/calcium-v0.5.toml [buffers.*]` (issue #32 Phase 2
# slice 9). See calibration-coupling-2026-05-25.qmd chain 4.

K_CALR   = dict(_KINETICS['buffers']['calr'])
N_CALR   = _KINETICS['buffers']['calr_pool']['n_total']

K_CALR_P = dict(_KINETICS['buffers']['calr_p'])
N_CALR_P = _KINETICS['buffers']['calr_p_pool']['n_total']


# ── HSP90B1 / GRP94 / endoplasmin — Phase 3 / issue #25 ──────────────────
# Second-most-abundant ER luminal Ca²⁺-binding chaperone after CALR. Argon
# & Simen 1999 report 15 Ca²⁺-binding sites per molecule, split into:
#   - 4 medium-affinity sites at Kd ~ 2 µM (matters during DTS depletion)
#   - 11 low-affinity sites at Kd ~ 600 µM (matters at resting [Ca²⁺]_DTS)
#
# Modelling these separately captures the dynamics: at rest both bind
# Ca²⁺; during the IP3R-driven transient, the low-affinity sites release
# fast (the bulk of the deliverable reserve), while the medium-affinity
# sites hold until free [Ca²⁺]_DTS drops below ~ 2 µM. The medium sites
# are what should produce a "floor" above zero in free [Ca²⁺]_DTS during
# peak — the headline biology improvement of this issue.
#
# Platelet copy number: ~10 000 (order-of-magnitude estimate; HSP90B1 is
# the chaperone-of-record alongside CALR and BiP in the ER, Burkhart 2012
# lists it among the top ER-associated proteins; precise count is in the
# supplementary table not directly extractable). Flagged in dissertation
# notes as a v0.3 stretch estimate.
K_HSP90B1_M = dict(_KINETICS['buffers']['hsp90b1_medium'])
K_HSP90B1_L = dict(_KINETICS['buffers']['hsp90b1_low'])

# Per-molecule × sites-per-molecule decomposition kept in TOML so a
# biologist can edit either independently.
_HSP90B1_POOL = _KINETICS['buffers']['hsp90b1_pool']
N_HSP90B1   = _HSP90B1_POOL['n_molecules']
N_HSP90B1_M = _HSP90B1_POOL['n_molecules'] * _HSP90B1_POOL['sites_medium']
N_HSP90B1_L = _HSP90B1_POOL['n_molecules'] * _HSP90B1_POOL['sites_low']


# ── BiP / HSPA5 / GRP78 — Phase 3 / issue #25 ────────────────────────────
# Most abundant ER chaperone. Primary role is protein folding, but
# Lièvremont 1997 demonstrated BiP contributes ~25 % of the ER Ca²⁺
# store, with a stoichiometry of 1–2 Ca²⁺ per molecule. Modelled as a
# single low-affinity pool with Kd ~ 500 µM (matches the Lièvremont
# "mM-range free [Ca²⁺]_ER" they describe).
#
# Platelet copy number: ~50 000 (order-of-magnitude estimate; BiP is the
# canonically most-abundant ER chaperone, more so than CALR or HSP90B1).
K_BIP = dict(_KINETICS['buffers']['bip'])
_BIP_POOL = _KINETICS['buffers']['bip_pool']
N_BIP = _BIP_POOL['n_molecules'] * _BIP_POOL['sites_per']


# ── CREC family pool — CALU + RCN1 + RCN2 lumped — Phase 3 / issue #25 ──
# The CREC family proteins are smaller multi-EF-hand low-affinity
# Ca²⁺ binders localised to the ER lumen and secretory pathway. Honoré
# & Vorum 2000 review: 6–7 EF-hands per molecule with Kd "up to mM",
# i.e. very low affinity. Vorum 1998 cloned and characterised calumenin.
# Aggregated here as one coarse-grained pool (split into individual
# proteins in v0.4 if granule secretion work requires CALU specifically).
#
# Combined platelet copy number estimate: ~15 000 (CALU ~5 k + RCN1 ~5 k
# + RCN2 ~5 k). Effective Ca²⁺-binding sites per molecule ~4 (most
# EF-hands have functional Ca²⁺ binding; some are structural).
K_CREC = dict(_KINETICS['buffers']['crec'])
_CREC_POOL = _KINETICS['buffers']['crec_pool']
N_CREC = _CREC_POOL['n_molecules'] * _CREC_POOL['sites_per']


# ── PMCA4b CaM-activated path (Caride 2007 Table 3 steps 8–12) ──────────
# Ca₄·CaM binds free PMCA (step 8), the complex binds and pumps Ca²⁺ with
# ~5× higher k_cat than basal (step 10 vs step 5). Step 11 dissociates 4
# Ca²⁺ from the active complex (Ca₄·CaM·PMCA → PMCA·CaM + 4 Ca²⁺_cyt) and
# step 12 is the slow CaM dissociation (PMCA·CaM → PMCA + CaM, k12 = 0.033
# s⁻¹, τ ~ 30 s). Both 11 and 12 are required for mass conservation: with
# 11 active but 12 absent, PMCA accumulates dead-end in PMCA·CaM (the bug
# previously worked around by omitting step 11 entirely; restored
# 2026-05-07 after Phase 0 audit found Caride k₁₂ missing).
# Values live in `reports/params/calcium-v0.5.toml [pmca.cam_activated]`
# (issue #32 Phase 2 slice 7). k12 in-vivo override commentary (1 s⁻¹
# vs Caride's in-vitro 0.033 s⁻¹) lives in the TOML section header.
K_CAM_PMCA = dict(_KINETICS['pmca']['cam_activated'])


# ── P2X1 ATP-gated cation channel (Phase 2.5, 2026-05-11) ────────────────
# P2X1 is the dominant *fast* Ca²⁺ entry pathway in activated platelets.
# It's a trimeric ATP-gated cation channel in the plasma membrane that
# opens within ms of extracellular-ATP exposure (released from dense
# granules during activation), conducts a mixed Na⁺/K⁺/Ca²⁺ current with
# ~5–10 % Ca²⁺ fraction, then desensitises within ~100 ms. Recovery
# from desensitisation is slow (~30 s).
#
# Source-info: Mahaut-Smith et al. 2000 / 2004 (platelet P2X1
# electrophysiology), Vial & Evans 2002 (P2X1 in platelets), Hechler
# et al. 2003 (P2X1 in thrombosis), Burnstock 2007 (P2X family review).
# Copy number: ~600–3000 channels per platelet (Mahaut-Smith); we use
# N_P2X1 = 1 000 as a mid-range estimate.
#
# Three-state coarse kinetic scheme:
#
#   Closed  --(k_act × [ATP_ex])-->  Open  --(k_des)-->  Desensitised
#     ^                               |                        |
#     |                            k_close                  k_rec
#     +-----(k_close)-----------------+   +--(k_rec)----------+
#
# Ca²⁺ flux through Open state: Nernst form, but with V_PM driving
# force (E_Ca,PM, not E_Ca,IM). Gated on CA_EX_UM > 0 — when
# extracellular Ca²⁺ is removed (EDTA condition), P2X1 contributes no
# Ca²⁺ even though the channel may still cycle.
#
# This is exactly what should close the SOCE-differential gap in
# Phase 3: the +Ca_ex peak now has a fast P2X1 contribution that −Ca_ex
# lacks (real biology: ~100 nM differential at peak).
# Values live in `reports/params/calcium-v0.5.toml [p2x1.kinetics]`
# (issue #32 Phase 2 slice 3).
K_P2X1 = dict(_KINETICS['p2x1']['kinetics'])

# Total P2X1 functional channels (trimers; mass per trimer = 3 × 45 kDa).
# Value lives in `[p2x1.channel] n_total`.
N_P2X1 = _KINETICS['p2x1']['channel']['n_total']

# P2X1 Ca²⁺-specific effective conductance.
# Single-channel current ~0.5–1 pA at -60 mV, Ca²⁺ fraction ~5–10 %,
# so effective Ca²⁺-specific γ ≈ 10–50 fS per channel. Starting from
# 0.01 pS — calibration anchor for Phase 3 SOCE-differential target
# (Dolan ~100 nM). See lab book. Value lives in `[p2x1.channel] gamma_s`.
GAMMA_P2X1_S = _KINETICS['p2x1']['channel']['gamma_s']


# ── Extracellular ATP forcing (drives P2X1) ───────────────────────────────
# ATP released from dense granules during activation reaches 1–10 µM near a
# forming thrombus, then is cleared by ectonucleotidases (CD39) over tens of
# seconds. CD39 keeps resting ATP_ex near zero — any small baseline leaks
# P2X1 over hundreds of seconds and overfills the DTS.
# Values live in `reports/params/calcium-v0.5.toml [agonists.atp_ex]`
# (issue #32 Phase 2 slice 2). The Python name `ATP_EX_TAU_RISE` (no
# `_S` suffix) is preserved verbatim from the pre-refactor module; the
# TOML uses the consistent `tau_rise_s` key.
ATP_EX_REST_UM   = _KINETICS['agonists']['atp_ex']['rest_uM']
ATP_EX_PEAK_UM   = _KINETICS['agonists']['atp_ex']['peak_uM']
ATP_EX_TAU_RISE  = _KINETICS['agonists']['atp_ex']['tau_rise_s']
ATP_EX_T_PEAK    = _KINETICS['agonists']['atp_ex']['t_peak_s']
ATP_EX_TAU_DECAY = _KINETICS['agonists']['atp_ex']['tau_decay_s']


def atp_ex_forcing_uM(t, delay=0.0, peak_uM=None):
	"""Plateau-decay approximation for extracellular ATP during activation.

	`peak_uM=None` (default) reads `ATP_EX_PEAK_UM` from the module at call
	time, so reassigning the module constant takes effect immediately.
	Passing `peak_uM=0` yields a flat resting curve.
	"""
	if peak_uM is None:
		peak_uM = ATP_EX_PEAK_UM
	t_eff = t - delay
	if t_eff < 0:
		return ATP_EX_REST_UM
	rise = 1.0 - np.exp(-t_eff / ATP_EX_TAU_RISE)
	decay = np.exp(-max(0.0, t_eff - ATP_EX_T_PEAK) / ATP_EX_TAU_DECAY)
	return ATP_EX_REST_UM + (peak_uM - ATP_EX_REST_UM) * rise * decay


# ── SOCE: Dolan 2014 MWC + STIM1 dimerisation (Hoover & Lewis 2011 frame) ─
# STIM1 cycle (mass-action) — keeps the dimer pool size as a state variable.
# `STIM1_dim` is counted in DIMER PARTICLES (matches Dolan Table S1
# "STIM1₂ (11)") — 1 dimerisation event consumes 2 free monomers and
# creates 1 dimer particle. Rate constants chosen so the Dolan 2014
# Table S1 resting IC (st_Ca=3805, st_free=438, st_dim=11) is at
# detailed balance.
# Values live in `reports/params/calcium-v0.5.toml [soce.stim]`
# (issue #32 Phase 2 slice 4). Detailed-balance derivation for
# `k_release_r` and `k_dim_f` is documented in the TOML section header.
K_STIM = dict(_KINETICS['soce']['stim'])

# Hoover & Lewis 2011 MWC parameters (Fig. 4 best-fit, verified PDF):
#   L  — intrinsic opening equilibrium constant (closed→open without STIM)
#   f  — opening cooperativity factor per bound STIM2 (each STIM2 stabilises
#        the open state by factor f via fL, f²L, ..., f⁴L)
#   a  — binding cooperativity factor (Hoover labels this `a`; <1 = negative
#        cooperativity for successive bindings)
#   Ka — STIM association constant. Hoover fits Ka=100 in HEK arbitrary units
#        (a.u.) where saturating STIM expression Stotal=3.2 a.u. To map onto
#        platelet dimer counts we rescale: in our model Sf ranges from ~0.1
#        dimers (rest) to ~170 dimers (full puncta entry at saturating Ca).
#        Setting Ka_platelet so that Ka_platelet × Sf_saturating ≈ Hoover's
#        Ka × Stotal = 320 gives Ka_platelet = 320/170 ≈ 1.9, which we round
#        to 2 (the MWC shape is insensitive to ~2× perturbations once at the
#        saturating end of the binding curve). f, a, and L are dimensionless
#        and transfer directly.
# Values live in `reports/params/calcium-v0.5.toml [soce.mwc]`
# (issue #32 Phase 2 slice 4). Ka rescaling from Hoover a.u. → platelet
# dimer counts is documented in the TOML section header.
K_MWC = dict(_KINETICS['soce']['mwc'])

# Dolan 2014 puncta entry (eq. 2): qp = α·[Ca]_cyt^n / (KM^n + [Ca]_cyt^n) + 0.01
#   qp gives the fraction of STIM2 dimers translocated into puncta where
#   they can engage Orai. α = 0.2 is the Dolan default. KM and n are the
#   two free parameters Dolan scans within homeostatic constraints.
# Values live in `reports/params/calcium-v0.5.toml [soce.puncta]`
# (issue #32 Phase 2 slice 4).
PUNCTA = dict(_KINETICS['soce']['puncta'])

# Orai single-channel Ca²⁺ conductance. The CRAC channel literature value is
# ~24 fS (Prakriya & Lewis 2002, Vig 2006), measured at saturating Po with
# patch-clamp in HEK cells. For the platelet model the *effective* γ_SOC is
# reduced by the integer-count realism of having <1 channel open at rest:
# Hoover's L=10⁻⁴ would give 0.04 fully-open channels with γ=24 fS, producing
# spurious µM/s leaks at rest. We calibrate γ_SOC analytically against the
# resting balance condition SOCE_rest ≈ PMCA_steady_rest ≈ 76 ions/s,
# which gives γ_SOC ≈ 0.3 fS at the Po(MWC, Sf_rest) ≈ 1.2×10⁻³ value our
# rescaled Ka produces. (Issue #46 — full single-channel current calibration.)
GAMMA_SOC_S = 0.3e-15            # 0.3 fS = effective single-channel conductance

# ── Basal plasma-membrane Ca²⁺ leak ──────────────────────────────────────
# A small constant cyt influx that compensates PMCA outflow at rest, keeping
# the cytosolic resting concentration at ~100 nM. Biologically this represents
# unidentified background Ca²⁺ entry pathways (TRPC, NCX reverse, residual
# constitutive permeability — Sage & Rink 1985–1990; Brandman & Liou 2010
# review). Calibrated against the steady-state PM balance condition
# J_SOCE + J_leak = J_PMCA at cyt=100 nM, DTS=250 µM:
#   PMCA quasi-eq outflow at cyt=100 nM ≈ k_cat · PMCA·Ca_eq ≈ 5.5 · 14 = 77
#   SOCE at full DTS / basal STIM1_dim ≈ 6
#   ⇒ leak ≈ 71 ions/s, rounded to 75
# This is the (ii) addition diagnosed in lab-book 2026-05-05; before this
# term the model had no PM-side cyt source large enough to balance PMCA.
J_PM_LEAK_IONS_S = 75.0          # ions/s, constant cyt influx


# ── NCX (Na⁺/Ca²⁺ exchanger) — v0.3.4 / second extrusion pathway ─────────
# Forward-mode Ca²⁺ extrusion: 3 Na⁺ in : 1 Ca²⁺ out, driven by the Na⁺
# gradient (no ATP). Provides a second Ca²⁺ extrusion pathway alongside
# PMCA, which is rate-limited at low cyt Ca²⁺. NCX has higher Vmax per
# transporter and a lower-affinity (higher K_m) substrate site, so it
# dominates at high cyt Ca²⁺ where PMCA saturates.
#
# Platelet NCX presence: NCX1 (SLC8A1) and NCX3 (SLC8A3) detected in
# Burkhart 2012 proteome. Functional contribution is contested
# (Sage & Rink 1985 reported limited activity; later work argues NCX
# contributes 10–30 % of platelet Ca²⁺ extrusion at peak). Modelled
# here as a plausible secondary extruder — see
# `dissertation-notes.md §7.3` for the uncertainty disclosure.
#
# Kinetic scheme: substrate Hill term × allosteric Ca²⁺-activation gate.
# The allosteric gate keeps NCX silent at rest (cyt = 100 nM) regardless
# of substrate kinetics — captures the regulatory Ca²⁺-binding site of
# real NCX. Only forward mode modelled (reverse mode requires cyt Na⁺
# state + membrane potential, both out of scope for v0.3).
# Values live in `reports/params/calcium-v0.5.toml [ncx.kinetics]`
# (issue #32 Phase 2 slice 4).
K_NCX = dict(_KINETICS['ncx']['kinetics'])


# ── GPCR cascade — P2Y1 + PAR1/4 → Gαq → PLCβ (issue #9) ─────────────────
# Explicit receptor → G-protein cascade. The model takes physiological
# agonist concentrations (thrombin in nM, ADP in µM) as input.
#
# Three Gq-coupled receptors in scope:
#   - P2Y1 (ADP, reversible binding; ~150 copies, Coller 1995)
#   - PAR1 (thrombin, high-affinity, irreversible proteolytic cleavage;
#     ~2500 copies)
#   - PAR4 (thrombin, low-affinity; ~500 copies)
#
# Total Gαq pool: 5 000 molecules per platelet (Mazet, Tindall, Gibbins
# & Fry 2020). Inactive Gαq is implicit (total - active).
#
# Out of scope: GPVI (collagen, PLCγ2 path — different from Gq cascade);
# P2Y12 (Gi-coupled, inhibitory — issue #10); receptor desensitisation /
# internalisation kinetics (lumped into a slow first-order decay).

# Values live in `reports/params/calcium-v0.5.toml [gpcr.*]`
# (issue #32 Phase 2 slice 5). Receptor/Gq calibration commentary
# (k_basal derivation from resting Gq fraction) lives in the TOML
# section header.
K_P2Y1 = dict(_KINETICS['gpcr']['p2y1'])
K_PAR1 = dict(_KINETICS['gpcr']['par1'])
K_PAR4 = dict(_KINETICS['gpcr']['par4'])
K_GQ   = dict(_KINETICS['gpcr']['gq'])

N_GQ_TOTAL = _KINETICS['gpcr']['gq_pool']['n_total']  # Mazet 2020 platelet Gαq


# ── Agonist forcing functions — stimulation inputs ───────────────────────
# Physiological agonist concentrations (thrombin in nM, ADP in µM, ATP_ex
# in µM). Each forcing function reads its peak from the module at call time
# (None sentinel pattern) so dose sweeps can reassign the constant or pass
# `peak_*=X` per call without monkey-patching captured defaults.

# Thrombin (PAR1/4 agonist). Standard stimulation = 1 nM peak (Dolan 2014).
# Clotting cascade gives a fast rise, sustained plateau, slow decay.
# Values live in `reports/params/calcium-v0.5.toml [agonists.thrombin]`
# (issue #32 Phase 2 slice 2).
THROMBIN_REST_NM     = _KINETICS['agonists']['thrombin']['rest_nM']
THROMBIN_PEAK_NM     = _KINETICS['agonists']['thrombin']['peak_nM']
THROMBIN_TAU_RISE_S  = _KINETICS['agonists']['thrombin']['tau_rise_s']
THROMBIN_T_PEAK_S    = _KINETICS['agonists']['thrombin']['t_peak_s']
THROMBIN_TAU_DECAY_S = _KINETICS['agonists']['thrombin']['tau_decay_s']

# ADP (P2Y1 agonist). Released from dense granules during activation,
# cleared by ectoNTPDases. Standard stimulation = 10 µM peak.
# Values live in `reports/params/calcium-v0.5.toml [agonists.adp]`
# (issue #32 Phase 2 slice 2).
ADP_REST_UM     = _KINETICS['agonists']['adp']['rest_uM']
ADP_PEAK_UM     = _KINETICS['agonists']['adp']['peak_uM']
ADP_TAU_RISE_S  = _KINETICS['agonists']['adp']['tau_rise_s']
ADP_T_PEAK_S    = _KINETICS['agonists']['adp']['t_peak_s']
ADP_TAU_DECAY_S = _KINETICS['agonists']['adp']['tau_decay_s']


def thrombin_nM(t, delay=0.0, peak_nM=None):
	"""Thrombin concentration timecourse (nM). `peak_nM=None` reads
	`THROMBIN_PEAK_NM` live from the module; `peak_nM=0` yields rest."""
	if peak_nM is None:
		peak_nM = THROMBIN_PEAK_NM
	t_eff = t - delay
	if t_eff <= 0:
		return THROMBIN_REST_NM
	rise  = 1.0 - np.exp(-t_eff / THROMBIN_TAU_RISE_S)
	decay = np.exp(-max(0.0, t_eff - THROMBIN_T_PEAK_S) / THROMBIN_TAU_DECAY_S)
	return THROMBIN_REST_NM + (peak_nM - THROMBIN_REST_NM) * rise * decay


def adp_uM(t, delay=0.0, peak_uM=None):
	"""ADP concentration timecourse (µM). `peak_uM=None` reads
	`ADP_PEAK_UM` live from the module; `peak_uM=0` yields rest."""
	if peak_uM is None:
		peak_uM = ADP_PEAK_UM
	t_eff = t - delay
	if t_eff <= 0:
		return ADP_REST_UM
	rise  = 1.0 - np.exp(-t_eff / ADP_TAU_RISE_S)
	decay = np.exp(-max(0.0, t_eff - ADP_T_PEAK_S) / ADP_TAU_DECAY_S)
	return ADP_REST_UM + (peak_uM - ADP_REST_UM) * rise * decay

# Number of monomers per Orai1 tetramer (CRAC channel pore-forming subunit).
ORAI_SUBUNITS_PER_CHANNEL = 4

# Number of monomers per STIM1 dimer (sensor unit that binds Orai). Used
# for total-monomer mass-balance accounting (free + STIM1_Ca + 2·STIM1_dim
# is conserved). The MWC and listener treat STIM1_dim as dimer particles.
STIM_MONOMERS_PER_DIMER = 2


# ── PI cycle / PLCβ — Phase 4 / issue #31 ────────────────────────────────
# IP3 is produced from PIP2 hydrolysis by Gq-activated PLCβ, following the
# framework of Mazet, Tindall, Gibbins & Fry 2020 *Sci. Rep.* 10:13889.
#
# Coarse-grained scheme (v0.3 scope; full PI/PI4P/PI45P2 chain is v0.4):
#
#   PLCb_i + Gq  ──(k_act)──►  PLCb_a                  (PLCβ activation)
#   PLCb_a      ──(k_inact)──►  PLCb_i                  (RGS / GTPase)
#   PLCb_a + PIP2  ──(k_cat)──►  PLCb_a + IP3 + DAG    (PIP2 hydrolysis)
#   PI pool     ──(k_resynth)──►  PIP2                  (lumped resynthesis)
#   IP3          ──(k_ip3_deg)──►  IP2/IP4 (out)        (3-kinase + 5-phosphatase)
#   DAG          ──(k_dag_deg)──►  PA      (out)        (DAG kinase)
#
# Mazet 2020 has 35 rate constants in the full PI cycle; here we use a
# reduced 5-rate-constant model with the rates *calibrated* to reproduce
# the Dolan Fig. S2 IP3 timecourse under standard Gq forcing. Direct
# transfer of Mazet's constants would require their full PI/PI4P chain,
# which is out of scope for v0.3.
# Values live in `reports/params/calcium-v0.5.toml [pi_cycle.*]`
# (issue #32 Phase 2 slice 6). Calibration commentary
# (PIP2-rescaling k_cat 2.26e-7 → 2.26e-8 from lab-book 2026-05-15;
# τ_IP3 ≈ 50 s anchoring to Dolan Fig. S2) lives in the TOML headers.
K_PLCB     = dict(_KINETICS['pi_cycle']['plcb'])
K_PI_CYCLE = dict(_KINETICS['pi_cycle']['metabolism'])


# ── Mitochondrial Ca²⁺ (MCU + NCLX) — issue #22 ──────────────────────────
# Real platelets have 2–5 mitochondria per cell. The mitochondrial Ca²⁺
# uniporter (MCU, inner-membrane channel) and the Na⁺/Ca²⁺ exchanger
# NCLX (efflux) together act as a fast cyt Ca²⁺ buffer during transients:
# mito takes up Ca²⁺ rapidly (Hill cooperativity) and releases it
# slowly (linear in [Ca²⁺]_mito over a ~minutes timescale).
#
# Sources (all in `source-info/calcium-papers/`):
#   - Ghatge et al. 2026 — MCU-/- platelets show elevated cyt Ca²⁺
#   - Ajanel et al. 2025 — MCU regulates ITAM-dependent activation
#   - Shehwar et al. 2025 — review of platelet mito-Ca²⁺ biology
#
# This bypasses the slow PMCA extrusion: Ca²⁺ goes cyt → mito (fast,
# during peak) → cyt (slow, over minutes) → DTS via SERCA → PMCA out
# (gradually). Without MCU, all the SOCE-imported Ca²⁺ must exit via
# the PMCA bottleneck.
# Values live in `reports/params/calcium-v0.5.toml [mito.kinetics]`
# (issue #32 Phase 2 slice 6).
K_MITO = dict(_KINETICS['mito']['kinetics'])



# ── MWC SOCE solver (Hoover & Lewis 2011 / Dolan 2014 eq. 3) ──────────────
# Statistical occupation factors for 4 binding sites and binomial coefficients
# c(4,i). Cooperativity factor a^(i(i-1)/2) gives the standard MWC form.
_BINOM4 = (1.0, 4.0, 6.0, 4.0, 1.0)


def _mwc_open_fraction(stim2_p, n_orai, max_iter=20, tol=1e-6):
	"""Solve the Hoover/Dolan MWC equilibrium for the Orai open fraction.

	Parameters
	----------
	stim2_p : float
		Total STIM2 dimers in the puncta region (count). Sum of free Sf and
		(STIM2)p bound to Orai across all CSi/OSi states.
	n_orai : float
		Total tetrameric Orai channels available for binding (count).

	Returns
	-------
	po : float
		Channel open probability (open channels / total channels).
	sf : float
		Free (unbound) STIM2 dimer count at MWC equilibrium.
	"""
	if stim2_p <= 0.0 or n_orai <= 0.0:
		return K_MWC['L'] / (1.0 + K_MWC['L']), max(stim2_p, 0.0)

	L = K_MWC['L']
	Ka = K_MWC['Ka']
	f = K_MWC['f']
	a = K_MWC['a']
	# a^(i(i-1)/2) for i=0..4 — cumulative cooperativity.
	a_pow = tuple(a ** (i * (i - 1) / 2.0) for i in range(5))
	f_pow = tuple(f ** i for i in range(5))

	def occupancy(sf):
		"""Returns (open_frac_per_C, total_per_C, bound_count_per_C)."""
		x = Ka * sf
		x_pow = (1.0, x, x * x, x ** 3, x ** 4)
		total = 0.0
		open_ = 0.0
		bound = 0.0
		for i in range(5):
			cs_i = _BINOM4[i] * a_pow[i] * x_pow[i]
			os_i = cs_i * L * f_pow[i]
			total += cs_i + os_i
			open_ += os_i
			bound += i * (cs_i + os_i)
		return open_, total, bound

	# Solve mass balance by Newton-like fixed-point iteration:
	#   bound_count(sf) = (bound_per_C / total_per_C) × n_orai
	#   sf = stim2_p − bound_count(sf)
	# We iterate sf, clamped to [0, stim2_p].
	sf = stim2_p
	for _ in range(max_iter):
		open_per_C, total_per_C, bound_per_C = occupancy(sf)
		bound_count = (bound_per_C / total_per_C) * n_orai
		sf_new = max(0.0, min(stim2_p, stim2_p - bound_count))
		if abs(sf_new - sf) <= tol * max(1.0, sf):
			sf = sf_new
			break
		# Damped update for stability when bound saturates.
		sf = 0.5 * (sf + sf_new)

	open_per_C, total_per_C, _ = occupancy(sf)
	po = open_per_C / total_per_C if total_per_C > 0.0 else 0.0
	return po, sf


def _ode_rhs(t, y, t_sim_start, agonist_delay=0.0,
		thrombin_peak_nM=None, adp_peak_uM=None, atp_ex_peak_uM=None):
	"""Right-hand side of the calcium ODE.

	`y` carries integer-equivalent counts (continuous floats during
	integration). `t` is the sub-step time *within* this 1-second outer
	step; `t_sim_start` is the wall-clock time at which the outer step
	began.

	Agonist forcing is controlled by the three `*_peak_*` kwargs. Each is
	the peak concentration of its species during the activation transient
	(thrombin in nM, ADP and ATP_ex in µM). `None` (default) reads the
	module-level default constant; `0` yields a resting / un-stimulated
	sim where the receptor sees only its REST level.
	"""
	dy = np.zeros(N_SPECIES)

	# Concentrations (µM) for Ca²⁺ and IP3. IP3 is produced endogenously by
	# the PI cycle (PLCβ-driven; Phase 4 / #31), with the upstream agonist
	# peak kwargs controlling receptor stimulation that drives Gαq → PLCβ.
	ca_cyt = max(y[_IDX['CA2_CYT[c]']], 0.0) * _UM_PER_COUNT_CYT
	ca_dts = max(y[_IDX['CA2_DTS[dts]']], 0.0) * _UM_PER_COUNT_DTS
	ip3    = max(y[_IDX['IP3[c]']],     0.0) * _UM_PER_COUNT_CYT

	# PI cycle state reads
	pip2_count = max(y[_IDX['PIP2[c]']], 0.0)
	dag_count  = max(y[_IDX['DAG[c]']],  0.0)
	plcb_i     = max(y[_IDX['PLCb_inactive[c]']], 0.0)
	plcb_a     = max(y[_IDX['PLCb_active[c]']],   0.0)

	# Mitochondrial Ca²⁺ state read
	ca_mito_count = max(y[_IDX['CA2_MITO[m]']], 0.0)

	# GPCR cascade state reads
	p2y1_i = max(y[_IDX['P2Y1_inactive[pl]']], 0.0)
	p2y1_a = max(y[_IDX['P2Y1_active[pl]']],   0.0)
	par1_i = max(y[_IDX['PAR1_inactive[pl]']], 0.0)
	par1_a = max(y[_IDX['PAR1_active[pl]']],   0.0)
	par4_i = max(y[_IDX['PAR4_inactive[pl]']], 0.0)
	par4_a = max(y[_IDX['PAR4_active[pl]']],   0.0)
	gq_a   = max(y[_IDX['Gq_active[c]']],      0.0)

	# IP3R inactivation variable (count of non-inhibited channels, 0–N_IP3R).
	ip3r_h_count = max(y[_IDX['IP3R_h[dts]']], 0.0)

	se1   = max(y[_IDX['SERCA_E1[dts]']],    0.0)
	se2   = max(y[_IDX['SERCA_E2[dts]']],    0.0)
	se1c  = max(y[_IDX['SERCA_E1Ca[dts]']],  0.0)
	se1pc = max(y[_IDX['SERCA_E1PCa[dts]']], 0.0)
	se2pc = max(y[_IDX['SERCA_E2PCa[dts]']], 0.0)
	se2p  = max(y[_IDX['SERCA_E2P[dts]']],   0.0)

	pmca   = max(y[_IDX['PMCA[pl]']],    0.0)
	pmcaca = max(y[_IDX['PMCA_Ca[pl]']], 0.0)

	cam_free         = max(y[_IDX['CaM_free[c]']],        0.0)
	ca2_cam          = max(y[_IDX['Ca2_CaM[c]']],         0.0)
	ca4_cam          = max(y[_IDX['Ca4_CaM[c]']],         0.0)
	ca4_cam_pmca     = max(y[_IDX['Ca4_CaM_PMCA[pl]']],   0.0)
	ca4_cam_pmca_ca  = max(y[_IDX['Ca4_CaM_PMCA_Ca[pl]']],0.0)
	pmca_cam         = max(y[_IDX['PMCA_CaM[pl]']],       0.0)

	gsn_free = max(y[_IDX['GSN_free[c]']], 0.0)
	gsn_ca   = max(y[_IDX['GSN_Ca[c]']],   0.0)

	calr_free = max(y[_IDX['CALR_free[dts]']], 0.0)
	calr_ca   = max(y[_IDX['CALR_Ca[dts]']],   0.0)
	calr_p_free = max(y[_IDX['CALR_P_free[dts]']], 0.0)
	calr_p_ca   = max(y[_IDX['CALR_P_Ca[dts]']],   0.0)

	hsp90_m_free = max(y[_IDX['HSP90B1_M_free[dts]']], 0.0)
	hsp90_m_ca   = max(y[_IDX['HSP90B1_M_Ca[dts]']],   0.0)
	hsp90_l_free = max(y[_IDX['HSP90B1_L_free[dts]']], 0.0)
	hsp90_l_ca   = max(y[_IDX['HSP90B1_L_Ca[dts]']],   0.0)
	bip_free     = max(y[_IDX['BiP_free[dts]']],       0.0)
	bip_ca       = max(y[_IDX['BiP_Ca[dts]']],         0.0)
	crec_free    = max(y[_IDX['CREC_free[dts]']],      0.0)
	crec_ca      = max(y[_IDX['CREC_Ca[dts]']],        0.0)

	p2x1_c = max(y[_IDX['P2X1[pl]']],   0.0)
	p2x1_o = max(y[_IDX['P2X1_O[pl]']], 0.0)
	p2x1_d = max(y[_IDX['P2X1_D[pl]']], 0.0)

	st_free = max(y[_IDX['STIM1_free[dts]']], 0.0)
	st_ca   = max(y[_IDX['STIM1_Ca[dts]']],   0.0)
	st_dim  = max(y[_IDX['STIM1_dim[dts]']],  0.0)

	orai_total = max(y[_IDX['ORAI1[pl]']], 0.0)

	# ── IP3R deYoung-Keizer 1992 / Li-Rinzel 1994 ─────────────────────────
	h = ip3r_h_count / N_IP3R

	# Quasi-steady activation and slow inactivation ODE.
	m_inf = (ip3 / (ip3 + K_DYK['d1'])) * (ca_cyt / (ca_cyt + K_DYK['d5']))
	po_channel = (m_inf ** 4) * h

	dh_dt = K_DYK['a2'] * (K_DYK['d2'] - (ca_cyt + K_DYK['d2']) * h)
	dy[_IDX['IP3R_h[dts]']] += dh_dt * N_IP3R

	# Ca²⁺ flux via the Nernst form (Dolan 2014 eq. 4 / Purvis eq. 13).
	# N_IP3R = 1 328 channels (Burkhart 2012 / Dolan Table S1 convention).
	# Flux is zero when either compartment is empty to avoid phantom currents
	# from a degenerate Nernst potential.
	if ca_cyt > 0.0 and ca_dts > 0.0:
		e_ca_im_v = RT_OVER_zF_V * math.log(ca_dts / ca_cyt)
		driving_v = V_IM_V - e_ca_im_v
		flux_ip3r_ions_s = (
			-GAMMA_IP3R_S * N_IP3R * po_channel * driving_v * NA_OVER_zF
		)
	else:
		flux_ip3r_ions_s = 0.0
	# Positive = into cytosol (DTS→cyt direction gives negative driving_v,
	# leading minus flips sign to positive).
	dy[_IDX['CA2_CYT[c]']]   += +flux_ip3r_ions_s
	dy[_IDX['CA2_DTS[dts]']] += -flux_ip3r_ions_s

	# ── SERCA cycle (mass-action, 2 Ca²⁺ per turnover) ───────────────
	v_shuttle = K_SERCA['k_shuttle_f'] * se2 - K_SERCA['k_shuttle_r'] * se1
	v_bind    = K_SERCA['k_bind_f'] * se1 * (ca_cyt ** 2) - K_SERCA['k_bind_r'] * se1c
	v_phos    = K_SERCA['k_phos_f'] * se1c - K_SERCA['k_phos_r'] * se1pc
	v_conf    = K_SERCA['k_conf_f'] * se1pc - K_SERCA['k_conf_r'] * se2pc
	v_release = K_SERCA['k_release_f'] * se2pc - K_SERCA['k_release_r'] * se2p * (ca_dts ** 2)
	v_dephos  = K_SERCA['k_dephos_f'] * se2p - K_SERCA['k_dephos_r'] * se2

	dy[_IDX['SERCA_E1[dts]']]    += +v_shuttle - v_bind
	dy[_IDX['SERCA_E2[dts]']]    += -v_shuttle + v_dephos
	dy[_IDX['SERCA_E1Ca[dts]']]  += +v_bind - v_phos
	dy[_IDX['SERCA_E1PCa[dts]']] += +v_phos - v_conf
	dy[_IDX['SERCA_E2PCa[dts]']] += +v_conf - v_release
	dy[_IDX['SERCA_E2P[dts]']]   += +v_release - v_dephos

	# Each SERCA turnover removes 2 Ca²⁺_cyt (at the bind step) and
	# delivers 2 Ca²⁺_dts (at the release step). Atoms are temporarily held
	# in the enzyme-bound states between bind and release.
	dy[_IDX['CA2_CYT[c]']]   += -2.0 * v_bind
	dy[_IDX['CA2_DTS[dts]']] += +2.0 * v_release

	# ── CaM Ca²⁺ binding (Caride 2007 steps 6–7) ────────────────────
	# Step 6: CaM + 2 Ca²⁺ ⇌ Ca₂·CaM  (slow N-lobe)
	v_cam_bind1 = K_CAM['k6'] * cam_free * (ca_cyt ** 2) - K_CAM['k6r'] * ca2_cam
	# Step 7: Ca₂·CaM + 2 Ca²⁺ ⇌ Ca₄·CaM  (fast C-lobe)
	v_cam_bind2 = K_CAM['k7'] * ca2_cam * (ca_cyt ** 2) - K_CAM['k7r'] * ca4_cam

	dy[_IDX['CaM_free[c]']] += -v_cam_bind1
	dy[_IDX['Ca2_CaM[c]']]  += +v_cam_bind1 - v_cam_bind2
	dy[_IDX['Ca4_CaM[c]']]  += +v_cam_bind2
	# CaM-bound Ca²⁺ is removed from the free cytosolic pool (buffering effect).
	dy[_IDX['CA2_CYT[c]']]  += -2.0 * v_cam_bind1 - 2.0 * v_cam_bind2

	# ── Coarse-grained cytosolic Ca²⁺ buffer (gelsolin proxy) ────────
	# See K_GSN comment block above for the biology / scope disclosure
	# and the Phase 2 calibration rationale.
	v_gsn = K_GSN['k_on'] * gsn_free * ca_cyt - K_GSN['k_off'] * gsn_ca
	dy[_IDX['GSN_free[c]']] += -v_gsn
	dy[_IDX['GSN_Ca[c]']]   += +v_gsn
	dy[_IDX['CA2_CYT[c]']]  += -v_gsn

	# ── Calreticulin DTS Ca²⁺ buffer (Phase 2 / #28) ─────────────────
	# Two sites per CALR molecule: 25 low-affinity C-domain (fast) +
	# 1 high-affinity P-domain (slow). See K_CALR / K_CALR_P blocks.
	v_calr = K_CALR['k_on'] * calr_free * ca_dts - K_CALR['k_off'] * calr_ca
	dy[_IDX['CALR_free[dts]']] += -v_calr
	dy[_IDX['CALR_Ca[dts]']]   += +v_calr
	dy[_IDX['CA2_DTS[dts]']]   += -v_calr

	v_calr_p = K_CALR_P['k_on'] * calr_p_free * ca_dts - K_CALR_P['k_off'] * calr_p_ca
	dy[_IDX['CALR_P_free[dts]']] += -v_calr_p
	dy[_IDX['CALR_P_Ca[dts]']]   += +v_calr_p
	dy[_IDX['CA2_DTS[dts]']]     += -v_calr_p

	# ── HSP90B1 — medium-affinity sites (Kd = 2 µM; slow-release floor) ──
	v_hsp90_m = K_HSP90B1_M['k_on'] * hsp90_m_free * ca_dts \
				- K_HSP90B1_M['k_off'] * hsp90_m_ca
	dy[_IDX['HSP90B1_M_free[dts]']] += -v_hsp90_m
	dy[_IDX['HSP90B1_M_Ca[dts]']]   += +v_hsp90_m
	dy[_IDX['CA2_DTS[dts]']]        += -v_hsp90_m

	# ── HSP90B1 — low-affinity sites (Kd = 600 µM; fast equilibrium) ─────
	v_hsp90_l = K_HSP90B1_L['k_on'] * hsp90_l_free * ca_dts \
				- K_HSP90B1_L['k_off'] * hsp90_l_ca
	dy[_IDX['HSP90B1_L_free[dts]']] += -v_hsp90_l
	dy[_IDX['HSP90B1_L_Ca[dts]']]   += +v_hsp90_l
	dy[_IDX['CA2_DTS[dts]']]        += -v_hsp90_l

	# ── BiP / HSPA5 — single low-affinity site (Kd = 500 µM) ─────────────
	v_bip = K_BIP['k_on'] * bip_free * ca_dts - K_BIP['k_off'] * bip_ca
	dy[_IDX['BiP_free[dts]']] += -v_bip
	dy[_IDX['BiP_Ca[dts]']]   += +v_bip
	dy[_IDX['CA2_DTS[dts]']]  += -v_bip

	# ── CREC pool (CALU + RCN1 + RCN2 lumped; Kd = 1 mM) ─────────────────
	v_crec = K_CREC['k_on'] * crec_free * ca_dts - K_CREC['k_off'] * crec_ca
	dy[_IDX['CREC_free[dts]']] += -v_crec
	dy[_IDX['CREC_Ca[dts]']]   += +v_crec
	dy[_IDX['CA2_DTS[dts]']]   += -v_crec

	# ── PMCA basal path (Caride 2007 steps 4–5) ──────────────────────
	v_pmca_bind = K_PMCA['k_on'] * pmca * ca_cyt - K_PMCA['k_off'] * pmcaca
	v_pmca_cat  = K_PMCA['k_cat'] * pmcaca

	dy[_IDX['PMCA[pl]']]    += -v_pmca_bind + v_pmca_cat
	dy[_IDX['PMCA_Ca[pl]']] += +v_pmca_bind - v_pmca_cat
	dy[_IDX['CA2_CYT[c]']]  += -v_pmca_bind
	# Catalysis ejects Ca²⁺ to the extracellular reservoir (not tracked).

	# ── PMCA CaM-activated path (Caride 2007 steps 8–12) ────────────
	# Full 5-state Caride scheme. Restored 2026-05-07 after Phase 0 audit
	# found step 12 (slow CaM dissociation, k₁₂ = 0.033 s⁻¹, τ ~ 30 s) was
	# missing and step 11 was therefore disabled to work around the
	# resulting accumulation bug. Now both are present.
	#
	# Mass balance (full cycle):
	#   d/dt(PMCA + Ca4_CaM_PMCA + Ca4_CaM_PMCA_Ca + PMCA_CaM) = 0  ✓
	#   d/dt(CaM_free + Ca4_CaM + Ca4_CaM_PMCA + Ca4_CaM_PMCA_Ca + PMCA_CaM)
	#       (plus Ca2_CaM via the upstream ladder) = 0  ✓
	#
	# Step 8: PMCA + Ca₄·CaM ⇌ Ca₄·CaM·PMCA
	# Ca₄·CaM in µM; PMCA as count — gives rate in count·s⁻¹.
	ca4_cam_uM = ca4_cam * _UM_PER_COUNT_CYT
	v_cam_bind_pmca = (
		K_CAM_PMCA['k8'] * pmca * ca4_cam_uM
		- K_CAM_PMCA['k8r'] * ca4_cam_pmca
	)
	# Step 9: Ca₄·CaM·PMCA + Ca²⁺ ⇌ Ca₄·CaM·PMCA·Ca
	v_cam_pmca_bind = (
		K_CAM_PMCA['k9'] * ca4_cam_pmca * ca_cyt
		- K_CAM_PMCA['k9r'] * ca4_cam_pmca_ca
	)
	# Step 10: Ca₄·CaM·PMCA·Ca → Ca₄·CaM·PMCA + Ca²⁺_ex  (activated extrusion)
	# Recycles the empty Ca₄·CaM·PMCA complex back to step 9 — PMCA is not consumed.
	v_cam_pmca_cat = K_CAM_PMCA['k10'] * ca4_cam_pmca_ca
	# Step 11: Ca₄·CaM·PMCA ⇌ PMCA·CaM + 4 Ca²⁺_cyt
	# k11r in µM⁻⁴·s⁻¹ — at rest cyt ≈ 0.1 µM the reverse is ~7e-8 s⁻¹
	# and effectively zero, so this is a near-irreversible Ca²⁺ release.
	v_cam_pmca_release = (
		K_CAM_PMCA['k11'] * ca4_cam_pmca
		- K_CAM_PMCA['k11r'] * pmca_cam * (ca_cyt ** 4)
	)
	# Step 12: PMCA·CaM → PMCA + CaM_free (slow CaM dissociation, no reverse)
	v_cam_dissoc = K_CAM_PMCA['k12'] * pmca_cam

	dy[_IDX['PMCA[pl]']]             += -v_cam_bind_pmca + v_cam_dissoc
	dy[_IDX['Ca4_CaM[c]']]           += -v_cam_bind_pmca  # CaM leaves free pool on step 8 fwd
	dy[_IDX['Ca4_CaM_PMCA[pl]']]     += +v_cam_bind_pmca - v_cam_pmca_bind + v_cam_pmca_cat - v_cam_pmca_release
	dy[_IDX['Ca4_CaM_PMCA_Ca[pl]']]  += +v_cam_pmca_bind - v_cam_pmca_cat
	dy[_IDX['PMCA_CaM[pl]']]         += +v_cam_pmca_release - v_cam_dissoc
	dy[_IDX['CaM_free[c]']]          += +v_cam_dissoc       # CaM returns to free pool on step 12
	# Step 9 fwd removes one Ca²⁺_cyt; step 10 ejects it extracellularly.
	# Step 11 fwd releases 4 Ca²⁺ to cyt (reverse consumes 4).
	dy[_IDX['CA2_CYT[c]']]           += -v_cam_pmca_bind + 4.0 * v_cam_pmca_release

	# ── STIM1 cycle (mass-action; calibrated to Dolan IC detailed balance) ─
	v_stim1_release = (
		K_STIM['k_release_f'] * st_ca
		- K_STIM['k_release_r'] * st_free * ca_dts
	)
	v_dim = K_STIM['k_dim_f'] * st_free * st_free - K_STIM['k_dim_r'] * st_dim

	dy[_IDX['STIM1_Ca[dts]']]   += -v_stim1_release
	dy[_IDX['STIM1_free[dts]']] += +v_stim1_release - 2.0 * v_dim
	dy[_IDX['STIM1_dim[dts]']]  += +v_dim
	# Ca²⁺ released from STIM EF-hand returns to the free DTS pool.
	dy[_IDX['CA2_DTS[dts]']]    += v_stim1_release

	# ── SOCE: Dolan eq. 2 (puncta entry) + MWC equilibrium + eq. 4 ──────
	# qp Hill function of [Ca²⁺]_cyt drives STIM2 dimers into puncta.
	if ca_cyt > 0.0:
		hill = (ca_cyt ** PUNCTA['n']
				/ (PUNCTA['KM_uM'] ** PUNCTA['n'] + ca_cyt ** PUNCTA['n']))
	else:
		hill = 0.0
	qp = PUNCTA['alpha'] * hill + PUNCTA['baseline']
	# STIM1_dim count in our state vector = dimer particles directly
	# (Dolan Table S1 lists 11 dimers; we now match that convention).
	stim2_p = qp * st_dim
	# Total Orai *channels* (tetramers).
	n_orai_channels = orai_total / ORAI_SUBUNITS_PER_CHANNEL
	# Solve MWC for channel-level open probability.
	po_orai, _sf = _mwc_open_fraction(stim2_p, n_orai_channels)

	# SOC current via Eq. 4: I = γ · N · Po · (ψ_PM − E_Ca,PM) · NA/(zF).
	# Only computed when there is extracellular Ca²⁺ to flow in. Under the
	# Dolan Fig. 4 EDTA / no-extracellular-Ca condition (CA_EX_UM = 0) both
	# the SOCE current and the basal PM leak are physically zero — both are
	# Ca²⁺ inflows from outside, and there is no outside Ca²⁺ to source.
	if CA_EX_UM > 0.0 and ca_cyt > 0.0:
		e_ca_pm_v = RT_OVER_zF_V * math.log(CA_EX_UM / ca_cyt)
		driving_pm_v = V_PM_V - e_ca_pm_v
		soce_ions_s = (
			-GAMMA_SOC_S * n_orai_channels * po_orai * driving_pm_v * NA_OVER_zF
		)
		dy[_IDX['CA2_CYT[c]']] += soce_ions_s

		# Basal plasma-membrane Ca²⁺ leak — compensates PMCA outflow at rest
		# so the resting cyt steady state sits at ~100 nM (lab-book 2026-05-05).
		dy[_IDX['CA2_CYT[c]']] += J_PM_LEAK_IONS_S

		# P2X1 Ca²⁺ entry — fast ionotropic channel gated on extracellular
		# ATP (forced curve), distinct from SOCE which is store-operated.
		# P2X1 is the +Ca_ex–specific entry pathway that distinguishes the
		# two Phase 3 conditions in real biology (Dolan ~100 nM SOCE
		# differential). See K_P2X1 / atp_ex_forcing_uM block.
		flux_p2x1_ions_s = (
			-GAMMA_P2X1_S * p2x1_o * driving_pm_v * NA_OVER_zF
		)
		dy[_IDX['CA2_CYT[c]']] += flux_p2x1_ions_s

		# NB: P2X1 *state* transitions run unconditionally below, since the
		# channel can still cycle C → O → D in response to ATP regardless
		# of whether external Ca²⁺ is available to flow through it.

		# NCX (Na⁺/Ca²⁺ exchanger) — Ca²⁺ extrusion gated on extracellular
		# Ca²⁺ availability (needs Na⁺ gradient + somewhere for Ca²⁺ to go).
		# See K_NCX block above for biology and the Hill formulation.
		g_act = (ca_cyt / K_NCX['K_a']) ** K_NCX['h']
		g_act /= (1.0 + g_act)
		v_ncx = K_NCX['V_max'] * g_act * ca_cyt / (K_NCX['K_m'] + ca_cyt)
		dy[_IDX['CA2_CYT[c]']] -= v_ncx
	# The extracellular reservoir is treated as infinite (no debit).

	# ── P2X1 state transitions (always run, even when CA_EX = 0) ────
	# Vial & Evans 2002: P2X1 desensitises whether Ca²⁺ is present or not.
	atp_ex_um = atp_ex_forcing_uM(
		t_sim_start + t, delay=agonist_delay, peak_uM=atp_ex_peak_uM)
	v_p2x1_act   = K_P2X1['k_act']   * p2x1_c * atp_ex_um
	v_p2x1_close = K_P2X1['k_close'] * p2x1_o
	v_p2x1_des   = K_P2X1['k_des']   * p2x1_o
	v_p2x1_rec   = K_P2X1['k_rec']   * p2x1_d
	dy[_IDX['P2X1[pl]']]   += -v_p2x1_act + v_p2x1_close + v_p2x1_rec
	dy[_IDX['P2X1_O[pl]']] += +v_p2x1_act - v_p2x1_close - v_p2x1_des
	dy[_IDX['P2X1_D[pl]']] +=                              +v_p2x1_des - v_p2x1_rec

	# ── PI cycle: PLCβ-driven IP3 production from PIP2 — Phase 4 / #31 ──
	# Gq concentration comes from the GPCR cascade state variable.
	# Conversion factor 5 was calibrated so resting Gq_active ≈ 100 maps
	# to gq_um ≈ 0.1 µM, which holds basal IP3 at its ~50 nM Purvis target.
	gq_um = (gq_a / N_GQ_TOTAL) * 5.0

	# PLCβ activation cycle.
	v_plcb_act   = K_PLCB['k_act']   * plcb_i * gq_um
	v_plcb_inact = K_PLCB['k_inact'] * plcb_a

	# PLCβ-catalysed PIP2 hydrolysis → IP3 + DAG (in count/s).
	v_plcb_cat = K_PLCB['k_cat'] * plcb_a * pip2_count

	# PIP2 resynthesis from the PI/PI4P chain (lumped order-zero rate).
	v_pip2_resynth = K_PI_CYCLE['k_resynth']

	# IP3 degradation (5-phosphatase + 3-kinase, lumped).
	v_ip3_deg = K_PI_CYCLE['k_ip3_deg'] * y[_IDX['IP3[c]']]
	# DAG kinase (DAG → PA).
	v_dag_deg = K_PI_CYCLE['k_dag_deg'] * dag_count

	dy[_IDX['PLCb_inactive[c]']] += -v_plcb_act + v_plcb_inact
	dy[_IDX['PLCb_active[c]']]   += +v_plcb_act - v_plcb_inact
	dy[_IDX['PIP2[c]']]          += -v_plcb_cat + v_pip2_resynth
	dy[_IDX['IP3[c]']]           += +v_plcb_cat - v_ip3_deg
	dy[_IDX['DAG[c]']]           += +v_plcb_cat - v_dag_deg

	# ── GPCR cascade — receptors + Gαq cycle (issue #9) ─────────────
	# Agonist concentrations → receptor activation → Gαq exchange →
	# PLCβ activation (in the PI cycle block above).
	adp_um_now = adp_uM(
		t_sim_start + t, delay=agonist_delay, peak_uM=adp_peak_uM)
	thr_nm_now = thrombin_nM(
		t_sim_start + t, delay=agonist_delay, peak_nM=thrombin_peak_nM)

	# P2Y1: reversible ADP binding
	v_p2y1_on  = K_P2Y1['k_on']  * p2y1_i * adp_um_now
	v_p2y1_off = K_P2Y1['k_off'] * p2y1_a
	dy[_IDX['P2Y1_inactive[pl]']] += -v_p2y1_on + v_p2y1_off
	dy[_IDX['P2Y1_active[pl]']]   += +v_p2y1_on - v_p2y1_off

	# PAR1: thrombin cleavage (essentially irreversible) + one-way
	# internalization (v0.4.1 fix). Previously internalization went back
	# to the inactive surface pool, allowing slow thrombin tails to
	# re-cleave the same receptors indefinitely; PAR1_active stayed
	# elevated for ~2000s after thrombin clears. Real biology has
	# internalized PAR1 either degraded or recycled on hours timescale
	# (out of our sim window). Modelled as a one-way sink.
	v_par1_cleave = K_PAR1['k_cleave'] * par1_i * thr_nm_now
	v_par1_int    = K_PAR1['k_internalize'] * par1_a
	dy[_IDX['PAR1_inactive[pl]']]     += -v_par1_cleave
	dy[_IDX['PAR1_active[pl]']]       += +v_par1_cleave - v_par1_int
	dy[_IDX['PAR1_internalized[pl]']] += +v_par1_int

	# PAR4: low-affinity thrombin receptor, sustained response.
	# Same one-way internalization as PAR1 (v0.4.1 fix).
	v_par4_cleave = K_PAR4['k_cleave'] * par4_i * thr_nm_now
	v_par4_int    = K_PAR4['k_internalize'] * par4_a
	dy[_IDX['PAR4_inactive[pl]']]     += -v_par4_cleave
	dy[_IDX['PAR4_active[pl]']]       += +v_par4_cleave - v_par4_int
	dy[_IDX['PAR4_internalized[pl]']] += +v_par4_int

	# Gαq cycle: all active receptors share the same Gq pool.
	# Inactive Gq is implicit (total - active).
	gq_i_count = max(N_GQ_TOTAL - gq_a, 0.0)
	total_active_R = p2y1_a + par1_a + par4_a
	v_gq_act = (K_GQ['k_act_per_R'] * total_active_R + K_GQ['k_basal']) * gq_i_count
	v_gq_rgs = K_GQ['k_rgs'] * gq_a
	dy[_IDX['Gq_active[c]']] += v_gq_act - v_gq_rgs

	# ── Mitochondrial Ca²⁺ (MCU + NCLX) — issue #22 ──────────────────
	# MCU uptake: cooperative Hill kinetics on cyt Ca²⁺.
	# NCLX efflux: linear in mito Ca²⁺ (slow release).
	cyt_n = ca_cyt ** K_MITO['n_MCU']
	km_n  = K_MITO['K_MCU'] ** K_MITO['n_MCU']
	v_mcu  = K_MITO['V_max_MCU'] * cyt_n / (km_n + cyt_n)
	v_nclx = K_MITO['k_NCLX'] * ca_mito_count
	dy[_IDX['CA2_CYT[c]']]  += -v_mcu + v_nclx
	dy[_IDX['CA2_MITO[m]']] += +v_mcu - v_nclx

	return dy


class CalciumSignalling:
	"""Parameters and ODE driver for the CalciumDynamics process.

	Public surface matches the wcEcoli convention used by other process
	dataclasses (see TwoComponentSystem):

	  * `molecule_names` — array of bulk-molecule IDs the process views.
	  * `molecules_to_next_time_step(counts, dt, t_sim, ...)` — integrate
	    the ODE for one outer timestep and return integer count deltas +
	    estimated ATP cost. Agonist forcing is controlled by the three
	    optional peak concentrations; passing 0 yields a resting sim.
	"""

	def __init__(self, sim_data):
		self.molecule_names = np.array(MOLECULE_NAMES, dtype='U30')
		self.n_species = N_SPECIES
		# Per-species fractional-delta residual carried over between
		# timesteps. Without it, any species whose |dy/dt| < 0.5 ions/s
		# is frozen at its current integer count (rounding artifact):
		# the ODE evolves with float precision but the commit to bulk
		# molecules rounds to int per-step. For slow approach to
		# equilibrium (e.g. IP3 near rest, PIP2 refill), this stranded
		# the model at a non-equilibrium fixed point. Tracking the
		# unused fractional part and adding it to the next step's delta
		# preserves the ODE dynamics on average. Diagnosed in
		# lab-book-2026-05-12-pi-cycle-design.md §"IP3 stuck-at-205".
		self._residual = np.zeros(N_SPECIES)

	def molecules_to_next_time_step(self, counts, dt, t_sim,
			agonist_delay=0.0, thrombin_peak_nM=None,
			adp_peak_uM=None, atp_ex_peak_uM=None):
		"""Run one outer timestep of the calcium ODE.

		Agonist forcing is controlled by the three `*_peak_*` kwargs;
		`None` (default) uses the module-level peak constants, `0` gives
		a resting / un-stimulated sim. See `_ode_rhs` docstring.

		Returns
		-------
		delta_counts : ndarray of int64, shape (N_SPECIES,)
			Net change to apply to bulk molecule counts.
		atp_cost : int
			Number of ATP molecules consumed by SERCA + PMCA pumping
			during this step.
		"""
		y0 = counts.astype(float)

		sol = solve_ivp(
			_ode_rhs, (0.0, dt), y0,
			method='BDF',
			args=(t_sim, agonist_delay,
				thrombin_peak_nM, adp_peak_uM, atp_ex_peak_uM),
			atol=1e-3,    # counts; ~0.001 molecule precision
			rtol=1e-6,
			max_step=dt,
		)
		if not sol.success:
			raise RuntimeError(
				f'CalciumSignalling ODE failed at t_sim={t_sim}: {sol.message}')

		y_final = np.maximum(sol.y[:, -1], 0.0)
		# Carry over fractional residuals across timesteps so that
		# sub-1-ion-per-second fluxes still integrate correctly over
		# many steps. The exact float delta is integer-rounded for
		# the commit, and the unused fraction is added to the next
		# step's accumulator. See `__init__` for the diagnosis.
		fractional_delta = y_final - y0 + self._residual
		delta = np.round(fractional_delta).astype(np.int64)
		self._residual = fractional_delta - delta

		# ATP cost: SERCA delivers 2 Ca²⁺ per ATP, PMCA delivers 1 per ATP.
		# Estimated from Ca²⁺ flux integrated over the step. Net Ca²⁺_cyt
		# loss to SERCA is paired with a net rise in CA2_DTS; net loss to
		# PMCA we infer from how many turnovers PMCA_Ca did (its catalytic
		# step).
		dts_gain = max(delta[_IDX['CA2_DTS[dts]']], 0)
		serca_atp = dts_gain // 2
		# PMCA ATP: 1 per turnover on both the basal path (step 5) and the
		# CaM-activated path (step 10).
		pmca_atp = max(int(
			K_PMCA['k_cat'] * counts[_IDX['PMCA_Ca[pl]']] * dt
			+ K_CAM_PMCA['k10'] * counts[_IDX['Ca4_CaM_PMCA_Ca[pl]']] * dt
		), 0)
		atp_cost = serca_atp + pmca_atp

		return delta, atp_cost
