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
  * IP3 forcing (Dolan 2014 Fig. S2 shape).

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
	# P2X1 ATP-gated cation channel (Phase 2.5; see K_P2X1 block).
	'P2X1[pl]',     # closed (resting)
	'P2X1_O[pl]',   # open (conducting)
	'P2X1_D[pl]',   # desensitised
)
# Index lookups for readability inside the rate function.
_IDX = {name: i for i, name in enumerate(MOLECULE_NAMES)}
N_SPECIES = len(MOLECULE_NAMES)


# ── Resting concentrations / IP3 forcing ──────────────────────────────────
# Source: Dolan & Diamond 2014 main-text + Fig. S2 fit.
CA_EX_UM = 1200.0       # extracellular Ca²⁺, fixed reservoir
IP3_REST_UM = 0.05      # cytosolic IP3 baseline (50 nM)

# IP3 forcing parameters (Dolan 2014 Fig. S2 shape; v0.3 replaces this with
# the upstream P2Y1/Gq/PLCβ cascade producing IP3 endogenously).
IP3_FOLD = 5.5
IP3_T_PEAK = 3.0
IP3_TAU_RISE = 3.0
IP3_TAU_DECAY = 60.0


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
K_DYK = {
	'd1':  0.13,     # IP3 activation half-saturation   (µM)      b₁/a₁ = 52/400
	'd2':  1.049,    # Ca²⁺ inhibition half-saturation  (µM)      b₄/a₄
	'd5':  0.08234,  # Ca²⁺ activation half-saturation  (µM)      b₅/a₅
	'a2':  0.2,      # Ca²⁺ inhibition on-rate          (µM⁻¹·s⁻¹) — sets τ_h
}

# Total IP3R channels: Burkhart 2012 ITPR2 count, Dolan 2014 Table S1.
N_IP3R = 1328

# IP3R Ca²⁺ flux: Nernst-based Purvis 2008 eq. 13 / Dolan 2014 eq. 4
#   I = γ · N · Po · (NA/(zF)) · (ψ_IM − E_Ca,IM)
#
# ⚠ CALIBRATION-COUPLED PARAMETER — read before changing.
# γ_IP3R is *not* an independently measured single-channel conductance.
# It is the value that balances `K_SERCA` (above) at the Dolan 2014
# platelet resting state (cyt = 100 nM, DTS = 250 µM). If you change
# any of `k_bind_f`, `k_phos_f`, `k_conf_f`, `k_release_f` or the SERCA
# total copy number, you MUST re-derive γ_IP3R from the new SERCA
# cycle flux. See `reports/dissertation-notes.md §3.1` for the
# coupling diagram and §3.2 for the open question of whether the
# Purvis 2008 SERCA rate constants themselves over-estimate the
# SERCA3b pump rate at low cyt Ca²⁺.
#
# Derivation: SERCA 6-state cycle steady-state flux at cyt = 100 nM,
# DTS = 250 µM, solved analytically as a linear system, gives J =
# 112 570 Ca²⁺ ions/s (56 285 cycles/s × 2 Ca²⁺/cycle). Setting IP3R
# resting flux equal to this and inverting the Nernst flux formula:
#   γ_required = J / (N · Po · |driving| · NA/zF)
#              = 112 570 / (1 328 · 4.91×10⁻⁴ · 0.1605 · 3.122×10¹⁸)
#              = 0.344 pS → rounded to 0.35 pS.
#
# Biological plausibility: Bezprozvanny 1991 and Mak & Foskett 1997
# measured effective IP3R Ca²⁺ conductance in cellular conditions at
# ~0.05–0.5 pS, with K⁺ carrying most of the unitary current. Our 0.35
# pS sits within that range. The historical 10 pS (Zschauer 1988, via
# Purvis 2008) was a bilayer measurement under symmetric high Ca²⁺
# where K⁺ contributes negligibly to current and is not transferable.
GAMMA_IP3R_S = 0.175e-12         # 0.175 pS = calibrated Ca²⁺ conductance, A/V


# ── SERCA cycle (Purvis 2008 Table 1, Dode 2002 isoform 3b kinetics) ──────
# Primary-source values restored. Earlier calibration reduced k_bind_f by
# ~470× to compensate for IP3R Po and flux bugs; with Po⁴ + Nernst the
# Purvis Vmax balances the corrected IP3R leak (~1.18×10⁵ ions/s) at rest.
#
# ⚠ CALIBRATION-COUPLED — any change to these rate constants requires
# re-deriving GAMMA_IP3R_S (above) to restore the resting-state flux
# balance. See `reports/dissertation-notes.md §3.1`.
# ⚠ OPEN BIOLOGY QUESTION — these constants imply ~4.7 cycles/s per
# pump at cyt = 100 nM, which is ~2–5× higher than the SERCA3b Vmax /
# Km literature (Inesi 1985; Nishi 1992; Dode 2002 itself: Vmax ~30–50
# cycles/s saturating, Km ~0.7–1.1 µM, so v/Vmax at 100 nM ≈ 2%). The
# Purvis 2008 rate constants appear to over-estimate the platelet SERCA
# pump rate at resting Ca²⁺. v0.3+ should re-derive from primary
# sources. See `reports/dissertation-notes.md §3.2`.
K_SERCA = {
	'k_shuttle_f':  600.0,    # E2 → E1                        (s⁻¹)
	'k_shuttle_r':  600.0,    # E1 → E2                        (s⁻¹)
	'k_bind_f':     500.0,    # E1 + 2 Ca²⁺_cyt → E1·Ca²⁺      (µM⁻²·s⁻¹)
	'k_bind_r':      10.0,    # reverse                        (s⁻¹)
	'k_phos_f':     700.0,    # E1·Ca → E1P·Ca                 (s⁻¹)
	'k_phos_r':       5.0,
	'k_conf_f':     600.0,    # E1P·Ca ⇌ E2P·Ca                (s⁻¹)
	'k_conf_r':      50.0,
	'k_release_f': 1000.0,    # E2P·Ca → E2P + 2 Ca²⁺_dts      (s⁻¹)
	'k_release_r':  4.0e-3,   # reverse (µM⁻²·s⁻¹; 4e9 M⁻²s⁻¹)
	'k_dephos_f':   500.0,    # E2P → E2                       (s⁻¹)
	'k_dephos_r':     1.0,
}


# ── PMCA4b basal path (Caride 2007 Table 3 steps 4–5) ────────────────────
# Steps 4–5 are unchanged; the CaM-activated path (steps 8–11) is below.
K_PMCA = {
	'k_on':   10.0,    # PMCA + Ca²⁺ ⇌ PMCA·Ca   (µM⁻¹·s⁻¹)  step 4 fwd
	'k_off':  50.0,    # reverse                  (s⁻¹)        step 4 rev
	'k_cat':   5.5,    # PMCA·Ca → PMCA + Ca²⁺_ex (s⁻¹)        step 5 (basal turnover)
}

# ── CaM Ca²⁺ binding (Caride 2007 Table 3 steps 6–7) ─────────────────────
# Two-lobe cooperative scheme: slow N-lobe (step 6) then fast C-lobe (step 7).
# Ca²⁺ concentrations in µM; rates in µM⁻²·s⁻¹ (forward) or s⁻¹ (reverse).
K_CAM = {
	'k6':    2.669,   # CaM + 2 Ca²⁺ → Ca₂·CaM  (µM⁻²·s⁻¹)  step 6 fwd
	'k6r':   2.682,   # reverse                  (s⁻¹)        step 6 rev
	'k7':  170.4,     # Ca₂·CaM + 2 Ca²⁺ → Ca₄·CaM (µM⁻²·s⁻¹) step 7 fwd
	'k7r':   1.551,   # reverse                  (s⁻¹)        step 7 rev
}

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
#
# See `reports/dissertation-notes.md §1.1` for the literature gap and
# v0.3+ plan to split this into explicit gelsolin / annexin / Ca-ATP.
K_GSN = {
	'k_on':  100.0,    # GSN_site + Ca²⁺ → GSN_site·Ca  (µM⁻¹·s⁻¹) — fast EF-hand binding
	'k_off': 100.0,    # reverse                        (s⁻¹)       — Kd = 1.0 µM
}

# Effective Ca²⁺-binding sites (Phase 2 cyt+DTS-coupled calibration).
# Biological gelsolin: ~100 000 copies × ~5 sites = 500 000 sites (Burkhart
# 2012; Yin & Stossel 1979). N_GSN below is calibrated against Phase 3
# peak heights with CALR active — see lab book for the iteration log.
N_GSN = 800_000


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
K_CALR = {
	'k_on':    1.0,    # CALR_site + Ca²⁺ → CALR_site·Ca  (µM⁻¹·s⁻¹) — fast equilibrium
	'k_off': 1000.0,   # reverse                          (s⁻¹)       — Kd = 1.0 mM
}

# Total CALR Ca²⁺-binding sites: 20 324 CALR × 25 C-domain sites.
N_CALR = 508_100

# CALR high-affinity P-domain: 1 site per CALR molecule, Kd ~ 1 µM, slow
# release kinetics (k_off ~ 1 s⁻¹). Source: Baksh & Michalak 1991,
# Vassilakos 1998. At resting DTS [Ca²⁺] = 250 µM this site is always
# saturated (occupancy > 99.6 %), so it doesn't matter much at rest — but
# during IP3R-driven DTS depletion, the slow release rate means these
# 20 324 Ca²⁺ ions take ~1 s to liberate after free [Ca²⁺]_DTS drops
# below the Kd. This adds a small "delayed reservoir" that smooths DTS
# recovery without preventing the transient depletion. See lab book.
K_CALR_P = {
	'k_on':    1.0,    # CALR_P + Ca²⁺ → CALR_P·Ca   (µM⁻¹·s⁻¹)
	'k_off':   1.0,    # reverse                     (s⁻¹)       — Kd = 1.0 µM
}
N_CALR_P = 20_324


# ── PMCA4b CaM-activated path (Caride 2007 Table 3 steps 8–12) ──────────
# Ca₄·CaM binds free PMCA (step 8), the complex binds and pumps Ca²⁺ with
# ~5× higher k_cat than basal (step 10 vs step 5). Step 11 dissociates 4
# Ca²⁺ from the active complex (Ca₄·CaM·PMCA → PMCA·CaM + 4 Ca²⁺_cyt) and
# step 12 is the slow CaM dissociation (PMCA·CaM → PMCA + CaM, k12 = 0.033
# s⁻¹, τ ~ 30 s). Both 11 and 12 are required for mass conservation: with
# 11 active but 12 absent, PMCA accumulates dead-end in PMCA·CaM (the bug
# previously worked around by omitting step 11 entirely; restored
# 2026-05-07 after Phase 0 audit found Caride k₁₂ missing).
K_CAM_PMCA = {
	'k8':   0.2,       # PMCA + Ca₄·CaM → Ca₄·CaM·PMCA  (µM⁻¹·s⁻¹) step 8 fwd
	'k8r':  8.0e-4,    # reverse                         (s⁻¹)       step 8 rev
	'k9':  50.0,       # Ca₄·CaM·PMCA + Ca²⁺ ⇌ Ca₄·CaM·PMCA·Ca (µM⁻¹·s⁻¹) step 9
	'k9r': 10.0,       # reverse                         (s⁻¹)
	'k10': 30.0,       # Ca₄·CaM·PMCA·Ca → Ca₄·CaM·PMCA + Ca²⁺_ex (s⁻¹) step 10
	'k11':  10.0,      # Ca₄·CaM·PMCA → PMCA·CaM + 4 Ca²⁺_cyt (s⁻¹) step 11 fwd
	'k11r':  7.332e-4, # reverse (µM⁻⁴·s⁻¹)                          step 11 rev
	'k12':  0.033,     # PMCA·CaM → PMCA + CaM_free (s⁻¹) step 12 (slow); no rev
}


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
K_P2X1 = {
	'k_act':    30.0,    # closed + ATP → open       (µM⁻¹·s⁻¹) — fast ATP binding
	'k_close':   5.0,    # open → closed             (s⁻¹)       — ATP unbinds
	'k_des':    10.0,    # open → desensitised       (s⁻¹)       — τ_des ≈ 100 ms
	'k_rec':     0.03,   # desensitised → closed     (s⁻¹)       — τ_rec ≈ 30 s
}

# Total P2X1 functional channels (trimers; mass per trimer = 3 × 45 kDa).
N_P2X1 = 1_000

# P2X1 Ca²⁺-specific effective conductance.
# Single-channel current ~0.5–1 pA at -60 mV, Ca²⁺ fraction ~5–10 %,
# so effective Ca²⁺-specific γ ≈ 10–50 fS per channel. Starting from
# 0.01 pS — calibration anchor for Phase 3 SOCE-differential target
# (Dolan ~100 nM). See lab book.
GAMMA_P2X1_S = 0.0006e-12   # 0.6 fS Ca²⁺-specific conductance, A/V — calibrated


# ── Extracellular ATP forcing (drives P2X1) ───────────────────────────────
# Parallels the IP3 forcing curve (Dolan Fig. S2). Real biology: ATP is
# released from dense granules during platelet activation, reaching
# 1–10 µM in the extracellular space near a forming thrombus, then
# cleared by ectonucleotidases (CD39) over tens of seconds.
#
# Same ip3_delay / ip3_forced gating as IP3 — when stimulus is off,
# ATP_ex stays at its near-zero resting level.
ATP_EX_REST_UM = 0.0          # exactly zero — CD39 ectonucleotidase keeps
                              # local extracellular ATP near zero at rest
                              # (any small baseline leaks P2X1 over hundreds
                              # of seconds, overfilling the DTS).
ATP_EX_PEAK_UM = 10.0         # 10 µM peak during activation
ATP_EX_TAU_RISE = 0.5         # s — fast rise (dense granule secretion)
ATP_EX_T_PEAK = 1.0           # s — peak time
ATP_EX_TAU_DECAY = 30.0       # s — ectonucleotidase clearance


def atp_ex_forcing_uM(t, delay=0.0):
	"""Plateau-decay approximation for extracellular ATP during activation.

	Rises fast (τ = 0.5 s) from ~1 nM baseline to 10 µM peak, then
	decays with τ = 30 s. `delay` shifts the stimulus onset so the
	curve is flat at baseline for t < delay.
	"""
	t_eff = t - delay
	if t_eff < 0:
		return ATP_EX_REST_UM
	rise = 1.0 - np.exp(-t_eff / ATP_EX_TAU_RISE)
	decay = np.exp(-max(0.0, t_eff - ATP_EX_T_PEAK) / ATP_EX_TAU_DECAY)
	return ATP_EX_REST_UM + (ATP_EX_PEAK_UM - ATP_EX_REST_UM) * rise * decay


# ── SOCE: Dolan 2014 MWC + STIM1 dimerisation (Hoover & Lewis 2011 frame) ─
# STIM1 cycle (mass-action) — keeps the dimer pool size as a state variable.
# `STIM1_dim` is counted in DIMER PARTICLES (matches Dolan Table S1
# "STIM1₂ (11)") — 1 dimerisation event consumes 2 free monomers and
# creates 1 dimer particle. Rate constants chosen so the Dolan 2014
# Table S1 resting IC (st_Ca=3805, st_free=438, st_dim=11) is at
# detailed balance.
K_STIM = {
	# STIM1·Ca²⁺_dts ↔ STIM1_free + Ca²⁺_dts (Ca²⁺ release from STIM EF-hand)
	'k_release_f':   0.1,      # forward (s⁻¹)
	# k_release_r derived from detailed balance at Dolan IC:
	#   k_release_r = k_release_f × st_Ca / (st_free × ca_dts)
	#               = 0.1 × 3805 / (438 × 250) = 3.475e-3 µM⁻¹·s⁻¹
	'k_release_r':   3.475e-3, # reverse (µM⁻¹·s⁻¹)
	# 2 STIM1_free ↔ STIM1_dim — diffusion-limited dimerisation.
	# k_dim_f from detailed balance at Dolan IC (dimer-particle count):
	#   k_dim_f = k_dim_r × st_dim / st_free² = 1.0 × 11 / 438² ≈ 5.73e-5
	'k_dim_f':      5.73e-5,   # forward (count⁻¹·s⁻¹)
	'k_dim_r':       1.0,      # reverse (s⁻¹)
}

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
K_MWC = {
	'L':   1.0e-4,      # opening equilibrium without STIM
	'Ka':  2.0,         # STIM2 association constant (rescaled from Hoover a.u.)
	'f':   14.2,        # opening cooperativity per bound STIM2
	'a':   0.5,         # binding cooperativity (negative)
}

# Dolan 2014 puncta entry (eq. 2): qp = α·[Ca]_cyt^n / (KM^n + [Ca]_cyt^n) + 0.01
#   qp gives the fraction of STIM2 dimers translocated into puncta where
#   they can engage Orai. α = 0.2 is the Dolan default. KM and n are the
#   two free parameters Dolan scans within homeostatic constraints.
PUNCTA = {
	'alpha':  0.2,      # max puncta fraction at saturating [Ca²⁺]_cyt
	'KM_uM':  0.5,      # half-activation [Ca²⁺]_cyt (chosen mid-range; Dolan-scanned)
	'n':      4.0,      # Hill coefficient (chosen mid-range; Dolan-scanned)
	'baseline': 0.01,   # constitutive puncta fraction at zero [Ca²⁺]_cyt
}

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

# Number of monomers per Orai1 tetramer (CRAC channel pore-forming subunit).
ORAI_SUBUNITS_PER_CHANNEL = 4

# Number of monomers per STIM1 dimer (sensor unit that binds Orai). Used
# for total-monomer mass-balance accounting (free + STIM1_Ca + 2·STIM1_dim
# is conserved). The MWC and listener treat STIM1_dim as dimer particles.
STIM_MONOMERS_PER_DIMER = 2


def ip3_forcing_uM(t, delay=0.0):
	"""Dolan 2014 Fig. S2 IP3 time curve, returning concentration in µM.

	Plateau-decay approximation: rises ~5.5× over 3 s, decays with τ=60 s.
	`delay` shifts the stimulus onset: the curve is flat at the resting
	baseline for t < delay and begins rising at t = delay.
	At t_eff <= 0 returns the resting baseline so the curve is well-defined
	for any sub-step the BDF solver evaluates.
	"""
	t_eff = t - delay
	if t_eff <= 0:
		return IP3_REST_UM
	rise = 1.0 - np.exp(-t_eff / IP3_TAU_RISE)
	decay = np.exp(-max(0.0, t_eff - IP3_T_PEAK) / IP3_TAU_DECAY)
	return IP3_REST_UM * (1.0 + (IP3_FOLD - 1.0) * rise * decay)



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


def _ode_rhs(t, y, t_sim_start, ip3_forced, ip3_delay=0.0):
	"""Right-hand side of the calcium ODE.

	`y` carries integer-equivalent counts (continuous floats during
	integration). `t` is the sub-step time *within* this 1-second outer step;
	`t_sim_start` is the wall-clock time at which the outer step began. IP3
	forcing is parameterised in absolute simulation time.
	"""
	dy = np.zeros(N_SPECIES)

	# Concentrations (µM) for Ca²⁺ and IP3.
	ca_cyt = max(y[_IDX['CA2_CYT[c]']], 0.0) * _UM_PER_COUNT_CYT
	ca_dts = max(y[_IDX['CA2_DTS[dts]']], 0.0) * _UM_PER_COUNT_DTS
	if ip3_forced:
		ip3 = ip3_forcing_uM(t_sim_start + t, delay=ip3_delay)
	else:
		ip3 = max(y[_IDX['IP3[c]']], 0.0) * _UM_PER_COUNT_CYT

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
	# Scaffold-only (N_GSN = 5 000); see K_GSN comment block for the
	# biology / scope disclosure.
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
	# The extracellular reservoir is treated as infinite (no debit).

	# ── P2X1 state transitions (always run, even when CA_EX = 0) ────
	# When CA_EX = 0 the channel can still cycle through C → O → D in
	# response to extracellular ATP — it just doesn't deliver any
	# Ca²⁺ (the flux block above is skipped). This matches the
	# Vial & Evans 2002 observation that P2X1 desensitises whether
	# Ca²⁺ is present or not.
	if ip3_forced:
		atp_ex_um = atp_ex_forcing_uM(t_sim_start + t, delay=ip3_delay)
	else:
		atp_ex_um = ATP_EX_REST_UM
	v_p2x1_act   = K_P2X1['k_act']   * p2x1_c * atp_ex_um
	v_p2x1_close = K_P2X1['k_close'] * p2x1_o
	v_p2x1_des   = K_P2X1['k_des']   * p2x1_o
	v_p2x1_rec   = K_P2X1['k_rec']   * p2x1_d
	dy[_IDX['P2X1[pl]']]   += -v_p2x1_act + v_p2x1_close + v_p2x1_rec
	dy[_IDX['P2X1_O[pl]']] += +v_p2x1_act - v_p2x1_close - v_p2x1_des
	dy[_IDX['P2X1_D[pl]']] +=                              +v_p2x1_des - v_p2x1_rec

	# IP3 is forced when `ip3_forced` is True; otherwise it free-floats
	# (decay/regeneration handled by upstream processes in v0.3+).
	if ip3_forced:
		# Drive IP3 count toward the curve without integrating the ODE on it.
		target_count = ip3_forcing_uM(t_sim_start + t, delay=ip3_delay) / _UM_PER_COUNT_CYT
		dy[_IDX['IP3[c]']] = (target_count - y[_IDX['IP3[c]']]) / 0.1

	return dy


class CalciumSignalling:
	"""Parameters and ODE driver for the CalciumDynamics process.

	Public surface matches the wcEcoli convention used by other process
	dataclasses (see TwoComponentSystem):

	  * `molecule_names` — array of bulk-molecule IDs the process views.
	  * `molecules_to_next_time_step(counts, dt, t_sim, ip3_forced)` —
	    integrate the ODE for one outer timestep and return integer count
	    deltas + estimated ATP cost.
	"""

	def __init__(self, sim_data):
		self.molecule_names = np.array(MOLECULE_NAMES, dtype='U30')
		self.n_species = N_SPECIES

	def molecules_to_next_time_step(self, counts, dt, t_sim, ip3_forced=False,
			ip3_delay=0.0):
		"""Run one outer timestep of the calcium ODE.

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
			args=(t_sim, ip3_forced, ip3_delay),
			atol=1e-3,    # counts; ~0.001 molecule precision
			rtol=1e-6,
			max_step=dt,
		)
		if not sol.success:
			raise RuntimeError(
				f'CalciumSignalling ODE failed at t_sim={t_sim}: {sol.message}')

		y_final = np.maximum(sol.y[:, -1], 0.0)
		delta = np.round(y_final - y0).astype(np.int64)

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
