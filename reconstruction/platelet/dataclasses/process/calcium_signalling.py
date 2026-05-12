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
	# v0.4 GPCR receptor cascade (issue #9) — replaces gq_signal_uM
	# forcing with explicit P2Y1 + PAR1 + PAR4 → Gαq cascade.
	'P2Y1_inactive[pl]', # P2Y1 ADP receptor (Gq-coupled), inactive
	'P2Y1_active[pl]',   # P2Y1 active (ADP-bound)
	'PAR1_inactive[pl]', # PAR1 thrombin receptor (high-affinity), inactive
	'PAR1_active[pl]',   # PAR1 cleaved/active (proteolytic, ~irreversible)
	'PAR4_inactive[pl]', # PAR4 thrombin receptor (low-affinity)
	'PAR4_active[pl]',
	'Gq_active[c]',      # Gαq-GTP (active); inactive Gq implicit (total - active)
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
GAMMA_IP3R_S = 0.075e-12         # 0.075 pS = calibrated Ca²⁺ conductance, A/V


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
	'k_bind_f':     210.0,    # E1 + 2 Ca²⁺_cyt → E1·Ca²⁺      (µM⁻²·s⁻¹)
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
N_GSN = 1_400_000


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
K_HSP90B1_M = {
	'k_on':    0.5,    # HSP90B1_M + Ca²⁺ → HSP90B1_M·Ca  (µM⁻¹·s⁻¹) — slow
	'k_off':   1.0,    # reverse                          (s⁻¹) — Kd = 2 µM
	                   # τ_release ≈ 1 s — matches transient timescale, so
	                   # these sites hold their Ca²⁺ during the ~1 s peak
	                   # and act as a "floor" keeping free DTS [Ca²⁺] > 0.
}
K_HSP90B1_L = {
	'k_on':    1.0,    # HSP90B1_L + Ca²⁺ → HSP90B1_L·Ca  (µM⁻¹·s⁻¹)
	'k_off': 600.0,    # reverse                          (s⁻¹) — Kd = 600 µM
}
N_HSP90B1 = 10_000              # molecules
N_HSP90B1_M = N_HSP90B1 * 4     # 40 000 medium-affinity sites
N_HSP90B1_L = N_HSP90B1 * 11    # 110 000 low-affinity sites


# ── BiP / HSPA5 / GRP78 — Phase 3 / issue #25 ────────────────────────────
# Most abundant ER chaperone. Primary role is protein folding, but
# Lièvremont 1997 demonstrated BiP contributes ~25 % of the ER Ca²⁺
# store, with a stoichiometry of 1–2 Ca²⁺ per molecule. Modelled as a
# single low-affinity pool with Kd ~ 500 µM (matches the Lièvremont
# "mM-range free [Ca²⁺]_ER" they describe).
#
# Platelet copy number: ~50 000 (order-of-magnitude estimate; BiP is the
# canonically most-abundant ER chaperone, more so than CALR or HSP90B1).
K_BIP = {
	'k_on':    2.0,    # BiP + Ca²⁺ → BiP·Ca   (µM⁻¹·s⁻¹)
	'k_off': 1000.0,   # reverse               (s⁻¹) — Kd = 500 µM
}
N_BIP = 50_000 * 1                 # 1 effective site per BiP molecule (mid of 1–2 range with the lower count used to stay biologically conservative — total sites 50 000)
# (For 1.5-site stoichiometry: N_BIP_TOTAL = 50 000 × 1.5 = 75 000; using
# 50 000 as the conservative estimate. Resulting bound at rest = 16 700
# instead of 25 000 — closer to the lower bound of Lièvremont's 25 % of
# store. Capacity can be scaled up in v0.4 with explicit two-site model.)


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
K_CREC = {
	'k_on':    0.5,    # CREC + Ca²⁺ → CREC·Ca  (µM⁻¹·s⁻¹)
	'k_off': 500.0,    # reverse                (s⁻¹) — Kd = 1 mM
}
N_CREC = 15_000 * 4                # 60 000 sites


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
	# k12 = CaM dissociation from PMCA. Caride 2007 measured 0.033 s⁻¹
	# (τ = 30 s) in purified in vitro preparation. In vivo, PIP2 binding
	# to the PMCA C-terminus competitively displaces CaM on a much faster
	# timescale (Penniston & Enyedi 1998 review; Mandal 2024). The in
	# vitro value traps PMCA in PMCA·CaM during sustained Ca²⁺ elevation
	# and prevents recovery; using an effective in vivo rate of 1 s⁻¹
	# (τ = 1 s, 30× faster) restores physiological PMCA Vmax during
	# transient recovery. v0.3.1 fix — see lab book 2026-05-12 (DTS
	# overshoot diagnosis).
	'k12':    1.0,     # PMCA·CaM → PMCA + CaM_free (s⁻¹) — in-vivo rate
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
GAMMA_P2X1_S = 0.0013e-12   # 1.3 fS Ca²⁺-specific conductance, A/V — calibrated


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
K_NCX = {
	'V_max':  5_000.0,   # ions/s — total per platelet; calibration anchor
	'K_m':    5.0,       # substrate Hill half-saturation (µM)
	'K_a':    0.2,       # allosteric activation half-point (µM) — slightly lower for more recovery-phase contribution
	'h':      4,         # allosteric Hill cooperativity (switch-like)
}


# ── v0.4 GPCR cascade — P2Y1 + PAR1/4 → Gαq → PLCβ (issue #9) ────────────
# Replaces the v0.3.x gq_signal_uM forcing with an explicit receptor →
# G-protein cascade. The model now takes physiological agonist
# concentrations (thrombin in nM, ADP in µM) as input rather than an
# abstract Gq curve.
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

K_P2Y1 = {
	'k_on':   1.0,    # P2Y1_i + ADP → P2Y1_a       (µM⁻¹·s⁻¹) — reversible
	'k_off':  0.5,    # P2Y1_a → P2Y1_i + ADP        (s⁻¹)
}

K_PAR1 = {
	'k_cleave':       2.0,    # PAR1_i + thrombin → PAR1_a   (nM⁻¹·s⁻¹) — fast cleavage
	'k_internalize':  0.02,   # PAR1_a → PAR1_i              (s⁻¹) — τ ~50s (lumped internalisation+recycling)
}

K_PAR4 = {
	'k_cleave':       0.2,    # PAR4_i + thrombin → PAR4_a   (nM⁻¹·s⁻¹) — slower than PAR1
	'k_internalize':  0.005,  # PAR4_a → PAR4_i              (s⁻¹) — τ ~200s (sustained response)
}

K_GQ = {
	'k_act_per_R':  0.001,    # active receptor catalyses Gq exchange  (s⁻¹·count⁻¹)
	'k_rgs':        0.033,    # Gq_active → Gq_inactive (RGS-accelerated GTPase)  (s⁻¹) — τ ~30s
	# k_basal calibrated so resting Gq_active ≈ 100 → gq_um = 0.1 µM
	# (matches the v0.3.x GQ_REST_UM that the PLCβ ODE expects).
	# Steady-state fraction = k_basal / (k_basal + k_rgs); with k_rgs =
	# 0.033 and target fraction 100/5000 = 0.02, k_basal = 0.033 × 0.02/0.98
	# = 6.7e-4 s⁻¹.
	'k_basal':      6.7e-4,
}

N_GQ_TOTAL = 5_000  # Mazet 2020 platelet Gαq molecules


# ── Agonist forcing functions — v0.4 stimulation inputs ──────────────────
# These replace the abstract gq_signal_uM. The model now takes
# physiological agonist concentrations (thrombin in nM, ADP in µM).

# Thrombin (PAR1/4 agonist).  Standard stimulation = 1 nM peak
# (matches the Dolan 2014 protocol). Rises within ms (clotting cascade
# is fast), plateaus, then is cleared.
THROMBIN_REST_NM     = 0.0      # no thrombin at rest
THROMBIN_PEAK_NM     = 1.0      # default 1 nM peak
THROMBIN_TAU_RISE_S  = 0.5      # fast rise
THROMBIN_T_PEAK_S    = 5.0      # sustained plateau
THROMBIN_TAU_DECAY_S = 120.0    # slow decay (thrombin lingers)

# ADP (P2Y1 agonist).  Released from dense granules during platelet
# activation; cleared by ectoNTPDases. Standard stimulation = 10 µM peak.
ADP_REST_UM     = 0.0           # no ADP at rest
ADP_PEAK_UM     = 10.0          # default 10 µM peak
ADP_TAU_RISE_S  = 1.0           # dense granule secretion timescale
ADP_T_PEAK_S    = 5.0
ADP_TAU_DECAY_S = 30.0          # ectoNTPDase clearance


def thrombin_nM(t, delay=0.0, peak_nM=THROMBIN_PEAK_NM):
	"""Thrombin concentration timecourse (nM)."""
	t_eff = t - delay
	if t_eff <= 0:
		return THROMBIN_REST_NM
	rise  = 1.0 - np.exp(-t_eff / THROMBIN_TAU_RISE_S)
	decay = np.exp(-max(0.0, t_eff - THROMBIN_T_PEAK_S) / THROMBIN_TAU_DECAY_S)
	return THROMBIN_REST_NM + (peak_nM - THROMBIN_REST_NM) * rise * decay


def adp_uM(t, delay=0.0, peak_uM=ADP_PEAK_UM):
	"""ADP concentration timecourse (µM)."""
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


def ip3_forcing_uM(t, delay=0.0):
	"""Dolan 2014 Fig. S2 IP3 time curve, returning concentration in µM.

	v0.2.x legacy: this was the *forced* IP3 timecourse driving the
	model. v0.3 (#31) replaces this with the PI cycle: IP3 becomes a
	model output, produced from PIP2 hydrolysis by Gq-activated PLCβ.
	The function is retained as a *calibration reference* — the PI
	cycle parameters are tuned so that the model-produced IP3
	timecourse approximates this curve under standard Gq forcing.
	"""
	t_eff = t - delay
	if t_eff <= 0:
		return IP3_REST_UM
	rise = 1.0 - np.exp(-t_eff / IP3_TAU_RISE)
	decay = np.exp(-max(0.0, t_eff - IP3_T_PEAK) / IP3_TAU_DECAY)
	return IP3_REST_UM * (1.0 + (IP3_FOLD - 1.0) * rise * decay)


# ── PI cycle / PLCβ — Phase 4 / issue #31 ────────────────────────────────
# Replaces the forced IP3 curve. IP3 is now produced from PIP2 hydrolysis
# by Gq-activated PLCβ, following the framework of Mazet, Tindall,
# Gibbins & Fry 2020 *Sci. Rep.* 10:13889.
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
K_PLCB = {
	'k_act':    0.5,      # PLCb_i + Gq → PLCb_a    (µM⁻¹·s⁻¹) — calibrated
	'k_inact':  0.3,      # PLCb_a    → PLCb_i      (s⁻¹)       — τ ~ 3 s
	'k_cat':    2.26e-7,  # PLCb_a + PIP2 → PLCb_a + IP3 + DAG  (count⁻¹·s⁻¹)
}                         # — calibrated against resting + peak IP3 targets

K_PI_CYCLE = {
	# PIP2 resynthesis — lumped PI → PI4P → PIP2 chain. Set equal to
	# basal hydrolysis rate so PIP2 sits at its resting value.
	'k_resynth':   3.62,   # PIP2 / s — calibrated to basal balance
	# IP3 degradation (5-phosphatase to IP2 + 3-kinase to IP4, lumped)
	# τ ~ 50 s matches Dolan Fig. S2 decay tail.
	'k_ip3_deg':   0.02,   # IP3 → IP2/IP4   (s⁻¹)
	# DAG kinase (DAG → PA)
	'k_dag_deg':   0.05,   # DAG → PA       (s⁻¹) — τ ~ 20 s
}


# ── Gq activity forcing — replaces ip3_forcing_uM ────────────────────────
# In v0.3 the stimulation is delivered as a Gαq activity signal that
# drives PLCβ activation (replacing the v0.2.x forced IP3 curve). When
# v0.4 receptor signalling lands (#9), this forcing is replaced by an
# explicit GPCR → Gαq cascade with receptor / agonist-specific dynamics.
GQ_REST_UM   = 0.1         # tonic basal activity — maintains resting IP3 ~ 50 nM
GQ_PEAK_UM   = 2.0         # peak Gq concentration during activation (~20× rest)
GQ_TAU_RISE  = 0.5         # s — fast onset (GPCR → Gαq exchange)
GQ_T_PEAK    = 1.0         # s
GQ_TAU_DECAY = 30.0        # s — RGS-mediated GTPase inactivation


def gq_signal_uM(t, delay=0.0):
	"""Coarse-grained Gαq active concentration; drives PLCβ activation.

	Replaces v0.2.x `ip3_forcing_uM` as the model's primary stimulation
	input. v0.4 (#9) will replace this with explicit receptor cascades.

	Returns `GQ_REST_UM` for `t_eff <= 0` (tonic basal activity that
	maintains the resting IP3 baseline at ~50 nM).
	"""
	t_eff = t - delay
	if t_eff <= 0:
		return GQ_REST_UM
	rise = 1.0 - np.exp(-t_eff / GQ_TAU_RISE)
	decay = np.exp(-max(0.0, t_eff - GQ_T_PEAK) / GQ_TAU_DECAY)
	return GQ_REST_UM + (GQ_PEAK_UM - GQ_REST_UM) * rise * decay


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
K_MITO = {
	# MCU uptake — cooperative Hill kinetics. n=4 gives a sharp switch
	# at K_MCU = 1 µM; at resting cyt (100 nM) MCU is effectively off
	# (only ~5 ions/s), preserving the cyt peak signal. During the
	# transient (cyt > 1 µM) MCU activates strongly.
	'V_max_MCU':  50_000.0,   # ions/s — total over all mitochondria; calibrated
	'K_MCU':      1.0,        # Hill half-saturation (µM); literature 0.5–10
	'n_MCU':      4,          # Hill cooperativity; literature 2–4 (we use upper end for switch-like)
	# NCLX efflux — linear in [Ca²⁺]_mito; slow release.
	'k_NCLX':     0.005,      # s⁻¹  (τ = 200 s ~ 3 min slow release)
}



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
	# IP3 is produced by the PI cycle (PLCβ-driven; Phase 4 / #31). No
	# longer forced — the `ip3_forced` flag now controls Gq forcing, not
	# direct IP3 substitution.
	ip3 = max(y[_IDX['IP3[c]']], 0.0) * _UM_PER_COUNT_CYT

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

		# NCX (Na⁺/Ca²⁺ exchanger) — Ca²⁺ extrusion gated on extracellular
		# Ca²⁺ availability (needs Na⁺ gradient + somewhere for Ca²⁺ to go).
		# See K_NCX block above for biology and the Hill formulation.
		g_act = (ca_cyt / K_NCX['K_a']) ** K_NCX['h']
		g_act /= (1.0 + g_act)
		v_ncx = K_NCX['V_max'] * g_act * ca_cyt / (K_NCX['K_m'] + ca_cyt)
		dy[_IDX['CA2_CYT[c]']] -= v_ncx
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

	# ── PI cycle: PLCβ-driven IP3 production from PIP2 — Phase 4 / #31 ──
	# v0.4 (#9): Gq concentration now comes from the explicit GPCR cascade
	# state variable rather than the forced gq_signal_uM curve. The
	# downstream PI cycle is unchanged; only the input source differs.
	# Conversion: 1 µM gq_um ≈ 1 000 Gq_active count (calibrated so that
	# resting Gq_active ~100 maps to gq_um ~0.1 µM, matching GQ_REST_UM).
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

	# ── v0.4 GPCR cascade — receptors + Gαq cycle (issue #9) ────────
	# Replaces the v0.3.x gq_signal_uM forcing. Agonist concentrations
	# from external forcing → receptor activation → Gαq exchange →
	# PLCβ activation (already in the PI cycle block above).
	if ip3_forced:
		adp_um_now    = adp_uM(t_sim_start + t, delay=ip3_delay)
		thr_nm_now    = thrombin_nM(t_sim_start + t, delay=ip3_delay)
	else:
		adp_um_now    = ADP_REST_UM
		thr_nm_now    = THROMBIN_REST_NM

	# P2Y1: reversible ADP binding
	v_p2y1_on  = K_P2Y1['k_on']  * p2y1_i * adp_um_now
	v_p2y1_off = K_P2Y1['k_off'] * p2y1_a
	dy[_IDX['P2Y1_inactive[pl]']] += -v_p2y1_on + v_p2y1_off
	dy[_IDX['P2Y1_active[pl]']]   += +v_p2y1_on - v_p2y1_off

	# PAR1: thrombin cleavage (essentially irreversible) + slow recycling
	v_par1_cleave = K_PAR1['k_cleave'] * par1_i * thr_nm_now
	v_par1_int    = K_PAR1['k_internalize'] * par1_a
	dy[_IDX['PAR1_inactive[pl]']] += -v_par1_cleave + v_par1_int
	dy[_IDX['PAR1_active[pl]']]   += +v_par1_cleave - v_par1_int

	# PAR4: low-affinity thrombin receptor, sustained response
	v_par4_cleave = K_PAR4['k_cleave'] * par4_i * thr_nm_now
	v_par4_int    = K_PAR4['k_internalize'] * par4_a
	dy[_IDX['PAR4_inactive[pl]']] += -v_par4_cleave + v_par4_int
	dy[_IDX['PAR4_active[pl]']]   += +v_par4_cleave - v_par4_int

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
	  * `molecules_to_next_time_step(counts, dt, t_sim, ip3_forced)` —
	    integrate the ODE for one outer timestep and return integer count
	    deltas + estimated ATP cost.
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
