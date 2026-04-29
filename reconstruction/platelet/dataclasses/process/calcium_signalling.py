"""
Calcium signalling dataclass for the v0.2 platelet whole-cell model.

Holds species ordering, rate constants, compartment volumes, and the
ODE right-hand-side used by the CalciumDynamics process. The ODE covers:

  * IP3R (6-state Markov, Sneyd & Dufour 2002 type 2 — mass-action form,
    Purvis 2008 Table 1)
  * SERCA cycle (E1/E2, Purvis 2008 Table 1 / Dode 2002)
  * PMCA (2-state Michaelis–Menten, Caride 2007 Table 3 basal)
  * SOCE (STIM1 free/Ca-bound/dimer, Dolan 2014 MWC reduction)
  * IP3 forcing (Dolan 2014 Fig. S2 shape)

State layout: integer counts per species, indexed by `MOLECULE_NAMES`.
ODE works in count units; rate laws convert to concentration internally
where needed (volumes for cytosol vs DTS).

Numerical regime: 1-second outer timestep, scipy.integrate.solve_ivp with
the BDF method (the system is stiff — SERCA cycle rates run up to 1000 s⁻¹).
"""

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
	'IP3R_n[dts]',
	'IP3R_o[dts]',
	'IP3R_a[dts]',
	'IP3R_i1[dts]',
	'IP3R_i2[dts]',
	'IP3R_s[dts]',
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


# ── Rate constants ────────────────────────────────────────────────────────
# IP3R Sneyd & Dufour 2002 type-2 (Purvis 2008 Table 1).
# Mass-action approximation: drops the φ-function L-modulation; this is
# adequate for the v0.2 resting-stability and small-transient targets.
K_IP3R = {
	'k1':  0.64,    # n + Ca²⁺ → i1   (µM⁻¹·s⁻¹)
	'km1': 0.04,    # i1 → n           (s⁻¹)
	'k2':  37.4,    # n + IP3 → o     (µM⁻¹·s⁻¹)
	'km2': 1.4,     # o → n            (s⁻¹)
	'k3':  11.0,    # o + Ca²⁺ → s    (µM⁻¹·s⁻¹)  (shut)
	'km3': 29.8,    # s → o            (s⁻¹)
	'k4':  4.0,     # o + Ca²⁺ → a    (µM⁻¹·s⁻¹)
	'km4': 0.54,    # a → o            (s⁻¹)
	'l2':  1.7,     # i ↔ s            (s⁻¹)
	'lm2': 0.8,     # s → i            (s⁻¹)
}

# IP3R Ca²⁺ flux: empirical conductance × open fraction × concentration
# gradient. Replaces the Nernst formulation in v0.2; the conductance value
# is calibrated so that resting J_IP3R ≈ 69 nM/s (balanceable by SERCA).
# Original value (0.30) gave 5143 nM/s — 74× too large, draining DTS in 1.6 s.
K_IP3R_FLUX = 0.004     # µM/s per (µM gradient × Po)

# SERCA cycle (Purvis 2008 Table 1, Dode 2002). Mass-action.
K_SERCA = {
	'k_shuttle_f':  600.0,    # E2 → E1                        (s⁻¹)
	'k_shuttle_r':  600.0,    # E1 → E2                        (s⁻¹)
	# k_bind_f calibrated so SERCA throughput ≈ 124 ct/s at rest (100 nM cytosol).
	# Original value (1e3 µM⁻²·s⁻¹) gave 59,140 ct/s — draining cytosol in ms.
	# Km kept at 0.1 µM: k_bind_r = k_bind_f × 0.01.
	'k_bind_f':     2.1101,   # E1 + 2 Ca²⁺_cyt → E1·Ca       (µM⁻²·s⁻¹)
	'k_bind_r':     0.021101, # E1·Ca → E1 + 2 Ca²⁺_cyt        (s⁻¹)
	'k_phos_f':     700.0,    # E1·Ca → E1P·Ca                 (s⁻¹)
	'k_phos_r':     5.0,
	'k_conf_f':     600.0,    # E1P·Ca ⇌ E2P·Ca                (s⁻¹)
	'k_conf_r':     50.0,
	'k_release_f':  1000.0,   # E2P·Ca → E2P + 2 Ca²⁺_dts      (s⁻¹)
	'k_release_r':  4.0e-3,   # reverse (µM⁻²·s⁻¹) [scaled]
	'k_dephos_f':   500.0,    # E2P → E2                       (s⁻¹)
	'k_dephos_r':   1.0,
}

# PMCA (Caride 2007 Table 3, basal).
K_PMCA = {
	'k_on':   10.0,    # PMCA + Ca²⁺ → PMCA·Ca   (µM⁻¹·s⁻¹)
	'k_off':  50.0,    # PMCA·Ca → PMCA          (s⁻¹)
	'k_cat':  5.5,     # PMCA·Ca → PMCA + Ca²⁺_ex (s⁻¹)
}

# SOCE (Dolan 2014, simplified MWC reduction).
# STIM1_Ca ⇌ STIM1_free is set by [Ca²⁺]_dts; STIM1_free ⇌ STIM1_dim is
# the activation step; STIM1_dim gates Orai1 to allow Ca²⁺_ex → Ca²⁺_cyt.
K_SOCE = {
	'k_release_f':  0.1,     # STIM1_Ca → STIM1_free (s⁻¹) — slow at high Ca²⁺_dts
	# k_release_r calibrated so v_STIM1 = 0 at Dolan 2014 Table S1 IC
	# (st_ca=3805, st_free=438, ca_dts=250 µM):
	#   k_release_r = k_release_f × st_ca / (st_free × ca_dts) = 3.475e-3
	# Original value (1e-3) was 3.5× too small, causing continuous Ca release
	# from DTS into the STIM1-Ca pool at rest.
	'k_release_r':  3.475e-3, # STIM1_free + Ca²⁺_dts → STIM1_Ca (µM⁻¹·s⁻¹)
	# k_dim_f: set so that the Dolan 2014 Table S1 resting initial conditions
	# (st_free=438, st_dim=22) are at equilibrium: k_dim_f = k_dim_r × st_dim / st_free².
	# The original value (0.05) was ~436× too large — an ad hoc estimate that caused
	# runaway STIM1 dimerisation in the first timestep.  See GitHub issue #46 for the
	# longer-term fix (implement the full Dolan 2014 MWC allosteric model).
	'k_dim_f':      1.15e-4, # 2 STIM1_free → STIM1_dim (count⁻¹·s⁻¹) [equilibrium-derived]
	'k_dim_r':      1.0,     # STIM1_dim → 2 STIM1_free  (s⁻¹)
	# k_orai calibrated so SOCE = PMCA efflux at rest (21 nM/s each) with
	# st_dim=22, Ca_ex−Ca_cyt ≈ 1200 µM:
	#   k_orai = PMCA_efflux_µMs / (st_dim × ΔCa) = 21.09e-3 / (22 × 1200) = 7.99e-7
	# Original value (0.001) gave SOCE = 26,400 nM/s — ~1260× too large.
	'k_orai':       7.99e-7, # SOCE flux per STIM1_dim per (Ca²⁺_ex − Ca²⁺_cyt) µM/s
}


def ip3_forcing_uM(t):
	"""Dolan 2014 Fig. S2 IP3 time curve, returning concentration in µM.

	Plateau-decay approximation: rises ~5.5× over 3 s, decays with τ=60 s.
	At t<=0 returns the resting baseline so the curve is well-defined for
	any sub-step the BDF solver evaluates.
	"""
	if t <= 0:
		return IP3_REST_UM
	rise = 1.0 - np.exp(-t / IP3_TAU_RISE)
	decay = np.exp(-max(0.0, t - IP3_T_PEAK) / IP3_TAU_DECAY)
	return IP3_REST_UM * (1.0 + (IP3_FOLD - 1.0) * rise * decay)


def _ode_rhs(t, y, t_sim_start, ip3_forced):
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
		ip3 = ip3_forcing_uM(t_sim_start + t)
	else:
		ip3 = max(y[_IDX['IP3[c]']], 0.0) * _UM_PER_COUNT_CYT

	# Counts (kept as-is; mass-action between protein states is volumeless).
	n  = y[_IDX['IP3R_n[dts]']]
	o  = y[_IDX['IP3R_o[dts]']]
	a  = y[_IDX['IP3R_a[dts]']]
	i1 = y[_IDX['IP3R_i1[dts]']]
	i2 = y[_IDX['IP3R_i2[dts]']]
	s  = y[_IDX['IP3R_s[dts]']]

	se1   = y[_IDX['SERCA_E1[dts]']]
	se2   = y[_IDX['SERCA_E2[dts]']]
	se1c  = y[_IDX['SERCA_E1Ca[dts]']]
	se1pc = y[_IDX['SERCA_E1PCa[dts]']]
	se2pc = y[_IDX['SERCA_E2PCa[dts]']]
	se2p  = y[_IDX['SERCA_E2P[dts]']]

	pmca   = y[_IDX['PMCA[pl]']]
	pmcaca = y[_IDX['PMCA_Ca[pl]']]

	st_free = y[_IDX['STIM1_free[dts]']]
	st_ca   = y[_IDX['STIM1_Ca[dts]']]
	st_dim  = y[_IDX['STIM1_dim[dts]']]

	# ── IP3R 6-state Markov (mass-action) ─────────────────────────────
	# n + Ca → i1, n + IP3 → o
	v_n_i1 = K_IP3R['k1']  * ca_cyt * n - K_IP3R['km1'] * i1
	v_n_o  = K_IP3R['k2']  * ip3    * n - K_IP3R['km2'] * o
	# o + Ca → a, o + Ca → s
	v_o_a  = K_IP3R['k4']  * ca_cyt * o - K_IP3R['km4'] * a
	v_o_s  = K_IP3R['k3']  * ca_cyt * o - K_IP3R['km3'] * s
	# a + Ca → i2 (same kinetics as n + Ca → i1)
	v_a_i2 = K_IP3R['k1']  * ca_cyt * a - K_IP3R['km1'] * i2
	# inhibited ↔ shut connector (via l2/lm2); routes inhibited states out
	v_i1_s = K_IP3R['l2'] * i1 - K_IP3R['lm2'] * s
	v_i2_s = K_IP3R['l2'] * i2 - K_IP3R['lm2'] * s

	dy[_IDX['IP3R_n[dts]']]  += -v_n_i1 - v_n_o
	dy[_IDX['IP3R_o[dts]']]  += +v_n_o - v_o_a - v_o_s
	dy[_IDX['IP3R_a[dts]']]  += +v_o_a - v_a_i2
	dy[_IDX['IP3R_i1[dts]']] += +v_n_i1 - v_i1_s
	dy[_IDX['IP3R_i2[dts]']] += +v_a_i2 - v_i2_s
	dy[_IDX['IP3R_s[dts]']]  += +v_o_s + v_i1_s + v_i2_s

	# IP3R Ca²⁺ flux: open and active states conduct (a fully, o at 10%).
	ip3r_total = max(n + o + a + i1 + i2 + s, 1.0)
	po = (a + 0.1 * o) / ip3r_total
	# µM/s in the cytosolic compartment, driven by the Ca²⁺ gradient
	flux_ip3r_uMs = K_IP3R_FLUX * po * (ca_dts - ca_cyt)
	# Convert µM/s → counts/s in each compartment (mass-conserving).
	# The same number of atoms leave DTS as enter CYT; dividing by the
	# DTS volume factor would create phantom Ca²⁺ (bug fix, issue #46).
	flux_ip3r_count_cyt = flux_ip3r_uMs / _UM_PER_COUNT_CYT
	flux_ip3r_count_dts = flux_ip3r_count_cyt

	dy[_IDX['CA2_CYT[c]']]   += +flux_ip3r_count_cyt
	dy[_IDX['CA2_DTS[dts]']] += -flux_ip3r_count_dts

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

	# Each SERCA turnover removes 2 Ca²⁺_cyt and adds 2 Ca²⁺_dts.
	dy[_IDX['CA2_CYT[c]']]   += -2.0 * v_bind
	dy[_IDX['CA2_DTS[dts]']] += +2.0 * v_release

	# ── PMCA (2-state MM; basal kinetics, no CaM coupling) ───────────
	v_pmca_bind = K_PMCA['k_on'] * pmca * ca_cyt - K_PMCA['k_off'] * pmcaca
	v_pmca_cat  = K_PMCA['k_cat'] * pmcaca

	dy[_IDX['PMCA[pl]']]    += -v_pmca_bind + v_pmca_cat
	dy[_IDX['PMCA_Ca[pl]']] += +v_pmca_bind - v_pmca_cat
	dy[_IDX['CA2_CYT[c]']]  += -v_pmca_bind  # one Ca²⁺ leaves cytosol per binding
	# Catalysis ejects Ca²⁺ to extracellular; the reservoir is not modelled.

	# ── SOCE (STIM1 dimerisation + Orai1 flux) ────────────────────────
	# Rename from v_release to avoid collision with SERCA's v_release above.
	v_stim1_release = (
		K_SOCE['k_release_f'] * st_ca
		- K_SOCE['k_release_r'] * st_free * ca_dts
	)
	# Dimerisation as a 2nd-order step in free monomers (count units).
	v_dim = K_SOCE['k_dim_f'] * st_free * st_free - K_SOCE['k_dim_r'] * st_dim

	dy[_IDX['STIM1_Ca[dts]']]   += -v_stim1_release
	dy[_IDX['STIM1_free[dts]']] += +v_stim1_release - 2.0 * v_dim
	dy[_IDX['STIM1_dim[dts]']]  += +v_dim
	# Ca²⁺ released from STIM1 binding returns to the free DTS pool (bug fix, issue #46).
	dy[_IDX['CA2_DTS[dts]']]    += v_stim1_release

	# Orai1 flux into cytosol, gated by STIM1_dim count.
	soce_flux_uMs = K_SOCE['k_orai'] * st_dim * (CA_EX_UM - ca_cyt)
	dy[_IDX['CA2_CYT[c]']] += soce_flux_uMs / _UM_PER_COUNT_CYT

	# IP3 is forced when `ip3_forced` is True; otherwise it free-floats
	# (decay/regeneration handled by upstream processes in v0.3+).
	if ip3_forced:
		# Drive IP3 count toward the curve without integrating the ODE on it.
		target_count = ip3_forcing_uM(t_sim_start + t) / _UM_PER_COUNT_CYT
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

	def molecules_to_next_time_step(self, counts, dt, t_sim, ip3_forced=False):
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
			args=(t_sim, ip3_forced),
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
		pmca_atp = max(int(K_PMCA['k_cat'] * counts[_IDX['PMCA_Ca[pl]']] * dt), 0)
		atp_cost = serca_atp + pmca_atp

		return delta, atp_cost
