"""Unit tests for the calcium_signalling ODE solver.

Focused on parameter-binding regressions that are easy to miss from
end-to-end sim tests — in particular, whether agonist-peak overrides
take effect at call time rather than being baked into function
defaults at import time.
"""

import unittest

import numpy as np

from reconstruction.platelet.dataclasses.process import calcium_signalling as cs
from reconstruction.platelet.run_config import RunConfig


class TestAgonistForcingPeaks(unittest.TestCase):
	"""Peak overrides must take effect at call time (issue #44)."""

	def test_default_thrombin_peak_matches_module_constant(self):
		"""At t = T_peak, the forcing curve should approach the configured peak."""
		t = cs.THROMBIN_T_PEAK_S  # plateau time
		val = cs.thrombin_nM(t)
		self.assertGreater(val, 0.5 * cs.THROMBIN_PEAK_NM,
			'Default thrombin forcing should be near peak at T_peak.')

	def test_thrombin_peak_kwarg_override_is_live(self):
		"""Passing peak_nM=X must change the trajectory immediately."""
		t = cs.THROMBIN_T_PEAK_S
		baseline = cs.thrombin_nM(t, peak_nM=1.0)
		overridden = cs.thrombin_nM(t, peak_nM=5.0)
		# Linear in peak (REST is 0) → 5× peak ≈ 5× value.
		self.assertAlmostEqual(overridden / baseline, 5.0, places=6)

	def test_thrombin_peak_zero_yields_rest(self):
		"""peak_nM=0 must give zero forcing at all times."""
		for t in (1.0, 5.0, 60.0):
			self.assertEqual(cs.thrombin_nM(t, peak_nM=0.0), 0.0)

	def test_adp_peak_kwarg_override_is_live(self):
		t = cs.ADP_T_PEAK_S
		baseline = cs.adp_uM(t, peak_uM=1.0)
		overridden = cs.adp_uM(t, peak_uM=10.0)
		self.assertAlmostEqual(overridden / baseline, 10.0, places=6)

	def test_atp_ex_peak_kwarg_override_is_live(self):
		t = cs.ATP_EX_T_PEAK
		baseline = cs.atp_ex_forcing_uM(t, peak_uM=1.0)
		overridden = cs.atp_ex_forcing_uM(t, peak_uM=20.0)
		self.assertAlmostEqual(overridden / baseline, 20.0, places=6)

	def test_module_constant_reassignment_takes_effect(self):
		"""Reassigning cs.ADP_PEAK_UM at runtime must change the call result.

		This is the dose-sweep gating bug from issue #44: defaulting via
		`peak_uM=ADP_PEAK_UM` in the signature baked the value at import.
		"""
		original = cs.ADP_PEAK_UM
		try:
			cs.ADP_PEAK_UM = 0.1
			val_low = cs.adp_uM(cs.ADP_T_PEAK_S)
			cs.ADP_PEAK_UM = 50.0
			val_high = cs.adp_uM(cs.ADP_T_PEAK_S)
		finally:
			cs.ADP_PEAK_UM = original
		self.assertAlmostEqual(val_high / val_low, 50.0 / 0.1, places=6)


class TestOdeRhsResting(unittest.TestCase):
	"""All-zero peaks should give a resting RHS (no agonist excitation)."""

	def test_rest_ode_keeps_receptors_inactive(self):
		"""With all peaks = 0, the GPCR receptors at rest receive no signal."""
		from reconstruction.platelet.simulation_data import SimulationDataPlatelet
		sim_data = SimulationDataPlatelet()
		all_ids = list(sim_data.internal_state.bulk_molecules.bulk_data['id'])
		all_counts = sim_data.internal_state.bulk_molecules.initial_counts
		y0 = np.zeros(cs.N_SPECIES)
		for i, name in enumerate(cs.MOLECULE_NAMES):
			y0[i] = float(all_counts[all_ids.index(name)])

		# t_sim_start = 100 (past the agonist t_peak); with peaks = 0 the
		# forcing functions still return REST (0), so PAR1/4 cleavage rates
		# from extracellular thrombin must be exactly zero. Empty step_inputs →
		# no autocrine contribution.
		config = RunConfig(
			thrombin_peak_nM=0.0, adp_peak_uM=0.0, atp_ex_peak_uM=0.0)
		dy_rest = cs._ode_rhs(0.0, y0, 100.0, config, {})
		par1_idx = list(cs.MOLECULE_NAMES).index('PAR1_active[pl]')
		par4_idx = list(cs.MOLECULE_NAMES).index('PAR4_active[pl]')

		# At rest with thrombin = 0, the only sink for PAR*_active is
		# internalization (negative). New cleavage must be 0, so total dy/dt
		# for active receptors is ≤ 0.
		self.assertLessEqual(dy_rest[par1_idx], 0.0,
			'PAR1_active should not accumulate at rest.')
		self.assertLessEqual(dy_rest[par4_idx], 0.0,
			'PAR4_active should not accumulate at rest.')


class TestMcuCoupling(unittest.TestCase):
	"""#76 Part 2 — the MCU → IP3R-relief gate: the helper's clamps/behaviour,
	and the ODE's *application* of it (not just the algebra)."""

	def test_ip3r_relief_factor_clamps_and_behaviour(self):
		f = cs.ip3r_relief_factor
		# WT (scale=1): full relief → factor 1, regardless of Ca²⁺ or gain.
		self.assertEqual(f(0.5, 1.0, 1.0), 1.0)
		# Over-expression (scale>1): clamped to 1 — cannot *raise* release.
		# (Regression for the runPerturbation mcu ×4 sweep, which sets scale=4.)
		self.assertEqual(f(0.5, 4.0, 1.0), 1.0)
		# Resting (low Ca²⁺), KO: activation gate ≈ off → factor ≈ 1 (rest spared).
		self.assertGreater(f(0.1, 0.0, 1.0), 0.97)
		# Evoked (high Ca²⁺), KO: relief lost → factor < 1.
		self.assertLess(f(0.5, 0.0, 1.0), 1.0)
		# No misconfiguration can drive it negative (would reverse the flux).
		self.assertGreaterEqual(f(0.5, 0.0, 100.0), 0.0)
		# gain = 0 disables the coupling entirely → factor 1 even at full KO.
		self.assertEqual(f(0.5, 0.0, 0.0), 1.0)

	def test_coupling_cuts_ip3r_release_in_ode_and_is_inert_at_wt(self):
		"""At a fixed evoked state, the gain toggle isolates the coupling: gain=1
		(KO) cuts IP3R release vs gain=0; at WT (scale=1) the gate is inert."""
		from reconstruction.platelet.simulation_data import SimulationDataPlatelet
		sim_data = SimulationDataPlatelet()
		all_ids = list(sim_data.internal_state.bulk_molecules.bulk_data['id'])
		all_counts = sim_data.internal_state.bulk_molecules.initial_counts
		names = list(cs.MOLECULE_NAMES)
		y = np.zeros(cs.N_SPECIES)
		for i, name in enumerate(names):
			y[i] = float(all_counts[all_ids.index(name)])
		# Evoked: high cyt Ca²⁺ and IP3 so the IP3R is actively releasing.
		y[names.index('CA2_CYT[c]')] = round(0.5 / cs._UM_PER_COUNT_CYT)
		y[names.index('IP3[c]')] = round(0.3 / cs._UM_PER_COUNT_CYT)
		cyt, dts = names.index('CA2_CYT[c]'), names.index('CA2_DTS[dts]')

		def dy(scale, gain):
			cfg = RunConfig(mcu_vmax_scale=scale, mito_coupling_gain=gain)
			return cs._ode_rhs(0.0, y, 5.0, cfg, {})

		ko_on, ko_off = dy(0.0, 1.0), dy(0.0, 0.0)
		self.assertLess(ko_on[cyt], ko_off[cyt],
			'KO coupling should cut Ca²⁺ release into the cytosol.')
		self.assertGreater(ko_on[dts], ko_off[dts],
			'KO coupling should leave the store fuller (less depletion).')
		wt_on, wt_off = dy(1.0, 1.0), dy(1.0, 0.0)
		self.assertAlmostEqual(wt_on[cyt], wt_off[cyt], places=9,
			msg='Coupling must be inert at wild type (scale=1).')


if __name__ == '__main__':
	unittest.main()
