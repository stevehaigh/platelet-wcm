"""Unit tests for the calcium_signalling ODE solver.

Focused on parameter-binding regressions that are easy to miss from
end-to-end sim tests — in particular, whether agonist-peak overrides
take effect at call time rather than being baked into function
defaults at import time.
"""

import unittest

import numpy as np

from reconstruction.platelet.dataclasses.process import calcium_signalling as cs


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
		# from extracellular thrombin must be exactly zero.
		dy_rest = cs._ode_rhs(0.0, y0, 100.0, 0.0, 0.0, 0.0, 0.0)
		par1_idx = list(cs.MOLECULE_NAMES).index('PAR1_active[pl]')
		par4_idx = list(cs.MOLECULE_NAMES).index('PAR4_active[pl]')

		# At rest with thrombin = 0, the only sink for PAR*_active is
		# internalization (negative). New cleavage must be 0, so total dy/dt
		# for active receptors is ≤ 0.
		self.assertLessEqual(dy_rest[par1_idx], 0.0,
			'PAR1_active should not accumulate at rest.')
		self.assertLessEqual(dy_rest[par4_idx], 0.0,
			'PAR4_active should not accumulate at rest.')


if __name__ == '__main__':
	unittest.main()
