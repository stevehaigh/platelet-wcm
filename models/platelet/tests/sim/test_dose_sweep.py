"""Smoke test for the 2-D dose-response sweep driver.

Runs a 2 × 2 / 10-s grid (4 cells, ~10 s wall-clock total) and asserts the
expected artefacts are written, all observables harvest, and the corner-
to-corner trend on peak Ca²⁺ goes the right way for a time-limited sim.
"""

import os
import tempfile
import unittest

import numpy as np
import pytest

from runscripts.manual.runDoseSweep import (
	OBSERVABLES,
	run_dose_sweep,
	write_all_outputs,
)


@pytest.mark.slow
class TestDoseSweepSmoke(unittest.TestCase):
	"""End-to-end on a 2x2 grid; CI-safe."""

	@classmethod
	def setUpClass(cls):
		cls._tmpdir = tempfile.TemporaryDirectory()
		cls.sweep = run_dose_sweep(
			cls._tmpdir.name,
			grid=2,
			length_sec=10,
			adp_range_uM=(0.1, 10.0),
			thr_range_nM=(0.01, 1.0),
			seed=0,
			keep_cell_output=False,
			log_to_shell=False,
		)
		write_all_outputs(cls.sweep, cls._tmpdir.name)

	@classmethod
	def tearDownClass(cls):
		cls._tmpdir.cleanup()

	def test_sweep_shape_and_grids_log_spaced(self):
		"""Grids should be log-spaced and every matrix shape (n_thr, n_adp)."""
		self.assertEqual(self.sweep.adp_grid.shape, (2,))
		self.assertEqual(self.sweep.thr_grid.shape, (2,))
		np.testing.assert_allclose(self.sweep.adp_grid, [0.1, 10.0])
		np.testing.assert_allclose(self.sweep.thr_grid, [0.01, 1.0])
		for key, *_ in OBSERVABLES:
			self.assertIn(key, self.sweep.matrices,
				f'observable {key} not harvested')
			self.assertEqual(self.sweep.matrices[key].shape, (2, 2),
				f'observable {key} has wrong shape')

	def test_all_artefacts_written(self):
		expected = ['sweep.csv', 'sweep.npz', 'sweep_summary.json',
			'sweep_panel.png', 'sweep_surface.png']
		expected += [f'sweep_{k}.png' for k, *_ in OBSERVABLES]
		for name in expected:
			path = os.path.join(self._tmpdir.name, name)
			self.assertTrue(os.path.isfile(path), f'missing artefact: {name}')
			self.assertGreater(os.path.getsize(path), 0, f'empty artefact: {name}')

	def test_peak_ca_positive(self):
		"""Every cell should produce a positive peak Ca²⁺."""
		self.assertTrue(np.all(self.sweep.matrices['peak_ca_nM'] > 0.0))

	def test_peak_ip3_positive(self):
		"""Every cell should produce a positive peak IP3."""
		self.assertTrue(np.all(self.sweep.matrices['peak_ip3_nM'] > 0.0))

	def test_corner_monotone_peak_ca(self):
		"""High-corner (high ADP, high thrombin) > low-corner peak Ca²⁺.

		Weak monotonicity only — strict cell-by-cell monotonicity is
		biology-fragile and the model saturates at long sim length. The
		smoke test uses 10 s where time-limited PAR1 cleavage gives a
		clear corner gradient.
		"""
		m = self.sweep.matrices['peak_ca_nM']
		self.assertGreater(m[-1, -1], m[0, 0],
			'High-agonist corner should produce a higher peak than the '
			'low-agonist corner; check the sweep harness.')

	def test_cell_output_pruned_by_default(self):
		"""With keep_cell_output=False, only CalciumTrace survives per cell."""
		for i in range(2):
			for j in range(2):
				cell_dir = os.path.join(self._tmpdir.name, 'cells', f'{i:02d}_{j:02d}')
				simout_dirs = [d for d, _, _ in os.walk(cell_dir)
					if os.path.basename(d) == 'simOut']
				self.assertEqual(len(simout_dirs), 1,
					f'expected exactly one simOut/ under {cell_dir}')
				entries = set(os.listdir(simout_dirs[0]))
				self.assertIn('CalciumTrace', entries)
				self.assertNotIn('BulkMolecules', entries)
				self.assertNotIn('Mass', entries)

	def test_npz_roundtrip(self):
		"""npz reload preserves grids and all observable matrices."""
		data = np.load(os.path.join(self._tmpdir.name, 'sweep.npz'))
		np.testing.assert_array_equal(data['adp_grid'], self.sweep.adp_grid)
		np.testing.assert_array_equal(data['thr_grid'], self.sweep.thr_grid)
		for key, *_ in OBSERVABLES:
			np.testing.assert_array_equal(data[key], self.sweep.matrices[key])


if __name__ == '__main__':
	unittest.main()
