"""Golden-run regression tests for the platelet simulation.

These tests run the simulation with a fixed seed and length and assert that
key metrics stay within tolerances of the baseline captured from the source
repo (stevehaigh/wcEcoli, platelet branch, 2026-05-05).  They protect the
migration by detecting any unexpected change in simulation output.

Baseline (seed=42, length_sec=60):
    ca_cyt[0]   ≈  99.9 nM   (resting, CaM pre-equilibrated)
    ca_cyt peak ≈ 279.8 nM   (step 1, IP3R burst)
    ca_dts[0]   ≈ 250.0 µM   (initial DTS store)
    ip3[0]      ≈  50.1 nM   (resting IP3)
    dry_mass     > 0 throughout

Note: DTS drains to zero by ~60 s — this is a known Phase 2 issue
(GAMMA_IP3R_S needs recalibration; tracked in issue #48).  The final Ca_dts
value is therefore NOT asserted here to avoid locking in the broken behaviour.
"""

import os
import tempfile
import unittest

import numpy as np

from runscripts.manual.runPlateletSim import run_platelet_sim
from wholecell.io.tablereader import TableReader

GOLDEN_SEED = 42
GOLDEN_LENGTH_SEC = 60


class TestPlateletGoldenRun(unittest.TestCase):
	"""Regression suite: key metrics must stay within tolerance of baseline."""

	@classmethod
	def setUpClass(cls):
		"""Run one simulation shared by all test methods in this class."""
		cls._tmpdir = tempfile.TemporaryDirectory()
		paths = run_platelet_sim(
			cls._tmpdir.name,
			length_sec=GOLDEN_LENGTH_SEC,
			seed=GOLDEN_SEED,
			log_to_shell=False,
		)
		sim_out = paths['sim_out_dir']

		ca_reader = TableReader(os.path.join(sim_out, 'CalciumTrace'))
		cls.ca_cyt_nM     = ca_reader.readColumn('ca_cyt_nM').flatten()
		cls.ca_dts_uM     = ca_reader.readColumn('ca_dts_uM').flatten()
		cls.ip3_nM        = ca_reader.readColumn('ip3_nM').flatten()
		cls.soce_flux_nMs = ca_reader.readColumn('soce_flux_nMs').flatten()

		mass_reader = TableReader(os.path.join(sim_out, 'Mass'))
		cls.dry_mass = mass_reader.readColumn('dryMass').flatten()

	@classmethod
	def tearDownClass(cls):
		cls._tmpdir.cleanup()

	# ── Resting-state initial conditions ──────────────────────────────────────

	def test_resting_ca_cyt_near_100_nM(self):
		"""Initial cytosolic Ca²⁺ should be ~100 nM (CaM pre-equilibrated ICs)."""
		self.assertGreater(self.ca_cyt_nM[0], 80.0,
			'Resting Ca²⁺ too low — check CaM ICs or PMCA balance')
		self.assertLess(self.ca_cyt_nM[0], 120.0,
			'Resting Ca²⁺ too high — check CaM ICs or PMCA balance')

	def test_resting_ca_dts_near_250_uM(self):
		"""Initial DTS Ca²⁺ should be ~250 µM (Purvis 2008 baseline)."""
		self.assertGreater(self.ca_dts_uM[0], 230.0,
			'Initial DTS store too low')
		self.assertLess(self.ca_dts_uM[0], 270.0,
			'Initial DTS store too high')

	def test_resting_ip3_near_50_nM(self):
		"""Resting IP3 should be ~50 nM (Purvis 2008 baseline)."""
		self.assertGreater(self.ip3_nM[0], 40.0,
			'Resting IP3 too low')
		self.assertLess(self.ip3_nM[0], 60.0,
			'Resting IP3 too high')

	# ── Peak calcium transient ─────────────────────────────────────────────────

	def test_peak_ca_cyt_in_physiological_range(self):
		"""Peak cytosolic Ca²⁺ should be 200–800 nM (Dolan 2014 Fig 4 target)."""
		peak = self.ca_cyt_nM.max()
		self.assertGreater(peak, 200.0,
			f'Peak Ca²⁺ {peak:.1f} nM below physiological range (200–800 nM)')
		self.assertLess(peak, 800.0,
			f'Peak Ca²⁺ {peak:.1f} nM above physiological range (200–800 nM)')

	def test_peak_ca_cyt_near_baseline(self):
		"""Peak cytosolic Ca²⁺ should be within 30% of the 280 nM baseline."""
		peak = self.ca_cyt_nM.max()
		self.assertGreater(peak, 196.0,   # 280 * 0.70
			f'Peak Ca²⁺ {peak:.1f} nM regressed below baseline −30%')
		self.assertLess(peak, 364.0,      # 280 * 1.30
			f'Peak Ca²⁺ {peak:.1f} nM regressed above baseline +30%')

	# ── Flux sanity checks ─────────────────────────────────────────────────────

	def test_soce_flux_non_negative(self):
		"""SOCE flux is an influx — must be ≥ 0 throughout."""
		self.assertTrue(
			np.all(self.soce_flux_nMs >= 0.0),
			'SOCE flux went negative — check SOCE rate equation sign')

	# ── Mass conservation ──────────────────────────────────────────────────────

	def test_dry_mass_always_positive(self):
		"""Dry mass must remain positive at every timestep."""
		self.assertTrue(
			np.all(self.dry_mass > 0),
			f'Dry mass reached zero or below (min={self.dry_mass.min():.6f})')

	def test_dry_mass_near_baseline(self):
		"""Initial dry mass should be within 1% of the 196.17 fg baseline."""
		initial = self.dry_mass[0]
		self.assertGreater(initial, 194.2,   # 196.17 * 0.99
			f'Initial dry mass {initial:.3f} below baseline −1%')
		self.assertLess(initial, 198.1,      # 196.17 * 1.01
			f'Initial dry mass {initial:.3f} above baseline +1%')
