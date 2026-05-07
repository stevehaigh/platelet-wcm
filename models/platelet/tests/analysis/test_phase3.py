"""Regression tests for the Phase 3 Dolan 2014 Fig. 4 validation.

Runs both conditions (with and without extracellular Ca²⁺), invokes the
plotting + summary code, and asserts the headline numbers stay within
tolerance of the current baseline. Locks Phase 3 output against future
drift — analogous to ``test_regression.py`` but for the two-condition
comparison.

Baseline (seed=0, length=200 s, both conditions, post-IP3R-subunit-/4
fix 2026-05-07):
  +Ca_ex peak Ca_cyt:  ~391 nM at t=1 s
  −Ca_ex peak Ca_cyt:  ~362 nM at t=1 s   (1% over Dolan ±30% ceiling)
  SOCE differential:   ~29 nM   (still FAILS Dolan ≥100 nM criterion)

Phase 3 acceptance: 3 of 5 (regressed from 4/5 with the IP3R subunit
correction). The −Ca_ex peak criterion flipped FAIL because the
slower drain after correcting the 4× IP3R-channel-count overstatement
lets cyt build higher in the EDTA condition. Diagnosis in
lab-book-2026-05-07-dts-drain-investigation.md; expected to recover
once downstream calibration (J_PM_LEAK, γ_SOC, k_dim_f) is rerun
against the corrected leak.
"""

import os
import tempfile
import unittest

from models.platelet.analysis.phase3_dolan_fig4 import make_phase3_plot
from runscripts.manual.runPlateletSim import run_platelet_sim


GOLDEN_SEED = 0
GOLDEN_LENGTH_SEC = 200


class TestPhase3DolanFig4(unittest.TestCase):
	"""Regression: Phase 3 two-condition pipeline must produce stable peaks."""

	@classmethod
	def setUpClass(cls):
		cls._tmpdir = tempfile.TemporaryDirectory()
		root = cls._tmpdir.name

		with_ca_dir = os.path.join(root, 'with_ca')
		no_ca_dir   = os.path.join(root, 'no_ca')
		os.makedirs(with_ca_dir)
		os.makedirs(no_ca_dir)

		paths_with = run_platelet_sim(
			with_ca_dir, length_sec=GOLDEN_LENGTH_SEC, seed=GOLDEN_SEED,
			log_to_shell=False, ca_ex_mM=1.2)
		paths_no = run_platelet_sim(
			no_ca_dir, length_sec=GOLDEN_LENGTH_SEC, seed=GOLDEN_SEED,
			log_to_shell=False, ca_ex_mM=0.0)

		plot_path = os.path.join(root, 'phase3_dolan_fig4.png')
		cls.summary = make_phase3_plot(
			paths_with['sim_out_dir'],
			paths_no['sim_out_dir'],
			plot_path,
		)
		cls.plot_path = plot_path

	@classmethod
	def tearDownClass(cls):
		cls._tmpdir.cleanup()

	# ── Plot artefacts ────────────────────────────────────────────────────────

	def test_plot_file_written(self):
		"""The PNG figure should be written to the requested path."""
		self.assertTrue(os.path.isfile(self.plot_path),
			f'Phase 3 plot not written at {self.plot_path}')
		self.assertGreater(os.path.getsize(self.plot_path), 50_000,
			'Phase 3 PNG suspiciously small; render likely failed')

	# ── Headline metrics ──────────────────────────────────────────────────────

	def test_with_ca_peak_in_range(self):
		"""+Ca_ex peak Ca²⁺ in the active band (≥ 200 nM) and within 30% of baseline."""
		peak = self.summary['with_extracellular_ca']['peak_cyt_nM']
		self.assertGreater(peak, 200.0,
			f'+Ca_ex peak {peak:.1f} nM below the active threshold')
		# ±30% of the 380 nM post-k12 baseline
		self.assertGreater(peak, 266.0,
			f'+Ca_ex peak {peak:.1f} regressed below baseline −30%')
		self.assertLess(peak, 494.0,
			f'+Ca_ex peak {peak:.1f} regressed above baseline +30%')

	def test_no_ca_peak_in_range(self):
		"""−Ca_ex peak Ca²⁺ in the active band and within 30% of baseline."""
		peak = self.summary['without_extracellular_ca']['peak_cyt_nM']
		self.assertGreater(peak, 200.0,
			f'−Ca_ex peak {peak:.1f} nM below the active threshold')
		# ±30% of the 325 nM post-k12 baseline
		self.assertGreater(peak, 228.0,
			f'−Ca_ex peak {peak:.1f} regressed below baseline −30%')
		self.assertLess(peak, 423.0,
			f'−Ca_ex peak {peak:.1f} regressed above baseline +30%')

	# ── Acceptance-criteria pass/fail count ───────────────────────────────────

	def test_criteria_pass_count(self):
		"""3 of 5 acceptance criteria pass after the IP3R subunit /4 fix.

		Two failing criteria:
		- SOCE differential (29 nM, target ≥ 100 nM) — depends on the
		  resting-state DTS depletion issue tracked in #24 and the
		  candidate-2 work (DTS Ca²⁺ buffers).
		- −Ca_ex peak in Dolan ±30% (362 nM, ceiling 358 nM) — flipped
		  FAIL when the 4× IP3R-leak overstatement was corrected; the
		  slower drain lets cyt build higher in EDTA. Recovery expected
		  once J_PM_LEAK and other calibrated knobs are re-derived
		  against the corrected leak.

		Asserting the count locks the current state. A future fix that
		recovers either criterion ticks this test up to 4/5 or 5/5.
		"""
		passed = sum(1 for c in self.summary['criteria'] if c['passed'])
		self.assertEqual(passed, 3,
			f'Expected 3/5 Phase 3 criteria at the post-IP3R-subunit-/4 '
			f'baseline; got {passed}. If a fix landed, update this '
			f'assertion and the lab book entry.')


if __name__ == '__main__':
	unittest.main()
