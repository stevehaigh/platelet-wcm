"""Regression tests for the Phase 3 Dolan 2014 Fig. 4 validation.

Runs both conditions (with and without extracellular Ca²⁺), invokes the
plotting + summary code, and asserts the headline numbers stay within
tolerance of the current baseline. Locks Phase 3 output against future
drift — analogous to ``test_regression.py`` but for the two-condition
comparison.

Baseline (seed=0, length=200 s, both conditions, post-Caride-k12-fix
2026-05-07):
  +Ca_ex peak Ca_cyt:  ~380 nM at t=1 s
  −Ca_ex peak Ca_cyt:  ~325 nM at t=1 s
  SOCE differential:   ~54 nM   (FAILS Dolan ≥100 nM criterion;
                                 documented as the remaining open issue
                                 19b / #22, not a regression)

Phase 3 acceptance now passes 4 of 5 criteria (was 3/5 before
Caride k12 fix; the +Ca_ex Dolan ±30% peak criterion flipped to PASS
when k11/k12 reactions were restored, releasing 4 Ca²⁺ during PMCA-CaM
deactivation and lifting the +Ca_ex peak from 299 → 380 nM).
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
		"""4 of 5 acceptance criteria pass at the post-Caride-k12 baseline.

		The remaining failing criterion — SOCE differential — traces to
		the resting-state DTS depletion issue (DTS empties before SOCE can
		establish a plateau). Documented in design doc §6.8 D7 and the
		2026-05-06/07 lab books; ticketed as #19b (close cyt/DTS biology
		gap) and #22 (MCU). Asserting the count locks the current state —
		a future change that fixes the SOCE differential would tick this
		test up to 5/5.
		"""
		passed = sum(1 for c in self.summary['criteria'] if c['passed'])
		self.assertEqual(passed, 4,
			f'Expected 4/5 Phase 3 criteria pass at post-k12 baseline; got '
			f'{passed}. If a fix landed, update this assertion and the '
			f'lab book entry.')


if __name__ == '__main__':
	unittest.main()
