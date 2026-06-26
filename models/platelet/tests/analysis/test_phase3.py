"""Regression tests for the Phase 3 Dolan 2014 Fig. 4 validation.

Runs both conditions (with and without extracellular Ca²⁺), invokes the
plotting + summary code, and asserts the headline numbers stay within
tolerance of the current baseline. Locks Phase 3 output against future
drift — analogous to ``test_regression.py`` but for the two-condition
comparison.

**Role of "Dolan 5/5" (re-reviewed 2026-06-19).** As the model has grown
beyond the Ca²⁺ transient (PI cycle, GPCR cascade, PKC feedback, secretion,
thromboxane, integrin, and the v0.7 inhibitory axis), this 5/5 is best read
as a **regression invariant on the Ca²⁺ core**, not the model's primary
validation. Because cytosolic Ca²⁺ is store-limited/SOCE-clamped (see
``project_pkc_ca_invisible`` / the validation-map doc), many parameterisations
pass it — it has low discriminating power for new biology, and each new layer
is deliberately normalised at rest to *keep* it passing. The model's evidence
now rests on a portfolio of subsystem targets (VASP/PRI, PAC-1, lumi-
aggregometry, drug dose-response) — see
``reports/design/validation-map-2026-06-19.qmd``. Keep this test as the cheap
guard that the Ca²⁺ core is not accidentally broken.

Current (seed 0, length 200 s, both conditions): +Ca_ex peak ≈ 521 nM,
−Ca_ex peak ≈ 296 nM, SOCE differential ≈ 225 nM.

The all-or-nothing "Dolan 5/5" gate was retired 2026-06-26: the behavioural
Dolan validation (paired peaks in band + the SOCE differential) now lives in
``test_acceptance.py::TestDolanTransient``. ``_eval_criteria`` still annotates
the plot, but is no longer a pass/fail gate. See
``docs/validation-and-regressions.md``.
"""

import os
import tempfile
import unittest

import pytest

from models.platelet.analysis.phase3_dolan_fig4 import make_phase3_plot
from runscripts.manual.runPlateletSim import run_platelet_sim


GOLDEN_SEED = 0
GOLDEN_LENGTH_SEC = 200


@pytest.mark.slow
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
		"""+Ca_ex peak Ca²⁺ within 30% of 530 nM (v0.3 PI-cycle baseline).

		Baseline raised 2026-05-12 (Phase 4 / #31) after the forced IP3
		curve was replaced by the Mazet 2020 PI cycle. The +Ca_ex peak
		now includes both IP3R-driven release (from PI cycle IP3) and
		P2X1 entry, sustained by the PLCβ-mediated IP3 timecourse.
		"""
		peak = self.summary['with_extracellular_ca']['peak_cyt_nM']
		self.assertGreater(peak, 200.0,
			f'+Ca_ex peak {peak:.1f} nM below the active threshold')
		# ±30% of the 530 nM PI-cycle baseline
		self.assertGreater(peak, 371.0,
			f'+Ca_ex peak {peak:.1f} regressed below baseline −30%')
		self.assertLess(peak, 689.0,
			f'+Ca_ex peak {peak:.1f} regressed above baseline +30%')

	def test_no_ca_peak_in_range(self):
		"""−Ca_ex peak Ca²⁺ within 30% of 325 nM baseline (IP3R-only response)."""
		peak = self.summary['without_extracellular_ca']['peak_cyt_nM']
		self.assertGreater(peak, 200.0,
			f'−Ca_ex peak {peak:.1f} nM below the active threshold')
		# ±30% of the 325 nM post-k12 baseline (unchanged since v0.2.5)
		self.assertGreater(peak, 228.0,
			f'−Ca_ex peak {peak:.1f} regressed below baseline −30%')
		self.assertLess(peak, 423.0,
			f'−Ca_ex peak {peak:.1f} regressed above baseline +30%')


if __name__ == '__main__':
	unittest.main()
