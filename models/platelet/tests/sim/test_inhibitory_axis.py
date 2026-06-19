"""
Behavioural tests for the inhibitory axis — P2Y12 / Gi / cAMP / PKA (issue #10).

These lock in the issue #10 acceptance criteria as integration tests:

  1. extracellular ADP (via P2Y12 → Gi) lowers cAMP;
  2. the resulting drop in PKA dis-inhibits the αIIbβ3 inside-out step, so the
     PAC-1 readout is larger with P2Y12 active than under P2Y12 blockade
     (the clopidogrel mechanism);
  3. blocking P2Y12 (``p2y12_block = 1``) keeps cAMP at its basal tone and
     reduces activation.

Plus the load-bearing invariant: the PKA brake is normalised to 1.0 at resting
cAMP, so the resting Ca²⁺ fixed point is unchanged (Dolan goldens preserved —
see test_validation_targets / test_byte_identical).

The clinical VASP/PRI direction (phospho-VASP falls under ADP, preserved under
blockade) is checked as the quantitative pharmacology readout.
"""

import os
import tempfile
import unittest

import pytest

from reconstruction.platelet.run_config import RunConfig
from runscripts.manual.runPlateletSim import run_platelet_sim
from wholecell.io.tablereader import TableReader


def _run(**cfg):
	"""Run a 150 s sim and harvest the inhibitory-axis traces + PAC-1."""
	rc = RunConfig(ca_ex_mM=1.2, **cfg)
	paths = run_platelet_sim(
		tempfile.mkdtemp(), length_sec=150, seed=0,
		log_to_shell=False, run_config=rc)
	sim_out = paths['sim_out_dir']
	ct = TableReader(os.path.join(sim_out, 'CalciumTrace'))
	rb = TableReader(os.path.join(sim_out, 'BulkMolecules'))
	ids = list(rb.readAttribute('objectNames'))
	counts = rb.readColumn('counts')

	def bcol(name):
		return counts[:, ids.index(name)].astype(float)

	act = bcol('aIIbb3_active[pl]')
	rest = bcol('aIIbb3_resting[pl]')
	total = act + rest
	return {
		'camp_uM': ct.readColumn('camp_uM').flatten(),
		'vasp_phos_frac': ct.readColumn('vasp_phos_frac').flatten(),
		'p2y12_active_frac': ct.readColumn('p2y12_active_frac').flatten(),
		'ca_cyt_nM': ct.readColumn('ca_cyt_nM').flatten(),
		'pac1': 100.0 * act / total,
	}


@pytest.mark.slow
class TestInhibitoryAxis(unittest.TestCase):
	"""P2Y12 / cAMP / PKA behavioural criteria (issue #10)."""

	@classmethod
	def setUpClass(cls):
		cls.rest = _run(thrombin_peak_nM=0, adp_peak_uM=0, atp_ex_peak_uM=0)
		cls.wt = _run()
		cls.clopidogrel = _run(p2y12_block=1.0)

	def test_resting_camp_at_basal(self):
		"""At rest (ADP = 0) P2Y12 is silent → cAMP stays at its 1 µM tone."""
		self.assertAlmostEqual(self.rest['camp_uM'][0], 1.0, places=2)
		self.assertGreater(self.rest['camp_uM'].min(), 0.98)
		# No P2Y12 occupancy at rest (the basis of the brake = 1 invariant).
		self.assertEqual(self.rest['p2y12_active_frac'].max(), 0.0)

	def test_agonist_lowers_camp(self):
		"""Criterion 1: ADP → P2Y12 → Gi lowers cAMP well below basal."""
		self.assertLess(self.wt['camp_uM'].min(), 0.6)
		self.assertGreater(self.wt['p2y12_active_frac'].max(), 0.3)

	def test_clopidogrel_keeps_camp_high(self):
		"""Criterion 3: P2Y12 blockade → no P2Y12 activation, cAMP at basal."""
		self.assertEqual(self.clopidogrel['p2y12_active_frac'].max(), 0.0)
		self.assertGreater(self.clopidogrel['camp_uM'].min(), 0.98)

	def test_clopidogrel_reduces_integrin_activation(self):
		"""Criterion 2 (functional output): less PAC-1 with P2Y12 blocked.

		The PKA brake bites in the rising phase, so compare at t≈100 s.
		"""
		self.assertGreater(
			self.wt['pac1'][100], self.clopidogrel['pac1'][100] + 5.0)

	def test_vasp_pri_readout_direction(self):
		"""Clinical VASP/PRI: phospho-VASP falls under ADP, preserved on block."""
		self.assertLess(self.wt['vasp_phos_frac'].min(), self.rest['vasp_phos_frac'][0])
		self.assertGreater(
			self.clopidogrel['vasp_phos_frac'].min(),
			self.wt['vasp_phos_frac'].min())

	def test_resting_fixed_point_preserved(self):
		"""Brake = 1 at rest → resting cytosolic Ca²⁺ fixed point unchanged."""
		self.assertGreater(self.rest['ca_cyt_nM'][0], 80.0)
		self.assertLess(self.rest['ca_cyt_nM'][0], 120.0)
