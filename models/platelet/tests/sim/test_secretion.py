"""Tests for the GranuleSecretion process + SecretionTrace listener (v0.61).

Covers the three invariants the design (pkc-downstream-effects-2026-06-12 §1)
calls for:
  * resting quiescence — no secretion when the platelet is un-stimulated;
  * activation release — PKC + Ca²⁺ drive graded dense/α-granule release;
  * cargo conservation — secretion relocates cargo, never creates/destroys it.
"""

import os
import tempfile

import numpy as np
import pytest

from wholecell.io.tablereader import TableReader
from runscripts.manual.runPlateletSim import run_platelet_sim

SEC_FRACTIONS = (
	'adp_released_frac', 'serotonin_released_frac',
	'fibrinogen_released_frac', 'pselectin_surface_frac',
)


def _run(out_dir, length, **sim_kwargs):
	paths = run_platelet_sim(out_dir, length_sec=length, seed=0,
		log_to_shell=False, ca_ex_mM=1.2, **sim_kwargs)
	return paths['sim_out_dir']


def _secretion(sim_out_dir, col):
	reader = TableReader(os.path.join(sim_out_dir, 'SecretionTrace'))
	return reader.readColumn(col).flatten()


def _bulk(sim_out_dir):
	reader = TableReader(os.path.join(sim_out_dir, 'BulkMolecules'))
	ids = list(reader.readAttribute('objectNames'))
	counts = reader.readColumn('counts')
	return ids, counts


@pytest.mark.slow
class TestGranuleSecretion:

	def test_resting_quiescence_no_secretion(self):
		"""Un-stimulated platelet (all agonists off): secretion stays exactly 0."""
		with tempfile.TemporaryDirectory() as td:
			sim_out = _run(td, 30, thrombin_peak_nM=0.0,
				adp_peak_uM=0.0, atp_ex_peak_uM=0.0)
			assert _secretion(sim_out, 'secretion_gate').max() == 0.0
			for col in SEC_FRACTIONS:
				assert _secretion(sim_out, col).max() == 0.0, col
			assert _secretion(sim_out, 'adp_e').max() == 0

	def test_activation_releases_cargo(self):
		"""Standard +Ca²⁺ transient drives graded dense + α-granule release."""
		with tempfile.TemporaryDirectory() as td:
			sim_out = _run(td, 60)
			gate = _secretion(sim_out, 'secretion_gate')
			# Gate is zero at t=0 (PKC* not yet activated) and rises.
			assert gate[0] == 0.0
			assert gate.max() > 0.1
			# Dense granules (ADP/5-HT) release faster than α-granules.
			adp = _secretion(sim_out, 'adp_released_frac')
			fga = _secretion(sim_out, 'fibrinogen_released_frac')
			psel = _secretion(sim_out, 'pselectin_surface_frac')
			assert adp[-1] > 0.5
			assert psel[-1] > 0.0
			assert adp[-1] > fga[-1]          # dense kinetics faster than α
			# Monotonic, non-decreasing release (cargo only leaves the granule).
			assert np.all(np.diff(adp) >= -1e-12)

	def test_cargo_count_conserved(self):
		"""Secretion relocates cargo: granule + secreted totals are constant."""
		with tempfile.TemporaryDirectory() as td:
			sim_out = _run(td, 60)
			ids, counts = _bulk(sim_out)
			idx = {n: i for i, n in enumerate(ids)}
			for lumen, dest in (
					('ADP[dg]', 'ADP[e]'),
					('5HT[dg]', '5HT[e]'),
					('FGA[ag]', 'FGA[e]'),
					('SELP[ag]', 'SELP_surface[pl]')):
				total = counts[:, idx[lumen]] + counts[:, idx[dest]]
				assert np.all(total == total[0]), f'{lumen}+{dest} not conserved'
				# Something was actually released by the end of activation.
				assert counts[-1, idx[dest]] > 0
