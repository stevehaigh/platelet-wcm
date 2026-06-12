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
		"""Secretion relocates cargo: granule + secreted (+ cleared) totals
		are constant. ADP additionally passes ADP[e] → AMP[e] via NTPDase."""
		with tempfile.TemporaryDirectory() as td:
			sim_out = _run(td, 60)
			ids, counts = _bulk(sim_out)
			idx = {n: i for i, n in enumerate(ids)}
			pools = {
				'ADP': ('ADP[dg]', 'ADP[e]', 'AMP[e]'),   # + NTPDase product
				'5HT': ('5HT[dg]', '5HT[e]'),
				'FGA': ('FGA[ag]', 'FGA[e]'),
				'SELP': ('SELP[ag]', 'SELP_surface[pl]'),
			}
			for name, members in pools.items():
				total = sum(counts[:, idx[m]] for m in members)
				assert np.all(total == total[0]), f'{name} pool not conserved'
				# Something actually left the lumen by the end of activation.
				assert counts[-1, idx[members[1]]] + (
					counts[-1, idx[members[2]]] if len(members) > 2 else 0) > 0

	def test_autocrine_adp_drives_p2y1(self):
		"""Thrombin-only (no exogenous ADP): secreted ADP closes the loop —
		it drives P2Y1 (→ desensitisation) and is self-limited by NTPDase."""
		with tempfile.TemporaryDirectory() as td:
			sim_out = _run(td, 150, adp_peak_uM=0.0, atp_ex_peak_uM=0.0)
			adp_e_uM = _secretion(sim_out, 'adp_e_uM')
			amp_e = _secretion(sim_out, 'amp_e')
			p2y1_des = TableReader(os.path.join(sim_out, 'CalciumTrace')
				).readColumn('p2y1_desensitised_frac').flatten()
			# Secreted ADP reaches a P2Y1-relevant pericellular level...
			assert adp_e_uM.max() > 1.0
			# ...and drives P2Y1 even with zero exogenous ADP (loop closed).
			assert p2y1_des.max() > 0.1
			# Self-limiting: ADP[e] peaks then clears (NTPDase → AMP rises).
			assert adp_e_uM[-1] < adp_e_uM.max()
			assert amp_e[-1] > 0
