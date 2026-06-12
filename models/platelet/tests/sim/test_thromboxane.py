"""Tests for ThromboxaneSynthesis + ThromboxaneTrace (v0.61 Slice A).

Production-only slice: TXA₂ is synthesised (Ca×PKC×COX-1 gated) and decays to
TXB₂, with no Gq feedback yet. Covers:
  * resting quiescence — no TXA₂ when un-stimulated;
  * activation — PKC + Ca²⁺ drive TXA₂ synthesis, TXB₂ accumulates;
  * aspirin — COX-1 knockout (cox1_factor = 0) abolishes thromboxane.
"""

import os
import tempfile

import numpy as np
import pytest

from wholecell.io.tablereader import TableReader
from runscripts.manual.runPlateletSim import run_platelet_sim
import reconstruction.platelet.dataclasses.process.thromboxane_synthesis as tx_mod


def _run(out_dir, length, **sim_kwargs):
	paths = run_platelet_sim(out_dir, length_sec=length, seed=0,
		log_to_shell=False, ca_ex_mM=1.2, **sim_kwargs)
	return paths['sim_out_dir']


def _tx(sim_out_dir, col):
	reader = TableReader(os.path.join(sim_out_dir, 'ThromboxaneTrace'))
	return reader.readColumn(col).flatten()


@pytest.mark.slow
class TestThromboxaneSynthesis:

	def test_resting_quiescence_no_txa2(self):
		"""Un-stimulated platelet: TXA₂ synthesis stays exactly zero."""
		with tempfile.TemporaryDirectory() as td:
			sim_out = _run(td, 30, thrombin_peak_nM=0.0,
				adp_peak_uM=0.0, atp_ex_peak_uM=0.0)
			assert _tx(sim_out, 'txa2_synth_gate').max() == 0.0
			assert _tx(sim_out, 'txa2_uM').max() == 0.0
			assert _tx(sim_out, 'txb2').max() == 0

	def test_activation_produces_thromboxane(self):
		"""Standard +Ca²⁺ transient drives TXA₂ synthesis; TXB₂ accumulates."""
		with tempfile.TemporaryDirectory() as td:
			sim_out = _run(td, 100)
			gate = _tx(sim_out, 'txa2_synth_gate')
			txa2 = _tx(sim_out, 'txa2_uM')
			txb2 = _tx(sim_out, 'txb2')
			assert gate[0] == 0.0           # zero until PKC activates
			assert txa2.max() > 0.3         # reaches a TP-relevant level
			assert txb2[-1] > 0
			# TXB₂ is the stable metabolite — monotonically non-decreasing.
			assert np.all(np.diff(txb2) >= 0)

	def test_aspirin_abolishes_thromboxane(self):
		"""COX-1 knockout (cox1_factor = 0) → no TXA₂ / TXB₂ at all."""
		try:
			tx_mod.COX1_FACTOR = 0.0
			with tempfile.TemporaryDirectory() as td:
				sim_out = _run(td, 100)
				assert _tx(sim_out, 'txa2_uM').max() == 0.0
				assert _tx(sim_out, 'txb2').max() == 0
		finally:
			tx_mod.COX1_FACTOR = 1.0
