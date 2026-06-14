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
from reconstruction.platelet.run_config import RunConfig


def _run(out_dir, length, **config_kwargs):
	"""Run one sim; config_kwargs are RunConfig fields (peaks, cox1_factor, …)."""
	run_config = RunConfig(ca_ex_mM=1.2, **config_kwargs)
	paths = run_platelet_sim(out_dir, length_sec=length, seed=0,
		log_to_shell=False, run_config=run_config)
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
		with tempfile.TemporaryDirectory() as td:
			sim_out = _run(td, 100, cox1_factor=0.0)
			assert _tx(sim_out, 'txa2_uM').max() == 0.0
			assert _tx(sim_out, 'txb2').max() == 0

	def test_txa2_drives_tp_and_amplifies_gq(self):
		"""Slice B: synthesised TXA₂ activates TP and amplifies the Gq cascade
		(IP₃); aspirin (no TXA₂) leaves TP inactive and the loop open."""
		from wholecell.io.tablereader import TableReader

		def ip3(sim_out):
			return TableReader(os.path.join(sim_out, 'CalciumTrace')
				).readColumn('ip3_nM').flatten()

		with tempfile.TemporaryDirectory() as td:
			on = _run(os.path.join(td, 'on'), 150)
			tp_on = _tx(on, 'tp_active_frac')
			ip3_on = ip3(on)
		with tempfile.TemporaryDirectory() as td:
			asp = _run(os.path.join(td, 'asp'), 150, cox1_factor=0.0)
			tp_asp = _tx(asp, 'tp_active_frac')
			ip3_asp = ip3(asp)

		# TXA₂ activates TP only when COX-1 is intact.
		assert tp_on.max() > 0.3
		assert tp_asp.max() == 0.0
		# Resting state is unperturbed by the loop (TP inactive at rest).
		assert ip3_on[0] == ip3_asp[0]
		# The closed loop adds Gq drive → IP₃ is ≥ the aspirin (open-loop) case,
		# and strictly greater once TXA₂ has accumulated.
		assert np.all(ip3_on >= ip3_asp - 1e-9)
		assert ip3_on[-1] > ip3_asp[-1]
