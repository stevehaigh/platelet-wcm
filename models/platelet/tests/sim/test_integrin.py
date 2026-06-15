"""Tests for IntegrinActivation + IntegrinTrace (v0.61 §3).

Minimal 2-state αIIbβ3 inside-out activation: a resting ⇌ active conformational
switch gated by PKC × Ca²⁺. The active fraction is the per-cell PAC-1 readout.
Covers:
  * resting quiescence — no activation when un-stimulated;
  * activation — PKC + Ca²⁺ drive a high active (PAC-1⁺) fraction;
  * conservation — the conformational switch is mass-neutral (count conserved);
  * Glanzmann / antagonist — integrin_act_scale = 0 abolishes activation.
"""

import os
import tempfile

import numpy as np
import pytest

from wholecell.io.tablereader import TableReader
from runscripts.manual.runPlateletSim import run_platelet_sim
from reconstruction.platelet.run_config import RunConfig

_TOTAL_INTEGRIN = 80000   # aIIbb3_resting[pl] initial count (species TSV)


def _run(out_dir, length, **config_kwargs):
	"""Run one sim; config_kwargs are RunConfig fields (peaks, scales, …)."""
	run_config = RunConfig(ca_ex_mM=1.2, **config_kwargs)
	paths = run_platelet_sim(out_dir, length_sec=length, seed=0,
		log_to_shell=False, run_config=run_config)
	return paths['sim_out_dir']


def _ig(sim_out_dir, col):
	reader = TableReader(os.path.join(sim_out_dir, 'IntegrinTrace'))
	return reader.readColumn(col).flatten()


@pytest.mark.slow
class TestIntegrinActivation:

	def test_resting_quiescence_no_activation(self):
		"""Un-stimulated platelet: αIIbβ3 stays in the resting conformation."""
		with tempfile.TemporaryDirectory() as td:
			sim_out = _run(td, 30, thrombin_peak_nM=0.0,
				adp_peak_uM=0.0, atp_ex_peak_uM=0.0)
			assert _ig(sim_out, 'integrin_gate').max() == 0.0
			assert _ig(sim_out, 'active_frac').max() == 0.0
			assert _ig(sim_out, 'aIIbb3_active').max() == 0

	def test_activation_drives_pac1(self):
		"""Standard +Ca²⁺ transient drives αIIbβ3 to the high-affinity state."""
		with tempfile.TemporaryDirectory() as td:
			sim_out = _run(td, 100)
			gate = _ig(sim_out, 'integrin_gate')
			frac = _ig(sim_out, 'active_frac')
			active = _ig(sim_out, 'aIIbb3_active')
			assert gate[0] == 0.0           # zero until PKC activates
			assert frac[0] == 0.0
			assert frac.max() > 0.3         # a clear PAC-1⁺ fraction
			assert frac.max() < 1.0         # reversible — not fully saturated
			# Activation is monotonic non-decreasing while the agonist drives it.
			assert active[-1] > active[0]

	def test_total_integrin_conserved(self):
		"""resting ⇌ active is a conformational switch — count is conserved."""
		with tempfile.TemporaryDirectory() as td:
			sim_out = _run(td, 100)
			total = (_ig(sim_out, 'aIIbb3_active')
				+ _ig(sim_out, 'aIIbb3_resting'))
			assert np.all(total == _TOTAL_INTEGRIN)

	def test_glanzmann_knockout_abolishes_activation(self):
		"""integrin_act_scale = 0 (αIIbβ3 antagonist / Glanzmann) → no active."""
		with tempfile.TemporaryDirectory() as td:
			sim_out = _run(td, 100, integrin_act_scale=0.0)
			assert _ig(sim_out, 'active_frac').max() == 0.0
			assert _ig(sim_out, 'aIIbb3_active').max() == 0
