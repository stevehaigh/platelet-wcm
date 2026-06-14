"""Smoke tests for the v0.61 §3 / recovery-phase figure scripts.

These reproduce the headline figures end-to-end (runscripts/manual/
plotIntegrin.py and plotRecoveryPhase.py): each must run, write its artefacts,
and show the core directional result the figure exists to make.
"""

import os

import pytest

from runscripts.manual.plotIntegrin import run_integrin_figure
from runscripts.manual.plotRecoveryPhase import run_recovery_phase
from runscripts.manual.runDoseResponse import run_dose_response, write_outputs


@pytest.mark.slow
class TestIntegrinFigure:
	def test_runs_and_shows_graded_activation(self, tmp_path):
		out = str(tmp_path / 'integrin')
		os.makedirs(out)
		results = run_integrin_figure(out, length=80, log_to_shell=False)

		# Artefacts written.
		for name in ('integrin_activation.png', 'integrin_activation.npz'):
			assert os.path.exists(os.path.join(out, name)), name

		# PAC-1 activation: the standard transient activates αIIbβ3, while the
		# Glanzmann / antagonist knockout and the resting cell stay at zero.
		assert results['standard']['active_frac'].max() > 0.2
		assert results['glanzmann']['active_frac'].max() == 0.0
		assert results['rest']['active_frac'].max() == 0.0
		# Graded by agonist strength: standard > thrombin-only.
		assert (results['standard']['active_frac'].max()
			> results['thrombin']['active_frac'].max())


@pytest.mark.slow
class TestRecoveryPhaseFigure:
	def test_runs_and_shows_loop_sustained_vs_recovery(self, tmp_path):
		out = str(tmp_path / 'recovery')
		os.makedirs(out)
		results = run_recovery_phase(out, length=250, log_to_shell=False)

		# Artefacts written.
		for name in ('recovery_phase_traces.png', 'recovery_phase.npz'):
			assert os.path.exists(os.path.join(out, name)), name

		# The full v0.61 model sustains IP3 (the autocrine second wave); the
		# Dolan-equivalent (loops off + transient stimulus) recovers it toward
		# the 50 nM baseline.
		full_ip3_end = results['full_ca']['ip3_nM'][-1]
		dolan_ip3_end = results['dolan_ca']['ip3_nM'][-1]
		assert full_ip3_end > 150.0
		assert dolan_ip3_end < full_ip3_end - 50.0


@pytest.mark.slow
class TestDoseResponse:
	def test_runs_and_shows_graded_refill_vs_commitment(self, tmp_path):
		out = str(tmp_path / 'dr')
		os.makedirs(out)
		payload = run_dose_response(out, agonist='adp', grid=3, length=300,
			log_to_shell=False)
		write_outputs(payload, out)

		# Artefacts written.
		for name in ('dose_response_adp.png', 'dose_response_adp.npz',
				'dose_response_adp.json'):
			assert os.path.exists(os.path.join(out, name)), name

		off = payload['results']['loops_off']['dts_end']
		full = payload['results']['full']['dts_end']
		# Loops off: the sustained store level grades with dose — it refills far
		# more at the lowest dose than the highest (graded, reversible).
		assert off[0] > off[-1] + 20.0
		# Full v0.61: the autocrine loops keep the store empty (committed) even at
		# the highest dose — the commitment switch.
		assert full[-1] < 30.0
