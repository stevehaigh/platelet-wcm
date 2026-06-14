"""Smoke test for the second-wave experiment (runSecondWave.py).

Mechanics + the core directional result: at a weak transient ADP stimulus the
closed-loop (v0.61) model sustains cytosolic Ca2+ above the open-loop (v0.6)
model — the autocrine second wave. Both loop knobs must be restored afterwards.
"""

import json
import os

import pytest

from runscripts.manual.runSecondWave import run_second_wave


@pytest.mark.slow
class TestSecondWave:
	def test_runs_restores_knobs_and_shows_second_wave(self, tmp_path):
		out = str(tmp_path / 'sw')
		os.makedirs(out)

		# The recovery-phase divergence establishes over ~200-300 s (the open
		# loop decays while the closed loop stays activated), so the test needs
		# a long enough window to see the second wave.
		results, scalars = run_second_wave(out, adp_uM=0.5, length=250,
			log_to_shell=False)

		# Same store-limited peak; the loops act in the recovery phase.
		assert abs(scalars['full']['peak_cyt_nM']
			- scalars['v06']['peak_cyt_nM']) < 5.0
		# The second wave: closed loop sustains Ca2+ well above the open loop.
		assert scalars['full']['end_cyt_nM'] > scalars['v06']['end_cyt_nM'] + 10.0
		# The sustained Ca2+ arm is the autocrine ADP loop (aspirin ~ full),
		# while TXA2 adds IP3 drive (full IP3 >= aspirin IP3).
		assert scalars['aspirin']['end_cyt_nM'] > scalars['v06']['end_cyt_nM']
		assert scalars['full']['ip3_end_nM'] >= scalars['aspirin']['ip3_end_nM']

		# Artefacts written.
		for name in ('second_wave_traces.png', 'second_wave.npz',
				'second_wave_summary.json'):
			assert os.path.exists(os.path.join(out, name)), name
		with open(os.path.join(out, 'second_wave_summary.json')) as f:
			payload = json.load(f)
		assert payload['sustained_cyt_gap_full_minus_v06_nM'] > 0
