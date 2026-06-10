"""Smoke tests for the PMCA / MCU perturbation runner (issue #53).

Mechanics only — that each experiment runs, restores the knob it mutates, and
writes its artefacts. The scientific results (monotonic recovery-tail AUC for
PMCA; monotonic cytosolic buffering for MCU) are length-sensitive and are
verified by the full-length reference run, not here.
"""

import json
import os

import pytest

from reconstruction.platelet.dataclasses.process import calcium_signalling as cs_mod
from runscripts.manual.runPerturbation import (
	EXPERIMENTS,
	run_perturbation,
	write_outputs,
	write_summary,
)


@pytest.mark.slow
class TestPerturbationSmoke:
	def test_pmca_runs_restores_knob_and_writes_outputs(self, tmp_path):
		out = str(tmp_path / 'pert')
		os.makedirs(out)
		baseline = cs_mod.K_PMCA['k_cat']

		scan = run_perturbation(out, 'pmca', length_override=15,
			log_to_shell=False)

		# The mutated knob must be restored to its baseline afterwards.
		assert cs_mod.K_PMCA['k_cat'] == baseline
		# One trace row per factor; cyt and DTS traces aligned in length.
		assert scan.cyt.shape[0] == len(EXPERIMENTS['pmca']['factors'])
		assert scan.cyt.shape[1] == scan.dts.shape[1] > 10
		assert all('recovery_tail_auc_nMs' in s for s in scan.scalars)

		write_outputs(scan, out)
		write_summary([scan], out)
		for name in ('pmca.npz', 'pmca_traces.png', 'perturbation_summary.json'):
			assert os.path.exists(os.path.join(out, name)), name
		with open(os.path.join(out, 'perturbation_summary.json')) as f:
			payload = json.load(f)
		assert 'pmca' in payload['experiments']

	def test_mcu_runs_and_restores_knob(self, tmp_path):
		out = str(tmp_path / 'pert_mcu')
		os.makedirs(out)
		baseline = cs_mod.K_MITO['V_max_MCU']

		scan = run_perturbation(out, 'mcu', length_override=15,
			log_to_shell=False)

		assert cs_mod.K_MITO['V_max_MCU'] == baseline
		# MCU knockout (factor 0) is the first condition.
		assert scan.factors[0] == 0.0
		write_outputs(scan, out)
		assert os.path.exists(os.path.join(out, 'mcu_traces.png'))
