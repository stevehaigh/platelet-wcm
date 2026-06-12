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

	def test_pkc_runs_restores_knob_and_harvests_aux(self, tmp_path):
		"""v0.6: PKC P2Y1-desensitisation knockout scan + aux trace + figure."""
		out = str(tmp_path / 'pert_pkc')
		os.makedirs(out)
		baseline = cs_mod.K_P2Y1_DES['k_des']

		scan = run_perturbation(out, 'pkc', length_override=20,
			log_to_shell=False)

		# Knob restored; knockout (k_des×0) is the first condition.
		assert cs_mod.K_P2Y1_DES['k_des'] == baseline
		assert scan.factors[0] == 0.0
		# The desensitised-fraction aux trace was harvested, aligned with cyt.
		assert scan.aux is not None
		frac = scan.aux['p2y1_desensitised_frac']
		assert frac.shape == scan.cyt.shape
		# Knockout must show no desensitisation; baseline must show some.
		assert frac[0].max() == 0.0
		assert frac[1].max() > 0.0
		assert all('p2y1_desensitised_frac_max' in s for s in scan.scalars)

		write_outputs(scan, out)
		assert os.path.exists(os.path.join(out, 'pkc_traces.png'))
		assert os.path.exists(os.path.join(out, 'pkc.npz'))

	def test_plcb_runs_restores_knob_and_harvests_dual_aux(self, tmp_path):
		"""v0.6 Slice 3: PKC→PLCβ-phosphorylation knockout + IP₃/phos aux."""
		out = str(tmp_path / 'pert_plcb')
		os.makedirs(out)
		baseline = cs_mod.K_PLCB_PHOS['k_plcb_phos']

		scan = run_perturbation(out, 'plcb', length_override=20,
			log_to_shell=False)

		assert cs_mod.K_PLCB_PHOS['k_plcb_phos'] == baseline
		assert scan.factors[0] == 0.0
		# Two aux columns harvested (IP₃ + phospho fraction), both aligned.
		assert scan.aux is not None
		assert set(scan.aux) == {'ip3_nM', 'plcb_phosphorylated_frac'}
		phos = scan.aux['plcb_phosphorylated_frac']
		assert phos.shape == scan.cyt.shape
		# Knockout: no PLCβ phosphorylation; baseline: some.
		assert phos[0].max() == 0.0
		assert phos[1].max() > 0.0

		write_outputs(scan, out)
		assert os.path.exists(os.path.join(out, 'plcb_traces.png'))
		assert os.path.exists(os.path.join(out, 'plcb.npz'))
