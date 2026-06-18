"""Tests for the JSON-config simulation runner (runscripts/manual/runFromConfig)."""

import json

from runscripts.manual import runFromConfig


def test_load_spec_builds_run_config_and_forces_live(tmp_path):
	spec = {
		'length_sec': 30,
		'seed': 7,
		'run_config': {'ca_ex_mM': 0.0, 'cox1_factor': 0.0, 'adp_peak_uM': 0.5},
	}
	path = tmp_path / 'run.json'
	path.write_text(json.dumps(spec))

	length_sec, seed, run_config = runFromConfig.load_spec(str(path))
	assert length_sec == 30 and seed == 7
	assert run_config.ca_ex_mM == 0.0
	assert run_config.cox1_factor == 0.0
	assert run_config.adp_peak_uM == 0.5
	assert run_config.live is True  # the on-demand runner always streams live
	# unspecified fields keep their RunConfig defaults
	assert run_config.pmca_kcat_scale == 1.0


def test_load_spec_defaults_when_minimal(tmp_path):
	path = tmp_path / 'min.json'
	path.write_text(json.dumps({'run_config': {}}))
	length_sec, seed, run_config = runFromConfig.load_spec(str(path))
	assert length_sec == 60 and seed == 0
	assert run_config.live is True
