"""Tests for the TUI form schema / RunConfig assembly (wholecell/tui/runspec)."""

import pytest

from reconstruction.platelet.run_config import RunConfig
from wholecell.tui import runspec


def _default_values():
	return {key: runspec.field_default(key) for key in runspec.ALL_INPUT_KEYS}


def test_schema_keys_match_run_config_fields():
	# every form knob (besides run params) is a real RunConfig field; `live` is
	# forced by the runner and `count_overrides` is driven by the knockout
	# entity map rather than a numeric input.
	fields = set(RunConfig().to_metadata())
	assert set(runspec.RUN_CONFIG_KEYS) | {'live', 'count_overrides'} == fields


def test_build_spec_defaults_round_trip_into_run_config():
	spec = runspec.build_spec(_default_values(), knockouts={}, at_rest=False)
	assert spec['length_sec'] == 60 and spec['seed'] == 0
	rc = RunConfig(**spec['run_config'])  # accepts every key, no extras
	assert rc.thrombin_peak_nM == 1.0 and rc.ca_ex_mM == 1.2
	assert rc.cox1_factor == 1.0 and rc.live is True


def test_build_spec_knockout_forces_zero():
	spec = runspec.build_spec(
		_default_values(), knockouts={'cox1_factor': True}, at_rest=False)
	assert spec['run_config']['cox1_factor'] == 0.0


def test_build_spec_entity_knockouts_become_count_overrides():
	spec = runspec.build_spec(
		_default_values(), knockouts={}, at_rest=False,
		entity_knockouts=['PAR1 (thrombin)'])
	co = spec['run_config']['count_overrides']
	assert co == {
		'PAR1_inactive[pl]': 0, 'PAR1_active[pl]': 0,
		'PAR1_internalized[pl]': 0}


def test_build_spec_no_entity_knockouts_empty_count_overrides():
	spec = runspec.build_spec(_default_values(), knockouts={}, at_rest=False)
	assert spec['run_config']['count_overrides'] == {}


def test_build_spec_at_rest_zeroes_agonists():
	spec = runspec.build_spec(_default_values(), knockouts={}, at_rest=True)
	for key in runspec.AGONIST_KEYS:
		assert spec['run_config'][key] == 0.0


def test_build_spec_rejects_bad_fields():
	bad = _default_values()
	bad['ca_ex_mM'] = ''
	with pytest.raises(ValueError, match='Ca2'):
		runspec.build_spec(bad, knockouts={}, at_rest=False)

	bad = _default_values()
	bad['adp_peak_uM'] = '-1'
	with pytest.raises(ValueError, match='ADP'):
		runspec.build_spec(bad, knockouts={}, at_rest=False)

	bad = _default_values()
	bad['seed'] = 'x'
	with pytest.raises(ValueError, match='Seed'):
		runspec.build_spec(bad, knockouts={}, at_rest=False)


def test_diff_from_defaults_lists_changed_knobs():
	values = _default_values()
	knockouts = {k: False for k in runspec.KNOCKOUT_KEYS}
	assert runspec.diff_from_defaults(values, knockouts, at_rest=False) == []

	values['adp_peak_uM'] = '0.5'
	knockouts['cox1_factor'] = True
	out = runspec.diff_from_defaults(values, knockouts, at_rest=True)
	assert 'ADP (uM)' in out
	assert any('COX-1' in s for s in out)
	assert 'at rest' in out


def test_write_spec_round_trips(tmp_path):
	spec = runspec.build_spec(_default_values(), knockouts={}, at_rest=False)
	path = tmp_path / 'run.json'
	runspec.write_spec(spec, str(path))
	import json
	assert json.loads(path.read_text()) == spec
