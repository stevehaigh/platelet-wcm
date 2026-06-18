"""Tests for TUI presets (wholecell/tui/presets)."""

from wholecell.tui import presets, runspec


def test_resolve_builtin_defaults():
	values, knockouts, at_rest = presets.resolve('Defaults')
	assert values['ca_ex_mM'] == '1.2'
	assert values['thrombin_peak_nM'] == '1.0'
	assert at_rest is False
	assert not any(knockouts.values())


def test_resolve_edta_zeroes_ca_and_sets_length():
	values, _ko, _ar = presets.resolve('EDTA (no Ca_ex)')
	assert values['ca_ex_mM'] == '0.0'
	assert values['length_sec'] == '200'


def test_resolve_resting_sets_at_rest():
	_v, _ko, at_rest = presets.resolve('Resting')
	assert at_rest is True


def test_resolve_aspirin_knocks_out_cox1():
	_v, knockouts, _ar = presets.resolve('Aspirin (COX-1 KO)')
	assert knockouts['cox1_factor'] is True
	assert not knockouts['mcu_vmax_scale']


def test_safe_name_sanitises_path_chars():
	assert presets._safe_name('a/b.c') == 'a_b_c'
	assert presets._safe_name('   ') == 'preset'


def test_save_and_reload_user_preset(tmp_path, monkeypatch):
	monkeypatch.setenv('PLATELET_TUI_PRESETS_DIR', str(tmp_path))
	values = {k: runspec.field_default(k) for k in runspec.ALL_INPUT_KEYS}
	values['adp_peak_uM'] = '0.5'
	knockouts = {k: False for k in runspec.KNOCKOUT_KEYS}
	knockouts['mcu_vmax_scale'] = True

	saved = presets.save('My Run', values, knockouts, at_rest=True)
	assert saved == 'My Run'
	assert saved in presets.all_names()

	v2, k2, ar2 = presets.resolve(saved)
	assert v2['adp_peak_uM'] == '0.5'
	assert k2['mcu_vmax_scale'] is True
	assert ar2 is True
