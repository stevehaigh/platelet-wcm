"""Tests for the platelet TUI experiment bench app (phases P0–P2).

These drive the Textual app, so they need `textual` / `textual-plotext`. CI's
minimal requirements omit front-end deps (as for the Dash webapp), so the whole
module skips there; the textual-free TUI tests (runspec / presets / runFromConfig)
still run in CI.
"""

import asyncio
import json
import os

import pytest

pytest.importorskip('textual')
pytest.importorskip('textual_plotext')

from wholecell.tui import app as tui


# ── pure helpers (fast) ────────────────────────────────────────────────────

def test_resolve_sim_path_relative_goes_under_out():
	assert tui.resolve_sim_path('myrun', '/repo') == '/repo/out/myrun'
	assert tui.resolve_sim_path('out/myrun', '/repo') == '/repo/out/myrun'
	assert tui.resolve_sim_path('/abs/run', '/repo') == '/abs/run'


def test_live_csv_path_mirrors_run_nesting():
	path = tui.live_csv_path('tui_run', 0, '/repo')
	assert path == (
		'/repo/out/tui_run/platelet_stub_000000/000000/'
		'generation_000000/000000/simOut/live.csv')


def test_read_live_csv_skips_partial_final_row(tmp_path):
	csv_path = tmp_path / 'live.csv'
	csv_path.write_text(
		'time,ca_cyt_nM,ca_dts_uM,ip3_nM,soce_flux_nMs\n'
		'0,100.0,250.0,50.0,0.0\n'
		'1,180.5,248.0,120.0,1.2\n'
		'2,partial')  # truncated row written mid-timestep
	cols = tui.read_live_csv(str(csv_path))
	assert cols['time'] == [0.0, 1.0]
	assert cols['ca_cyt_nM'] == [100.0, 180.5]


def test_read_live_csv_missing_file_is_empty():
	cols = tui.read_live_csv('/no/such/live.csv')
	assert cols['time'] == []


# ── headless app (fast: no simulation) ─────────────────────────────────────

def test_app_builds_with_grouped_form_and_knockouts():
	from textual.widgets import Input, Switch, Checkbox
	from wholecell.tui import runspec

	async def scenario():
		app = tui.PlateletBenchApp()
		async with app.run_test(size=(120, 50)) as pilot:
			await pilot.pause()
			# results panes + a numeric input for every RunConfig/run field
			assert app.query_one('#run-btn')
			assert app.query_one('#plot-cyt') and app.query_one('#plot-dts')
			for key in runspec.ALL_INPUT_KEYS:
				assert app.query_one('#in-' + key, Input)
			# "run at rest" disables the three agonist inputs
			app.query_one('#sw-at-rest', Switch).value = True
			await pilot.pause()
			for key in runspec.AGONIST_KEYS:
				assert app.query_one('#in-' + key, Input).disabled
			# a KO checkbox disables (greys out) its numeric input
			app.query_one('#ko-cox1_factor', Checkbox).value = True
			await pilot.pause()
			assert app.query_one('#in-cox1_factor', Input).disabled
			# an expression-knockout checkbox feeds _entity_knockouts()
			assert app._entity_knockouts() == []
			app.query_one('#ent-1', Checkbox).value = True  # PAR1
			await pilot.pause()
			assert len(app._entity_knockouts()) == 1
			assert 'KO' in app._defaults_diff_text

	asyncio.run(scenario())


def test_baseline_overlay_set_guard_draw_and_clear():
	async def scenario():
		app = tui.PlateletBenchApp()
		async with app.run_test(size=(120, 50)) as pilot:
			await pilot.pause()
			# set-baseline before any run is a guarded no-op
			app.action_set_baseline()
			await pilot.pause()
			assert app._baseline is None
			# a synthetic baseline redraws both panels without error
			app._baseline = {
				'time': [0.0, 1.0, 2.0], 'ca_cyt_nM': [100.0, 150.0, 200.0],
				'ca_dts_uM': [250.0, 240.0, 230.0], 'ip3_nM': [50.0, 60.0, 70.0],
				'soce_flux_nMs': [0.0, 1.0, 1.0]}
			app._reset_plots()
			await pilot.pause()
			app.action_clear_baseline()
			await pilot.pause()
			assert app._baseline is None

	asyncio.run(scenario())


def test_modified_indicators_track_edits_and_presets(tmp_path, monkeypatch):
	monkeypatch.setenv('PLATELET_TUI_PRESETS_DIR', str(tmp_path))
	from textual.widgets import Input

	async def scenario():
		app = tui.PlateletBenchApp()
		async with app.run_test(size=(120, 50)) as pilot:
			await pilot.pause()
			# fresh app matches Defaults, nothing changed
			assert app._preset_status_text == '● Defaults'
			assert app._defaults_diff_text == 'Δ from defaults: none'
			# editing a field flags the stale preset + the defaults diff
			app.query_one('#in-adp_peak_uM', Input).value = '0.5'
			await pilot.pause()
			assert 'edited' in app._preset_status_text
			assert 'ADP' in app._defaults_diff_text
			# applying a preset re-syncs the marker (no longer edited)
			app._apply_preset('EDTA (no Ca_ex)')
			await pilot.pause()
			assert app._preset_status_text == '● EDTA (no Ca_ex)'
			# EDTA differs from defaults (Ca_ex + length), so the diff persists
			assert app._defaults_diff_text != 'Δ from defaults: none'

	asyncio.run(scenario())


def test_figure_action_guarded_without_a_run():
	async def scenario():
		app = tui.PlateletBenchApp()
		async with app.run_test(size=(120, 50)) as pilot:
			await pilot.pause()
			app.action_figure()  # nothing has run yet
			await pilot.pause()
			assert app._figure_ok is None
			assert 'Run a simulation first' in app._last_status

	asyncio.run(scenario())


@pytest.mark.slow
def test_figure_handoff_renders_pdf(monkeypatch):
	monkeypatch.setenv('PLATELET_TUI_NO_OPEN', '1')  # don't launch a viewer
	from textual.widgets import Input

	async def scenario():
		app = tui.PlateletBenchApp()
		async with app.run_test(size=(120, 50)) as pilot:
			await pilot.pause()
			app.query_one('#in-length_sec', Input).value = '5'
			await pilot.pause()
			await pilot.press('r')
			for _ in range(80):
				await pilot.pause(0.25)
				if app._sim_ok is not None:
					break
			assert app._sim_ok is True
			app.action_figure()
			for _ in range(80):
				await pilot.pause(0.25)
				if app._figure_ok is not None:
					break
			assert app._figure_ok is True
			root = tui.fp.ROOT_PATH
			sim_out = os.path.dirname(tui.live_csv_path('tui_run', 0, root))
			pdf = os.path.join(
				sim_out.replace('simOut', 'plotOut'), 'calcium_trace.pdf')
			assert os.path.exists(pdf)

	asyncio.run(scenario())


def test_presets_apply_to_form_and_save_round_trips(tmp_path, monkeypatch):
	monkeypatch.setenv('PLATELET_TUI_PRESETS_DIR', str(tmp_path))
	from textual.widgets import Input, Switch, Select
	from wholecell.tui import presets

	async def scenario():
		app = tui.PlateletBenchApp()
		async with app.run_test(size=(120, 50)) as pilot:
			await pilot.pause()
			# applying the EDTA preset zeroes Ca_ex and stretches the run
			app._apply_preset('EDTA (no Ca_ex)')
			await pilot.pause()
			assert app.query_one('#in-ca_ex_mM', Input).value == '0.0'
			assert app.query_one('#in-length_sec', Input).value == '200'
			# applying Resting flips the at-rest switch on
			app._apply_preset('Resting')
			await pilot.pause()
			assert app.query_one('#sw-at-rest', Switch).value is True
			# saving the current form adds a user preset to the dropdown
			app.query_one('#preset-name', Input).value = 'MyPreset'
			app._save_preset()
			await pilot.pause()
			assert 'MyPreset' in presets.all_names()
			assert app.query_one('#preset-select', Select).value == 'MyPreset'

	asyncio.run(scenario())


# ── end-to-end: a short sim driven through the TUI worker (slow) ───────────

@pytest.mark.slow
def test_short_run_through_tui_streams_live_csv():
	from textual.widgets import Input, Checkbox

	async def scenario():
		app = tui.PlateletBenchApp()
		async with app.run_test(size=(120, 50)) as pilot:
			await pilot.pause()
			app.query_one('#in-length_sec', Input).value = '5'
			app.query_one('#ko-cox1_factor', Checkbox).value = True  # aspirin
			app.query_one('#ent-1', Checkbox).value = True  # PAR1 expression KO
			await pilot.pause()
			# 'r' is a priority binding, so it fires even with a field focused
			await pilot.press('r')
			# a length=5 run takes ~1 s; poll for the worker's outcome (bounded)
			for _ in range(80):
				await pilot.pause(0.25)
				if app._sim_ok is not None:
					break
			assert app._sim_ok is True
			root = tui.fp.ROOT_PATH
			# the sim wrote live.csv where the TUI computed it
			cols = tui.read_live_csv(tui.live_csv_path('tui_run', 0, root))
			assert len(cols['time']) >= 1
			# the aspirin knockout flowed through to the written run spec
			spec_path = os.path.join(
				tui.resolve_sim_path('tui_run', root), 'run_config.json')
			with open(spec_path) as handle:
				spec = json.load(handle)
			assert spec['run_config']['cox1_factor'] == 0.0
			# the PAR1 expression knockout reached count_overrides
			assert spec['run_config']['count_overrides']['PAR1_active[pl]'] == 0

	asyncio.run(scenario())
