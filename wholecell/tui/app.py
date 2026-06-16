"""
Platelet WCM — terminal experiment bench (Textual TUI), phase P1.

Builds on the P0 skeleton: the full `RunConfig` is now surfaced, grouped by
subsystem (Stimulus / Feedback loops / Pumps & brakes), with inline "KO"
checkboxes for the curated knockouts (each forces a knob to 0). Editing the
form builds a JSON run spec that is run on demand by
`runscripts/manual/runFromConfig.py` as a subprocess (decision §12), with the
cytosolic / DTS Ca²⁺ trace streamed live from the run's `live.csv`.

Schema and config assembly live in `wholecell/tui/runspec.py`. Still to come
(P1): presets (save/load) and compare-to-baseline. See
`reports/design/tui-tinkering-dashboard-2026-06-15.qmd`.
"""

from __future__ import annotations

import csv
import os
import subprocess
import sys
from typing import Dict, List, Optional

from textual import work
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical, VerticalScroll
from textual.widgets import (
	Button, Checkbox, Collapsible, Footer, Header, Input, Label, Rule, Select,
	Static, Switch)
from textual_plotext import PlotextPlot

from wholecell.tui import presets, runspec
from wholecell.utils import filepath as fp

# Columns the CalciumTrace listener writes to live.csv (one row per timestep).
_LIVE_COLUMNS = ('time', 'ca_cyt_nM', 'ca_dts_uM', 'ip3_nM', 'soce_flux_nMs')

# Fixed output directory for TUI runs (overwritten each run).
_OUTDIR = 'tui_run'

# Plot panels: (widget id, y-column, axis label, line colour).
_PANELS = (
	('plot-cyt', 'ca_cyt_nM', 'Ca2+_cyt (nM)', 'blue'),
	('plot-dts', 'ca_dts_uM', 'Ca2+_dts (uM)', 'green'),
)


def resolve_sim_path(outdir: str, root: str) -> str:
	"""Mirror runPlateletSim.resolve_sim_path: relative dirs live under out/."""
	if os.path.isabs(outdir):
		return outdir
	if outdir.startswith('out/'):
		return os.path.join(root, outdir)
	return os.path.join(root, 'out', outdir)


def live_csv_path(outdir: str, seed: int, root: str) -> str:
	"""Path to a run's live.csv (mirrors run_platelet_sim's output nesting)."""
	sim_path = resolve_sim_path(outdir, root)
	return os.path.join(
		sim_path, f'platelet_stub_{seed:06d}', f'{seed:06d}',
		'generation_000000', '000000', 'simOut', 'live.csv')


def read_live_csv(path: str) -> Dict[str, List[float]]:
	"""Read live.csv into column lists, tolerating a partial final row."""
	cols: Dict[str, List[float]] = {k: [] for k in _LIVE_COLUMNS}
	try:
		with open(path, newline='') as handle:
			for row in csv.DictReader(handle):
				try:
					vals = [float(row[k]) for k in _LIVE_COLUMNS]
				except (ValueError, KeyError, TypeError):
					continue  # last row may be written mid-timestep
				for key, val in zip(_LIVE_COLUMNS, vals):
					cols[key].append(val)
	except FileNotFoundError:
		pass
	return cols


class PlateletBenchApp(App):
	"""Edit the full RunConfig, knock out pathways, run, watch Ca2+ live."""

	CSS = """
	#body { height: 1fr; }
	#form {
		width: 42;
		border: round $primary;
		padding: 0 1;
	}
	#form Label { color: $text-muted; }
	#preset-status { height: 1; color: $text-muted; }
	#defaults-diff { height: 1; color: $text-muted; margin: 1 0 0 0; }
	.edited { color: $warning; }
	.heading { text-style: bold; color: $text; margin: 1 0 0 0; }
	.field-row { height: 1; }
	.field-label { width: 18; height: 1; content-align: left middle; }
	.field-row Input {
		border: none;
		height: 1;
		padding: 0 1;
		width: 1fr;
		background: $boost;
	}
	.field-row Checkbox { border: none; height: 1; width: auto; margin: 0 0 0 1; }
	.save-row { height: 3; margin: 1 0 0 0; }
	.save-row Input { width: 1fr; }
	.save-row Button { width: auto; min-width: 9; margin: 0 0 0 1; }
	#rest-row { height: 1; }
	#rest-row Label { width: 18; height: 1; content-align: left middle; }
	#run-btn { margin: 1 0 1 0; width: 100%; }
	#results { width: 1fr; }
	PlotextPlot { height: 1fr; border: round $accent; }
	#baseline-row { height: 3; }
	#baseline-row Button { width: auto; min-width: 11; margin: 0 1 0 0; }
	#status {
		height: 3;
		padding: 0 1;
		content-align: left middle;
		color: $text-muted;
	}
	"""

	# priority=True so the keys fire even while a numeric field has focus
	# (the run-condition inputs never accept letters, so this steals nothing).
	BINDINGS = [
		Binding('r', 'run', 'Run', priority=True),
		Binding('b', 'set_baseline', 'Baseline', priority=True),
		Binding('f', 'figure', 'Figure', priority=True),
		Binding('q', 'quit', 'Quit', priority=True),
	]

	def __init__(self) -> None:
		super().__init__()
		self._sim_running = False
		self._sim_timer = None
		self._sim_live_path: Optional[str] = None
		self._sim_length = 0
		self._sim_ok: Optional[bool] = None  # last run outcome (for tests)
		self._baseline: Optional[Dict[str, List[float]]] = None  # overlay trace
		self._figure_ok: Optional[bool] = None  # last figure render (for tests)
		self._last_status = ''  # mirror of the status line (for tests)
		self._preset_status_text = ''  # mirror of the stale-preset marker
		self._defaults_diff_text = ''  # mirror of the modified-from-defaults line
		self._ui_ready = False  # True once mounted (guards early change events)
		# form state of the last applied preset, for the stale-preset marker
		self._applied_state = presets.resolve('Defaults')
		self._applied_preset_name = 'Defaults'

	def compose(self) -> ComposeResult:
		yield Header()
		with Horizontal(id='body'):
			with VerticalScroll(id='form'):
				yield Label('Preset', classes='heading')
				yield Select(
					[(n, n) for n in presets.all_names()],
					id='preset-select', value='Defaults',
					allow_blank=False, compact=True)
				yield Static('● Defaults', id='preset-status')
				yield Horizontal(
					Input(placeholder='save as…', id='preset-name'),
					Button('Save', id='save-btn'),
					classes='save-row')
				yield Rule()
				yield Label('RUN', classes='heading')
				for key, label, default, kind in runspec.RUN_PARAMS:
					yield from self._field_row(key, label, default, kind, False)
				for title, fields in runspec.GROUPS:
					with Collapsible(
							title=title, collapsed=(title != 'Stimulus')):
						if title == 'Stimulus':
							yield Horizontal(
								Label('Run at rest'), Switch(id='sw-at-rest'),
								id='rest-row')
						for key, label, default, kind, ko in fields:
							yield from self._field_row(
								key, label, default, kind, ko)
				with Collapsible(
						title='Knockouts (remove protein)', collapsed=True):
					for i, name in enumerate(runspec.KNOCKOUT_ENTITIES):
						yield Checkbox(name, id='ent-' + str(i))
				yield Static('Δ from defaults: none', id='defaults-diff')
				yield Button('Run ▶', id='run-btn', variant='success')
			with Vertical(id='results'):
				yield PlotextPlot(id='plot-cyt')
				yield PlotextPlot(id='plot-dts')
				yield Horizontal(
					Button('Set baseline', id='baseline-btn'),
					Button('Clear', id='clear-baseline-btn'),
					Button('Figure (5-panel)', id='figure-btn'),
					id='baseline-row')
				yield Static('Ready.', id='status')
		yield Footer()

	def _field_row(self, key: str, label: str, default: str, kind: str,
			ko: bool) -> ComposeResult:
		"""Yield one compact row: label + inline input (+ KO checkbox)."""
		itype = 'integer' if kind == 'int' else 'number'
		children = [
			Label(label, classes='field-label'),
			Input(value=default, id='in-' + key, type=itype),
		]
		if ko:
			children.append(Checkbox('KO', id='ko-' + key))
		yield Horizontal(*children, classes='field-row')

	def on_mount(self) -> None:
		self.title = 'Platelet WCM — experiment bench'
		self.sub_title = 'P2'
		self._reset_plots()
		self._ui_ready = True
		self._update_modified()

	# ── actions / events ────────────────────────────────────────────────────

	def action_run(self) -> None:
		self._start_run()

	def on_button_pressed(self, event: Button.Pressed) -> None:
		if event.button.id == 'run-btn':
			self._start_run()
		elif event.button.id == 'save-btn':
			self._save_preset()
		elif event.button.id == 'baseline-btn':
			self.action_set_baseline()
		elif event.button.id == 'clear-baseline-btn':
			self.action_clear_baseline()
		elif event.button.id == 'figure-btn':
			self.action_figure()

	def action_set_baseline(self) -> None:
		"""Snapshot the last run's trace as the overlay baseline."""
		cols = (read_live_csv(self._sim_live_path)
			if self._sim_live_path else {})
		if not cols.get('time'):
			self._set_status('⚠ Run a simulation first, then set it as baseline')
			return
		self._baseline = cols
		self._poll_or_redraw()
		self._set_status(
			f'Baseline set ({len(cols["time"])} steps) — new runs overlay it '
			f'in grey.')

	def action_clear_baseline(self) -> None:
		self._baseline = None
		self._poll_or_redraw()
		self._set_status('Baseline cleared.')

	def _poll_or_redraw(self) -> None:
		"""Redraw the plots (live data if a run path exists, else just baseline)."""
		if self._sim_live_path:
			self._poll()
		else:
			self._reset_plots()

	# ── full-figure handoff ─────────────────────────────────────────────────

	def action_figure(self) -> None:
		"""Render the 5-panel matplotlib calcium figure and open it externally."""
		if self._sim_running:
			self._set_status('⚠ Wait for the run to finish, then render the figure')
			return
		if not (self._sim_ok and self._sim_live_path):
			self._set_status('⚠ Run a simulation first, then render its figure')
			return
		root = fp.ROOT_PATH
		sim_out_dir = os.path.dirname(self._sim_live_path)
		fig_path = os.path.join(
			sim_out_dir.replace('simOut', 'plotOut'), 'calcium_trace.pdf')
		env = dict(os.environ)
		env['PYTHONPATH'] = root
		env['OPENBLAS_NUM_THREADS'] = '1'
		env['MPLBACKEND'] = 'Agg'  # render to file, no GUI backend
		self._figure_ok = None
		self._set_status('Rendering 5-panel figure…')
		self._figure_worker(root, env, fig_path)

	@work(thread=True, group='figure', exclusive=True)
	def _figure_worker(self, cwd: str, env: dict, fig_path: str) -> None:
		cmd = [
			sys.executable, 'runscripts/manual/analysisPlatelet.py', _OUTDIR,
			'--plot', 'calcium_trace']
		try:
			proc = subprocess.run(
				cmd, cwd=cwd, env=env, capture_output=True, text=True)
			ok = proc.returncode == 0 and os.path.exists(fig_path)
			err = '' if ok else (
				proc.stderr or proc.stdout or 'figure not produced').strip()
		except Exception as exc:  # pragma: no cover - defensive
			ok, err = False, str(exc)
		self.call_from_thread(self._figure_finished, ok, err, fig_path)

	def _figure_finished(self, ok: bool, err: str, fig_path: str) -> None:
		self._figure_ok = ok
		if not ok:
			tail = err.splitlines()[-1] if err else 'unknown error'
			self._set_status(f'✗ Figure failed: {tail}')
			return
		# PLATELET_TUI_NO_OPEN lets tests/headless runs skip launching a viewer
		if os.environ.get('PLATELET_TUI_NO_OPEN'):
			self._set_status(f'✓ Figure rendered: {fig_path}')
			return
		opener = {'darwin': 'open', 'win32': 'start'}.get(
			sys.platform, 'xdg-open')
		try:
			subprocess.Popen([opener, fig_path])
		except Exception as exc:  # pragma: no cover - defensive
			self._set_status(f'✓ Figure at {fig_path} (open failed: {exc})')
			return
		self._set_status(f'✓ Figure rendered & opened: calcium_trace.pdf')

	def on_select_changed(self, event: Select.Changed) -> None:
		if event.select.id == 'preset-select' and event.value is not Select.BLANK:
			self._apply_preset(str(event.value))

	def on_switch_changed(self, event: Switch.Changed) -> None:
		if event.switch.id == 'sw-at-rest':
			for key in runspec.AGONIST_KEYS:
				self.query_one('#in-' + key, Input).disabled = event.value
		self._update_modified()

	def on_checkbox_changed(self, event: Checkbox.Changed) -> None:
		cid = event.checkbox.id or ''
		if cid.startswith('ko-'):
			# a knocked-out knob is forced to 0, so grey out its numeric input
			self.query_one('#in-' + cid[3:], Input).disabled = event.value
		self._update_modified()

	def on_input_changed(self, event: Input.Changed) -> None:
		if event.input.id != 'preset-name':  # ignore the save-name field
			self._update_modified()

	def _update_modified(self) -> None:
		"""Refresh the stale-preset marker and the modified-from-defaults line."""
		if not self._ui_ready:
			return
		values, knockouts, at_rest = self._form_state()
		entity_kos = self._entity_knockouts()
		# presets carry no entity knockouts, so any ticked entity is an edit
		edited = ((values, knockouts, at_rest) != self._applied_state
			or bool(entity_kos))
		self._preset_status_text = (
			f'✎ edited (≠ {self._applied_preset_name})' if edited
			else f'● {self._applied_preset_name}')
		status = self.query_one('#preset-status', Static)
		status.set_class(edited, 'edited')
		status.update(self._preset_status_text)

		changed = runspec.diff_from_defaults(values, knockouts, at_rest)
		changed += [name + ' KO' for name in entity_kos]
		if not changed:
			self._defaults_diff_text = 'Δ from defaults: none'
		else:
			shown = ', '.join(changed[:3])
			more = f' +{len(changed) - 3}' if len(changed) > 3 else ''
			self._defaults_diff_text = (
				f'Δ from defaults ({len(changed)}): {shown}{more}')
		diff = self.query_one('#defaults-diff', Static)
		diff.set_class(bool(changed), 'edited')
		diff.update(self._defaults_diff_text)

	# ── presets ─────────────────────────────────────────────────────────────

	def _form_state(self):
		"""Current form state: (values, knockouts, at_rest)."""
		values = {key: self.query_one('#in-' + key, Input).value
			for key in runspec.ALL_INPUT_KEYS}
		knockouts = {key: self.query_one('#ko-' + key, Checkbox).value
			for key in runspec.KNOCKOUT_KEYS}
		at_rest = self.query_one('#sw-at-rest', Switch).value
		return values, knockouts, at_rest

	def _entity_knockouts(self) -> List[str]:
		"""Logical entities whose expression-knockout checkbox is ticked."""
		return [name for i, name in enumerate(runspec.KNOCKOUT_ENTITIES)
			if self.query_one('#ent-' + str(i), Checkbox).value]

	def _apply_preset(self, name: str) -> None:
		state = presets.resolve(name)
		values, knockouts, at_rest = state
		for key, val in values.items():
			self.query_one('#in-' + key, Input).value = val
		for key, on in knockouts.items():
			self.query_one('#ko-' + key, Checkbox).value = on
		self.query_one('#sw-at-rest', Switch).value = at_rest
		self._applied_state = state
		self._applied_preset_name = name
		self._update_modified()

	def _save_preset(self) -> None:
		name = self.query_one('#preset-name', Input).value.strip()
		if not name:
			self._set_status('⚠ Enter a name to save the preset')
			return
		values, knockouts, at_rest = self._form_state()
		saved = presets.save(name, values, knockouts, at_rest)
		select = self.query_one('#preset-select', Select)
		select.set_options([(n, n) for n in presets.all_names()])
		select.value = saved
		self.query_one('#preset-name', Input).value = ''
		# the current form now *is* this preset → no longer "edited"
		self._applied_state = (values, knockouts, at_rest)
		self._applied_preset_name = saved
		self._update_modified()
		self._set_status(f'✓ Saved preset “{saved}”')

	# ── run lifecycle ───────────────────────────────────────────────────────

	def _start_run(self) -> None:
		if self._sim_running:
			return
		values, knockouts, at_rest = self._form_state()
		try:
			spec = runspec.build_spec(
				values, knockouts, at_rest, self._entity_knockouts())
		except ValueError as exc:
			self._set_status(f'⚠ {exc}')
			return

		root = fp.ROOT_PATH
		self._sim_length = spec['length_sec']
		self._sim_live_path = live_csv_path(_OUTDIR, spec['seed'], root)
		sim_path = resolve_sim_path(_OUTDIR, root)
		os.makedirs(sim_path, exist_ok=True)
		cfg_path = os.path.join(sim_path, 'run_config.json')
		runspec.write_spec(spec, cfg_path)
		try:
			os.remove(self._sim_live_path)  # drop a stale trace from a prior run
		except OSError:
			pass
		self._reset_plots()

		cmd = [
			sys.executable, 'runscripts/manual/runFromConfig.py',
			'--config', cfg_path, '--out', _OUTDIR]
		env = dict(os.environ)
		env['PYTHONPATH'] = root
		env['OPENBLAS_NUM_THREADS'] = '1'

		self._sim_running = True
		self._sim_ok = None
		self.query_one('#run-btn', Button).disabled = True
		self._set_status(f'Running… 0 / {self._sim_length} s')
		self._sim_timer = self.set_interval(0.5, self._poll)
		self._run_worker(cmd, root, env)

	@work(thread=True, exclusive=True)
	def _run_worker(self, cmd: List[str], cwd: str, env: dict) -> None:
		try:
			proc = subprocess.run(
				cmd, cwd=cwd, env=env, capture_output=True, text=True)
			ok = proc.returncode == 0
			err = '' if ok else (proc.stderr or proc.stdout or '').strip()
		except Exception as exc:  # pragma: no cover - defensive
			ok, err = False, str(exc)
		self.call_from_thread(self._run_finished, ok, err)

	def _run_finished(self, ok: bool, err: str) -> None:
		if self._sim_timer is not None:
			self._sim_timer.stop()
			self._sim_timer = None
		self._poll()  # final read of the completed trace
		self._sim_running = False
		self._sim_ok = ok
		self.query_one('#run-btn', Button).disabled = False
		if ok:
			self._set_status(
				f'✓ Done — {self._sim_length} s. Edit and run again, '
				f'or q to quit.')
		else:
			tail = err.splitlines()[-1] if err else 'unknown error'
			self._set_status(f'✗ Run failed: {tail}')

	def _poll(self) -> None:
		if not self._sim_live_path:
			return
		cols = read_live_csv(self._sim_live_path)
		time = cols['time']
		for pid, ycol, ylabel, colour in _PANELS:
			self._draw(pid, ycol, time, cols[ycol], ylabel, colour)
		if self._sim_running and time:
			self._set_status(
				f'Running… {time[-1]:.0f} / {self._sim_length} s  '
				f'({len(time)} steps)')

	# ── plotting ────────────────────────────────────────────────────────────

	def _draw(self, pid: str, ycol: str, time: List[float],
			yvals: List[float], ylabel: str, colour: str) -> None:
		widget = self.query_one('#' + pid, PlotextPlot)
		plt = widget.plt
		plt.clear_data()
		if self._baseline:
			base_t, base_y = self._baseline.get('time', []), \
				self._baseline.get(ycol, [])
			if base_t and base_y:
				plt.plot(base_t, base_y, color='gray')  # reference, behind
		title = ylabel
		if time and yvals:
			plt.plot(time, yvals, color=colour)
			title = f'{ylabel}   peak {max(yvals):.1f}'
		if self._baseline:
			title += '  (grey = baseline)'
		plt.title(title)
		plt.xlabel('time (s)')
		widget.refresh()

	def _reset_plots(self) -> None:
		for pid, ycol, ylabel, colour in _PANELS:
			self._draw(pid, ycol, [], [], ylabel, colour)

	# ── helpers ─────────────────────────────────────────────────────────────

	def _set_status(self, text: str) -> None:
		self._last_status = text
		self.query_one('#status', Static).update(text)
