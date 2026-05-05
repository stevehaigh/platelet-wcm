"""Inspect tab — interactive listener data browser."""

from __future__ import annotations

import json
import os
import shutil
import sys
from typing import List, Tuple

import dash
from dash import dcc, html
from dash.dependencies import ALL, Input, Output, State
import numpy as np
import plotly.graph_objs as go

from wholecell.webapp import results


def make_run_options(out_path: str) -> List[dict]:
	"""Build dropdown options from available simulation directories."""

	options = []
	for sim_dir in results.find_sim_dirs(out_path):
		ts = results.dir_timestamp(sim_dir)
		for variant in results.find_variants(sim_dir):
			cells = results.find_cells(sim_dir, variant)
			if cells:
				cell_count = f"{len(cells)} cell{'s' if len(cells) > 1 else ''}"
				label = f"{variant} ({cell_count}) {ts}" if ts else f"{variant} ({cell_count})"
				value = f"{sim_dir}|{variant}"
				options.append({'label': label, 'value': value})
	return options


def parse_run_value(value: str) -> Tuple[str, str]:
	"""Parse a run dropdown value into (sim_dir, variant)."""
	parts = value.split('|', 1)
	if len(parts) != 2:
		raise ValueError(f'Invalid run value (expected sim_dir|variant): {value!r}')
	return parts[0], parts[1]


def _trace_label(trace_info: dict) -> str:
	"""Human-readable label for an overlay trace."""
	_, variant = parse_run_value(trace_info['run'])
	return f"{variant} / {trace_info['listener']} / {trace_info['column']}"


def layout(out_path: str) -> html.Div:
	"""Create the Inspect tab layout."""

	run_options = make_run_options(out_path)

	return html.Div(children=[
		dcc.ConfirmDialog(id='inspect-delete-confirm', message=''),
		dcc.Store(id='inspect-delete-pending', data=None),

		html.Div(className='grid-3', style={'marginBottom': '15px'}, children=[
			html.Div([
				html.Label('Run'),
				html.Div(style={'display': 'flex', 'gap': '6px'}, children=[
					dcc.Dropdown(
						id='inspect-run',
						options=run_options,
						value=run_options[0]['value'] if run_options else None,
						style={'flex': '1', 'minWidth': '0'},
					),
					html.Button('Delete', id='inspect-delete-run', n_clicks=0,
						style={'padding': '6px 12px', 'cursor': 'pointer',
							'color': '#c00', 'border': '1px solid #c00',
							'borderRadius': '4px', 'background': '#fff',
							'whiteSpace': 'nowrap'}),
				]),
			]),
			html.Div([
				html.Label('Listener'),
				dcc.Dropdown(id='inspect-listener'),
			]),
			html.Div([
				html.Label('Column'),
				dcc.Dropdown(id='inspect-column'),
			]),
		]),

		html.Div(style={'marginBottom': '15px', 'display': 'flex', 'gap': '20px', 'alignItems': 'center'}, children=[
			html.Label('Transform:'),
			dcc.Checklist(
				id='inspect-transform',
				options=[
					{'label': ' Normalize to t=0', 'value': 'normalize'},
					{'label': ' Log scale', 'value': 'log'},
				],
				value=[],
				inline=True,
				style={'display': 'flex', 'gap': '15px'},
			),
		]),

		dcc.Graph(
			id='inspect-graph',
			style={'height': '550px'},
			config={'displayModeBar': True, 'scrollZoom': True},
		),

		html.Div(style={'marginTop': '15px', 'padding': '10px', 'background': '#f0f0f0', 'borderRadius': '5px'}, children=[
			html.Label('Overlay traces:', style={'fontWeight': 'bold', 'marginBottom': '5px', 'display': 'block'}),
			html.Div(style={'display': 'grid', 'gridTemplateColumns': '1fr 1fr 1fr auto', 'gap': '10px', 'alignItems': 'end'}, children=[
				dcc.Dropdown(
					id='inspect-overlay-run',
					options=run_options,
					placeholder='Select run...',
				),
				dcc.Dropdown(id='inspect-overlay-listener', placeholder='Listener...'),
				dcc.Dropdown(id='inspect-overlay-column', placeholder='Column...'),
				html.Button('Add Trace', id='inspect-add-trace', n_clicks=0,
					style={'padding': '8px 16px', 'cursor': 'pointer'}),
			]),
			html.Div(id='inspect-overlay-trace-list', style={'marginTop': '8px'}),
		]),

		# Hidden store for overlay traces
		dcc.Store(id='inspect-overlay-traces', data=[]),
	])


def register_callbacks(app: dash.Dash, out_path: str) -> None:
	"""Register all Inspect tab callbacks."""

	@app.callback(
		Output('inspect-listener', 'options'),
		Output('inspect-listener', 'value'),
		Input('inspect-run', 'value'),
	)
	def update_listeners(run_value):
		if not run_value:
			return [], None
		sim_dir, variant = parse_run_value(run_value)
		cells = results.find_cells(sim_dir, variant)
		if not cells:
			return [], None
		listeners = results.find_listeners(cells[0]['simout_path'])
		options = [{'label': l, 'value': l} for l in listeners]
		# Prefer CalciumTrace (the headline platelet listener), then Mass
		if 'CalciumTrace' in listeners:
			default = 'CalciumTrace'
		elif 'Mass' in listeners:
			default = 'Mass'
		else:
			default = listeners[0] if listeners else None
		return options, default

	@app.callback(
		Output('inspect-column', 'options'),
		Output('inspect-column', 'value'),
		Input('inspect-listener', 'value'),
		State('inspect-run', 'value'),
	)
	def update_columns(listener, run_value):
		if not listener or not run_value:
			return [], None
		sim_dir, variant = parse_run_value(run_value)
		cells = results.find_cells(sim_dir, variant)
		if not cells:
			return [], None
		columns = results.find_columns(cells[0]['simout_path'], listener)
		options = [{'label': c, 'value': c} for c in columns]
		# Sensible defaults per listener
		if listener == 'CalciumTrace' and 'ca_cyt_nM' in columns:
			default = 'ca_cyt_nM'
		elif listener == 'Mass' and 'dryMass' in columns:
			default = 'dryMass'
		elif listener == 'Main' and 'time' in columns:
			default = 'time'
		else:
			default = columns[0] if columns else None
		return options, default

	# Same cascading logic for overlay dropdowns
	@app.callback(
		Output('inspect-overlay-listener', 'options'),
		Output('inspect-overlay-listener', 'value'),
		Input('inspect-overlay-run', 'value'),
	)
	def update_overlay_listeners(run_value):
		if not run_value:
			return [], None
		sim_dir, variant = parse_run_value(run_value)
		cells = results.find_cells(sim_dir, variant)
		if not cells:
			return [], None
		listeners = results.find_listeners(cells[0]['simout_path'])
		options = [{'label': l, 'value': l} for l in listeners]
		return options, None

	@app.callback(
		Output('inspect-overlay-column', 'options'),
		Output('inspect-overlay-column', 'value'),
		Input('inspect-overlay-listener', 'value'),
		State('inspect-overlay-run', 'value'),
	)
	def update_overlay_columns(listener, run_value):
		if not listener or not run_value:
			return [], None
		sim_dir, variant = parse_run_value(run_value)
		cells = results.find_cells(sim_dir, variant)
		if not cells:
			return [], None
		columns = results.find_columns(cells[0]['simout_path'], listener)
		options = [{'label': c, 'value': c} for c in columns]
		return options, None

	@app.callback(
		Output('inspect-delete-confirm', 'displayed'),
		Output('inspect-delete-confirm', 'message'),
		Output('inspect-delete-pending', 'data'),
		Input('inspect-delete-run', 'n_clicks'),
		State('inspect-run', 'value'),
		State('inspect-run', 'options'),
		prevent_initial_call=True,
	)
	def show_delete_confirm(n_clicks, run_value, run_options):
		if not n_clicks or not run_value:
			return False, '', None
		label = next((o['label'] for o in (run_options or []) if o['value'] == run_value), run_value)
		sim_dir, variant = parse_run_value(run_value)
		variant_path = os.path.join(sim_dir, variant)
		msg = (f"Permanently delete run '{label}'?\n\n"
			f"This will remove all simulation output at:\n{variant_path}")
		return True, msg, run_value

	@app.callback(
		Output('inspect-run', 'options'),
		Output('inspect-run', 'value'),
		Output('inspect-overlay-run', 'options'),
		Output('explore-left-run', 'options'),
		Output('explore-left-run', 'value'),
		Output('explore-right-run', 'options'),
		Output('explore-right-run', 'value'),
		Output('inspect-overlay-traces', 'data', allow_duplicate=True),
		Input('inspect-delete-confirm', 'submit_n_clicks'),
		State('inspect-delete-pending', 'data'),
		State('inspect-run', 'value'),
		State('explore-left-run', 'value'),
		State('explore-right-run', 'value'),
		State('inspect-overlay-traces', 'data'),
		prevent_initial_call=True,
	)
	def execute_delete(submit_clicks, pending_run, current_run, explore_left, explore_right, overlay_traces):
		if not submit_clicks or not pending_run:
			raise dash.exceptions.PreventUpdate
		sim_dir, variant = parse_run_value(pending_run)
		variant_path = os.path.join(sim_dir, variant)
		if os.path.isdir(variant_path):
			try:
				shutil.rmtree(variant_path)
			except Exception as e:
				print(f'Error deleting {variant_path}: {e}', file=sys.stderr)
				raise dash.exceptions.PreventUpdate
		new_options = make_run_options(out_path)
		explore_options = results.explore_run_options(out_path)
		valid_inspect = {o['value'] for o in new_options}
		valid_explore = {o['value'] for o in explore_options}
		new_value = current_run if current_run in valid_inspect else (new_options[0]['value'] if new_options else None)
		new_left = explore_left if explore_left in valid_explore else (explore_options[0]['value'] if explore_options else None)
		new_right = explore_right if explore_right in valid_explore else None
		new_traces = [t for t in (overlay_traces or []) if t['run'] != pending_run]
		return new_options, new_value, new_options, explore_options, new_left, explore_options, new_right, new_traces

	@app.callback(
		Output('inspect-overlay-traces', 'data'),
		Input('inspect-add-trace', 'n_clicks'),
		Input({'type': 'remove-overlay-trace', 'index': ALL}, 'n_clicks'),
		State('inspect-overlay-run', 'value'),
		State('inspect-overlay-listener', 'value'),
		State('inspect-overlay-column', 'value'),
		State('inspect-overlay-traces', 'data'),
	)
	def update_overlay_traces(add_clicks, remove_clicks, run_value, listener, column, existing_traces):
		traces = list(existing_traces or [])
		ctx = dash.callback_context
		if not ctx.triggered:
			return traces

		trigger_id = ctx.triggered[0]['prop_id']

		if trigger_id == 'inspect-add-trace.n_clicks':
			if add_clicks and all([run_value, listener, column]):
				new_trace = {'run': run_value, 'listener': listener, 'column': column}
				if new_trace not in traces:
					traces.append(new_trace)
		elif '"type":"remove-overlay-trace"' in trigger_id:
			# Extract the index from the triggered component id
			prop_id = trigger_id.rsplit('.', 1)[0]
			try:
				component_id = json.loads(prop_id)
				idx = component_id['index']
				if ctx.triggered[0]['value'] and 0 <= idx < len(traces):
					traces.pop(idx)
			except (json.JSONDecodeError, KeyError, IndexError):
				pass

		return traces

	@app.callback(
		Output('inspect-overlay-trace-list', 'children'),
		Input('inspect-overlay-traces', 'data'),
	)
	def render_overlay_trace_list(traces):
		if not traces:
			return []
		rows = []
		for i, trace_info in enumerate(traces):
			label = _trace_label(trace_info)
			rows.append(html.Div(
				style={'display': 'flex', 'alignItems': 'center', 'gap': '8px', 'marginTop': '4px'},
				children=[
					html.Span(label, style={'flex': '1', 'fontSize': '0.9em'}),
					html.Button('×', id={'type': 'remove-overlay-trace', 'index': i},
						n_clicks=0,
						style={
							'padding': '2px 8px',
							'cursor': 'pointer',
							'border': '1px solid #ccc',
							'borderRadius': '3px',
							'background': '#fff',
							'lineHeight': '1',
						}),
				],
			))
		return rows

	@app.callback(
		Output('inspect-graph', 'figure'),
		Input('inspect-column', 'value'),
		Input('inspect-transform', 'value'),
		Input('inspect-overlay-traces', 'data'),
		State('inspect-run', 'value'),
		State('inspect-listener', 'value'),
	)
	def update_graph(column, transforms, overlay_traces, run_value, listener):
		if not all([column, run_value, listener]):
			return go.Figure()

		transforms = transforms or []
		fig = go.Figure()

		# Load primary trace
		_add_traces_to_fig(fig, run_value, listener, column, transforms,
			is_primary=True)

		# Load overlay traces
		for trace_info in (overlay_traces or []):
			_add_traces_to_fig(
				fig,
				trace_info['run'],
				trace_info['listener'],
				trace_info['column'],
				transforms,
				is_primary=False,
			)

		fig.update_layout(
			xaxis_title='Time (s)',
			yaxis_title=f'{listener} / {column}',
			hovermode='closest',
			template='plotly_white',
			legend=dict(orientation='h', yanchor='bottom', y=1.02),
		)
		if 'log' in transforms:
			fig.update_yaxis(type='log')

		return fig


def _add_traces_to_fig(
		fig: go.Figure,
		run_value: str,
		listener: str,
		column: str,
		transforms: list,
		is_primary: bool = True,
		) -> None:
	"""Load data and add traces to a Plotly figure."""

	try:
		sim_dir, variant = parse_run_value(run_value)
		cells = results.find_cells(sim_dir, variant)
		if not cells:
			return

		simout = cells[0]['simout_path']
		time = results.load_time(simout)
		data, labels = results.load_column(simout, listener, column)
	except Exception:
		return

	if 'normalize' in transforms and data.shape[0] > 0:
		norm_row = data[0, :]
		norm_row = np.where(norm_row == 0, 1, norm_row)
		data = data / norm_row

	prefix = '' if is_primary else f'{variant} · '
	dash_style = None if is_primary else 'dash'
	x_vals = time if time is not None else np.arange(data.shape[0])
	n_series = min(len(labels), data.shape[1])

	# For columns with many series (>20), just plot the sum
	if n_series > 20:
		fig.add_trace(go.Scatter(
			x=x_vals,
			y=data.sum(axis=1),
			name=f'{prefix}{column} (sum of {n_series})',
			line=dict(dash=dash_style),
			mode='lines+markers',
		))
	else:
		for i in range(n_series):
			name = f'{prefix}{labels[i]}' if n_series > 1 else f'{prefix}{column}'
			fig.add_trace(go.Scatter(
				x=x_vals,
				y=data[:, i],
				name=name,
				line=dict(dash=dash_style),
				mode='lines+markers',
			))
