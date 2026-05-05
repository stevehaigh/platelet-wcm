"""Configure tab — platelet simulation setup form."""

from __future__ import annotations

import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State


# Preset simulation configurations.
# Each preset fills the form below with biologically interesting settings.
PRESETS = [
	{
		'id': 'preset-smoke',
		'label': '⚡ Smoke test (60 s)',
		'description': (
			'60-second platelet simulation — fastest check that the engine runs '
			'end-to-end and writes listener data.'
		),
		'config': {
			'length_sec': 60,
			'seed': 0,
			'description': 'Smoke test — 60 s',
		},
	},
	{
		'id': 'preset-transient',
		'label': '🩸 IP3 Ca²⁺ transient (200 s)',
		'description': (
			'200-second simulation with IP3 forcing on (Dolan 2014 Fig. S2 curve). '
			'Reproduces the Phase 1 Ca²⁺ transient — see calcium_trace plot.'
		),
		'config': {
			'length_sec': 200,
			'seed': 0,
			'description': 'IP3 Ca²⁺ transient — 200 s',
		},
	},
	{
		'id': 'preset-resting',
		'label': '🛌 Resting (300 s)',
		'description': (
			'300-second resting simulation — useful for verifying steady state at '
			'the current resting initial conditions.'
		),
		'config': {
			'length_sec': 300,
			'seed': 0,
			'description': 'Resting — 300 s',
		},
	},
]


def layout() -> html.Div:
	"""Create the Configure tab layout."""

	preset_buttons = [
		html.Button(
			p['label'],
			id=p['id'],
			n_clicks=0,
			className='btn-preset',
			title=p['description'],
		)
		for p in PRESETS
	]

	return html.Div(style={'maxWidth': '700px'}, children=[

		# Presets section
		html.Div(style={'marginBottom': '24px'}, children=[
			html.Label('Quick start presets'),
			html.P(
				'Click a preset to fill in the form below, then adjust as needed '
				'and click Run Simulation.',
				style={'color': '#57606a', 'fontSize': '13px', 'margin': '4px 0 10px 0'},
			),
			html.Div(preset_buttons, style={'display': 'flex', 'flexWrap': 'wrap', 'gap': '8px'}),
			html.Div(id='preset-description', style={
				'marginTop': '8px', 'fontSize': '13px', 'color': '#57606a',
				'fontStyle': 'italic', 'minHeight': '18px',
			}),
		]),

		html.Hr(style={'border': 'none', 'borderTop': '1px solid #d0d7de', 'marginBottom': '20px'}),

		html.Div(className='grid-2', style={'marginBottom': '20px'}, children=[
			html.Div([
				html.Label('Length (seconds)'),
				dcc.Input(id='config-length-sec', type='number', value=200, min=1,
					style={'width': '100%'}),
				html.Div(
					'Simulated wall-clock seconds. Typical: 60 (smoke), 200 (transient), 300+ (resting).',
					style={'color': '#57606a', 'fontSize': '12px', 'marginTop': '4px'},
				),
			]),
			html.Div([
				html.Label('Random seed'),
				dcc.Input(id='config-seed', type='number', value=0, min=0,
					style={'width': '100%'}),
				html.Div(
					'Currently no stochastic processes use the seed; included for future use.',
					style={'color': '#57606a', 'fontSize': '12px', 'marginTop': '4px'},
				),
			]),
		]),

		html.Div(style={'marginBottom': '20px'}, children=[
			html.Label('Description'),
			dcc.Input(
				id='config-description', type='text',
				placeholder='Brief description of this run...',
				style={'width': '100%'},
			),
		]),

		html.Button(
			'Run Simulation',
			id='config-run-button',
			n_clicks=0,
			className='btn-primary',
		),
		html.Div(id='config-status', style={'marginTop': '10px'}),
	])


def register_callbacks(app: dash.Dash, on_submit) -> None:
	"""Register Configure tab callbacks.

	Args:
		on_submit: callable(config_dict) that submits a job.
			Returns a status message string.
	"""

	preset_ids = [p['id'] for p in PRESETS]
	preset_lookup = {p['id']: p for p in PRESETS}

	@app.callback(
		Output('config-length-sec', 'value'),
		Output('config-seed', 'value'),
		Output('config-description', 'value'),
		Output('preset-description', 'children'),
		[Input(pid, 'n_clicks') for pid in preset_ids],
		prevent_initial_call=True,
	)
	def apply_preset(*args):
		from dash import ctx
		triggered = ctx.triggered_id
		if not triggered or triggered not in preset_lookup:
			raise dash.exceptions.PreventUpdate
		p = preset_lookup[triggered]['config']
		desc = preset_lookup[triggered]['description']
		return p['length_sec'], p['seed'], p['description'], desc

	@app.callback(
		Output('config-status', 'children'),
		Input('config-run-button', 'n_clicks'),
		State('config-length-sec', 'value'),
		State('config-seed', 'value'),
		State('config-description', 'value'),
		prevent_initial_call=True,
	)
	def submit_run(n_clicks, length_sec, seed, description):
		if not n_clicks:
			return ''

		config = {
			'length_sec': int(length_sec or 60),
			'seed': int(seed or 0),
			'description': description or '',
		}

		msg = on_submit(config)
		return html.Div(msg, style={'color': '#2ea44f', 'fontWeight': 'bold'})
