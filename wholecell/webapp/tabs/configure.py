"""Configure tab — platelet simulation setup form."""

from __future__ import annotations

import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State


# Preset simulation configurations.
# Each preset fills the form below with biologically interesting settings.
# A preset is genuinely different from another iff at least one of
# ca_ex_mM / at_rest / agonist_delay_s differs (these three set the
# biology); length_sec just controls how much of the response is
# observed and is allowed to vary independently.
PRESETS = [
	{
		'id': 'preset-transient',
		'label': 'Agonist Ca²⁺ transient (200 s, +Ca²⁺)',
		'description': (
			'Phase 1 reproduction: 200 s with the default thrombin + ADP '
			'time courses and 1.2 mM extracellular Ca²⁺. Headline '
			'transient — see the 5-panel calcium_trace plot.'
		),
		'config': {
			'length_sec': 200,
			'seed': 0,
			'ca_ex_mM': 1.2,
			'at_rest': False,
			'agonist_delay_s': 0,
			'description': 'Agonist transient (+Ca²⁺) — 200 s',
		},
	},
	{
		'id': 'preset-delayed-transient',
		'label': 'Agonist transient (60 s settle + 200 s, +Ca²⁺)',
		'description': (
			'Model settles at its natural fixed point for 60 s, then the '
			'thrombin / ADP stimulus is applied. Supervisor-suggested '
			'approach: ignore the start-up transient and read the Ca²⁺ '
			'response from a settled baseline.'
		),
		'config': {
			'length_sec': 260,
			'seed': 0,
			'ca_ex_mM': 1.2,
			'at_rest': False,
			'agonist_delay_s': 60,
			'description': 'Agonist transient (60 s settle, +Ca²⁺) — 260 s',
		},
	},
	{
		'id': 'preset-edta',
		'label': 'EDTA transient (200 s, no Ca²⁺_ex)',
		'description': (
			'Phase 3 EDTA condition: 200 s with default agonist stimulation '
			'but extracellular Ca²⁺ = 0. SOCE is correctly inactive; compare '
			'against the +Ca²⁺ transient to see the SOCE-dependent shape.'
		),
		'config': {
			'length_sec': 200,
			'seed': 0,
			'ca_ex_mM': 0.0,
			'at_rest': False,
			'agonist_delay_s': 0,
			'description': 'Agonist transient (EDTA) — 200 s',
		},
	},
	{
		'id': 'preset-resting',
		'label': 'Resting (300 s, no stimulus)',
		'description': (
			'300 s at rest — all agonist peaks zero, no extracellular '
			'stimulus. Useful for inspecting the model at rest and '
			'verifying the resting fixed point of the ODE.'
		),
		'config': {
			'length_sec': 300,
			'seed': 0,
			'ca_ex_mM': 1.2,
			'at_rest': True,
			'agonist_delay_s': 0,
			'description': 'Resting (no stimulus) — 300 s',
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
					'Simulated wall-clock seconds. Typical: 200 (transient response) or 300+ (resting / steady-state inspection).',
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

		html.Div(className='grid-2', style={'marginBottom': '20px'}, children=[
			html.Div([
				html.Label('Extracellular Ca²⁺ (mM)'),
				dcc.Input(id='config-ca-ex-mM', type='number', value=1.2, min=0,
					step=0.1, style={'width': '100%'}),
				html.Div(
					'Dolan 2014 nominal = 1.2. Set to 0 for the EDTA / no-extracellular-Ca²⁺ condition (SOCE inactive).',
					style={'color': '#57606a', 'fontSize': '12px', 'marginTop': '4px'},
				),
			]),
			html.Div([
				html.Label('Agonist stimulus delay (seconds)'),
				dcc.Input(id='config-agonist-delay-s', type='number', value=0, min=0,
					step=10, style={'width': '100%'}),
				html.Div(
					'Settling time before the thrombin / ADP / ext-ATP time courses start. 0 = immediate. 60 = let the model reach its natural fixed point first.',
					style={'color': '#57606a', 'fontSize': '12px', 'marginTop': '4px'},
				),
			]),
		]),

		html.Div(className='grid-2', style={'marginBottom': '20px'}, children=[
			html.Div([
				html.Label('Stimulation'),
				dcc.Checklist(
					id='config-at-rest',
					options=[{'label': ' Run at rest (zero all agonist peaks)',
						'value': 'on'}],
					value=[],
				),
				html.Div(
					'Default OFF — applies the standard thrombin / ADP / ext-ATP time courses via the GPCR cascade. Check for a resting / un-stimulated sim (all agonist peaks set to 0).',
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
		Output('config-ca-ex-mM', 'value'),
		Output('config-agonist-delay-s', 'value'),
		Output('config-at-rest', 'value'),
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
		at_rest_value = ['on'] if p['at_rest'] else []
		return (p['length_sec'], p['seed'], p['ca_ex_mM'],
			p.get('agonist_delay_s', 0), at_rest_value, p['description'], desc)

	@app.callback(
		Output('config-status', 'children'),
		Input('config-run-button', 'n_clicks'),
		State('config-length-sec', 'value'),
		State('config-seed', 'value'),
		State('config-ca-ex-mM', 'value'),
		State('config-agonist-delay-s', 'value'),
		State('config-at-rest', 'value'),
		State('config-description', 'value'),
		prevent_initial_call=True,
	)
	def submit_run(n_clicks, length_sec, seed, ca_ex_mM, agonist_delay_s,
			at_rest_list, description):
		if not n_clicks:
			return ''

		config = {
			'length_sec': int(length_sec or 60),
			'seed': int(seed or 0),
			'ca_ex_mM': float(ca_ex_mM if ca_ex_mM is not None else 1.2),
			'agonist_delay_s': float(agonist_delay_s if agonist_delay_s is not None else 0),
			'at_rest': 'on' in (at_rest_list or []),
			'description': description or '',
		}

		msg = on_submit(config)
		return html.Div(msg, style={'color': '#2ea44f', 'fontWeight': 'bold'})
