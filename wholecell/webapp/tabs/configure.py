"""Configure tab — simulation setup form."""

from __future__ import annotations

import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State

# Hardcoded to avoid importing the heavy model stack.
# Keep in sync with models/ecoli/sim/variants/__init__.py
VARIANT_NAMES = [
	'wildtype',
	'aa_synthesis_ko',
	'aa_synthesis_ko_shift',
	'aa_synthesis_sensitivity',
	'aa_uptake_sensitivity',
	'add_one_aa',
	'add_one_aa_shift',
	'condition',
	'gene_knockout',
	'mene_params',
	'metabolism_kinetic_objective_weight',
	'metabolism_secretion_penalty',
	'new_gene_internal_shift',
	'param_sensitivity',
	'ppgpp_conc',
	'ppgpp_limitations',
	'ppgpp_limitations_ribosome',
	'remove_aa_inhibition',
	'remove_aas_shift',
	'remove_one_aa',
	'remove_one_aa_shift',
	'rrna_operon_knockout',
	'rrna_location',
	'rrna_orientation',
	'tf_activity',
	'time_step',
	'timelines',
]


# Key simulation toggles (name, label, default)
TOGGLES = [
	('ppgpp_regulation', 'ppGpp regulation', True),
	('trna_charging', 'tRNA charging', True),
	('d_period_division', 'D-period division', True),
	('translation_supply', 'Translation supply', True),
	('superhelical_density', 'Superhelical density', False),
	('variable_elongation_transcription', 'Variable elongation (transcription)', False),
	('variable_elongation_translation', 'Variable elongation (translation)', False),
	('mechanistic_translation_supply', 'Mechanistic translation supply', False),
	('mechanistic_aa_transport', 'Mechanistic AA transport', False),
	('trna_attenuation', 'tRNA attenuation', False),
]


def layout() -> html.Div:
	"""Create the Configure tab layout."""

	# Build toggle checkboxes with defaults pre-checked
	default_toggles = [name for name, _, default in TOGGLES if default]

	return html.Div(style={'padding': '20px', 'maxWidth': '600px'}, children=[
		html.Div(style={'display': 'grid', 'gridTemplateColumns': '1fr 1fr', 'gap': '15px', 'marginBottom': '20px'}, children=[
			html.Div([
				html.Label('Variant type', style={'fontWeight': 'bold'}),
				dcc.Dropdown(
					id='config-variant',
					options=[{'label': v, 'value': v} for v in VARIANT_NAMES],
					value='wildtype',
				),
			]),
			html.Div([
				html.Label('Variant index range', style={'fontWeight': 'bold'}),
				html.Div(style={'display': 'flex', 'gap': '5px', 'alignItems': 'center'}, children=[
					dcc.Input(id='config-variant-first', type='number', value=0, min=0,
						style={'width': '60px'}),
					html.Span('to'),
					dcc.Input(id='config-variant-last', type='number', value=0, min=0,
						style={'width': '60px'}),
				]),
			]),
		]),

		html.Div(style={'display': 'grid', 'gridTemplateColumns': '1fr 1fr 1fr', 'gap': '15px', 'marginBottom': '20px'}, children=[
			html.Div([
				html.Label('Generations', style={'fontWeight': 'bold'}),
				dcc.Input(id='config-generations', type='number', value=1, min=1,
					style={'width': '100%'}),
			]),
			html.Div([
				html.Label('Seeds', style={'fontWeight': 'bold'}),
				dcc.Input(id='config-seeds', type='number', value=1, min=1,
					style={'width': '100%'}),
			]),
			html.Div([
				html.Label('Seed start', style={'fontWeight': 'bold'}),
				dcc.Input(id='config-seed-start', type='number', value=0, min=0,
					style={'width': '100%'}),
			]),
		]),

		html.Div(style={'marginBottom': '20px'}, children=[
			html.Label('Regulation toggles', style={'fontWeight': 'bold', 'marginBottom': '8px', 'display': 'block'}),
			dcc.Checklist(
				id='config-toggles',
				options=[{'label': f' {label}', 'value': name} for name, label, _ in TOGGLES],
				value=default_toggles,
				style={'display': 'grid', 'gridTemplateColumns': '1fr 1fr', 'gap': '4px'},
			),
		]),

		html.Div(style={'marginBottom': '20px'}, children=[
			html.Label('Description', style={'fontWeight': 'bold'}),
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
			style={
				'padding': '10px 24px', 'fontSize': '16px',
				'background': '#2ea44f', 'color': 'white', 'border': 'none',
				'borderRadius': '6px', 'cursor': 'pointer',
			},
		),
		html.Div(id='config-status', style={'marginTop': '10px'}),
	])


def register_callbacks(app: dash.Dash, on_submit) -> None:
	"""Register Configure tab callbacks.

	Args:
		on_submit: callable(config_dict) that submits a job.
			Returns a status message string.
	"""

	@app.callback(
		Output('config-status', 'children'),
		Input('config-run-button', 'n_clicks'),
		State('config-variant', 'value'),
		State('config-variant-first', 'value'),
		State('config-variant-last', 'value'),
		State('config-generations', 'value'),
		State('config-seeds', 'value'),
		State('config-seed-start', 'value'),
		State('config-toggles', 'value'),
		State('config-description', 'value'),
		prevent_initial_call=True,
	)
	def submit_run(n_clicks, variant, first_idx, last_idx, generations,
			seeds, seed_start, toggles, description):
		if not n_clicks:
			return ''

		config = {
			'variant': variant,
			'first_variant_index': first_idx or 0,
			'last_variant_index': last_idx or 0,
			'generations': generations or 1,
			'init_sims': seeds or 1,
			'seed': seed_start or 0,
			'description': description or '',
			'toggles': {name: (name in (toggles or [])) for name, _, _ in TOGGLES},
		}

		msg = on_submit(config)
		return html.Div(msg, style={'color': '#2ea44f', 'fontWeight': 'bold'})
