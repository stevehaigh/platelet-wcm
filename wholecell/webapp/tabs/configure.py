"""Configure tab — simulation setup form."""

from __future__ import annotations

import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State

# Hardcoded to avoid importing the heavy model stack.
# Keep in sync with models/ecoli/sim/variants/__init__.py
VARIANT_NAMES = [
	# Platelet model (no ParCa required; "variant index" = length in days)
	'platelet',
	# E. coli model variants
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


# Key simulation toggles (name, label, default, hint)
TOGGLES = [
	('ppgpp_regulation', 'ppGpp regulation', True,
		'ppGpp is the stringent-response alarmone: high levels slow ribosome synthesis when amino acids are scarce. '
		'Disable to use a simplified constant ribosome allocation.'),
	('trna_charging', 'tRNA charging', True,
		'Models the enzymatic charging of tRNA with amino acids before each translation step. '
		'Disable for a simplified translation model that ignores charging kinetics.'),
	('d_period_division', 'D-period division', True,
		'Enforces the D-period: a fixed delay between finishing DNA replication and cell division. '
		'Disable to allow division immediately after replication completes.'),
	('translation_supply', 'Translation supply', True,
		'Links ribosome elongation rate to the instantaneous amino acid supply. '
		'Disable to use a fixed average elongation rate regardless of nutrient availability.'),
	('superhelical_density', 'Superhelical density', False,
		'Models the effect of DNA supercoiling on transcription rates. '
		'Enable for a more mechanistic treatment of transcription; increases simulation time.'),
	('variable_elongation_transcription', 'Variable elongation (transcription)', False,
		'Allows RNAP elongation speed to vary per gene rather than using a single average rate. '
		'Enables more realistic transcript length distributions.'),
	('variable_elongation_translation', 'Variable elongation (translation)', False,
		'Allows ribosome elongation speed to vary per mRNA rather than using a single average rate. '
		'Enables more realistic protein synthesis dynamics.'),
	('mechanistic_translation_supply', 'Mechanistic translation supply', False,
		'Uses a detailed kinetic model of amino acid supply to ribosomes instead of a lumped flux. '
		'More accurate but significantly slower to simulate.'),
	('mechanistic_aa_transport', 'Mechanistic AA transport', False,
		'Models amino acid import and export across the inner membrane mechanistically. '
		'Enable to study transporter saturation effects during nutrient shifts.'),
	('trna_attenuation', 'tRNA attenuation', False,
		'Models tRNA-mediated transcriptional attenuation of amino acid biosynthesis operons (e.g. trp, his). '
		'Enable for a more detailed treatment of biosynthetic gene regulation.'),
]


# Preset simulation configurations.
# Each preset fills the form with biologically interesting settings.
# variant_index: the variant index to use (first=last for single index).
# Timelines index reference (from timelines_def.tsv):
#   2 = add_aa (minimal → minimal+AA at t=1200s)
#   25 = cut_aa (minimal+AA → minimal at t=1200s)
#   1 = cut_glucose (minimal → no_glucose at t=1200s)
#   18 = cut_oxygen (minimal → anaerobic at t=1200s)
# Condition index reference (from condition_defs.tsv):
#   0 = minimal media (control)
#   1 = with amino acids (faster growth)
#   2 = acetate (slower growth)
PRESETS = [
	{
		'id': 'preset-platelet',
		'label': '🩸 Platelet resting state (1 day)',
		'description': (
			'Platelet in resting state for 1 day — no ParCa required. '
			'The "Variant index" field controls duration in days. '
			'Inspect results in the Inspect Data tab (Mass → dryMass).'
		),
		'config': {
			'variant': 'platelet',
			'variant_first': 1,
			'variant_last': 1,
			'generations': 1,
			'seeds': 1,
			'seed_start': 0,
			'description': 'Platelet resting state — 1 day',
		},
	},
	{
		'id': 'preset-wildtype',
		'label': '🧫 Wildtype baseline',
		'description': 'Single wildtype cell in minimal media — the standard reference simulation.',
		'config': {
			'variant': 'wildtype',
			'variant_first': 0,
			'variant_last': 0,
			'generations': 1,
			'seeds': 1,
			'seed_start': 0,
			'description': 'Wildtype baseline — minimal media',
		},
	},
	{
		'id': 'preset-aa-shift',
		'label': '📈 Nutrient upshift (add amino acids)',
		'description': 'Cell growing in minimal media, then amino acids added at 20 min. Watch growth rate accelerate.',
		'config': {
			'variant': 'timelines',
			'variant_first': 2,
			'variant_last': 2,
			'generations': 1,
			'seeds': 1,
			'seed_start': 0,
			'description': 'Nutrient upshift: minimal → minimal+AA at t=1200s',
		},
	},
	{
		'id': 'preset-aa-downshift',
		'label': '📉 Nutrient downshift (remove amino acids)',
		'description': 'Cell growing in rich media (with amino acids), then amino acids removed at 20 min. Watch ppGpp stress response.',
		'config': {
			'variant': 'timelines',
			'variant_first': 25,
			'variant_last': 25,
			'generations': 1,
			'seeds': 1,
			'seed_start': 0,
			'description': 'Nutrient downshift: minimal+AA → minimal at t=1200s',
		},
	},
	{
		'id': 'preset-anaerobic',
		'label': '🔴 Switch to anaerobic',
		'description': 'Cell growing aerobically, then oxygen removed at 20 min. Watch metabolism shift.',
		'config': {
			'variant': 'timelines',
			'variant_first': 18,
			'variant_last': 18,
			'generations': 1,
			'seeds': 1,
			'seed_start': 0,
			'description': 'Anaerobic shift: minimal → minimal-oxygen at t=1200s',
		},
	},
	{
		'id': 'preset-rich-media',
		'label': '⚡ Rich media (fast growth)',
		'description': 'Wildtype cell growing in minimal+amino acids media. Compare growth rate and ribosome allocation to minimal.',
		'config': {
			'variant': 'condition',
			'variant_first': 1,
			'variant_last': 1,
			'generations': 1,
			'seeds': 1,
			'seed_start': 0,
			'description': 'Rich media: minimal + amino acids',
		},
	},
	{
		'id': 'preset-multi-seed',
		'label': '🎲 Wildtype × 3 seeds',
		'description': 'Three independent wildtype cells with different random seeds — shows cell-to-cell variability.',
		'config': {
			'variant': 'wildtype',
			'variant_first': 0,
			'variant_last': 0,
			'generations': 1,
			'seeds': 3,
			'seed_start': 0,
			'description': 'Wildtype × 3 seeds — cell-to-cell variability',
		},
	},
	{
		'id': 'preset-ppgpp',
		'label': '⚠️ ppGpp sweep (starvation signal)',
		'description': 'Three ppGpp levels: low (0.4×), control (1×), high (1.6×). High ppGpp mimics amino acid starvation — watch ribosome counts drop.',
		'config': {
			'variant': 'ppgpp_conc',
			'variant_first': 1,
			'variant_last': 3,
			'generations': 1,
			'seeds': 1,
			'seed_start': 0,
			'description': 'ppGpp sweep: low/control/high — stringent response',
		},
	},
	{
		'id': 'preset-acetate',
		'label': '🐢 Slow growth (acetate carbon source)',
		'description': 'Acetate as sole carbon source — much slower growth than glucose. Compare mass accumulation and ribosome allocation.',
		'config': {
			'variant': 'condition',
			'variant_first': 2,
			'variant_last': 2,
			'generations': 1,
			'seeds': 1,
			'seed_start': 0,
			'description': 'Slow growth: acetate carbon source (condition index 2)',
		},
	},
]


def layout() -> html.Div:
	"""Create the Configure tab layout."""

	# Build toggle checkboxes with defaults pre-checked
	default_toggles = [name for name, _, default, _ in TOGGLES if default]

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
				'Click a preset to fill in the form below, then adjust as needed and click Run Simulation.',
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
				html.Label('Variant type'),
				dcc.Dropdown(
					id='config-variant',
					options=[{'label': v, 'value': v} for v in VARIANT_NAMES],
					value='wildtype',
				),
			]),
			html.Div([
				html.Label('Variant index range'),
				html.Div(style={'display': 'flex', 'gap': '5px', 'alignItems': 'center'}, children=[
					dcc.Input(id='config-variant-first', type='number', value=0, min=0,
						style={'width': '60px'}),
					html.Span('to'),
					dcc.Input(id='config-variant-last', type='number', value=0, min=0,
						style={'width': '60px'}),
				]),
			]),
		]),

		html.Div(className='grid-3', style={'marginBottom': '20px'}, children=[
			html.Div([
				html.Label('Generations'),
				dcc.Input(id='config-generations', type='number', value=1, min=1,
					style={'width': '100%'}),
			]),
			html.Div([
				html.Label('Seeds'),
				dcc.Input(id='config-seeds', type='number', value=1, min=1,
					style={'width': '100%'}),
			]),
			html.Div([
				html.Label('Seed start'),
				dcc.Input(id='config-seed-start', type='number', value=0, min=0,
					style={'width': '100%'}),
			]),
		]),

		html.Div(style={'marginBottom': '20px'}, children=[
			html.Label('Regulation toggles', style={'marginBottom': '8px'}),
			dcc.Checklist(
				id='config-toggles',
				options=[
					{
						'label': html.Span([
							f' {label} ',
							html.Span('ℹ', className='toggle-hint', title=hint),
						]),
						'value': name,
					}
					for name, label, _, hint in TOGGLES
				],
				value=default_toggles,
				className='grid-2',
				style={'gap': '4px'},
			),
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

	# One combined callback for all preset buttons + form submission
	preset_ids = [p['id'] for p in PRESETS]
	preset_lookup = {p['id']: p for p in PRESETS}

	@app.callback(
		Output('config-variant', 'value'),
		Output('config-variant-first', 'value'),
		Output('config-variant-last', 'value'),
		Output('config-generations', 'value'),
		Output('config-seeds', 'value'),
		Output('config-seed-start', 'value'),
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
		return (
			p['variant'], p['variant_first'], p['variant_last'],
			p['generations'], p['seeds'], p['seed_start'],
			p['description'], desc,
		)

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
			'toggles': {name: (name in (toggles or [])) for name, _, _, _ in TOGGLES},
		}

		msg = on_submit(config)
		return html.Div(msg, style={'color': '#2ea44f', 'fontWeight': 'bold'})
