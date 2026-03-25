"""wcEcoli Web UI — Dash application factory."""

from __future__ import annotations

import os

import dash
from dash import dcc, html
from dash.dependencies import Input, Output

from wholecell.webapp.jobs import DB_FILENAME, JobManager


def create_app(out_path: str = None, wcecoli_root: str = None) -> dash.Dash:
	"""Create and configure the Dash application.

	Args:
		out_path: Path to the 'out/' directory with simulation results.
		wcecoli_root: Path to the wcEcoli repository root.
	"""

	if wcecoli_root is None:
		wcecoli_root = os.path.dirname(os.path.dirname(
			os.path.dirname(os.path.abspath(__file__))))
	if out_path is None:
		out_path = os.path.join(wcecoli_root, 'out')

	# Job manager
	db_path = os.path.join(wcecoli_root, DB_FILENAME)
	job_manager = JobManager(db_path, wcecoli_root)

	app = dash.Dash(
		__name__,
		suppress_callback_exceptions=True,
		title='wcEcoli Web UI',
	)

	# Import tabs
	from wholecell.webapp.tabs import configure, explore, inspect_data, runs

	# Build all tab layouts upfront so Dash registers every component ID
	# in the initial layout. Tabs are shown/hidden via CSS display property.
	tab_ids = ['inspect', 'explore', 'configure', 'runs']
	tab_contents = {
		'inspect': inspect_data.layout(out_path),
		'explore': explore.layout(out_path),
		'configure': configure.layout(),
		'runs': runs.layout(),
	}

	app.layout = html.Div([
		html.Div(
			style={
				'background': '#24292f', 'color': 'white',
				'padding': '12px 24px', 'display': 'flex',
				'alignItems': 'center', 'gap': '20px',
			},
			children=[
				html.H2('wcEcoli', style={'margin': 0, 'fontWeight': '300'}),
				html.Span('Whole-Cell E. coli Simulation',
					style={'opacity': '0.7', 'fontSize': '14px'}),
			],
		),
		dcc.Tabs(
			id='main-tabs',
			value='inspect',
			style={'borderBottom': '1px solid #ddd'},
			children=[
				dcc.Tab(label='Inspect Data', value='inspect'),
				dcc.Tab(label='Explore Plots', value='explore'),
				dcc.Tab(label='Configure', value='configure'),
				dcc.Tab(label='Run Status', value='runs'),
			],
		),
		# All tabs rendered, only one visible at a time
		*[html.Div(
			tab_contents[tid],
			id=f'tab-panel-{tid}',
			style={'display': 'block' if tid == 'inspect' else 'none'},
		) for tid in tab_ids],
	])

	# Register tab callbacks
	inspect_data.register_callbacks(app, out_path)
	explore.register_callbacks(app, out_path)
	configure.register_callbacks(app, on_submit=lambda cfg: _submit_job(job_manager, cfg))
	runs.register_callbacks(app, job_manager)

	@app.callback(
		[Output(f'tab-panel-{tid}', 'style') for tid in tab_ids],
		Input('main-tabs', 'value'),
	)
	def render_tab(tab):
		return [
			{'display': 'block'} if tid == tab else {'display': 'none'}
			for tid in tab_ids
		]

	return app


def _submit_job(job_manager: JobManager, config: dict) -> str:
	"""Submit a job and return a status message."""

	job_id = job_manager.submit(config)
	return f'Job #{job_id} submitted! Switch to Run Status tab to monitor progress.'
