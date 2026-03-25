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
		html.Div(id='tab-content'),
	])

	# Lazy-import tabs to avoid circular imports
	from wholecell.webapp.tabs import configure, explore, inspect_data, runs

	# Register tab callbacks
	inspect_data.register_callbacks(app, out_path)
	explore.register_callbacks(app, out_path)
	configure.register_callbacks(app, on_submit=lambda cfg: _submit_job(job_manager, cfg))
	runs.register_callbacks(app, job_manager)

	@app.callback(
		Output('tab-content', 'children'),
		Input('main-tabs', 'value'),
	)
	def render_tab(tab):
		if tab == 'inspect':
			return inspect_data.layout(out_path)
		elif tab == 'explore':
			return explore.layout(out_path)
		elif tab == 'configure':
			return configure.layout()
		elif tab == 'runs':
			return runs.layout()
		return html.P('Unknown tab')

	return app


def _submit_job(job_manager: JobManager, config: dict) -> str:
	"""Submit a job and return a status message."""

	job_id = job_manager.submit(config)
	return f'Job #{job_id} submitted! Switch to Run Status tab to monitor progress.'
