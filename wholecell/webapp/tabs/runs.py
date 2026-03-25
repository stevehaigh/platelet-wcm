"""Run Status tab — monitor submitted jobs."""

from __future__ import annotations

import json

import dash
from dash import dcc, html
from dash.dependencies import Input, Output

from wholecell.webapp.jobs import PHASE_DURATIONS, PHASES, JobManager


def layout() -> html.Div:
	"""Create the Run Status tab layout."""

	return html.Div(style={'padding': '20px'}, children=[
		html.Div(style={'display': 'flex', 'justifyContent': 'space-between', 'alignItems': 'center', 'marginBottom': '15px'}, children=[
			html.H3('Job Queue', style={'margin': 0}),
		]),
		dcc.Interval(id='runs-interval', interval=5000, n_intervals=0),
		html.Div(id='runs-table'),
	])


def register_callbacks(app: dash.Dash, job_manager: JobManager) -> None:
	"""Register Run Status tab callbacks."""

	@app.callback(
		Output('runs-table', 'children'),
		Input('runs-interval', 'n_intervals'),
	)
	def update_table(n):
		jobs = job_manager.list_jobs()
		if not jobs:
			return html.P('No jobs yet. Go to Configure to submit one.',
				style={'color': '#888', 'fontStyle': 'italic'})

		rows = []
		for job in jobs:
			status = job['status']
			status_color = {
				'queued': '#6e6e6e',
				'parca': '#0969da',
				'simulating': '#bf8700',
				'analyzing': '#8250df',
				'done': '#1a7f37',
				'failed': '#cf222e',
			}.get(status, '#333')

			config = json.loads(job.get('config_json', '{}'))
			variant = config.get('variant', '?')
			gens = config.get('generations', '?')
			seeds = config.get('init_sims', '?')

			rows.append(html.Tr([
				html.Td(f'#{job["id"]}'),
				html.Td(html.Span(status, style={
					'color': 'white', 'background': status_color,
					'padding': '2px 8px', 'borderRadius': '12px',
					'fontSize': '12px', 'fontWeight': 'bold',
				})),
				html.Td(variant),
				html.Td(f'{gens}G × {seeds}S'),
				html.Td(job.get('description', '')),
				html.Td(_format_time(job.get('started_at', ''))),
				html.Td(_format_time(job.get('finished_at', ''))),
				html.Td(
					html.Details([
						html.Summary('Error'),
						html.Pre(job.get('error_message', ''),
							style={'fontSize': '11px', 'maxHeight': '200px', 'overflow': 'auto'}),
					]) if job.get('error_message') else '',
				),
			]))

		return html.Table(
			style={'width': '100%', 'borderCollapse': 'collapse'},
			children=[
				html.Thead(html.Tr([
					html.Th('ID'), html.Th('Status'), html.Th('Variant'),
					html.Th('Scale'), html.Th('Description'),
					html.Th('Started'), html.Th('Finished'), html.Th(''),
				], style={'borderBottom': '2px solid #ddd', 'textAlign': 'left'})),
				html.Tbody(rows),
			],
		)


def _format_time(iso_str: str) -> str:
	"""Format an ISO timestamp for display."""
	if not iso_str:
		return '—'
	# Show just the time portion
	try:
		return iso_str.split('T')[1][:8]
	except (IndexError, AttributeError):
		return iso_str
