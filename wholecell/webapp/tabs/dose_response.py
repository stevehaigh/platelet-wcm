"""Dose Response tab — browse 2-D dose-sweep results.

Lists every sweep directory under `out/` (those containing both
`sweep_summary.json` and `sweep.npz`, written by `runDoseSweep.py`) in a
dropdown, plus an observable selector. Renders an interactive 3-D Plotly
surface and a 2-D heatmap for the selected sweep + observable.

The surface-building helper is imported from `plotDoseSweepInteractive.py`
to keep the rendering identical to the standalone HTML preview.
"""

from __future__ import annotations

import json
import os

import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import numpy as np
import plotly.graph_objects as go

from runscripts.manual.plotDoseSweepInteractive import (
	make_surface, mathtext_to_html,
)
from runscripts.manual.runDoseSweep import OBSERVABLES
from wholecell.webapp import results


_OBS_MAP: dict[str, tuple[str, str, str]] = {
	k: (label, unit, cmap) for k, label, unit, cmap in OBSERVABLES
}
_PLOTLY_CMAPS: dict[str, str] = {
	'viridis': 'Viridis', 'plasma': 'Plasma',
	'cividis': 'Cividis', 'magma': 'Magma',
}


def _observable_dropdown_options() -> list[dict]:
	return [
		{'label': f'{mathtext_to_html(label)} ({mathtext_to_html(unit)})',
		 'value': key}
		for key, label, unit, _ in OBSERVABLES
	]


def layout(out_path: str) -> html.Div:
	"""Create the Dose Response tab layout."""
	sweep_options = results.dose_sweep_options(out_path)
	default_sweep = sweep_options[0]['value'] if sweep_options else None

	return html.Div(children=[
		html.Div(className='grid-2', style={'marginBottom': '15px'}, children=[
			html.Div([
				html.Label('Sweep'),
				dcc.Dropdown(
					id='dose-sweep-select',
					options=sweep_options,
					value=default_sweep,
					placeholder='Select a sweep…',
				),
				html.Div(
					'Sweeps are produced by '
					'runscripts/manual/runDoseSweep.py and listed newest-first.',
					style={'color': '#57606a', 'fontSize': '12px',
						'marginTop': '4px'},
				),
			]),
			html.Div([
				html.Label('Observable'),
				dcc.Dropdown(
					id='dose-observable-select',
					options=_observable_dropdown_options(),
					value='peak_ip3_nM',
				),
				html.Div(
					'Peak IP3 gives the cleanest dose-response surface; '
					'AUC Ca²⁺ is the most discriminating integrated metric.',
					style={'color': '#57606a', 'fontSize': '12px',
						'marginTop': '4px'},
				),
			]),
		]),

		html.Div(id='dose-sweep-meta', style={
			'color': '#57606a', 'fontSize': '13px',
			'marginBottom': '10px', 'fontStyle': 'italic',
		}),

		html.Div(className='grid-2', style={'gap': '10px'}, children=[
			dcc.Graph(id='dose-surface-graph', style={'height': '560px'}),
			dcc.Graph(id='dose-heatmap-graph', style={'height': '560px'}),
		]),

		# Live refresh: re-scan out/ every 5 s so newly-completed sweeps
		# appear in the dropdown without restarting the server.
		dcc.Interval(id='dose-refresh-interval', interval=5000, n_intervals=0),
	])


def _empty_figure(message: str) -> go.Figure:
	"""Render an empty figure with a centered message — used when no sweep
	is selected or the sweep is missing on disk."""
	fig = go.Figure()
	fig.add_annotation(
		text=message, x=0.5, y=0.5, xref='paper', yref='paper',
		showarrow=False, font=dict(size=14, color='#57606a'),
	)
	fig.update_layout(
		xaxis=dict(visible=False), yaxis=dict(visible=False),
		margin=dict(l=0, r=0, t=0, b=0),
	)
	return fig


def _make_heatmap(adp_grid: np.ndarray, thr_grid: np.ndarray,
		matrix: np.ndarray, label: str, unit: str, colorscale: str) -> go.Figure:
	"""2-D heatmap with log-axis ticks back in original units."""
	log_adp = np.log10(adp_grid)
	log_thr = np.log10(thr_grid)
	label_html = mathtext_to_html(label)
	unit_html = mathtext_to_html(unit)
	fig = go.Figure(data=[go.Heatmap(
		x=log_adp, y=log_thr, z=matrix, colorscale=colorscale,
		colorbar=dict(title=f'{label_html}<br>({unit_html})'),
		hovertemplate=(
			'ADP: 10<sup>%{x:.2f}</sup> µM<br>'
			'Thrombin: 10<sup>%{y:.2f}</sup> nM<br>'
			f'{label_html}: %{{z:.1f}} {unit_html}<extra></extra>'
		),
	)])
	fig.update_layout(
		title=dict(text=f'{label_html} — 2-D heatmap', x=0.5, xanchor='center'),
		xaxis=dict(
			title='ADP peak (µM)',
			tickmode='array',
			tickvals=log_adp.tolist(),
			ticktext=[f'{x:g}' for x in adp_grid],
		),
		yaxis=dict(
			title='Thrombin peak (nM)',
			tickmode='array',
			tickvals=log_thr.tolist(),
			ticktext=[f'{y:g}' for y in thr_grid],
		),
		margin=dict(l=40, r=0, t=50, b=40),
	)
	return fig


def register_callbacks(app: dash.Dash, out_path: str) -> None:
	"""Register Dose Response tab callbacks."""

	@app.callback(
		Output('dose-sweep-select', 'options'),
		Input('dose-refresh-interval', 'n_intervals'),
		prevent_initial_call=True,
	)
	def refresh_options(_n):
		return results.dose_sweep_options(out_path)

	@app.callback(
		Output('dose-surface-graph', 'figure'),
		Output('dose-heatmap-graph', 'figure'),
		Output('dose-sweep-meta', 'children'),
		Input('dose-sweep-select', 'value'),
		Input('dose-observable-select', 'value'),
	)
	def update_plots(sweep_dir, observable):
		if not sweep_dir or not observable:
			empty = _empty_figure('Select a sweep and observable above.')
			return empty, empty, ''
		npz_path = os.path.join(sweep_dir, 'sweep.npz')
		summary_path = os.path.join(sweep_dir, 'sweep_summary.json')
		if not (os.path.isfile(npz_path) and os.path.isfile(summary_path)):
			empty = _empty_figure(
				f'Sweep no longer on disk: {os.path.basename(sweep_dir)}')
			return empty, empty, ''

		data = np.load(npz_path)
		if observable not in data.files:
			empty = _empty_figure(
				f'Observable {observable!r} not in this sweep (older format?).')
			return empty, empty, ''

		adp_grid = data['adp_grid']
		thr_grid = data['thr_grid']
		matrix = data[observable]
		label, unit, mpl_cmap = _OBS_MAP[observable]
		plotly_cmap = _PLOTLY_CMAPS.get(mpl_cmap, 'Viridis')

		surface_fig = make_surface(adp_grid, thr_grid, matrix,
			label=label, unit=unit, colorscale=plotly_cmap)
		heatmap_fig = _make_heatmap(adp_grid, thr_grid, matrix,
			label=label, unit=unit, colorscale=plotly_cmap)

		with open(summary_path) as f:
			summary = json.load(f)
		meta = (
			f'Grid: {summary["grid_thr"]} × {summary["grid_adp"]} · '
			f'length {summary["length_sec"]} s · seed {summary["seed"]} · '
			f'ADP {summary["adp_range_uM"][0]:g}–{summary["adp_range_uM"][1]:g} µM · '
			f'thrombin {summary["thr_range_nM"][0]:g}–{summary["thr_range_nM"][1]:g} nM'
		)
		return surface_fig, heatmap_fig, meta
