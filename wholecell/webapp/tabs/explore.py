"""Explore tab — browse and compare analysis plots."""

from __future__ import annotations

import base64
import os

import dash
from dash import dcc, html
from dash.dependencies import Input, Output

from wholecell.webapp import results


def layout(out_path: str) -> html.Div:
	"""Create the Explore tab layout."""

	run_options = results.explore_run_options(out_path)

	return html.Div(children=[
		html.Div(className='grid-2', style={'marginBottom': '15px'}, children=[
			html.Div([
				html.Label('Left run'),
				dcc.Dropdown(
					id='explore-left-run',
					options=run_options,
					value=run_options[0]['value'] if run_options else None,
				),
			]),
			html.Div([
				html.Label('Right run (for comparison)'),
				dcc.Dropdown(
					id='explore-right-run',
					options=run_options,
					placeholder='Select to compare...',
				),
			]),
		]),

		html.Div(id='explore-plots-container'),

		# Periodic refresh of run-options dropdowns so newly-completed
		# webapp jobs appear without restarting the server.
		dcc.Interval(id='explore-refresh-interval', interval=5000, n_intervals=0),
	])


def register_callbacks(app: dash.Dash, out_path: str) -> None:
	"""Register Explore tab callbacks."""

	@app.callback(
		Output('explore-left-run', 'options', allow_duplicate=True),
		Output('explore-right-run', 'options', allow_duplicate=True),
		Input('explore-refresh-interval', 'n_intervals'),
		prevent_initial_call=True,
	)
	def refresh_run_options(n_intervals):
		"""Re-scan out/ for new runs and update both Explore dropdowns.

		Polls every 5 s. User's current selection is preserved by Dash
		so long as the same value is still present in the refreshed list.
		"""
		new_options = results.explore_run_options(out_path)
		return new_options, new_options

	@app.callback(
		Output('explore-plots-container', 'children'),
		Input('explore-left-run', 'value'),
		Input('explore-right-run', 'value'),
	)
	def update_plots(left_value, right_value):
		if not left_value:
			return html.P('Select a run to view plots.',
				style={'color': '#888', 'fontStyle': 'italic'})

		if '|' not in left_value:
			return html.P('Invalid run selection.', style={'color': '#c00'})
		left_dir, left_variant = left_value.split('|', 1)
		left_images = results.find_plot_images(left_dir, left_variant)

		right_images = {}
		if right_value and '|' in right_value:
			right_dir, right_variant = right_value.split('|', 1)
			for img in results.find_plot_images(right_dir, right_variant):
				right_images[img['name']] = img

		if not left_images:
			return html.P(
				'No analysis plots found. Run analysis first '
				'(e.g., python runscripts/manual/analysisPlatelet.py).',
				style={'color': '#888', 'fontStyle': 'italic'})

		comparing = bool(right_value)
		children = []

		for img in left_images:
			left_encoded = _encode_image(img['path'])
			if left_encoded is None:
				continue

			row_children = [
				html.Div(
					style={'flex': '1', 'textAlign': 'center'},
					children=[html.Img(
						src=left_encoded,
						style={'maxWidth': '100%', 'border': '1px solid #ddd'},
					)],
				),
			]

			if comparing:
				right_img = right_images.get(img['name'])
				if right_img:
					right_encoded = _encode_image(right_img['path'])
					row_children.append(html.Div(
						style={'flex': '1', 'textAlign': 'center'},
						children=[html.Img(
							src=right_encoded,
							style={'maxWidth': '100%', 'border': '1px solid #ddd'},
						)] if right_encoded else [
							html.P('—', style={'color': '#ccc'})
						],
					))
				else:
					row_children.append(html.Div(
						style={'flex': '1', 'textAlign': 'center'},
						children=[html.P('No matching plot',
							style={'color': '#ccc', 'fontStyle': 'italic'})],
					))

			children.append(html.Div(
				style={'marginBottom': '25px'},
				children=[
					html.H4(img['name'],
						style={'marginBottom': '5px', 'color': '#333'}),
					html.Div(
						style={'display': 'flex', 'gap': '10px'},
						children=row_children,
					),
				],
			))

		return children


MAX_IMAGE_SIZE = 10 * 1024 * 1024  # 10 MB


def _encode_image(path: str) -> str | None:
	"""Read an image file and return a base64-encoded data URI."""

	if not os.path.isfile(path):
		return None
	if os.path.getsize(path) > MAX_IMAGE_SIZE:
		return None
	ext = os.path.splitext(path)[1].lower()
	mime = {'png': 'image/png', 'svg': 'image/svg+xml'}.get(ext.lstrip('.'), 'image/png')
	try:
		with open(path, 'rb') as f:
			encoded = base64.b64encode(f.read()).decode('utf-8')
		return f'data:{mime};base64,{encoded}'
	except (OSError, ValueError):
		return None
