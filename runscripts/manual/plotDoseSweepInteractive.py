"""One-shot prototype: interactive 3-D Plotly surface from a dose-sweep NPZ.

Reads `sweep.npz` from a runDoseSweep.py output directory and emits an
HTML file with a rotatable 3-D surface for each observable. Intended as
a quick preview of what the deferred webapp Dose Response tab (#46) will
look like; the embedded plotly.js makes the HTML standalone (no server,
just open in a browser).

Usage:
    PYTHONPATH=$PWD python runscripts/manual/plotDoseSweepInteractive.py \\
        <sweep_dir> [--observable peak_ip3_nM] [--out path.html]

Default observable is peak_ip3_nM (the cleanest dose-response surface);
pass `--observable all` to render every observable as a separate
<div> in the same HTML file.
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import plotly.graph_objects as go

from runscripts.manual.runDoseSweep import OBSERVABLES
from runscripts.manual.runPlateletSim import resolve_sim_path


def _resolve_sweep_dir(arg: str) -> str:
	"""Allow bare names ('dose_sweep_9x9_focus') or absolute/relative paths."""
	if os.path.isabs(arg) or os.path.isdir(arg):
		return arg
	return resolve_sim_path(arg)


def _make_surface(adp_grid: np.ndarray, thr_grid: np.ndarray,
		matrix: np.ndarray, label: str, unit: str, colorscale: str) -> go.Figure:
	"""Build a Plotly Surface figure for one observable."""
	# X = log10(ADP), Y = log10(thrombin). Plotly's 3-D log scale is
	# unreliable, so we plot log10-transformed coords on linear axes and
	# customise the tick labels back to original units.
	log_adp = np.log10(adp_grid)
	log_thr = np.log10(thr_grid)

	fig = go.Figure(data=[go.Surface(
		x=log_adp, y=log_thr, z=matrix,
		colorscale=colorscale,
		colorbar=dict(title=f'{label}<br>({unit})'),
		hovertemplate=(
			'ADP: 10^%{x:.2f} µM<br>'
			'Thrombin: 10^%{y:.2f} nM<br>'
			f'{label}: %{{z:.1f}} {unit}<extra></extra>'
		),
	)])

	fig.update_layout(
		title=dict(
			text=f'{label} — interactive dose-response surface',
			x=0.5, xanchor='center',
		),
		scene=dict(
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
			zaxis=dict(title=f'{label} ({unit})'),
			camera=dict(eye=dict(x=1.8, y=-1.6, z=0.9)),
		),
		margin=dict(l=0, r=0, t=60, b=0),
		height=620,
	)
	return fig


def render_html(sweep_dir: str, observable: str = 'peak_ip3_nM',
		out_path: str | None = None) -> str:
	"""Render one or more 3-D surfaces from a sweep NPZ into a single HTML."""
	npz_path = os.path.join(sweep_dir, 'sweep.npz')
	data = np.load(npz_path)
	adp_grid = data['adp_grid']
	thr_grid = data['thr_grid']

	obs_map = {k: (label, unit, cmap) for k, label, unit, cmap in OBSERVABLES}
	if observable == 'all':
		keys = [k for k, *_ in OBSERVABLES if k in data.files]
	else:
		if observable not in obs_map:
			raise ValueError(
				f'unknown observable {observable!r}; '
				f'choices: {", ".join(obs_map)} or "all"')
		keys = [observable]

	# Plotly cmap names differ from matplotlib's; map ours over.
	cmap_translate = {
		'viridis': 'Viridis',
		'plasma':  'Plasma',
		'cividis': 'Cividis',
		'magma':   'Magma',
	}

	html_parts = [
		'<!DOCTYPE html>',
		'<html lang="en"><head><meta charset="utf-8">',
		f'<title>Dose-response surfaces — {os.path.basename(sweep_dir)}</title>',
		'<style>body{font-family:sans-serif;margin:24px;color:#24292f}'
		'h1{font-size:20px}h2{font-size:14px;color:#57606a;margin-top:32px}'
		'.note{color:#57606a;font-size:13px;max-width:780px}</style>',
		'</head><body>',
		f'<h1>Dose-response surfaces — {os.path.basename(sweep_dir)}</h1>',
		'<p class="note">Drag to rotate, scroll to zoom. '
		'Surfaces are plotted on log10-transformed axes; '
		'tick labels are in original units.</p>',
	]

	include_plotlyjs: str | bool = 'cdn'  # first plot embeds, rest reuse
	for key in keys:
		label, unit, cmap = obs_map[key]
		fig = _make_surface(adp_grid, thr_grid, data[key],
			label=label, unit=unit,
			colorscale=cmap_translate.get(cmap, 'Viridis'))
		html_parts.append(f'<h2>{label} ({unit}) — key=<code>{key}</code></h2>')
		html_parts.append(fig.to_html(
			full_html=False, include_plotlyjs=include_plotlyjs,
			div_id=f'plot-{key}'))
		include_plotlyjs = False  # only embed once
	html_parts.append('</body></html>')

	out_path = out_path or os.path.join(sweep_dir, 'sweep_interactive.html')
	with open(out_path, 'w') as f:
		f.write('\n'.join(html_parts))
	return out_path


def main(argv: list[str] | None = None) -> None:
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('sweep_dir',
		help='Sweep output directory (under out/ or absolute).')
	parser.add_argument('--observable', default='peak_ip3_nM',
		help='Observable key, or "all" for one surface per observable. '
		'Default = peak_ip3_nM.')
	parser.add_argument('--out', dest='out_path', default=None,
		help='Output HTML path. Default = <sweep_dir>/sweep_interactive.html.')
	args = parser.parse_args(argv)

	sweep_dir = _resolve_sweep_dir(args.sweep_dir)
	out = render_html(sweep_dir, observable=args.observable,
		out_path=args.out_path)
	print(f'Wrote {out}')


if __name__ == '__main__':
	main()
