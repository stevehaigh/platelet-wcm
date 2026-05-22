"""2-D dose-response sweep over ADP × thrombin agonist peaks.

For each cell of an N × N log-spaced grid over (ADP_peak_uM, thrombin_peak_nM)
this script runs a platelet sim with the rest of the model at default,
harvests several observables from the CalciumTrace, and renders one heatmap
per observable plus a combined panel and a 3-D surface (peak Ca²⁺,
supplemental). Semantically maps "lower agonist peak ≡ stronger competitive
antagonism" — no antagonist species, no Kᵢ values.

Observables harvested per cell:
    peak_ca_nM         peak cytosolic Ca²⁺ (DTS-reservoir-dominated → ~binary)
    peak_ip3_nM        peak cytosolic IP3 (upstream cascade → graded)
    t_to_thresh_ca_s   time to first crossing of 250 nM cyt Ca²⁺ (rise time;
                       cells that never crossed are right-censored at the
                       sim length, indicating they did not fire)
    t_to_thresh_ip3_s  time to first crossing of 150 nM cyt IP3 (cascade
                       rise time, same censoring convention)
    auc_ca_nMs         Ca²⁺ AUC above resting baseline (integrated response)

Issue #45. Unblocked by #44 — agonist peaks are now live-readable via per-call
kwargs.

Usage:
    PYTHONPATH=$PWD python runscripts/manual/runDoseSweep.py [sim_outdir] \\
        [--length 200] [--grid 5] \\
        [--adp-min 0.1 --adp-max 10] \\
        [--thr-min 0.01 --thr-max 1] \\
        [--seed 0] [--keep-cell-output] [--no-log-to-shell] \\
        [--replot-only EXISTING_SWEEP_DIR]

`--replot-only` re-harvests + re-plots an existing sweep without rerunning
sims. Useful for iterating on the plot code (the cell sims may have been
pruned to CalciumTrace only — that is enough to recover every observable).

Outputs (under out/<sim_outdir>/):
    sweep.csv              one row per (i, j) cell × all observables
    sweep.npz              grids + matrices for replot
    sweep_summary.json     grid metadata + per-observable ranges
    sweep_<observable>.png one heatmap per observable
    sweep_panel.png        5-panel combined figure (headline)
    sweep_surface.png      3-D surface for peak Ca²⁺ (supplemental)
    cells/<i>_<j>/         per-cell sim output (CalciumTrace only by default)
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import textwrap
from dataclasses import dataclass, field
from datetime import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 — registers 3D projection
import numpy as np

from runscripts.manual.runPlateletSim import (
	resolve_sim_path,
	run_platelet_sim,
)
import wholecell.utils.filepath as fp
from wholecell.io.tablereader import TableReader


DEFAULT_GRID = 5
DEFAULT_LENGTH_SEC = 200
DEFAULT_ADP_RANGE_UM = (0.1, 10.0)
DEFAULT_THR_RANGE_NM = (0.01, 1.0)


# ── Observable registry ───────────────────────────────────────────────────
# (key, human label, unit, cmap). Order controls panel layout.
OBSERVABLES: tuple[tuple[str, str, str, str], ...] = (
	('peak_ca_nM',         r'Peak cytosolic Ca$^{2+}$',         'nM',           'viridis'),
	('peak_ip3_nM',        r'Peak cytosolic IP3',               'nM',           'plasma'),
	('t_to_thresh_ca_s',   r'Time to 250 nM cyt Ca$^{2+}$',     's',            'cividis'),
	('t_to_thresh_ip3_s',  r'Time to 150 nM cyt IP3',           's',            'cividis'),
	('auc_ca_nMs',         r'Ca$^{2+}$ AUC above rest',         r'nM$\cdot$s',  'magma'),
)

# Thresholds for the rise-time observables. Chosen to sit cleanly between
# the resting baseline and the saturated plateau, so the metric is well-
# defined: cells that fire cross during their rise; cells that don't fire
# are right-censored at the sim length (i.e. value = sim length means
# "never crossed", not "crossed at the very end").
_CA_RISE_THRESH_NM  = 250.0   # rest ~100 nM, saturated ~436 nM
_IP3_RISE_THRESH_NM = 150.0   # rest ~50 nM, saturated ~300 nM


def _first_crossing_s(y: np.ndarray, t: np.ndarray, threshold: float) -> float:
	"""First time `y` reaches `threshold`. If never, return the last
	sample time (right-censoring convention — caller is expected to know
	the sim length and treat boundary values as 'did not cross')."""
	above = y >= threshold
	if not above.any():
		return float(t[-1])
	return float(t[int(above.argmax())])


def harvest_cell(simout_dir: str) -> dict[str, float]:
	"""Read CalciumTrace and reduce to the 5 scalar observables."""
	reader = TableReader(os.path.join(simout_dir, 'CalciumTrace'))
	ca = reader.readColumn('ca_cyt_nM').flatten()
	ip3 = reader.readColumn('ip3_nM').flatten()
	t = np.arange(len(ca), dtype=float)   # 1-s timesteps
	ca_rest = float(ca[0])
	return {
		'peak_ca_nM':         float(ca.max()),
		'peak_ip3_nM':        float(ip3.max()),
		't_to_thresh_ca_s':   _first_crossing_s(ca,  t, _CA_RISE_THRESH_NM),
		't_to_thresh_ip3_s':  _first_crossing_s(ip3, t, _IP3_RISE_THRESH_NM),
		'auc_ca_nMs':         float(np.trapz(np.maximum(ca - ca_rest, 0.0), t)),
	}


@dataclass
class DoseSweep:
	"""Outcome of a 2-D dose-response sweep.

	`matrices[key][i, j]` is the value of observable `key` for the cell
	with thrombin = thr_grid[i] and ADP = adp_grid[j]. Rows = thrombin,
	columns = ADP — chosen so heatmaps with `pcolormesh(adp, thr, m)` put
	ADP on X and thrombin on Y.
	"""
	adp_grid: np.ndarray
	thr_grid: np.ndarray
	matrices: dict[str, np.ndarray] = field(default_factory=dict)
	length_sec: int = DEFAULT_LENGTH_SEC
	seed: int = 0

	def to_csv(self, path: str) -> None:
		keys = [k for k, *_ in OBSERVABLES if k in self.matrices]
		with open(path, 'w') as f:
			header = ['thr_idx', 'adp_idx', 'thrombin_peak_nM', 'adp_peak_uM'] + keys
			f.write(','.join(header) + '\n')
			for i, thr in enumerate(self.thr_grid):
				for j, adp in enumerate(self.adp_grid):
					row = [str(i), str(j), f'{float(thr):.6g}', f'{float(adp):.6g}']
					row.extend(f'{float(self.matrices[k][i, j]):.4f}' for k in keys)
					f.write(','.join(row) + '\n')

	def to_npz(self, path: str) -> None:
		np.savez(path,
			adp_grid=self.adp_grid, thr_grid=self.thr_grid,
			length_sec=self.length_sec, seed=self.seed,
			**self.matrices)

	def to_summary_json(self, path: str) -> None:
		obs_summary: dict[str, dict[str, object]] = {}
		for key in self.matrices:
			m = self.matrices[key]
			obs_summary[key] = {
				'min': float(m.min()),
				'max': float(m.max()),
				'argmax': [int(x) for x in np.unravel_index(int(np.argmax(m)), m.shape)],
				'argmin': [int(x) for x in np.unravel_index(int(np.argmin(m)), m.shape)],
			}
		payload = {
			'adp_range_uM': [float(self.adp_grid[0]), float(self.adp_grid[-1])],
			'thr_range_nM': [float(self.thr_grid[0]), float(self.thr_grid[-1])],
			'grid_adp': int(len(self.adp_grid)),
			'grid_thr': int(len(self.thr_grid)),
			'length_sec': int(self.length_sec),
			'seed': int(self.seed),
			'observables': obs_summary,
		}
		with open(path, 'w') as f:
			json.dump(payload, f, indent='\t')


# ── Driver ────────────────────────────────────────────────────────────────

def run_dose_sweep(out_path: str, grid: int = DEFAULT_GRID,
		length_sec: int = DEFAULT_LENGTH_SEC,
		adp_range_uM: tuple[float, float] = DEFAULT_ADP_RANGE_UM,
		thr_range_nM: tuple[float, float] = DEFAULT_THR_RANGE_NM,
		seed: int = 0, keep_cell_output: bool = False,
		log_to_shell: bool = True) -> DoseSweep:
	"""Run an N × N log-spaced sweep, return a DoseSweep with all matrices."""
	adp_grid = np.logspace(
		np.log10(adp_range_uM[0]), np.log10(adp_range_uM[1]), grid)
	thr_grid = np.logspace(
		np.log10(thr_range_nM[0]), np.log10(thr_range_nM[1]), grid)
	keys = [k for k, *_ in OBSERVABLES]
	matrices = {k: np.zeros((grid, grid)) for k in keys}

	cells_root = os.path.join(out_path, 'cells')
	fp.makedirs(cells_root)

	n_total = grid * grid
	for i, thr in enumerate(thr_grid):
		for j, adp in enumerate(adp_grid):
			cell_dir = os.path.join(cells_root, f'{i:02d}_{j:02d}')
			fp.makedirs(cell_dir)
			paths = run_platelet_sim(
				cell_dir,
				length_sec=length_sec,
				seed=seed,
				log_to_shell=False,
				adp_peak_uM=float(adp),
				thrombin_peak_nM=float(thr),
			)
			obs = harvest_cell(paths['sim_out_dir'])
			for k, v in obs.items():
				matrices[k][i, j] = v

			if log_to_shell:
				idx = i * grid + j + 1
				print(f'  [{idx:>3}/{n_total}] '
					f'adp={float(adp):>8.4f} µM, thr={float(thr):>8.4f} nM '
					f'→ peak Ca {obs["peak_ca_nM"]:>6.1f} nM, '
					f'peak IP3 {obs["peak_ip3_nM"]:>6.1f} nM')

			if not keep_cell_output:
				_prune_cell(cell_dir)

	return DoseSweep(adp_grid=adp_grid, thr_grid=thr_grid, matrices=matrices,
		length_sec=length_sec, seed=seed)


def replot_from_dir(sweep_dir: str) -> DoseSweep:
	"""Re-harvest + re-build a DoseSweep from an existing sweep_dir.

	Walks `sweep_dir/cells/<i>_<j>/` for surviving CalciumTrace data
	and reads the grid metadata from `sweep_summary.json`.
	"""
	summary_path = os.path.join(sweep_dir, 'sweep_summary.json')
	with open(summary_path) as f:
		summary = json.load(f)
	n_thr = int(summary['grid_thr'])
	n_adp = int(summary['grid_adp'])
	adp_grid = np.logspace(np.log10(summary['adp_range_uM'][0]),
		np.log10(summary['adp_range_uM'][1]), n_adp)
	thr_grid = np.logspace(np.log10(summary['thr_range_nM'][0]),
		np.log10(summary['thr_range_nM'][1]), n_thr)

	keys = [k for k, *_ in OBSERVABLES]
	matrices = {k: np.zeros((n_thr, n_adp)) for k in keys}

	for i in range(n_thr):
		for j in range(n_adp):
			cell_dir = os.path.join(sweep_dir, 'cells', f'{i:02d}_{j:02d}')
			simout_candidates = [d for d, _, _ in os.walk(cell_dir)
				if os.path.basename(d) == 'simOut']
			if not simout_candidates:
				raise FileNotFoundError(
					f'no simOut/ found under {cell_dir} — re-harvest impossible')
			obs = harvest_cell(simout_candidates[0])
			for k, v in obs.items():
				matrices[k][i, j] = v

	return DoseSweep(adp_grid=adp_grid, thr_grid=thr_grid, matrices=matrices,
		length_sec=int(summary['length_sec']), seed=int(summary['seed']))


def _prune_cell(cell_dir: str) -> None:
	"""Delete everything in a cell directory except simOut/CalciumTrace/."""
	for root, _, _ in os.walk(cell_dir):
		if os.path.basename(root) == 'simOut':
			for entry in os.listdir(root):
				full = os.path.join(root, entry)
				if entry == 'CalciumTrace':
					continue
				if os.path.isdir(full):
					shutil.rmtree(full)
				else:
					os.remove(full)
	for d in ('kb', 'metadata'):
		full = os.path.join(cell_dir, d)
		if os.path.isdir(full):
			shutil.rmtree(full)


# ── Plot helpers ──────────────────────────────────────────────────────────

def _log_edges(centers: np.ndarray) -> np.ndarray:
	"""Cell edges in log space so pcolormesh quads are centered on grid points."""
	log_c = np.log10(centers)
	mid = 0.5 * (log_c[:-1] + log_c[1:])
	edges = np.concatenate([
		[2.0 * log_c[0] - mid[0]],
		mid,
		[2.0 * log_c[-1] - mid[-1]],
	])
	return 10.0 ** edges


def _draw_heatmap(ax: plt.Axes, sweep: DoseSweep, key: str,
		label: str, unit: str, cmap: str) -> matplotlib.collections.QuadMesh:
	"""Render a heatmap for one observable into the given axes."""
	matrix = sweep.matrices[key]
	adp_edges = _log_edges(sweep.adp_grid)
	thr_edges = _log_edges(sweep.thr_grid)
	pcm = ax.pcolormesh(adp_edges, thr_edges, matrix,
		cmap=cmap, shading='flat')
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel('ADP peak (µM)')
	ax.set_ylabel('Thrombin peak (nM)')
	ax.set_title(f'{label} ({unit})', fontsize=10)

	# Annotate cells. White text on dark cells, black on light.
	threshold = 0.5 * (matrix.min() + matrix.max())
	max_abs = max(abs(matrix.min()), abs(matrix.max()), 1.0)
	for i, thr in enumerate(sweep.thr_grid):
		for j, adp in enumerate(sweep.adp_grid):
			val = matrix[i, j]
			color = 'white' if val < threshold else 'black'
			fmt = '{:.0f}' if max_abs >= 10 else '{:.2f}'
			ax.text(adp, thr, fmt.format(val),
				ha='center', va='center', fontsize=7, color=color)
	return pcm


def plot_observable_heatmap(sweep: DoseSweep, key: str, label: str,
		unit: str, cmap: str, png_path: str,
		caption: str | None = None) -> None:
	"""Single-observable heatmap PNG.

	`caption`, if provided, is rendered as a wrapped block below the
	heatmap — use it to explain what the observable means for readers
	who open the standalone figure without the surrounding prose.
	"""
	if caption:
		# Pre-wrap the caption so matplotlib doesn't have to guess. Lines
		# are sized to fit in the 7.5"-wide figure at 8 pt font.
		caption_wrapped = textwrap.fill(caption, width=100)
		n_lines = caption_wrapped.count('\n') + 1
		fig_h = 5.5 + 0.22 * n_lines + 0.3   # axes + caption rows + padding
		bottom_frac = (0.22 * n_lines + 0.3) / fig_h
	else:
		fig_h = 5.5
		bottom_frac = 0.0

	fig, ax = plt.subplots(figsize=(7.5, fig_h))
	pcm = _draw_heatmap(ax, sweep, key, label, unit, cmap)
	fig.colorbar(pcm, ax=ax, label=f'{label} ({unit})')
	fig.suptitle(
		f'Dose-response: {label} vs ADP × thrombin '
		f'(length {sweep.length_sec} s, seed {sweep.seed})', fontsize=11)
	# Tighten first so the axes sit cleanly, then place the caption in
	# the reserved bottom band as a separate figure-level text box.
	fig.tight_layout(rect=(0, bottom_frac, 1, 1))
	if caption:
		fig.text(0.5, bottom_frac * 0.55, caption_wrapped,
			ha='center', va='center', fontsize=8.5, color='#24292f',
			bbox=dict(boxstyle='round,pad=0.5',
				facecolor='#f6f8fa', edgecolor='#d0d7de'))
	fig.savefig(png_path, dpi=140, bbox_inches='tight')
	plt.close(fig)


def plot_panel(sweep: DoseSweep, png_path: str) -> None:
	"""5-panel combined figure: one heatmap per observable."""
	# 5 observables → 3×2 layout, leave the 6th slot empty.
	fig, axes = plt.subplots(3, 2, figsize=(13, 14))
	flat_axes = axes.flatten()
	for ax, (key, label, unit, cmap) in zip(flat_axes, OBSERVABLES):
		pcm = _draw_heatmap(ax, sweep, key, label, unit, cmap)
		fig.colorbar(pcm, ax=ax, label=unit)
	# Hide the trailing unused axes.
	for ax in flat_axes[len(OBSERVABLES):]:
		ax.set_visible(False)
	fig.suptitle(
		f'Dose-response panel — {sweep.length_sec} s sims, seed {sweep.seed}, '
		f'grid {len(sweep.thr_grid)}×{len(sweep.adp_grid)}', fontsize=12)
	fig.tight_layout()
	fig.savefig(png_path, dpi=140, bbox_inches='tight')
	plt.close(fig)


def plot_surface(sweep: DoseSweep, key: str, label: str, unit: str,
		png_path: str) -> None:
	"""3-D surface of one observable. Coords are log10-transformed."""
	matrix = sweep.matrices[key]
	A, T = np.meshgrid(np.log10(sweep.adp_grid), np.log10(sweep.thr_grid))
	fig = plt.figure(figsize=(8.5, 6))
	ax = fig.add_subplot(111, projection='3d')
	surf = ax.plot_surface(A, T, matrix,
		cmap='viridis', edgecolor='gray', linewidth=0.3, antialiased=True)

	ax.set_xticks(np.log10(sweep.adp_grid))
	ax.set_xticklabels([f'{x:g}' for x in sweep.adp_grid], fontsize=8)
	ax.set_yticks(np.log10(sweep.thr_grid))
	ax.set_yticklabels([f'{y:g}' for y in sweep.thr_grid], fontsize=8)
	ax.set_xlabel('ADP peak (µM)')
	ax.set_ylabel('Thrombin peak (nM)')
	ax.set_zlabel(f'{label} ({unit})')
	ax.set_title(f'Dose-response surface — {label}', fontsize=11)
	fig.colorbar(surf, shrink=0.55, label=unit, pad=0.1)
	fig.tight_layout()
	fig.savefig(png_path, dpi=140, bbox_inches='tight')
	plt.close(fig)


def write_all_outputs(sweep: DoseSweep, out_path: str) -> None:
	"""Write CSV, NPZ, summary JSON, all heatmaps, panel, and Ca surface."""
	sweep.to_csv(os.path.join(out_path, 'sweep.csv'))
	sweep.to_npz(os.path.join(out_path, 'sweep.npz'))
	sweep.to_summary_json(os.path.join(out_path, 'sweep_summary.json'))
	for key, label, unit, cmap in OBSERVABLES:
		plot_observable_heatmap(sweep, key, label, unit, cmap,
			os.path.join(out_path, f'sweep_{key}.png'),
			caption=_OBSERVABLE_CAPTIONS.get(key))
	plot_panel(sweep, os.path.join(out_path, 'sweep_panel.png'))
	# 3-D surface for the canonical Ca²⁺ observable only.
	plot_surface(sweep, 'peak_ca_nM', r'Peak cytosolic Ca$^{2+}$', 'nM',
		os.path.join(out_path, 'sweep_surface.png'))


# Per-observable captions rendered below each standalone heatmap. Keep
# concise — they share figure space with the colourbar. The panel figure
# omits these (axes are already tight) and relies on the prose caption
# the user provides at the use site.
_OBSERVABLE_CAPTIONS: dict[str, str] = {
	'peak_ca_nM': (
		r'Peak cytosolic [Ca$^{2+}$] over the 200 s window. The '
		r'DTS-reservoir release through IP3R is essentially '
		r'all-or-nothing: once the cell crosses threshold, the peak '
		r'locks at $\sim$436 nM regardless of stimulus strength. '
		r'Sub-threshold cells (bottom-left) sit near the 100 nM resting '
		r'baseline.'
	),
	'auc_ca_nMs': (
		r'Time-integrated cytosolic [Ca$^{2+}$] above the resting '
		r'baseline: $\int \max(\mathrm{Ca}_{cyt}(t) - \mathrm{Ca}_{rest}, 0)\,dt$ '
		r'over the 200 s window, units nM$\cdot$s. This is the closest '
		r'single-number proxy for how much downstream Ca$^{2+}$-dependent '
		r'machinery (PKC, CaM-kinases, granule-fusion machinery) was '
		r'activated $-$ sustained elevation matters more than peak '
		r'height. The cell makes a binary spike decision (see peak '
		r'Ca$^{2+}$) but the integrated dose is graded with agonist '
		r'input via spike duration.'
	),
	't_to_thresh_ca_s': (
		r'First time cyt [Ca$^{2+}$] crosses 250 nM (threshold sits '
		r'between $\sim$100 nM rest and the $\sim$436 nM saturated '
		r'plateau). Captures rise time without the flat-top artefact of '
		r'an argmax-based metric. Cells that never fired are right-'
		r'censored at the 200 s sim length, i.e. a value of 200 means '
		r'"did not cross", not "crossed at the very end".'
	),
	't_to_thresh_ip3_s': (
		r'First time cyt [IP3] crosses 150 nM (between $\sim$50 nM rest '
		r'and the $\sim$300 nM saturated plateau). Same censoring '
		r'convention as the Ca$^{2+}$ rise-time observable above.'
	),
}


# ── CLI ───────────────────────────────────────────────────────────────────

def _build_parser() -> argparse.ArgumentParser:
	parser = argparse.ArgumentParser(
		description='2-D dose-response sweep over ADP × thrombin agonist peaks.')
	parser.add_argument('sim_outdir', nargs='?', default=None,
		help='Output dir under out/. Default = dose_sweep_<timestamp>.')
	parser.add_argument('--length', '--length-sec', dest='length_sec',
		type=int, default=DEFAULT_LENGTH_SEC,
		help=f'Sim length per cell (s). Default = {DEFAULT_LENGTH_SEC}.')
	parser.add_argument('--grid', type=int, default=DEFAULT_GRID,
		help=f'N × N grid size. Default = {DEFAULT_GRID}.')
	parser.add_argument('--adp-min', type=float, default=DEFAULT_ADP_RANGE_UM[0],
		help=f'Min ADP peak (µM). Default = {DEFAULT_ADP_RANGE_UM[0]}.')
	parser.add_argument('--adp-max', type=float, default=DEFAULT_ADP_RANGE_UM[1],
		help=f'Max ADP peak (µM). Default = {DEFAULT_ADP_RANGE_UM[1]}.')
	parser.add_argument('--thr-min', type=float, default=DEFAULT_THR_RANGE_NM[0],
		help=f'Min thrombin peak (nM). Default = {DEFAULT_THR_RANGE_NM[0]}.')
	parser.add_argument('--thr-max', type=float, default=DEFAULT_THR_RANGE_NM[1],
		help=f'Max thrombin peak (nM). Default = {DEFAULT_THR_RANGE_NM[1]}.')
	parser.add_argument('--seed', type=int, default=0,
		help='RNG seed (shared across cells). Default = 0.')
	parser.add_argument('--keep-cell-output', action='store_true',
		help='Keep the full simOut for each cell (default: prune to '
		'CalciumTrace only — ~50 KB/cell instead of ~5 MB/cell).')
	parser.add_argument('--no-log-to-shell', dest='log_to_shell',
		action='store_false', help='Suppress per-cell progress prints.')
	parser.add_argument('--replot-only', dest='replot_only', default=None,
		help='Re-harvest + re-plot an existing sweep directory (under out/) '
		'without rerunning the sims. Cell CalciumTrace data must still exist.')
	parser.set_defaults(log_to_shell=True)
	return parser


def main(argv: list[str] | None = None) -> None:
	args = _build_parser().parse_args(argv)

	if args.replot_only:
		sweep_dir = resolve_sim_path(args.replot_only)
		print(f'Re-harvesting + re-plotting existing sweep at {sweep_dir}')
		sweep = replot_from_dir(sweep_dir)
		write_all_outputs(sweep, sweep_dir)
		_print_summary(sweep, sweep_dir)
		return

	sim_outdir = args.sim_outdir or 'dose_sweep_{}'.format(
		datetime.now().strftime('%Y%m%d.%H%M%S'))
	sim_path = resolve_sim_path(sim_outdir)
	fp.makedirs(sim_path)

	print(f'Dose sweep — {args.grid}×{args.grid} = {args.grid**2} cells '
		f'@ {args.length_sec} s each')
	print(f'  ADP      range: {args.adp_min} … {args.adp_max} µM (log-spaced)')
	print(f'  Thrombin range: {args.thr_min} … {args.thr_max} nM (log-spaced)')
	print(f'  Output: {sim_path}')
	print()

	sweep = run_dose_sweep(
		sim_path,
		grid=args.grid,
		length_sec=args.length_sec,
		adp_range_uM=(args.adp_min, args.adp_max),
		thr_range_nM=(args.thr_min, args.thr_max),
		seed=args.seed,
		keep_cell_output=args.keep_cell_output,
		log_to_shell=args.log_to_shell,
	)
	write_all_outputs(sweep, sim_path)
	_print_summary(sweep, sim_path)


def _print_summary(sweep: DoseSweep, out_path: str) -> None:
	print()
	print('Observable ranges:')
	for key, label, unit, _ in OBSERVABLES:
		m = sweep.matrices[key]
		print(f'  {label:<24} {m.min():>9.2f} … {m.max():>9.2f} {unit}')
	print()
	print(f'Outputs in: {out_path}')


if __name__ == '__main__':
	main()
