"""Sketch wireframes for the five graphical-abstract concepts.

These are layout references for the BioRender drawing, not finished art.
Boxes/arrows/labels stand in for the platelet cartoons; embedded real
plots (from reports/figures/v0.5/dose-sweep-9x9-transition.npz) show where
the actual heatmaps go and at what proportion. Output goes to
`reports/figures/v0.5/graphical-abstract-concept-{1..5}.png`.

Usage:
    PYTHONPATH=$PWD python runscripts/manual/sketchGraphicalAbstracts.py
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
FIGS_DIR = REPO_ROOT / 'reports' / 'figures'
NPZ_PATH = FIGS_DIR / 'dose-sweep-9x9-transition.npz'


# ── Helpers ───────────────────────────────────────────────────────────────

def _box(ax, x, y, w, h, label, facecolor='#f0f4f8', edgecolor='#3d4147',
		fontsize=9, lw=1.4, alpha=0.95):
	"""Draw a rounded rectangle with a centred label."""
	rect = patches.FancyBboxPatch(
		(x, y), w, h, boxstyle='round,pad=0.02',
		linewidth=lw, edgecolor=edgecolor, facecolor=facecolor, alpha=alpha)
	ax.add_patch(rect)
	ax.text(x + w / 2, y + h / 2, label,
		ha='center', va='center', fontsize=fontsize)


def _arrow(ax, x1, y1, x2, y2, label='', color='#3d4147', lw=1.4,
		fontsize=8, label_offset=(0, 0.015)):
	ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
		arrowprops=dict(arrowstyle='->', color=color, lw=lw))
	if label:
		ax.text((x1 + x2) / 2 + label_offset[0],
			(y1 + y2) / 2 + label_offset[1],
			label, ha='center', va='center',
			fontsize=fontsize, color=color, style='italic')


def _platelet_cartoon(ax, cx, cy, r=0.18, label='platelet',
		dts_color='#fde68a', cyto_color='#f0f4f8', edge='#3d4147'):
	"""Draw a stylised platelet cross-section. Returns the patch bounds
	so caller can layer receptors / arrows on top."""
	# Outer plasma membrane
	pm = patches.Ellipse((cx, cy), 2 * r, 1.4 * r,
		linewidth=1.6, edgecolor=edge, facecolor=cyto_color)
	ax.add_patch(pm)
	# Inner DTS reservoir (eccentric blob)
	dts = patches.Ellipse((cx - 0.02 * r, cy + 0.03 * r),
		0.95 * r, 0.55 * r,
		linewidth=1.2, edgecolor='#b45309', facecolor=dts_color)
	ax.add_patch(dts)
	ax.text(cx - 0.02 * r, cy + 0.03 * r, 'DTS',
		ha='center', va='center', fontsize=8, color='#b45309')
	# Mitochondrion (small purple bean)
	mito = patches.Ellipse((cx + 0.55 * r, cy - 0.3 * r),
		0.32 * r, 0.18 * r,
		linewidth=1.0, edgecolor='#7c3aed', facecolor='#ede9fe')
	ax.add_patch(mito)
	ax.text(cx + 0.55 * r, cy - 0.3 * r, 'mito',
		ha='center', va='center', fontsize=6.5, color='#7c3aed')
	# Label
	ax.text(cx, cy - 0.78 * r, label,
		ha='center', va='top', fontsize=9, fontweight='medium')
	return cx - r, cy - 0.7 * r, 2 * r, 1.4 * r


def _embed_heatmap(fig, sweep, key, x, y, w, h, title=None, cmap='viridis'):
	"""Embed a real heatmap from the sweep NPZ at figure-fraction coords.

	Uses imshow (not pcolormesh) since the wireframes only need the
	visual gist — the precise log-axis is irrelevant on a thumbnail.
	`origin='lower'` matches the row-thrombin convention used by the
	main heatmap renderer.
	"""
	ax = fig.add_axes((x, y, w, h))
	ax.imshow(sweep[key], cmap=cmap, aspect='auto',
		origin='lower', interpolation='nearest')
	ax.set_xticks([])
	ax.set_yticks([])
	if title:
		ax.set_title(title, fontsize=8)
	for spine in ax.spines.values():
		spine.set_linewidth(0.6)
	return ax


def _setup_axes(ax):
	ax.set_xlim(0, 1)
	ax.set_ylim(0, 1)
	ax.set_aspect('auto')
	ax.set_xticks([])
	ax.set_yticks([])
	for spine in ax.spines.values():
		spine.set_visible(False)


def _wireframe_note(fig, text):
	fig.text(0.5, 0.01,
		f'(wireframe — boxes replace BioRender cartoons; '
		f'live plots show placement & proportion)  · {text}',
		ha='center', va='bottom', fontsize=7,
		color='#57606a', fontstyle='italic')


def _load_sweep() -> dict[str, np.ndarray]:
	if not NPZ_PATH.is_file():
		raise FileNotFoundError(
			f'Need {NPZ_PATH.relative_to(REPO_ROOT)} from a prior dose sweep.')
	return dict(np.load(NPZ_PATH))


# ── Concept 1 — Three-layer signal shaping ────────────────────────────────

def concept1(sweep, out: Path) -> None:
	fig = plt.figure(figsize=(10, 9))
	ax = fig.add_axes((0.04, 0.06, 0.92, 0.88))
	_setup_axes(ax)

	fig.suptitle('Concept 1 — Three-layer signal shaping',
		fontsize=13, fontweight='bold', y=0.97)

	# Layer titles down the left
	layers = [
		(0.78, 'GRADED INPUT',  'P2Y1 / PAR1 / PAR4 saturate smoothly'),
		(0.49, 'BINARY SWITCH', 'IP3R / DTS releases all-or-nothing'),
		(0.20, 'GRADED OUTPUT', 'Integrated Ca²⁺ (AUC) scales with duration'),
	]
	for y, title, sub in layers:
		ax.text(0.025, y + 0.08, title, ha='left', va='center',
			fontsize=10, fontweight='bold', color='#3d4147')
		ax.text(0.025, y + 0.045, sub, ha='left', va='center',
			fontsize=8, color='#57606a', style='italic')

	# Top row — graded input
	_box(ax, 0.27, 0.74, 0.18, 0.10, 'agonist\ngradient\n(BioRender)',
		facecolor='#fef3c7')
	_arrow(ax, 0.46, 0.79, 0.55, 0.79)
	_box(ax, 0.55, 0.74, 0.18, 0.10, 'platelet +\nGPCR receptors\n(BioRender)',
		facecolor='#dbeafe')

	# Middle row — IP3R / DTS
	_box(ax, 0.27, 0.45, 0.18, 0.10, 'IP3 above\nthreshold ?\n(BioRender)',
		facecolor='#dbeafe')
	_arrow(ax, 0.46, 0.50, 0.55, 0.50)
	_box(ax, 0.55, 0.45, 0.18, 0.10,
		'DTS released\n(empty reservoir)\n(BioRender)',
		facecolor='#fde68a')

	# Bottom row — integrated effect
	_box(ax, 0.27, 0.16, 0.18, 0.10, 'cyt Ca²⁺ trace\nfor 200 s\n(BioRender)',
		facecolor='#dbeafe')
	_arrow(ax, 0.46, 0.21, 0.55, 0.21)
	_box(ax, 0.55, 0.16, 0.18, 0.10,
		'effectors\n(PKC, granule\nfusion, ...) ',
		facecolor='#bbf7d0')

	# Right column — actual heatmaps from the sweep
	_embed_heatmap(fig, sweep, 'peak_ip3_nM',
		0.76, 0.72, 0.18, 0.16, title='peak IP3', cmap='plasma')
	_embed_heatmap(fig, sweep, 'peak_ca_nM',
		0.76, 0.43, 0.18, 0.16, title='peak Ca²⁺', cmap='viridis')
	_embed_heatmap(fig, sweep, 'auc_ca_nMs',
		0.76, 0.14, 0.18, 0.16, title='Ca²⁺ AUC', cmap='magma')

	# Vertical rhythm arrows on the far left
	for y_from, y_to in ((0.72, 0.57), (0.43, 0.28)):
		_arrow(ax, 0.18, y_from, 0.18, y_to, lw=1.6, color='#9aa1a8')

	_wireframe_note(fig,
		'visual rhythm: ramp → cliff → ramp matches the panel-row heatmaps')
	fig.savefig(out, dpi=140, bbox_inches='tight')
	plt.close(fig)


# ── Concept 2 — Platelet cross-section with dose-surface inset ────────────

def concept2(sweep, out: Path) -> None:
	fig = plt.figure(figsize=(10, 7))
	ax = fig.add_axes((0.04, 0.06, 0.92, 0.88))
	_setup_axes(ax)

	fig.suptitle('Concept 2 — Single-cell biology constrains population dose-response',
		fontsize=13, fontweight='bold', y=0.97)

	# Big platelet on the left
	pl_cx, pl_cy = 0.30, 0.50
	_platelet_cartoon(ax, pl_cx, pl_cy, r=0.22,
		label='platelet cross-section (BioRender)')

	# Inflow / outflow arrows
	_arrow(ax, 0.05, 0.70, 0.16, 0.62,
		label='Ca²⁺ in (SOCE)', color='#0369a1', lw=1.6)
	_arrow(ax, 0.16, 0.36, 0.05, 0.30,
		label='Ca²⁺ out (PMCA)', color='#b45309', lw=1.6)
	_arrow(ax, 0.30, 0.60, 0.30, 0.50,
		label='IP3R release', color='#7c3aed', lw=1.6)
	_arrow(ax, 0.34, 0.50, 0.34, 0.60,
		label='SERCA refill', color='#15803d', lw=1.4)

	# Receptor counts callout
	ax.text(pl_cx, pl_cy - 0.32,
		'PAR1: 2 500  ·  PAR4: 500\n'
		'P2Y1: 150  ·  P2X1: 1 000\n'
		'DTS: 4.3 % of cytosol',
		ha='center', va='top', fontsize=8, color='#3d4147',
		bbox=dict(boxstyle='round,pad=0.4', facecolor='#f6f8fa',
			edgecolor='#d0d7de'))

	# Dose surface inset (real data, top right)
	_embed_heatmap(fig, sweep, 'auc_ca_nMs',
		0.62, 0.50, 0.30, 0.32, title='Ca²⁺ AUC vs ADP × thrombin',
		cmap='magma')

	# Caption tying the two halves
	ax.text(0.77, 0.32,
		'Single-cell ion machinery (left)\n'
		'constrains the integrated dose-\n'
		'response across the agonist plane\n'
		'(right). DTS reservoir size sets the\n'
		'peak; sustained IP3 sets the duration.',
		ha='center', va='center', fontsize=9, color='#3d4147',
		bbox=dict(boxstyle='round,pad=0.5',
			facecolor='#fff', edgecolor='#d0d7de'))

	_arrow(ax, 0.52, 0.55, 0.62, 0.65, color='#9aa1a8', lw=1.6)

	_wireframe_note(fig,
		'one platelet + one inset = ~30 min BioRender; tells the story in one image')
	fig.savefig(out, dpi=140, bbox_inches='tight')
	plt.close(fig)


# ── Concept 3 — Whole-cell-model funnel ───────────────────────────────────

def concept3(sweep, out: Path) -> None:
	fig = plt.figure(figsize=(10, 8))
	ax = fig.add_axes((0.04, 0.06, 0.92, 0.88))
	_setup_axes(ax)

	fig.suptitle('Concept 3 — Whole-cell-model funnel (12 mechanisms end-to-end)',
		fontsize=13, fontweight='bold', y=0.97)

	# Three-tier funnel
	tiers = [
		(0.78, 0.08, 0.78, 0.16,
			'INPUTS', '#fef3c7',
			'thrombin (nM)  ·  ADP (µM)  ·  extracellular ATP (µM)'),
		(0.50, 0.20, 0.54, 0.16,
			'CASCADE', '#dbeafe',
			'P2Y1 / PAR1 / PAR4 / P2X1  →  Gαq  →  PLCβ  →  PIP2 / IP3 / DAG'),
		(0.30, 0.36, 0.22, 0.16,
			'STORE-OPERATED Ca²⁺', '#fde68a',
			'IP3R / SERCA / DTS / SOCE / PMCA / NCX / MCU / NCLX / CaM / buffers'),
		(0.10, 0.55, 0.10, 0.16,
			'EFFECTORS', '#bbf7d0',
			'granule release  ·  integrin αIIbβ3  ·  shape change  ·  PS exposure'),
	]
	for cy, h, box_w_half, box_h, label, color, contents in tiers:
		x = 0.5 - box_w_half
		_box(ax, x, cy, 2 * box_w_half, box_h,
			f'{label}\n\n{contents}', facecolor=color, fontsize=8.5)

	# Down arrows between tiers
	for y_from in (0.78, 0.50, 0.30):
		_arrow(ax, 0.5, y_from, 0.5, y_from - 0.06, color='#3d4147', lw=2)

	# Dose-sweep callout in the bottom right
	ax.text(0.78, 0.20,
		'This whole pathway,\nswept across two agonist axes:',
		ha='center', va='center', fontsize=9, color='#3d4147')
	_embed_heatmap(fig, sweep, 'auc_ca_nMs',
		0.66, 0.04, 0.22, 0.13, title='AUC Ca²⁺', cmap='magma')

	_wireframe_note(fig,
		'positions the model as a comprehensive platform, not a single insight')
	fig.savefig(out, dpi=140, bbox_inches='tight')
	plt.close(fig)


# ── Concept 4 — Peak vs AUC pharmacology angle ────────────────────────────

def concept4(sweep, out: Path) -> None:
	fig = plt.figure(figsize=(11, 6))
	ax = fig.add_axes((0.04, 0.06, 0.92, 0.86))
	_setup_axes(ax)

	fig.suptitle('Concept 4 — Why peak Ca²⁺ hides the dose-graded response',
		fontsize=13, fontweight='bold', y=0.97)

	# Centre platelet
	_platelet_cartoon(ax, 0.50, 0.55, r=0.10,
		label='platelet (BioRender)')

	# Thought bubbles — two questions
	ax.text(0.40, 0.30,
		'"Did I fire?"',
		ha='center', va='center', fontsize=11, fontstyle='italic',
		color='#3d4147',
		bbox=dict(boxstyle='round,pad=0.4',
			facecolor='#fef3c7', edgecolor='#3d4147'))
	ax.text(0.60, 0.30,
		'"How hard?"',
		ha='center', va='center', fontsize=11, fontstyle='italic',
		color='#3d4147',
		bbox=dict(boxstyle='round,pad=0.4',
			facecolor='#bbf7d0', edgecolor='#3d4147'))

	# Arrows from platelet to left/right thought
	_arrow(ax, 0.46, 0.50, 0.42, 0.36, color='#9aa1a8', lw=1.6)
	_arrow(ax, 0.54, 0.50, 0.58, 0.36, color='#9aa1a8', lw=1.6)

	# Left: peak Ca²⁺ heatmap (binary)
	_embed_heatmap(fig, sweep, 'peak_ca_nM',
		0.04, 0.20, 0.30, 0.50,
		title='peak Ca²⁺ — saturates ≈ 436 nM', cmap='viridis')

	# Right: AUC heatmap (graded)
	_embed_heatmap(fig, sweep, 'auc_ca_nMs',
		0.66, 0.20, 0.30, 0.50,
		title='Ca²⁺ AUC — 8× dynamic range', cmap='magma')

	# Caption strap below
	ax.text(0.50, 0.08,
		'Same sweep, two different stories. Peak amplitude is reservoir-determined and binary; '
		'time-integrated exposure (AUC) captures the dose-graded biology that downstream '
		'effectors actually integrate. ',
		ha='center', va='center', fontsize=8.5, color='#3d4147',
		wrap=True)

	_wireframe_note(fig,
		'pharmacology / drug-screening framing — methodologically pointed')
	fig.savefig(out, dpi=140, bbox_inches='tight')
	plt.close(fig)


# ── Concept 5 — Before / after architecture ───────────────────────────────

def concept5(sweep, out: Path) -> None:
	fig = plt.figure(figsize=(11, 6.5))
	ax = fig.add_axes((0.04, 0.06, 0.92, 0.86))
	_setup_axes(ax)

	fig.suptitle('Concept 5 — Prior single-pathway models  →  whole-cell platform',
		fontsize=13, fontweight='bold', y=0.97)

	# Centre divider
	ax.plot([0.50, 0.50], [0.10, 0.85],
		color='#d0d7de', linewidth=1.5, linestyle=(0, (4, 4)))

	# LEFT — sparse prior model
	ax.text(0.25, 0.86, 'PRIOR MODELS',
		ha='center', va='center', fontsize=10, fontweight='bold',
		color='#57606a')
	_box(ax, 0.15, 0.55, 0.20, 0.16,
		'single ODE\n(IP3 → IP3R → cyt)',
		facecolor='#f6f8fa')
	_arrow(ax, 0.05, 0.63, 0.15, 0.63, label='IP3 in')
	_arrow(ax, 0.35, 0.63, 0.45, 0.63, label='Ca²⁺ trace out')
	ax.text(0.25, 0.30,
		'one input, one output,\nno cross-pathway coupling',
		ha='center', va='center', fontsize=8.5, color='#57606a',
		fontstyle='italic')

	# RIGHT — whole-cell platform
	ax.text(0.75, 0.86, 'PLATELET WCM (this work)',
		ha='center', va='center', fontsize=10, fontweight='bold',
		color='#3d4147')

	# Stack of mechanism boxes
	mech_y = 0.74
	mechanisms = [
		('GPCR cascade', '#dbeafe'),
		('PI cycle', '#dbeafe'),
		('IP3R / SERCA / DTS', '#fde68a'),
		('PMCA / SOCE / NCX', '#fde68a'),
		('MCU / NCLX / CaM / buffers', '#fde68a'),
	]
	for label, color in mechanisms:
		_box(ax, 0.55, mech_y, 0.30, 0.06, label,
			facecolor=color, fontsize=7.5)
		mech_y -= 0.08

	# Output: dose surface
	_embed_heatmap(fig, sweep, 'auc_ca_nMs',
		0.86, 0.28, 0.12, 0.34,
		title='dose surface\n(2 inputs)', cmap='magma')
	_arrow(ax, 0.85, 0.45, 0.86, 0.45, color='#3d4147', lw=1.6)

	ax.text(0.70, 0.16,
		'12 mechanisms end-to-end\ncoupled in a single ODE system\n'
		'→ 5 observables × 2-D dose plane',
		ha='center', va='center', fontsize=8.5, color='#3d4147',
		fontstyle='italic')

	_wireframe_note(fig,
		'methods-paper framing — "from single-pathway to whole-cell"')
	fig.savefig(out, dpi=140, bbox_inches='tight')
	plt.close(fig)


# ── Driver ───────────────────────────────────────────────────────────────

def main() -> int:
	FIGS_DIR.mkdir(parents=True, exist_ok=True)
	sweep = _load_sweep()
	concepts = [
		(1, concept1),
		(2, concept2),
		(3, concept3),
		(4, concept4),
		(5, concept5),
	]
	for n, fn in concepts:
		out = FIGS_DIR / f'graphical-abstract-concept-{n}.png'
		print(f'  rendering {out.relative_to(REPO_ROOT)}')
		fn(sweep, out)
	print(f'\nWrote {len(concepts)} wireframes into {FIGS_DIR.relative_to(REPO_ROOT)}/')
	return 0


if __name__ == '__main__':
	raise SystemExit(main())
