"""Animated terminal replay of a completed platelet simulation.

Reads a `simOut/` directory and replays the time series through a
Textual-based ASCII schematic of the platelet: receptors, GPCR cascade,
Ca²⁺ pools, PMCA and SOCE pumps, plus a scrolling sparkline of cytosolic
Ca²⁺. Designed as a viva-demo prop and teaching aid — kinetics that read
as a step-function in static figures show their actual shape when
animated in slow-mo.

Usage:
    PYTHONPATH=$PWD python runscripts/manual/replayTui.py <sim_outdir>
                                            [--speed 0.2] [--start 0]

    `sim_outdir` may be:
      - the run dir (e.g. out/phase3_xxx/with_ca/...), OR
      - the inner simOut/ path directly.

    --speed 1.0 is real-time (1 second of wall clock = 1 sim-second).
              0.2 is the recommended slow-mo for kinetics intuition.
              5.0 fast-forwards to skim the recovery.

Controls (any time):
    q         quit
    space     pause / resume
    +  /  -   speed up / slow down by 1.5×
    ←  /  →   step backward / forward 1 s (works any time; pauses if running)
    r         restart from t = 0

Issue #47. Replayer (option B) — reads finished simOut data; doesn't
attach to a live process.

Implementation note: this used to use rich.Live in alt-screen mode but
the cell-region rendering went blank on three different terminals
(Terminal.app, iTerm2, VSCode integrated terminal) despite a 4-mode
diagnostic confirming rich.Live itself was fine. Switched to Textual
(same author, framework-level TUI) which handles keyboard events,
event loop, and screen rendering natively — drops ~400 LoC of
hand-rolled termios + signal + select plumbing in the process.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from rich.box import SIMPLE
from rich.console import Console, Group
from rich.panel import Panel
from rich.text import Text

from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Container
from textual.reactive import reactive
from textual.timer import Timer
from textual.widgets import Footer, Static

from reconstruction.platelet.dataclasses.process import calcium_signalling as cs
from wholecell.io.tablereader import TableReader


# ── Playback constants ────────────────────────────────────────────────────
# Speed control (each +/- keypress multiplies/divides current speed by STEP).
# Textual's set_interval handles the actual tick scheduling.
SPEED_STEP            = 1.5
SPEED_MIN             = 0.01      # 1/100 real-time floor
SPEED_MAX             = 100.0     # 100× real-time ceiling

# Sparkline window — how many recent samples we show in the history strip.
SPARKLINE_WIDTH       = 70


# ── Snapshot type ─────────────────────────────────────────────────────────

@dataclass
class Snapshot:
	"""All values needed to render one frame.

	Indices map onto a single timestep `t` of the source CalciumTrace +
	BulkMolecules listeners; agonist concentrations are recomputed from
	the metadata-stored peaks/delay via the cs_mod forcing functions.
	"""
	t: float
	# Ca²⁺ pools
	ca_cyt_nM: float
	ca_dts_uM: float
	ca_mito_count: int        # raw model count (mito volume isn't tracked separately)
	ip3_nM: float
	# Agonist concentrations at this t (computed)
	thrombin_nM: float
	adp_uM: float
	atp_ex_uM: float
	# Receptors (counts)
	par1_inactive: int
	par1_active: int
	par1_internalized: int
	par4_inactive: int
	par4_active: int
	par4_internalized: int
	p2y1_inactive: int
	p2y1_active: int
	p2x1_closed: int
	p2x1_open: int
	p2x1_desens: int
	# Cascade
	gq_active: int
	plcb_inactive: int
	plcb_active: int
	pip2: int
	dag: int
	# DTS buffers (informational)
	calr_free: int
	calr_ca: int
	# Pumps / flux
	pmca_free: int
	pmca_ca: int
	pmca_cam: int
	soce_flux_nMs: float
	# Derived metrics
	peak_ca_so_far: float
	auc_above_rest: float   # ∫ max(ca_cyt − ca_rest, 0) dt, up to and incl. t


# Total counts (denominators for bars) — pulled at load time from the
# initial-state counts so the bars reflect each cell's per-receptor pool.
@dataclass
class Totals:
	par1: int = 0
	par4: int = 0
	p2y1: int = 0
	p2x1: int = 0
	gq: int = 0
	plcb: int = 0
	pmca: int = 0
	calr: int = 0


# ── Loading ───────────────────────────────────────────────────────────────

def _resolve_simout(path_str: str) -> Path:
	"""Accept either a run dir or the inner simOut/ path.

	A simOut/ must contain both `CalciumTrace/` and `BulkMolecules/`
	to be considered valid. If `path_str` itself satisfies that, it's
	returned unchanged. Otherwise we rglob for matching subdirs:
	exactly one match → use it; zero or multiple → raise (Phase 3
	driver outputs e.g. `with_ca/` and `no_ca/` siblings; the user
	must disambiguate).
	"""
	p = Path(path_str).expanduser().resolve()

	def _is_valid(d: Path) -> bool:
		return d.is_dir() and (d / 'CalciumTrace').is_dir() \
			and (d / 'BulkMolecules').is_dir()

	if _is_valid(p):
		return p
	# Nested layout: <run>/.../platelet_stub_*/<seed>/generation_*/<cell>/simOut
	candidates = [c for c in sorted(p.rglob('simOut')) if _is_valid(c)]
	if not candidates:
		raise FileNotFoundError(
			f'No simOut/ containing CalciumTrace + BulkMolecules found '
			f'under {p}. Pass the directory containing those listener dirs.')
	if len(candidates) > 1:
		listing = '\n  '.join(str(c.relative_to(p)) for c in candidates)
		raise FileNotFoundError(
			f'Ambiguous: {len(candidates)} simOut/ candidates under {p}.\n'
			f'Pass one explicitly:\n  {listing}')
	return candidates[0]


def _metadata_for(simout: Path) -> dict:
	"""Walk up from simOut/ to find the run-level metadata.json."""
	for ancestor in (simout, *simout.parents[:6]):
		candidate = ancestor / 'metadata' / 'metadata.json'
		if candidate.is_file():
			try:
				return json.loads(candidate.read_text())
			except (OSError, json.JSONDecodeError):
				continue
	return {}


def load_snapshots(simout: Path) -> tuple[list[Snapshot], Totals, dict]:
	"""Pre-load every timestep into a list of Snapshot. Sims are short
	(30-300 s), so loading everything upfront is fine and lets us
	step forward/backward instantly."""
	ca_reader = TableReader(str(simout / 'CalciumTrace'))
	bm_reader = TableReader(str(simout / 'BulkMolecules'))
	meta = _metadata_for(simout)

	# CalciumTrace cols
	t = ca_reader.readColumn('time').flatten()
	ca_cyt = ca_reader.readColumn('ca_cyt_nM').flatten()
	ca_dts = ca_reader.readColumn('ca_dts_uM').flatten()
	ip3 = ca_reader.readColumn('ip3_nM').flatten()
	soce = ca_reader.readColumn('soce_flux_nMs').flatten()
	pmca_free_col = ca_reader.readColumn('pmca_free').flatten()
	pmca_ca_col = ca_reader.readColumn('pmca_ca').flatten()
	pmca_cam_col = ca_reader.readColumn('pmca_cam').flatten()

	# BulkMolecules: counts matrix (T, N)
	counts = bm_reader.readColumn('counts')
	bm_ids = list(bm_reader.readAttribute('objectNames'))
	idx = {name: i for i, name in enumerate(bm_ids)}

	def col(name: str) -> np.ndarray:
		return counts[:, idx[name]] if name in idx else np.zeros(len(t))

	# Agonist peaks (recompute the forcing curves)
	thr_peak = meta.get('thrombin_peak_nM')
	adp_peak = meta.get('adp_peak_uM')
	atp_peak = meta.get('atp_ex_peak_uM')
	delay = float(meta.get('agonist_delay') or 0.0)

	thrombin_curve = np.array(
		[cs.thrombin_nM(float(ti), delay=delay, peak_nM=thr_peak) for ti in t])
	adp_curve = np.array(
		[cs.adp_uM(float(ti), delay=delay, peak_uM=adp_peak) for ti in t])
	atp_curve = np.array(
		[cs.atp_ex_forcing_uM(float(ti), delay=delay, peak_uM=atp_peak) for ti in t])

	# Totals (sum of all sub-states at t=0)
	totals = Totals(
		par1=int(col('PAR1_inactive[pl]')[0] + col('PAR1_active[pl]')[0]
			+ col('PAR1_internalized[pl]')[0]),
		par4=int(col('PAR4_inactive[pl]')[0] + col('PAR4_active[pl]')[0]
			+ col('PAR4_internalized[pl]')[0]),
		p2y1=int(col('P2Y1_inactive[pl]')[0] + col('P2Y1_active[pl]')[0]),
		p2x1=int(col('P2X1[pl]')[0] + col('P2X1_O[pl]')[0] + col('P2X1_D[pl]')[0]),
		gq=cs.N_GQ_TOTAL,
		plcb=int(col('PLCb_inactive[c]')[0] + col('PLCb_active[c]')[0]),
		pmca=int(pmca_free_col[0] + pmca_ca_col[0] + pmca_cam_col[0]),
		calr=cs.N_CALR,
	)

	# Running peak + AUC of cyt Ca²⁺ above rest
	rest = float(ca_cyt[0])
	running_peak = np.maximum.accumulate(ca_cyt)
	# Trapezoidal AUC at each timestep (above rest, clipped to 0)
	above = np.maximum(ca_cyt - rest, 0.0)
	dt = np.diff(t, prepend=t[0])
	auc_above = np.cumsum(0.5 * (above + np.roll(above, 1)) * dt)
	auc_above[0] = 0.0

	snaps: list[Snapshot] = []
	for i in range(len(t)):
		snaps.append(Snapshot(
			t=float(t[i]),
			ca_cyt_nM=float(ca_cyt[i]),
			ca_dts_uM=float(ca_dts[i]),
			# Mito Ca²⁺ stays as a raw count — mito volume isn't tracked
			# separately by the engine, so a concentration would be a
			# misleading derivation.
			ca_mito_count=int(col('CA2_MITO[m]')[i]),
			ip3_nM=float(ip3[i]),
			thrombin_nM=float(thrombin_curve[i]),
			adp_uM=float(adp_curve[i]),
			atp_ex_uM=float(atp_curve[i]),
			par1_inactive=int(col('PAR1_inactive[pl]')[i]),
			par1_active=int(col('PAR1_active[pl]')[i]),
			par1_internalized=int(col('PAR1_internalized[pl]')[i]),
			par4_inactive=int(col('PAR4_inactive[pl]')[i]),
			par4_active=int(col('PAR4_active[pl]')[i]),
			par4_internalized=int(col('PAR4_internalized[pl]')[i]),
			p2y1_inactive=int(col('P2Y1_inactive[pl]')[i]),
			p2y1_active=int(col('P2Y1_active[pl]')[i]),
			p2x1_closed=int(col('P2X1[pl]')[i]),
			p2x1_open=int(col('P2X1_O[pl]')[i]),
			p2x1_desens=int(col('P2X1_D[pl]')[i]),
			gq_active=int(col('Gq_active[c]')[i]),
			plcb_inactive=int(col('PLCb_inactive[c]')[i]),
			plcb_active=int(col('PLCb_active[c]')[i]),
			pip2=int(col('PIP2[c]')[i]),
			dag=int(col('DAG[c]')[i]),
			calr_free=int(col('CALR_free[dts]')[i]),
			calr_ca=int(col('CALR_Ca[dts]')[i]),
			pmca_free=int(pmca_free_col[i]),
			pmca_ca=int(pmca_ca_col[i]),
			pmca_cam=int(pmca_cam_col[i]),
			soce_flux_nMs=float(soce[i]),
			peak_ca_so_far=float(running_peak[i]),
			auc_above_rest=float(auc_above[i]),
		))
	return snaps, totals, meta


# ── Rendering helpers ─────────────────────────────────────────────────────

_BAR_FULL = '▰'
_BAR_EMPTY = '▱'
_SPARK_CHARS = ' ▁▂▃▄▅▆▇█'


def _bar(filled: float, width: int = 12) -> str:
	"""filled is 0.0–1.0; return a ▰▱ progress bar of `width` segments."""
	n = max(0, min(width, int(round(filled * width))))
	return _BAR_FULL * n + _BAR_EMPTY * (width - n)


_SPARKLINE_FLAT_TOL = 1.0e-9


def _sparkline(values: np.ndarray, width: int = 60) -> str:
	"""U+2580-block sparkline of recent history.

	Always returns exactly `width` characters, right-aligned: spaces on
	the left, glyphs on the right (so the newest sample sits at the
	right edge and the panel doesn't jitter while history fills up).
	"""
	if len(values) == 0:
		return ' ' * width
	tail = values[-width:]
	v_min, v_max = float(np.min(tail)), float(np.max(tail))
	if v_max - v_min < _SPARKLINE_FLAT_TOL:
		# Constant series — render as a flat bar of the smallest block.
		return (_SPARK_CHARS[1] * len(tail)).rjust(width)
	normalised = (tail - v_min) / (v_max - v_min)
	indices = np.clip((normalised * (len(_SPARK_CHARS) - 1)).astype(int),
		0, len(_SPARK_CHARS) - 1)
	return ''.join(_SPARK_CHARS[i] for i in indices).rjust(width)


def _ca_colour(ca_nM: float) -> str:
	"""rich colour name for a Ca²⁺ value — green at rest, red at peak."""
	if ca_nM < 150:
		return 'green'
	if ca_nM < 250:
		return 'yellow'
	if ca_nM < 400:
		return 'orange1'
	return 'red'


def _trend_arrow(curr: float, prev: float, tol: float = 0.0) -> str:
	if curr > prev + tol:
		return '[green]▲[/green]'
	if curr < prev - tol:
		return '[red]▼[/red]'
	return ' '


# ── Layout builders ───────────────────────────────────────────────────────

# Cell-schematic dimensions. Each top/bottom membrane carries three
# fixed-width receptor tags. Inner DTS box sits on the cytosolic left,
# Mito box on the right. The whole thing is one big rich.Text block.
_CELL_INNER_W = 96   # interior of the cell, excluding the side ║ walls
_RX_TAG_W     = 28   # width of one "[ LABEL bar pct% ]" block on the membrane


def _receptor_tag(label: str, frac: float) -> str:
	"""Fixed-width receptor tag for embedding in a membrane line.

	Natural width = 22 chars: "[ LABEL ▰..▱ XX% ]". No internal padding
	— spacing between tags is handled by `_membrane_line` distributing
	fill characters.
	"""
	bar = _bar(frac, width=10)
	pct = max(0, min(99, int(round(frac * 100))))
	return f'[ {label:<4} {bar} {pct:>2}% ]'


def _p2x1_tag(s: Snapshot) -> str:
	"""P2X1 is ionotropic — show open/desens counts not a fraction bar.
	22 chars to match `_receptor_tag` width."""
	return f'[ P2X1  O {s.p2x1_open:>3} D {s.p2x1_desens:>3}  ]'


def _flux_tag(label: str, value_str: str) -> str:
	"""A non-receptor membrane-edge tag (SOCE flux, PM leak, etc.).
	22 chars to match `_receptor_tag` width."""
	return f'[ {label:<6} {value_str:<10} ]'


def _membrane_line(corner_left: str, corner_right: str, fill: str,
		segments: list[str]) -> str:
	"""Build a membrane line by joining segments with fill chars.

	The total width is _CELL_INNER_W + 2 (for the two corner glyphs).
	Spreads any leftover fill space evenly between segments.
	"""
	n = len(segments)
	used = sum(len(s) for s in segments)
	# Leftover fill chars to distribute: cell_inner_width − sum(segments)
	leftover = _CELL_INNER_W - used
	# n+1 gaps (before first segment, between each pair, after last segment).
	if leftover < 0:
		leftover = 0
	gap = leftover // (n + 1)
	extra = leftover - gap * (n + 1)
	parts = [fill * gap]
	for i, seg in enumerate(segments):
		parts.append(seg)
		parts.append(fill * (gap + (1 if i < extra else 0)))
	body = ''.join(parts)
	return f'{corner_left}{body}{corner_right}'


def _inside_line(content_markup: str) -> str:
	"""Wrap a line in the cell's side walls, padding to inner width."""
	# Note: this measures the *visible* width of `content_markup` only if
	# the markup has no bracketed tags. For our use-case we pre-pad the
	# content ourselves before calling this, so we just append ║ on each side.
	return f'║{content_markup}║'


def _cell_schematic(s: Snapshot, totals: Totals, ca_ex_uM: float) -> Panel:
	"""The headline view: a stylised ASCII platelet cross-section.

	Cell outer border = plasma membrane (double-line box). Receptors
	sit *on* the membrane (top + bottom). DTS and Mito are nested
	single-line boxes inside the cytosol. Ca²⁺ readouts, cascade bars,
	and pump fluxes fill the cytosolic space.

	Returns a `Panel` with `box=SIMPLE` and no title — Live mode tends
	to silently drop bare `Text` from a Layout region (the failure
	mode is the whole region rendering blank, observed by users on
	wide terminals with plenty of vertical room). Wrapping in a Panel
	with a minimal box style keeps the schematic visible without
	adding redundant border weight on top of the cell's own ╔═╗ chars.

	All width arithmetic is hand-calibrated for `_CELL_INNER_W = 96`.
	If you change that constant, walk through each line below and
	confirm widths still add up.
	"""
	# ── Convenience metrics ───────────────────────────────────────
	par1_frac = s.par1_active / max(totals.par1, 1)
	par4_frac = s.par4_active / max(totals.par4, 1)
	p2y1_frac = s.p2y1_active / max(totals.p2y1, 1)
	gq_frac = s.gq_active / max(totals.gq, 1)
	plcb_frac = s.plcb_active / max(totals.plcb, 1)
	calr_frac = s.calr_ca / max(totals.calr, 1)
	ca_colour = _ca_colour(s.ca_cyt_nM)
	ip3_arrow = _trend_arrow(s.ip3_nM, 50.0, tol=2.0)
	dts_arrow = _trend_arrow(s.ca_dts_uM, 250.0, tol=2.0)

	# ── Membrane lines (top = receptors, bottom = receptor + fluxes) ──
	top_membrane = _membrane_line('╔', '╗', '═', [
		_receptor_tag('PAR1', par1_frac),
		_p2x1_tag(s),
		_receptor_tag('P2Y1', p2y1_frac),
	])
	bottom_membrane = _membrane_line('╚', '╝', '═', [
		_receptor_tag('PAR4', par4_frac),
		_flux_tag('SOCE', f'{s.soce_flux_nMs:>+6.2f} nM/s'),
		_flux_tag('PMleak', '75 ions/s'),
	])

	# ── Inside the cell ───────────────────────────────────────────
	# Helper to right-pad an inside line to exact width (after stripping
	# rich markup tags so the visible width still lines up).
	def inside(line_markup: str, visible_len: int) -> str:
		pad = max(0, _CELL_INNER_W - visible_len)
		return f'║{line_markup}{" " * pad}║'

	inside_lines: list[str] = []
	# Compact: no padding row after the membrane line.

	# Cytosol header
	hdr = '  [b cyan]Cytosol[/b cyan]'
	inside_lines.append(inside(hdr, 9))  # "  Cytosol" = 9 visible chars

	# Headline Ca²⁺ / IP3 / cascade readouts
	# All these are right-padded inside the inner cell width.
	ca_line = (
		f'  Ca²⁺  [{ca_colour}]{s.ca_cyt_nM:>6.1f}[/{ca_colour}] nM   '
		f'peak {s.peak_ca_so_far:>5.0f} nM   '
		f'AUC {s.auc_above_rest:>9,.0f} nM·s'
	)
	# Compute visible len: prefix + numbers + suffix (markup tags don't render)
	ca_visible = len(
		f'  Ca²⁺  {s.ca_cyt_nM:>6.1f} nM   '
		f'peak {s.peak_ca_so_far:>5.0f} nM   '
		f'AUC {s.auc_above_rest:>9,.0f} nM·s')
	inside_lines.append(inside(ca_line, ca_visible))

	ip3_line = (
		f'  IP3   [magenta]{s.ip3_nM:>6.1f}[/magenta] nM   {ip3_arrow}   '
		f'PIP2 {s.pip2:>10,}   DAG {s.dag:>6,}'
	)
	ip3_visible = len(
		f'  IP3   {s.ip3_nM:>6.1f} nM       '  # the arrow markup ▲ / ▼ / ' ' is 1 visible
		f'PIP2 {s.pip2:>10,}   DAG {s.dag:>6,}')
	inside_lines.append(inside(ip3_line, ip3_visible))

	# Gq + PLCβ activation bars
	cas_line = (
		f'  [b]Gαq[/b]  {_bar(gq_frac, 18)}  {s.gq_active:>5}/{totals.gq:<5}   '
		f'[b]PLCβ[/b] {_bar(plcb_frac, 18)}  {s.plcb_active:>4}/{totals.plcb:<4}'
	)
	cas_visible = len(
		f'  Gαq  {"▰" * 18}  {s.gq_active:>5}/{totals.gq:<5}   '
		f'PLCβ {"▰" * 18}  {s.plcb_active:>4}/{totals.plcb:<4}')
	inside_lines.append(inside(cas_line, cas_visible))

	# ── DTS box (left) and Mito box (right), drawn side by side ──
	#
	# DTS is 52 cols wide; Mito is 28 cols wide; gap between = 6 cols;
	# leading indent inside the cell = 4 cols; trailing pad = 6 cols.
	# 4 + 52 + 6 + 28 + 6 = 96 = _CELL_INNER_W. ✓
	DTS_W, GAP_W, MITO_W = 52, 6, 28
	INDENT_W = 4
	TRAIL_W = _CELL_INNER_W - INDENT_W - DTS_W - GAP_W - MITO_W
	indent = ' ' * INDENT_W
	gap_str = ' ' * GAP_W
	trail = ' ' * TRAIL_W

	def organelle_row(dts_text: str, dts_visible: int,
			mito_text: str = '', mito_visible: int = 0) -> str:
		# dts_text and mito_text are the *content* lines (no borders).
		dts_pad = max(0, DTS_W - 2 - dts_visible)   # -2 for the │ │ borders
		mito_pad = max(0, MITO_W - 2 - mito_visible)
		dts_box = f'│{dts_text}{" " * dts_pad}│'
		if mito_text:
			mito_box = f'│{mito_text}{" " * mito_pad}│'
		else:
			mito_box = ' ' * MITO_W
		return f'║{indent}{dts_box}{gap_str}{mito_box}{trail}║'

	# Top borders of DTS + Mito
	dts_top = '┌' + '─' * (DTS_W - 2) + '┐'
	mito_top = '┌' + '─' * (MITO_W - 2) + '┐'
	inside_lines.append(f'║{indent}{dts_top}{gap_str}{mito_top}{trail}║')

	# DTS row 1 + Mito row 1 — titles
	inside_lines.append(organelle_row(
		' [b yellow]DTS[/b yellow]                                              ',
		3 + 1 + 46,    # " " + "DTS" + spaces
		' [b magenta]Mito[/b magenta]                ',
		4 + 1 + 16))

	# DTS row 2 — Ca²⁺ + arrow ; Mito row 2 — Ca²⁺ count (not concentration)
	dts_ca_line = f' Ca²⁺  [b yellow]{s.ca_dts_uM:>6.1f}[/b yellow] µM   {dts_arrow}'
	dts_ca_visible = len(f' Ca²⁺  {s.ca_dts_uM:>6.1f} µM    ')
	mito_ca_line = f' Ca²⁺  [b magenta]{s.ca_mito_count:>7,}[/b magenta] ions'
	mito_ca_visible = len(f' Ca²⁺  {s.ca_mito_count:>7,} ions')
	inside_lines.append(organelle_row(dts_ca_line, dts_ca_visible,
		mito_ca_line, mito_ca_visible))

	# DTS row 3 — CALR buffer state ; Mito row 3 — pump activity label
	calr_line = (
		f' CALR-bound  {s.calr_ca:>7,} / {totals.calr:>7,}  '
		f'({100 * calr_frac:>4.1f} %)'
	)
	calr_visible = len(
		f' CALR-bound  {s.calr_ca:>7,} / {totals.calr:>7,}  '
		f'({100 * calr_frac:>4.1f} %)')
	mito_pump_line = ' [dim]MCU + NCLX active[/dim]'
	mito_pump_visible = len(' MCU + NCLX active')
	inside_lines.append(organelle_row(calr_line, calr_visible,
		mito_pump_line, mito_pump_visible))

	# Mito bottom border on the same row as DTS-only IP3R/SERCA label,
	# then DTS continues solo for two more rows + its own bottom border.
	mito_bottom = '└' + '─' * (MITO_W - 2) + '┘'
	dts_ip3r = ' [b yellow]IP3R[/b yellow]                       [b green]SERCA[/b green]'
	dts_ip3r_visible = len(' IP3R                       SERCA')
	dts_ip3r_pad = max(0, DTS_W - 2 - dts_ip3r_visible)
	inside_lines.append(
		f'║{indent}│{dts_ip3r}{" " * dts_ip3r_pad}│{gap_str}{mito_bottom}{trail}║')

	# DTS row 5 — flux arrows (DTS box only now)
	dts_arrows = '  [yellow]▲[/yellow]                          [green]▼[/green]'
	dts_arrows_visible = len('  ▲                          ▼')
	dts_arrows_pad = max(0, DTS_W - 2 - dts_arrows_visible)
	inside_lines.append(
		f'║{indent}│{dts_arrows}{" " * dts_arrows_pad}│{gap_str}{" " * MITO_W}{trail}║')

	# DTS bottom border (Mito already closed; trail with empty space)
	dts_bottom = '└' + '─' * (DTS_W - 2) + '┘'
	inside_lines.append(f'║{indent}{dts_bottom}{gap_str}{" " * MITO_W}{trail}║')

	# PMCA / SOCE pump details — bottom of cytosol space
	pmca_line = (
		f'  [b]PMCA[/b]  basal {s.pmca_ca:>4}   CaM-active {s.pmca_cam:>4}   '
		f'free {s.pmca_free:>4}'
	)
	pmca_visible = len(
		f'  PMCA  basal {s.pmca_ca:>4}   CaM-active {s.pmca_cam:>4}   '
		f'free {s.pmca_free:>4}')
	inside_lines.append(inside(pmca_line, pmca_visible))

	# Extracellular reservoir label — sits on the same row as the top
	# membrane in spirit, but we put it on its own line above so the
	# membrane border stays unbroken.
	ex_line = (
		f' [b cyan]Extracellular[/b cyan]   Ca²⁺ {ca_ex_uM/1000:>4.1f} mM   '
		f'thrombin {s.thrombin_nM:>6.3f} nM   '
		f'ADP {s.adp_uM:>6.3f} µM   '
		f'ATP {s.atp_ex_uM:>6.3f} µM'
	)

	# Assemble: extracellular line, then cell box.
	#
	# IMPORTANT: build as a Group of individual Text-per-line objects
	# rather than a single multi-line Text. rich.Live's rendering of a
	# Layout cell silently drops a multi-line Text whose lines exceed
	# the auto-wrap heuristics — observed empirically on a 232×73
	# terminal with plenty of room. One-Text-per-line + no_wrap=True
	# bypasses the wrap logic entirely and is what Live actually
	# wants. The outer Panel wrap was a red herring; the real issue
	# was the multi-line Text.
	all_lines = [
		ex_line,
		top_membrane,
		*inside_lines,
		bottom_membrane,
	]
	text_items: list[Text] = []
	for line in all_lines:
		t = Text.from_markup(line)
		t.no_wrap = True
		t.overflow = 'crop'
		text_items.append(t)
	return Panel(Group(*text_items),
		box=SIMPLE, border_style='dim', padding=(0, 0))


def _sparkline_panel(history_ca: np.ndarray, history_ip3: np.ndarray) -> Panel:
	width = 70
	ca_spark = _sparkline(history_ca, width=width)
	ip3_spark = _sparkline(history_ip3, width=width)
	current_ca = history_ca[-1] if len(history_ca) else 0.0
	ca_label = (
		f'[{_ca_colour(current_ca)}]cyt Ca²⁺[/{_ca_colour(current_ca)}]'
	)
	lines = [
		f'{ca_label:<22} {ca_spark}   '
		f'{history_ca.min():>6.0f}…{history_ca.max():>6.0f} nM',
		f'[magenta]IP3[/magenta]               '
		f'   {ip3_spark}   '
		f'{history_ip3.min():>6.0f}…{history_ip3.max():>6.0f} nM',
	]
	return Panel(Text.from_markup('\n'.join(lines)),
		title='[b]History (last 70 s)[/b]', border_style='dim')



def _header(snap: Snapshot, frame: int, n_snapshots: int,
		speed: float, paused: bool, extra: str) -> Panel:
	"""Header progress is computed from the frame index (0..n-1) rather
	than from sim time, so it behaves correctly when timesteps aren't
	exactly 1 s and shows 0 % at the first frame / 100 % at the last."""
	denom = max(n_snapshots - 1, 1)
	bar_frac = max(0.0, min(1.0, frame / denom))
	progress_bar = _bar(bar_frac, 30)
	state = '[red]⏸ PAUSED[/red]' if paused else f'⏵ {speed:.2f}×'
	t_end = float(n_snapshots - 1)
	body = (
		f'Platelet WCM replay   '
		f't = [bold]{snap.t:>6.1f}[/bold] / {t_end:>4.0f} s   '
		f'{progress_bar}   {state}'
	)
	if extra:
		body += f'   [dim]{extra}[/dim]'
	return Panel(Text.from_markup(body), border_style='white')


# ── Textual app ───────────────────────────────────────────────────────────

class PlateletReplayApp(App):
	"""Textual app — replay a finished platelet sim as an animated schematic.

	Three Static widgets stacked vertically (header / cell / sparkline) plus
	a Footer for the keybindings. Reactive `frame`, `speed`, `paused`
	drive a Textual interval timer that advances the frame and triggers
	re-rendering.

	All keyboard handling, event-loop scheduling, signal handling, and
	terminal-state restoration is delegated to Textual — the previous
	rich.Live implementation hand-rolled all of this and silently dropped
	the cell region on three different terminals despite the underlying
	rich library working fine in a 4-mode diagnostic.
	"""

	CSS = """
	Screen { layout: vertical; }
	#header { dock: top;    height: 3; padding: 0 1; }
	#cell   { height: 1fr;            padding: 0 0; content-align: left top; }
	#spark  { dock: bottom; height: 5; padding: 0 1; }
	"""

	BINDINGS = [
		Binding('q',            'quit',           'Quit'),
		Binding('space',        'toggle_pause',   'Pause/Resume'),
		Binding('plus',         'speed_up',       'Speed ×1.5'),
		Binding('equals_sign',  'speed_up',       show=False),
		Binding('minus',        'speed_down',     'Speed ÷1.5'),
		Binding('right',        'step_forward',   'Step +1 s'),
		Binding('left',         'step_back',      'Step −1 s'),
		Binding('r',            'restart',        'Restart'),
	]

	frame  = reactive(0)
	speed  = reactive(1.0)
	paused = reactive(False)

	def __init__(self, snapshots: list[Snapshot], totals: Totals,
			meta: dict, initial_speed: float = 1.0, start_frame: int = 0):
		super().__init__()
		self.snapshots = snapshots
		self.totals    = totals
		self.meta      = meta
		self._n        = len(snapshots)
		# Pre-compute the full sparkline history once; we slice per frame.
		self._ca_hist  = np.array([s.ca_cyt_nM for s in snapshots])
		self._ip3_hist = np.array([s.ip3_nM    for s in snapshots])
		# Set initial reactive values via __dict__ to avoid triggering
		# watchers before compose() has built the widgets.
		self.set_reactive(PlateletReplayApp.speed, max(SPEED_MIN,
			min(initial_speed, SPEED_MAX)))
		self.set_reactive(PlateletReplayApp.frame, max(0,
			min(start_frame, self._n - 1)))
		self._tick_handle: Timer | None = None

	def compose(self) -> ComposeResult:
		yield Static(id='header')
		yield Container(Static(id='cell'))
		yield Static(id='spark')
		yield Footer()

	def on_mount(self) -> None:
		self._schedule_tick()
		self._render_frame()

	def _schedule_tick(self) -> None:
		"""Re-arm the periodic tick at the current speed. Called on mount
		and whenever `speed` changes."""
		if self._tick_handle is not None:
			self._tick_handle.stop()
		# 1 sim-second per (1/speed) wall-seconds; pause is reactive so
		# the interval keeps firing but `_tick` is a no-op while paused.
		self._tick_handle = self.set_interval(1.0 / self.speed, self._tick)

	def _tick(self) -> None:
		if self.paused:
			return
		if self.frame < self._n - 1:
			self.frame += 1
		else:
			# Hold at the last frame so the user can read the final state.
			self.paused = True

	# Reactive watchers — trigger a re-render on any state change.
	def watch_frame(self, _old: int, _new: int) -> None:
		self._render_frame()

	def watch_paused(self, _old: bool, _new: bool) -> None:
		self._render_frame()

	def watch_speed(self, _old: float, _new: float) -> None:
		self._schedule_tick()
		self._render_frame()

	def _render_frame(self) -> None:
		if not self.is_running:
			return
		snap = self.snapshots[self.frame]
		ca_ex_mM = float(self.meta.get('ca_ex_mM') or 1.2)
		ca_ex_uM = ca_ex_mM * 1000.0
		extra = f'seed {self.meta.get("seed", "?")} · Ca_ex {ca_ex_mM} mM'

		self.query_one('#header', Static).update(
			_header(snap, self.frame, self._n, self.speed, self.paused, extra))
		self.query_one('#cell', Static).update(
			_cell_schematic(snap, self.totals, ca_ex_uM))
		w_start = max(0, self.frame - (SPARKLINE_WIDTH - 1))
		self.query_one('#spark', Static).update(_sparkline_panel(
			self._ca_hist[w_start:self.frame + 1],
			self._ip3_hist[w_start:self.frame + 1]))

	# ── Action handlers (bound in BINDINGS above) ───────────────────────

	def action_toggle_pause(self) -> None:
		self.paused = not self.paused

	def action_speed_up(self) -> None:
		self.speed = min(self.speed * SPEED_STEP, SPEED_MAX)

	def action_speed_down(self) -> None:
		self.speed = max(self.speed / SPEED_STEP, SPEED_MIN)

	def action_step_forward(self) -> None:
		self.paused = True
		if self.frame < self._n - 1:
			self.frame += 1

	def action_step_back(self) -> None:
		self.paused = True
		if self.frame > 0:
			self.frame -= 1

	def action_restart(self) -> None:
		self.frame = 0
		self.paused = False


def replay(simout: Path, initial_speed: float, start_t: int) -> int:
	"""Load the sim's snapshots and run the Textual app to completion."""
	snapshots, totals, meta = load_snapshots(simout)
	if not snapshots:
		print('No snapshots loaded — empty CalciumTrace?', file=sys.stderr)
		return 2
	app = PlateletReplayApp(snapshots, totals, meta,
		initial_speed=initial_speed, start_frame=start_t)
	app.run()
	return 0


# ── CLI ───────────────────────────────────────────────────────────────────

def main(argv: list[str] | None = None) -> int:
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('sim_outdir',
		help='Path to a run dir, or directly to a simOut/ directory.')
	parser.add_argument('--speed', type=float, default=1.0,
		help='Playback speed (1.0 = real-time, 0.2 = slow-mo). Default 1.0.')
	parser.add_argument('--start', type=int, default=0,
		help='Initial timestep index. Default 0.')
	parser.add_argument('--dump-frame', dest='dump_frame', type=int, default=None,
		help='Render frame N once to stdout via rich Console.print (no '
		'Textual, no keyboard) and exit. Useful for CI snapshot tests and '
		'for triaging "blank screen" reports — if --dump-frame N looks right '
		'but interactive mode does not, the bug is in the Textual layer.')
	args = parser.parse_args(argv)

	try:
		simout = _resolve_simout(args.sim_outdir)
	except FileNotFoundError as e:
		print(f'error: {e}', file=sys.stderr)
		return 1

	if args.dump_frame is not None:
		return _dump_one_frame(simout, args.dump_frame)

	if not sys.stdin.isatty() or not sys.stdout.isatty():
		print('replayTui needs a real terminal (tty) for the animated mode; '
			'pass --dump-frame N for a one-shot render to stdout instead.',
			file=sys.stderr)
		return 2

	return replay(simout, args.speed, args.start)


def _dump_one_frame(simout: Path, frame_idx: int) -> int:
	"""Render a single frame via Console.print and exit. Diagnostic mode.

	Prints layout cell dimensions + schematic line count alongside the
	rendered frame so anyone reporting "blank screen" can paste the
	output for triage.
	"""
	snapshots, totals, meta = load_snapshots(simout)
	if not snapshots:
		print('no snapshots loaded', file=sys.stderr)
		return 2
	idx = max(0, min(frame_idx, len(snapshots) - 1))
	snap = snapshots[idx]
	console = Console()
	ca_ex_mM = float(meta.get('ca_ex_mM') or 1.2)
	ca_ex_uM = ca_ex_mM * 1000.0

	# Diagnostic header
	print(f'== replayTui dump diagnostics ==')
	print(f'simout       : {simout}')
	print(f'frame        : {idx} of {len(snapshots) - 1}  (t = {snap.t:.1f} s)')
	print(f'terminal size: {console.size.width} cols x {console.size.height} rows')
	schematic = _cell_schematic(snap, totals, ca_ex_uM)
	# `schematic` is a Panel wrapping a Group of Texts (one per line).
	inner = schematic.renderable
	if isinstance(inner, Group):
		schematic_lines = [t.plain if hasattr(t, 'plain') else str(t)
			for t in inner.renderables]
	elif hasattr(inner, 'plain'):
		schematic_lines = inner.plain.split('\n')
	else:
		schematic_lines = [str(inner)]
	max_width = max(len(line) for line in schematic_lines) if schematic_lines else 0
	print(f'schematic    : {len(schematic_lines)} lines, max width {max_width} cols')
	print(f'rich version : {__import__("rich").__version__ if hasattr(__import__("rich"), "__version__") else "(no __version__)"}')
	print()

	# Render header + schematic + sparkline to stdout, one after another.
	# No Layout / Live — the Textual app handles those; this path is a
	# plain Console.print sequence for diagnostics + CI snapshot tests.
	ca_hist = np.array([s.ca_cyt_nM for s in snapshots[:idx + 1]])
	ip3_hist = np.array([s.ip3_nM for s in snapshots[:idx + 1]])
	console.print(_header(snap, idx, len(snapshots), 1.0, False,
		f'seed {meta.get("seed", "?")} · dump'))
	console.print(schematic)
	console.print(_sparkline_panel(ca_hist, ip3_hist))
	return 0


if __name__ == '__main__':
	sys.exit(main())
