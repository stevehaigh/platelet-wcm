"""
Phase 3 validation plot — Dolan 2014 Fig. 4 reproduction.

Compares the platelet-wcm v0.2 model's Ca²⁺ transient under two
conditions:

  +Ca_ex   :  extracellular Ca²⁺ = 1.2 mM (Dolan nominal)
  −Ca_ex   :  extracellular Ca²⁺ = 0     (Dolan EDTA condition)

Reference data: ``reports/data/dolan-2014-fig4-reference.json``
(literature values from Dolan & Diamond 2014 Fig. 4 + Fig. 3B
filtering criteria — not pixel-digitised).

Output: a 3-panel figure
  Panel 1 — [Ca²⁺]_cyt (nM) vs time, both conditions + Dolan
            +Ca peak / plateau / −Ca peak reference bands.
  Panel 2 — [Ca²⁺]_dts (µM) vs time, both conditions + Dolan
            DTS minimum reference bands.
  Panel 3 — Acceptance-criteria table (PASS / FAIL annotations
            against the four Dolan filtering criteria).

This module is not in ``analysis/single/`` because it operates on
two sim directories (one per condition), not one. Invoke via
``runscripts/manual/runPhase3.py``, which runs both sims and calls
``make_phase3_plot`` directly.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import Optional

import matplotlib as mp
mp.use('Agg')
from matplotlib import pyplot as plt
import numpy as np

from wholecell.io.tablereader import TableReader


REFERENCE_JSON = os.path.join(
	os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
		os.path.realpath(__file__))))),
	'reports', 'data', 'dolan-2014-fig4-reference.json',
)


@dataclass
class _ConditionTrace:
	"""One sim run's reduced state for the Phase 3 figure."""
	label: str
	colour: str
	t: np.ndarray
	ca_cyt_nM: np.ndarray
	ca_dts_uM: np.ndarray
	stim1_dim: np.ndarray
	soce_flux_nMs: np.ndarray
	peak_cyt_nM: float
	peak_t_s: float
	dts_min_uM: float

	@classmethod
	def from_simout(cls, simout_dir: str, label: str, colour: str
			) -> '_ConditionTrace':
		ct = TableReader(os.path.join(simout_dir, 'CalciumTrace'))
		t          = ct.readColumn('time').flatten()
		ca_cyt     = ct.readColumn('ca_cyt_nM').flatten()
		ca_dts     = ct.readColumn('ca_dts_uM').flatten()
		stim1_dim  = ct.readColumn('stim1_dim').flatten()
		soce       = ct.readColumn('soce_flux_nMs').flatten()
		peak_idx = int(np.argmax(ca_cyt))
		return cls(
			label=label, colour=colour, t=t,
			ca_cyt_nM=ca_cyt, ca_dts_uM=ca_dts,
			stim1_dim=stim1_dim, soce_flux_nMs=soce,
			peak_cyt_nM=float(ca_cyt[peak_idx]),
			peak_t_s=float(t[peak_idx]),
			dts_min_uM=float(np.min(ca_dts)),
		)


def _load_reference(path: str = REFERENCE_JSON) -> dict:
	with open(path) as f:
		return json.load(f)


def _eval_criteria(with_ca: _ConditionTrace, no_ca: _ConditionTrace,
		ref: dict) -> list[tuple[str, str, str, bool]]:
	"""Evaluate the Dolan acceptance criteria; return rows for the table panel."""
	rows = []
	# Active +Ca
	with_pass = with_ca.peak_cyt_nM > 200.0
	rows.append((
		'Active (+Ca_ex)',
		'peak Ca_cyt > 200 nM',
		f'{with_ca.peak_cyt_nM:.0f} nM',
		with_pass,
	))
	# Active −Ca
	no_pass = no_ca.peak_cyt_nM > 200.0
	rows.append((
		'Active (−Ca_ex)',
		'peak Ca_cyt > 200 nM',
		f'{no_ca.peak_cyt_nM:.0f} nM',
		no_pass,
	))
	# SOCE differential
	soce_diff = abs(with_ca.peak_cyt_nM - no_ca.peak_cyt_nM)
	rows.append((
		'SOCE differential',
		'|peak(+) − peak(−)| ≥ 100 nM',
		f'{soce_diff:.0f} nM',
		soce_diff >= 100.0,
	))
	# Peak amplitude bands (lab-book ±30% of Dolan)
	with_lo = ref['with_extracellular_ca']['peak_cyt_nM']['value'] * 0.7
	with_hi = ref['with_extracellular_ca']['peak_cyt_nM']['value'] * 1.3
	rows.append((
		'Peak in Dolan ±30% (+Ca_ex)',
		f'{with_lo:.0f}–{with_hi:.0f} nM',
		f'{with_ca.peak_cyt_nM:.0f} nM',
		with_lo <= with_ca.peak_cyt_nM <= with_hi,
	))
	no_lo = ref['without_extracellular_ca']['peak_cyt_nM']['value'] * 0.7
	no_hi = ref['without_extracellular_ca']['peak_cyt_nM']['value'] * 1.3
	rows.append((
		'Peak in Dolan ±30% (−Ca_ex)',
		f'{no_lo:.0f}–{no_hi:.0f} nM',
		f'{no_ca.peak_cyt_nM:.0f} nM',
		no_lo <= no_ca.peak_cyt_nM <= no_hi,
	))
	return rows


def _draw_cyt_panel(ax, with_ca: _ConditionTrace, no_ca: _ConditionTrace,
		ref: dict) -> None:
	# Dolan +Ca peak band
	w = ref['with_extracellular_ca']
	ax.axhspan(w['peak_cyt_nM']['range_lo'], w['peak_cyt_nM']['range_hi'],
		alpha=0.10, color='tab:blue',
		label=f'Dolan +Ca_ex peak ({w["peak_cyt_nM"]["range_lo"]}–{w["peak_cyt_nM"]["range_hi"]} nM)')
	if w['plateau_cyt_nM']['value'] is not None:
		ax.axhspan(w['plateau_cyt_nM']['range_lo'], w['plateau_cyt_nM']['range_hi'],
			alpha=0.07, color='tab:cyan',
			label=f'Dolan +Ca_ex plateau ({w["plateau_cyt_nM"]["range_lo"]}–{w["plateau_cyt_nM"]["range_hi"]} nM)')
	# Dolan −Ca peak band
	n = ref['without_extracellular_ca']
	ax.axhspan(n['peak_cyt_nM']['range_lo'], n['peak_cyt_nM']['range_hi'],
		alpha=0.10, color='tab:red',
		label=f'Dolan −Ca_ex peak ({n["peak_cyt_nM"]["range_lo"]}–{n["peak_cyt_nM"]["range_hi"]} nM)')

	ax.plot(with_ca.t, with_ca.ca_cyt_nM, color=with_ca.colour, linewidth=2.0,
		label=f'sim {with_ca.label}')
	ax.plot(no_ca.t, no_ca.ca_cyt_nM, color=no_ca.colour, linewidth=2.0,
		linestyle='--', label=f'sim {no_ca.label}')
	ax.axhline(200.0, color='grey', linewidth=0.8, linestyle=':',
		label='Active threshold (200 nM)')

	ax.set_ylabel('[Ca²⁺]_cyt (nM)', fontsize=10)
	ax.set_title('Panel 1 — Cytosolic Ca²⁺ transient: with vs without extracellular Ca²⁺  '
		f'(Dolan 2014 Fig. 4)', fontsize=10, fontweight='bold')
	ax.set_ylim(bottom=0)
	ax.legend(loc='upper right', fontsize=7)
	ax.grid(True, alpha=0.3)


def _draw_dts_panel(ax, with_ca: _ConditionTrace, no_ca: _ConditionTrace,
		ref: dict) -> None:
	w = ref['with_extracellular_ca']
	n = ref['without_extracellular_ca']
	ax.axhspan(w['dts_min_uM']['range_lo'], w['dts_min_uM']['range_hi'],
		alpha=0.10, color='tab:blue',
		label=f'Dolan +Ca_ex DTS min ({w["dts_min_uM"]["range_lo"]}–{w["dts_min_uM"]["range_hi"]} µM)')
	ax.axhspan(n['dts_min_uM']['range_lo'], n['dts_min_uM']['range_hi'],
		alpha=0.10, color='tab:red',
		label=f'Dolan −Ca_ex DTS min ({n["dts_min_uM"]["range_lo"]}–{n["dts_min_uM"]["range_hi"]} µM)')

	ax.plot(with_ca.t, with_ca.ca_dts_uM, color=with_ca.colour, linewidth=2.0,
		label=f'sim {with_ca.label}')
	ax.plot(no_ca.t, no_ca.ca_dts_uM, color=no_ca.colour, linewidth=2.0,
		linestyle='--', label=f'sim {no_ca.label}')
	ax.axhline(250.0, color='grey', linewidth=0.8, linestyle=':',
		label='Dolan resting (250 µM)')

	ax.set_ylabel('[Ca²⁺]_dts (µM)', fontsize=10)
	ax.set_title('Panel 2 — DTS store Ca²⁺  (depletion drives SOCE)',
		fontsize=10, fontweight='bold')
	ax.set_ylim(bottom=0)
	ax.legend(loc='upper right', fontsize=7)
	ax.grid(True, alpha=0.3)


def _draw_criteria_panel(ax, rows: list[tuple[str, str, str, bool]]) -> None:
	ax.axis('off')
	ax.set_title('Panel 3 — Phase 3 acceptance criteria  '
		'(Dolan 2014 Fig. 3B + lab-book Phase 3)',
		fontsize=10, fontweight='bold', loc='left')

	col_labels = ['Criterion', 'Rule', 'Measured', 'Result']
	cell_text = []
	cell_colours = []
	for (criterion, rule, measured, passed) in rows:
		cell_text.append([criterion, rule, measured,
			'PASS' if passed else 'FAIL'])
		row_colour = '#dfead6' if passed else '#f6d6d2'
		cell_colours.append([row_colour] * len(col_labels))

	tbl = ax.table(
		cellText=cell_text,
		colLabels=col_labels,
		cellColours=cell_colours,
		colWidths=[0.30, 0.32, 0.18, 0.10],
		loc='upper center',
		cellLoc='left',
	)
	tbl.auto_set_font_size(False)
	tbl.set_fontsize(9)
	tbl.scale(1.0, 1.6)


def make_phase3_plot(
		with_ca_simout: str,
		no_ca_simout: str,
		plot_out_path: str,
		reference_json: Optional[str] = None,
		) -> dict:
	"""Build and save the Phase 3 validation figure.

	Args:
		with_ca_simout: simOut/ directory of the +Ca_ex run.
		no_ca_simout:   simOut/ directory of the −Ca_ex run.
		plot_out_path:  output file path (PNG / PDF). Parent dir is created.
		reference_json: optional override for the Dolan reference JSON path.

	Returns:
		Dict summary of measured peaks and acceptance-criteria pass/fail —
		suitable for direct serialisation alongside the figure.
	"""
	ref = _load_reference(reference_json or REFERENCE_JSON)
	with_ca = _ConditionTrace.from_simout(with_ca_simout, '+Ca_ex (1.2 mM)', 'tab:blue')
	no_ca   = _ConditionTrace.from_simout(no_ca_simout,   '−Ca_ex (EDTA)',   'tab:red')
	rows    = _eval_criteria(with_ca, no_ca, ref)

	fig, axes = plt.subplots(3, 1, figsize=(11, 13),
		gridspec_kw={'height_ratios': [3, 3, 2]})
	fig.subplots_adjust(hspace=0.40)

	_draw_cyt_panel(axes[0], with_ca, no_ca, ref)
	_draw_dts_panel(axes[1], with_ca, no_ca, ref)
	_draw_criteria_panel(axes[2], rows)
	axes[1].set_xlabel('time (s)', fontsize=10)

	fig.suptitle(
		'Phase 3 validation — platelet-wcm v0.2 vs Dolan & Diamond 2014 Fig. 4',
		fontsize=12, fontweight='bold', y=0.995)

	os.makedirs(os.path.dirname(plot_out_path), exist_ok=True)
	fig.savefig(plot_out_path, dpi=150, bbox_inches='tight')
	plt.close(fig)

	summary = {
		'with_extracellular_ca': {
			'peak_cyt_nM': with_ca.peak_cyt_nM,
			'peak_t_s':    with_ca.peak_t_s,
			'dts_min_uM':  with_ca.dts_min_uM,
		},
		'without_extracellular_ca': {
			'peak_cyt_nM': no_ca.peak_cyt_nM,
			'peak_t_s':    no_ca.peak_t_s,
			'dts_min_uM':  no_ca.dts_min_uM,
		},
		'soce_differential_nM': abs(with_ca.peak_cyt_nM - no_ca.peak_cyt_nM),
		'criteria': [
			{'criterion': c, 'rule': r, 'measured': m, 'passed': p}
			for (c, r, m, p) in rows
		],
		'plot_path': plot_out_path,
	}
	return summary
