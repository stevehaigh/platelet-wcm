"""
Calcium demo figure (TUI) — cytosolic Ca²⁺, DTS store + SOCE, IP₃.

Focused figure for the calcium-centred TUI demos: the ±extracellular-Ca²⁺ /
EDTA contrast (demo 1) and the MCU-knockout result (demo 2). Grey baseline
overlay via PLATELET_BASELINE_SIMOUT (see _demo_common).

Usage: analysisPlatelet.py --plot demo_calcium [sim_dir]
"""

import matplotlib as mp
mp.use('Agg')
from matplotlib import pyplot as plt

from models.platelet.analysis import singleAnalysisPlot
from models.platelet.analysis.single._demo_common import (
	combined_legend, make_draw, provenance_footnote, resolve_baseline, suptitle)
from wholecell.analysis.analysis_tools import exportFigure


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	"""Cytosolic Ca²⁺, DTS store + SOCE, and IP₃ with baseline overlay."""

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		del simDataFile, validationDataFile
		base = resolve_baseline(simOutDir)
		draw = make_draw(simOutDir, base)

		fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
		fig.subplots_adjust(wspace=0.34, top=0.84, bottom=0.20)

		# Cytosolic Ca²⁺
		ax = axes[0]
		draw(ax, 'CalciumTrace', 'ca_cyt_nM', 'tab:blue',
			r'cytosolic $[\mathrm{Ca}^{2+}]$')
		ax.axhline(100.0, color='grey', lw=0.8, ls=':', label='100 nM rest')
		ax.set_ylabel(r'cytosolic $[\mathrm{Ca}^{2+}]$  (nM)')
		ax.set_ylim(bottom=0)
		ax.set_title('Cytosolic $\\mathrm{Ca}^{2+}$',
			fontsize=10, fontweight='bold')
		ax.legend(loc='best', fontsize=8)

		# DTS store + store-operated entry (SOCE on the right axis; baseline
		# overlay stays on the DTS trace only — two grey lines on two scales
		# would be ambiguous).
		ax = axes[1]
		ax_r = ax.twinx()
		draw(ax, 'CalciumTrace', 'ca_dts_uM', 'tab:red', 'DTS store')
		ax.set_ylabel(r'DTS store $[\mathrm{Ca}^{2+}]$  ($\mu$M)',
			color='tab:red')
		ax.tick_params(axis='y', labelcolor='tab:red')
		ax.set_ylim(bottom=0)
		draw(ax_r, 'CalciumTrace', 'soce_flux_nMs', 'darkorange',
			'SOCE influx', lw=1.4, ls='--', baseline=False)
		ax_r.set_ylabel('SOCE influx  (nM/s)', color='darkorange')
		ax_r.tick_params(axis='y', labelcolor='darkorange')
		ax.set_title('DTS store + store-operated entry  ($\\pm$Ca$_{ex}$ / EDTA)',
			fontsize=10, fontweight='bold')
		combined_legend(ax, ax_r)

		# IP₃
		ax = axes[2]
		draw(ax, 'CalciumTrace', 'ip3_nM', 'tab:green', r'$[\mathrm{IP}_3]$')
		ax.set_ylabel(r'$[\mathrm{IP}_3]$  (nM)')
		ax.set_ylim(bottom=0)
		ax.set_title(r'$\mathrm{IP}_3$  (upstream signal)',
			fontsize=10, fontweight='bold')
		ax.legend(loc='best', fontsize=8)

		for a in axes:
			a.set_xlabel('time (s)')
			a.grid(True, alpha=0.3)
		suptitle(fig, 'Calcium — cytosol, store / SOCE, IP$_3$', base)
		provenance_footnote(fig, 'demo_calcium', metadata)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
