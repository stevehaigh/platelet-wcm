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
	combined_legend, make_draw, resolve_baseline, suptitle)
from wholecell.analysis.analysis_tools import exportFigure


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	"""Cytosolic Ca²⁺, DTS store + SOCE, and IP₃ with baseline overlay."""

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		del simDataFile, validationDataFile
		base = resolve_baseline(simOutDir)
		draw = make_draw(simOutDir, base)

		fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
		fig.subplots_adjust(wspace=0.34, top=0.84, bottom=0.16)

		# Cytosolic Ca²⁺
		ax = axes[0]
		draw(ax, 'CalciumTrace', 'ca_cyt_nM', 'tab:blue',
			r'$[\mathrm{Ca}^{2+}]_\mathrm{cyt}$')
		ax.axhline(100.0, color='grey', lw=0.8, ls=':', label='100 nM rest')
		ax.set_ylabel(r'$[\mathrm{Ca}^{2+}]_\mathrm{cyt}$ (nM)')
		ax.set_ylim(bottom=0)
		ax.set_title('Cytosolic $\\mathrm{Ca}^{2+}$  (MCU KO raises it)',
			fontsize=10, fontweight='bold')
		ax.legend(loc='best', fontsize=8)

		# DTS store + SOCE flux
		ax = axes[1]
		ax_r = ax.twinx()
		draw(ax, 'CalciumTrace', 'ca_dts_uM', 'tab:red', 'DTS store')
		ax.set_ylabel(r'$[\mathrm{Ca}^{2+}]_\mathrm{DTS}$ ($\mu$M)',
			color='tab:red')
		ax.tick_params(axis='y', labelcolor='tab:red')
		ax.set_ylim(bottom=0)
		draw(ax_r, 'CalciumTrace', 'soce_flux_nMs', 'darkorange', 'SOCE flux',
			lw=1.4, ls='--')
		ax_r.set_ylabel('SOCE flux (nM/s)', color='darkorange')
		ax_r.tick_params(axis='y', labelcolor='darkorange')
		ax.set_title('DTS store + SOCE  ($\\pm$Ca$_{ex}$ / EDTA)',
			fontsize=10, fontweight='bold')
		combined_legend(ax, ax_r)

		# IP₃
		ax = axes[2]
		draw(ax, 'CalciumTrace', 'ip3_nM', 'tab:green', r'$\mathrm{IP}_3$')
		ax.set_ylabel(r'$\mathrm{IP}_3$ (nM)')
		ax.set_ylim(bottom=0)
		ax.set_title(r'$\mathrm{IP}_3$  (upstream readout)',
			fontsize=10, fontweight='bold')
		ax.legend(loc='best', fontsize=8)

		for a in axes:
			a.set_xlabel('Time (s)')
			a.grid(True, alpha=0.3)
		suptitle(fig, 'Calcium — cytosol, store / SOCE, IP$_3$', base)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
