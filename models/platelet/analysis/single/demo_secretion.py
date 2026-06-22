"""
Secretion demo figure (TUI) — autocrine ADP[e] and released cargo fractions.

Focused figure for the second-wave demo (demo 4): granule release and the
autocrine ADP[e] that feeds back onto P2Y1/P2Y12. Left axis: secreted ADP[e]
(µM, baseline overlay); right axis: released fractions (ADP, 5-HT, P-selectin,
current run only). Grey baseline overlay via PLATELET_BASELINE_SIMOUT.

Usage: analysisPlatelet.py --plot demo_secretion [sim_dir]
"""

import matplotlib as mp
mp.use('Agg')
from matplotlib import pyplot as plt

from models.platelet.analysis import singleAnalysisPlot
from models.platelet.analysis.single._demo_common import (
	combined_legend, make_draw, resolve_baseline, suptitle)
from wholecell.analysis.analysis_tools import exportFigure


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	"""Autocrine ADP[e] (left) and released fractions (right) with overlay."""

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		del simDataFile, validationDataFile
		base = resolve_baseline(simOutDir)
		draw = make_draw(simOutDir, base)

		fig, ax = plt.subplots(figsize=(9, 5.5))
		fig.subplots_adjust(top=0.9, bottom=0.12, right=0.86)

		draw(ax, 'SecretionTrace', 'adp_e_uM', 'tab:cyan',
			r'autocrine $[\mathrm{ADP}]_e$')
		ax.set_ylabel(r'secreted $[\mathrm{ADP}]_e$ ($\mu$M)', color='tab:cyan')
		ax.tick_params(axis='y', labelcolor='tab:cyan')
		ax.set_ylim(bottom=0)
		ax.set_xlabel('Time (s)')
		ax.grid(True, alpha=0.3)

		ax_r = ax.twinx()
		# Released fractions (current run only) share the right axis (0–1); the
		# baseline overlay stays on the primary ADP[e] trace for legibility.
		for col, c, lab in (
				('adp_released_frac', 'tab:red', 'ADP rel.'),
				('serotonin_released_frac', 'tab:orange', '5-HT rel.'),
				('pselectin_surface_frac', 'tab:green', 'P-selectin')):
			draw(ax_r, 'SecretionTrace', col, c, lab, lw=1.2, ls='--',
				baseline=False)
		ax_r.set_ylabel('released fraction')
		ax_r.set_ylim(0, 1)

		ax.set_title('Granule secretion  (2nd wave; autocrine ADP)',
			fontsize=11, fontweight='bold')
		combined_legend(ax, ax_r, loc='upper left')
		suptitle(fig, 'Secretion — autocrine ADP + cargo release', base)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
