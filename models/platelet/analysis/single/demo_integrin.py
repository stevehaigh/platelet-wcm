"""
Integrin demo figure (TUI) — PAC-1 activation, PKC×Ca gate, PKA brake.

Focused figure for the integrin demo (demo 3, clopidogrel) and the Glanzmann
foil. PAC-1 active fraction is the per-cell flow-cytometry readout; the PKA
brake (right axis) shows the cAMP-driven mechanism by which clopidogrel lowers
activation with the integrin intact. Grey baseline overlay via
PLATELET_BASELINE_SIMOUT.

Usage: analysisPlatelet.py --plot demo_integrin [sim_dir]
"""

import matplotlib as mp
mp.use('Agg')
from matplotlib import pyplot as plt

from models.platelet.analysis import singleAnalysisPlot
from models.platelet.analysis.single._demo_common import (
	combined_legend, make_draw, resolve_baseline, suptitle)
from wholecell.analysis.analysis_tools import exportFigure


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	"""PAC-1 active fraction (+ gate, PKA brake) with baseline overlay."""

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		del simDataFile, validationDataFile
		base = resolve_baseline(simOutDir)
		draw = make_draw(simOutDir, base)

		fig, ax = plt.subplots(figsize=(9, 5.5))
		fig.subplots_adjust(top=0.9, bottom=0.12, right=0.85)

		draw(ax, 'IntegrinTrace', 'active_frac', 'tab:purple', 'PAC-1 active')
		# Gate is current-run context only (no baseline) to keep the panel clean.
		draw(ax, 'IntegrinTrace', 'integrin_gate', 'tab:purple', 'PKC×Ca gate',
			lw=1.0, ls=':', baseline=False)
		ax.set_ylabel('fraction')
		ax.set_ylim(0, 1)
		ax.set_xlabel('Time (s)')
		ax.grid(True, alpha=0.3)

		ax_r = ax.twinx()
		# Brake gets the baseline overlay too: the control's dis-inhibited brake
		# (>1) vs the drug's re-engaged brake (≈1) IS the clopidogrel mechanism.
		draw(ax_r, 'IntegrinTrace', 'pka_brake', 'tab:brown',
			'PKA brake (≥1)', lw=1.4, ls='--')
		ax_r.set_ylabel('PKA brake (≥1; higher = more activation)',
			color='tab:brown')
		ax_r.tick_params(axis='y', labelcolor='tab:brown')

		ax.set_title(
			'PAC-1 integrin activation  (clopidogrel lowers; Glanzmann $\\to$ 0)',
			fontsize=11, fontweight='bold')
		combined_legend(ax, ax_r)
		suptitle(fig, 'Integrin — PAC-1 activation', base)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
