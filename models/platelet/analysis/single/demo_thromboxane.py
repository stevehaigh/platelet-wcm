"""
Thromboxane demo figure (TUI) — TXA₂ and TXB₂.

Focused figure for the aspirin demo (demo 4): COX-1 knockout (aspirin) abolishes
TXA₂ synthesis; TXA₂ decays to the stable TXB₂ ELISA metabolite. Grey baseline
overlay via PLATELET_BASELINE_SIMOUT (e.g. full-loop run vs aspirin).

Usage: analysisPlatelet.py --plot demo_thromboxane [sim_dir]
"""

import matplotlib as mp
mp.use('Agg')
from matplotlib import pyplot as plt

from models.platelet.analysis import singleAnalysisPlot
from models.platelet.analysis.single._demo_common import (
	combined_legend, make_draw, provenance_footnote, resolve_baseline, suptitle)
from wholecell.analysis.analysis_tools import exportFigure


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	"""TXA₂ (left) and TXB₂ (right) with baseline overlay on TXA₂."""

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		del simDataFile, validationDataFile
		base = resolve_baseline(simOutDir)
		draw = make_draw(simOutDir, base)

		fig, ax = plt.subplots(figsize=(9, 5.5))
		fig.subplots_adjust(top=0.9, bottom=0.16, right=0.86)
		ax_r = ax.twinx()

		# Active TXA₂ (left) carries the baseline overlay; the stable TXB₂
		# metabolite (right) is current-run only — a second grey line on the
		# other axis would be ambiguous.
		draw(ax, 'ThromboxaneTrace', 'txa2_uM', 'tab:brown',
			r'$[\mathrm{TXA}_2]$ (active)')
		ax.set_ylabel(r'$[\mathrm{TXA}_2]$  ($\mu$M)', color='tab:brown')
		ax.tick_params(axis='y', labelcolor='tab:brown')
		ax.set_ylim(bottom=0)
		ax.set_xlabel('time (s)')
		ax.grid(True, alpha=0.3)

		draw(ax_r, 'ThromboxaneTrace', 'txb2', 'tab:olive',
			r'$\mathrm{TXB}_2$ (stable metabolite)', lw=1.4, ls='--',
			baseline=False)
		ax_r.set_ylabel(r'$\mathrm{TXB}_2$  (molecules)', color='tab:olive')
		ax_r.tick_params(axis='y', labelcolor='tab:olive')
		ax_r.set_ylim(bottom=0)  # molecule counts are ≥ 0

		ax.set_title('Thromboxane synthesis  (aspirin / COX-1 KO $\\to$ 0)',
			fontsize=11, fontweight='bold')
		combined_legend(ax, ax_r)
		suptitle(fig, 'Thromboxane — $\\mathrm{TXA}_2$ synthesis and '
			'$\\mathrm{TXB}_2$ metabolite', base)
		provenance_footnote(fig, 'demo_thromboxane', metadata)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
